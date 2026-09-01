# -*- coding: utf-8 -*-
"""
60_outcome_discovery.py - 器官孟德尔化: 结局 GWAS 筛选 (P1-7)

对一组血细胞/器官功能性状, 从 GWAS Catalog REST API 找候选研究,
在 EBI FTP harmonised 目录定位全量汇总统计文件, 要求带 .tbi 索引
(远程 tabix 可查), 优先欧洲大样本。输出候选表供人工确认。
"""
import json
import re
import subprocess
import time
from collections import defaultdict

BASE = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"


def curl_json(url, timeout=90):
    for i in range(3):
        try:
            out = subprocess.run(["curl", "-s", "--ssl-no-revoke", "-m", str(timeout), url],
                                 capture_output=True, text=True, timeout=timeout + 10).stdout
            return json.loads(out)
        except Exception:
            time.sleep(3)
    return None


def head_status(url, timeout=60):
    for i in range(3):
        try:
            out = subprocess.run(["curl", "-sI", "--ssl-no-revoke", "-m", str(timeout), "-o", "/dev/null",
                                  "-w", "%{http_code}", url],
                                 capture_output=True, text=True, timeout=timeout + 10).stdout.strip()
            if out:
                return out
        except Exception:
            time.sleep(2)
    return "ERR"


# harmonised 清单 (本地)
harm = {}
with open("results/_tmp_organ/harmonised_list.txt") as f:
    for line in f:
        line = line.strip()
        if not line.startswith("./"):
            continue
        acc = line.split("/")[2]
        harm[acc] = line[2:]

# 性状 -> 候选 diseaseTrait 名称 (GWAS Catalog 中常用写法)
TRAIT_FAMILIES = {
    "WBC":        ["White blood cell count", "White blood cell count measurement"],
    "Neutrophil": ["Neutrophil count", "Neutrophil count measurement", "Absolute neutrophil count"],
    "Lymphocyte": ["Lymphocyte count", "Lymphocyte count measurement", "Absolute lymphocyte count"],
    "Monocyte":   ["Monocyte count", "Monocyte count measurement", "Absolute monocyte count"],
    "Eosinophil": ["Eosinophil count", "Eosinophil count measurement", "Absolute eosinophil count"],
    "Basophil":   ["Basophil count", "Basophil count measurement", "Absolute basophil count"],
    "RBC":        ["Red blood cell count", "Red blood cell count measurement", "Erythrocyte count"],
    "Hemoglobin": ["Hemoglobin measurement", "Hemoglobin level", "Hemoglobin concentration"],
    "Hematocrit": ["Hematocrit", "Hematocrit measurement"],
    "MCV":        ["Mean corpuscular volume"],
    "MCH":        ["Mean corpuscular hemoglobin"],
    "MCHC":       ["Mean corpuscular hemoglobin concentration"],
    "Platelet":   ["Platelet count", "Platelet count measurement"],
    "MPV":        ["Mean platelet volume"],
    "Reticulocyte": ["Reticulocyte count", "Reticulocyte percentage"],
    "eGFR":       ["Estimated glomerular filtration rate", "Glomerular filtration rate estimated",
                   "Estimated glomerular filtration rate measurement"],
    "Creatinine": ["Creatinine measurement", "Serum creatinine measurement", "Creatinine level",
                   "Serum creatinine level"],
    "Urate":      ["Urate level", "Serum urate level", "Uric acid level", "Uric acid measurement"],
    "Urea":       ["Urea level", "Blood urea nitrogen level", "Blood urea nitrogen measurement"],
    "ALT":        ["Alanine aminotransferase measurement", "Alanine aminotransferase level",
                   "Alanine aminotransferase activity"],
    "AST":        ["Aspartate aminotransferase measurement", "Aspartate aminotransferase level",
                   "Aspartate aminotransferase activity"],
    "GGT":        ["Gamma-glutamyl transferase level", "Gamma-glutamyl transferase measurement",
                   "Gamma glutamyl transferase level"],
    "Albumin":    ["Albumin measurement", "Serum albumin level", "Serum total protein measurement"],
    "CRP":        ["C-reactive protein measurement", "C-reactive protein level",
                   "High-sensitivity C-reactive protein measurement", "C reactive protein measurement"],
    "SBP":        ["Systolic blood pressure", "Systolic blood pressure measurement"],
    "DBP":        ["Diastolic blood pressure", "Diastolic blood pressure measurement"],
    "PP":         ["Pulse pressure", "Mean arterial pressure"],
    "IgG":        ["IgG measurement", "Immunoglobulin G measurement"],
    "TotalIgE":   ["Total IgE measurement", "Total immunoglobulin E measurement"],
}

seen = set()
results = []
for fam, names in TRAIT_FAMILIES.items():
    fam_best = []
    for nm in names:
        page = 0
        while True:
            d = curl_json("https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByDiseaseTrait"
                          "?diseaseTrait=%s&size=100&page=%d" % (nm.replace(" ", "%20"), page))
            if not d or not d.get("_embedded"):
                break
            for s in d["_embedded"]["studies"]:
                acc = s["accessionId"]
                if acc in seen:
                    continue
                seen.add(acc)
                init = s.get("initialSampleSize", "")
                rep = s.get("diseaseTrait", {}).get("trait", "")
                # 解析最大样本数和是否欧洲
                nums = [int(x.replace(",", "")) for x in re.findall(r"[\d,]{3,}", init)]
                nmax = max(nums) if nums else 0
                eur = "European" in init or "British" in init or "Iceland" in init
                if acc not in harm:
                    continue  # 无 harmonised 全量文件
                fam_best.append({"family": fam, "accession": acc, "trait": rep,
                                 "n": nmax, "eur": eur, "sample": init[:70],
                                 "harmonised": harm[acc]})
            tp = d.get("page", {})
            page += 1
            if page >= (tp.get("totalPages") or 1):
                break
    # 只保留欧洲 & n>=100k, 检查 tbi
    fam_best.sort(key=lambda x: -x["n"])
    checked = 0
    for cand in fam_best:
        if not cand["eur"] or cand["n"] < 100000:
            continue
        url = "%s/%s" % (BASE, cand["harmonised"] + ".tbi")
        code = head_status(url)
        cand["tbi"] = (code == "200")
        if cand["tbi"]:
            print("  [OK ] %-14s %-12s n=%-8d %s" % (fam, cand["accession"], cand["n"], cand["trait"]))
            results.append(cand)
            break   # 每性状取 tbi 可用的最大欧洲研究
        else:
            print("  [no ] %-14s %-12s n=%-8d (无 tbi)" % (fam, cand["accession"], cand["n"]))
            checked += 1
            if checked >= 8:
                break

with open("results/_tmp_organ/outcome_studies.tsv", "w") as f:
    f.write("family\taccession\ttrait\tn_eur\ttrait_reported\tsample\tftp_path\n")
    for r in results:
        f.write("%s\t%s\t%s\t%d\t%s\t%s\t%s\n" % (r["family"], r["accession"], r["trait"], r["n"], r["trait"], r["sample"], r["harmonised"]))
print("\n已选 %d 个性状, 见 results/_tmp_organ/outcome_studies.tsv" % len(results))
