# -*- coding: utf-8 -*-
"""
61_extract_outcome_effects.py - 器官孟德尔化: 远程 tabix 提取 IV 结局效应 (P1-7)

对 outcome_studies.tsv 中每个结局 GWAS, 用远程 tabix 读取 EBI FTP harmonised
文件, 提取 160 个 cis-meQTL IV 位点的等位基因匹配效应量。

输出: results/_tmp_organ/organ_outcome_effects.csv
"""
import csv
import json
import os
import re
import subprocess
import sys
import time
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from rtabix import RemoteTabix

os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

WALD = "results/GSE282512_mr_wald_snps.csv"
MQTL = "results/_tmp_smr/mqtl_full7.csv"
STUDIES = "results/_tmp_organ/outcome_studies.tsv"
OUT = "results/_tmp_organ/organ_outcome_effects.csv"
BASE = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"


def log(m):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), m), flush=True)


# ---------------------------------------------------------------- 1. IV 表 (含等位基因)
log("载入 IV 与等位基因")
# wald: 每 (cpg,pheno,rsid) 一行; 部分行缺 rsid/pos38 (Tyrmi 按位置匹配), 取并集
iv_merge = {}   # (cpg, pos37) -> dict
with open(WALD) as f:
    for r in csv.DictReader(f):
        if not r["pos37"]:
            continue
        key = (r["cpg"], int(r["pos37"]))
        d = iv_merge.setdefault(key, {"cpg": r["cpg"], "pos37": int(r["pos37"]),
                                      "rsid": "", "pos38": None, "region_id": r["region_id"],
                                      "beta_m": float(r["beta_m"]), "se_m": float(r["se_m"]),
                                      "p_m": float(r["p_m"])})
        if r["rsid"] and not d["rsid"]:
            d["rsid"] = r["rsid"]
        if r["pos38"] and not d["pos38"]:
            d["pos38"] = int(r["pos38"])
log("  IV (cpg,pos37): %d" % len(iv_merge))

# mqtl_full7: cpg,"chr:pos:SNP",beta,se,p,allele1,allele2,freq_a1,cistrans,clumped
alleles = {}   # (cpg, chr:pos37) -> (a1, a2, freq)
with open(MQTL) as f:
    for row in csv.reader(f):
        if len(row) >= 9 and row[0].startswith("cg"):
            loc = row[1].replace('"', "").replace(":SNP", "")
            try:
                alleles[(row[0], loc)] = (row[5], row[6], float(row[7]) if row[7] else None)
            except ValueError:
                pass

iv_list = []
for (cpg, pos37), d in iv_merge.items():
    iv_list.append({"cpg": cpg, "rsid": d["rsid"], "pos37": pos37,
                    "pos38": d["pos38"] or None, "region_id": d["region_id"],
                    "beta_m": d["beta_m"], "se_m": d["se_m"], "p_m": d["p_m"], "chr": ""})
# 染色体: 按 region_id 映射 (106528=chr3, 12268=chr6, 84604=chr17)
REG_CHR = {"106528": "3", "12268": "6", "84604": "17"}
for iv in iv_list:
    iv["chr"] = REG_CHR.get(iv["region_id"], "")
# 补 pos38: 从 58a lift 表
lift = {}
with open("results/_tmp_smr/lift_positions.csv") as f:
    for r in csv.DictReader(f):
        lift[(r["chr"], int(r["pos"]))] = int(r["pos38"])
n_fill = 0
for iv in iv_list:
    if iv["pos38"] is None and (iv["chr"], iv["pos37"]) in lift:
        iv["pos38"] = lift[(iv["chr"], iv["pos37"])]
        n_fill += 1
    if iv["pos38"] is None:
        iv["pos38"] = iv["pos37"]   # 最后回退 (仅当 lift 缺失)
log("  pos38 从 lift 表补齐: %d" % n_fill)
# 等位基因
n_al = 0
for iv in iv_list:
    k = (iv["cpg"], "chr%s:%d" % (iv["chr"], iv["pos37"]))
    if k in alleles:
        iv["a1"], iv["a2"], iv["freq_a1"] = alleles[k]
        n_al += 1
    else:
        iv["a1"], iv["a2"], iv["freq_a1"] = "", "", None
log("  等位基因匹配: %d/%d" % (n_al, len(iv_list)))

# 按染色体聚类区域 (gap>200kb 分簇)
regions = defaultdict(list)   # chr -> list of (beg,end)
by_chr = defaultdict(list)
for iv in iv_list:
    by_chr[iv["chr"]].append(iv)
region_list = []
for ch, items in by_chr.items():
    items.sort(key=lambda x: x["pos38"])
    cl = []
    for iv in items:
        p37, p38 = iv["pos37"], iv["pos38"]
        if cl and (p38 - cl[-1]["pos38"]) > 200000:
            region_list.append((ch, cl)); cl = []
        cl.append(iv)
    if cl:
        region_list.append((ch, cl))
log("  IV 区域簇: %d" % len(region_list))

# ---------------------------------------------------------------- 2. 逐结局提取
studies = []
with open(STUDIES) as f:
    rd = csv.DictReader(f, delimiter="\t")
    for r in rd:
        studies.append(r)
log("结局研究: %d" % len(studies))

out_rows = []
for st in studies:
    fam, acc, trait = st["family"], st["accession"], st["trait"]
    url = "%s/%s" % (BASE, st["ftp_path"])
    # meta yaml: build + 样本
    build = "GRCh38"
    try:
        y = subprocess.run(["curl", "-s", "--ssl-no-revoke", "-m", "60", url + "-meta.yaml"],
                           capture_output=True, text=True, timeout=70).stdout
        m = re.search(r"genome_assembly:\s*(\S+)", y)
        if m:
            build = m.group(1)
        eur_n = sum(int(x) for x in re.findall(r"sample_size:\s*(\d+)", y)[:2])  # 粗略
    except Exception:
        eur_n = int(st["n_eur"]) if st.get("n_eur") else 0
    try:
        tb = RemoteTabix(url, cache_dir="results/_tmp_organ/tbi_cache")
    except Exception as e:
        log("  [%s %s] tbi 失败: %s" % (fam, acc, e))
        continue
    hdr = tb.header()
    if not hdr:
        log("  [%s %s] header 失败" % (fam, acc))
        continue
    cols = hdr.split("\t")
    def cidx(*names):
        for nm in names:
            if nm in cols:
                return cols.index(nm)
        return None
    i_ch = cidx("chromosome", "chr", "CHR")
    i_pos = cidx("base_pair_location", "pos", "position", "BP")
    i_ea = cidx("effect_allele", "EA", "alleles") 
    i_oa = cidx("other_allele", "NEA", "OA")
    i_beta = cidx("beta", "BETA", "b")
    i_se = cidx("standard_error", "SE")
    i_p = cidx("p_value", "P", "pvalue", "PVALUE")
    i_af = cidx("effect_allele_frequency", "EAF", "af", "AF", "frequency")
    i_rsid = cidx("variant_id", "rsid", "rs_id", "RSID", "snp")
    if i_ch is None or i_pos is None or i_beta is None:
        log("  [%s %s] 列缺失: %s" % (fam, acc, hdr[:150]))
        continue

    # rsid 索引
    rows_by_rsid = {}
    rows_by_pos = {}
    nrow = 0
    for ch, cl in region_list:
        pmin = min(x["pos38"] if build == "GRCh38" else x["pos37"] for x in cl) - 2000
        pmax = max(x["pos38"] if build == "GRCh38" else x["pos37"] for x in cl) + 2000
        try:
            lines = tb.fetch(ch, pmin - 1, pmax)
        except Exception as e:
            log("  [%s %s] fetch %s 失败: %s" % (fam, acc, ch, e))
            continue
        for ln in lines:
            fld = ln.split("\t")
            if len(fld) <= max(x for x in [i_ch, i_pos, i_beta, i_se, i_p, i_af, i_rsid, i_ea, i_oa] if x is not None):
                continue
            try:
                pos = int(fld[i_pos])
            except ValueError:
                continue
            rs = fld[i_rsid] if i_rsid is not None else "NA"
            rec = {"pos": pos, "rsid": rs,
                   "ea": fld[i_ea] if i_ea is not None else "",
                   "oa": fld[i_oa] if i_oa is not None else "",
                   "beta": fld[i_beta],
                   "se": fld[i_se] if i_se is not None else "",
                   "p": fld[i_p] if i_p is not None else "",
                   "af": fld[i_af] if i_af is not None else ""}
            nrow += 1
            if rs and rs != "NA":
                rows_by_rsid[rs] = rec
            rows_by_pos[(ch, pos)] = rec
    # 匹配 IV
    n_match = 0
    for iv in iv_list:
        if iv["chr"] == "":
            continue
        rec = rows_by_rsid.get(iv["rsid"])
        if rec is None:
            pos_key = (iv["chr"], iv["pos38"] if build == "GRCh38" else iv["pos37"])
            rec = rows_by_pos.get(pos_key)
        if rec is None:
            continue
        try:
            beta = float(rec["beta"])
        except ValueError:
            continue
        se = float(rec["se"]) if rec["se"] else None
        out_rows.append({
            "family": fam, "accession": acc, "trait": trait, "build": build,
            "cpg": iv["cpg"], "region_id": iv["region_id"], "rsid": iv["rsid"],
            "pos": rec["pos"], "ea": rec["ea"], "oa": rec["oa"], "af": rec["af"],
            "beta_o": beta, "se_o": se, "p_o": rec["p"],
            "a1": iv["a1"], "a2": iv["a2"], "freq_a1": iv["freq_a1"],
            "beta_m": iv["beta_m"], "se_m": iv["se_m"], "p_m": iv["p_m"]})
        n_match += 1
    log("  [%-14s %s] %s | 区域行 %d | IV 匹配 %d/%d" % (fam, acc, build, nrow, n_match, len(iv_list)))

with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
    w.writeheader()
    w.writerows(out_rows)
log("完成: %d 行 -> %s" % (len(out_rows), OUT))
