# -*- coding: utf-8 -*-
"""
59_godmc_annotation.py - GoDMC mQTL 深度注释 (P1-6)

组件:
 A. 166 DMR (1030 CpG) 的 cis-mQTL 覆盖度 vs HM450 全背景 (Fisher 富集)
    - GoDMC cis 文件只含显著对 (p<~1e-5), "有 mQTL" = 在文件中出现
 B. 7 个 MR CpG 的 cis-mQTL 结构表 (对数/lead/p/效应)
 C. 160 个 IV SNP 的 Ensembl VEP 功能注释 (GRCh38)
 D. cfDNA 观测方向 (DMR hyper/hypo) x MR 因果方向 (OR) 三角验证表

输出:
 results/GSE282512_godmc_dmr_mqtl.csv          DMR 级汇总
 results/GSE282512_godmc_mr_cpg_annotation.csv 7 个 MR CpG 结构表 (含 VEP)
 results/GSE282512_godmc_iv_vep.csv            IV SNP VEP 注释
 results/GSE282512_godmc_direction_triangulation.csv 方向三角验证
 results/GSE282512_godmc_annotation_summary.txt 汇总
"""
import gzip
import csv
import json
import os
import subprocess
import sys
import time
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

GODMC = "data/mqtl/godmc_cis.csv.gz"
MANIFEST = "data/annot/HM450.hg38.manifest.tsv.gz"
DMR = "results/GSE282512_dmr_final_v2.csv"
WALD = "results/GSE282512_mr_wald_snps.csv"
IVW = "results/GSE282512_mr_ivw_results.csv"

MR_CPGS = ["cg18129748", "cg24308560", "cg27360567",   # 106528 chr3
           "cg00387872", "cg19490609", "cg25753631",   # 12268 SLC17A1 chr6
           "cg15445000"]                               # 84604 chr17


def log(msg):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), msg), flush=True)


# ---------------------------------------------------------------- 1. 载入
log("载入 HM450 manifest")
cpg_pos = {}   # cpg -> (chr, pos38)
with gzip.open(MANIFEST, "rt") as f:
    rd = csv.DictReader(f, delimiter="\t")
    for r in rd:
        try:
            cpg_pos[r["Probe_ID"]] = (r["CpG_chrm"].replace("chr", ""), int(r["CpG_beg"]))
        except Exception:
            pass
log("  manifest CpG: %d" % len(cpg_pos))

log("载入 166 DMR")
dmrs = []
with open(DMR) as f:
    for r in csv.DictReader(f):
        if r.get("final_call", "") in ("core_tube_robust", "tube_robust_fdr", "core", "") or True:
            dmrs.append(r)
# 只保留 166 个 final DMR (dmr_final_v2 已是 final 表)
log("  DMR: %d" % len(dmrs))

# DMR CpG 归属 (hg38)
dmr_by_chr = defaultdict(list)
for r in dmrs:
    dmr_by_chr[r["chr"].replace("chr", "")].append((int(r["start"]), int(r["end"]), r))
dmr_cpg_map = {}  # cpg -> dmr row
for cpg, (ch, p) in cpg_pos.items():
    for s, e, r in dmr_by_chr.get(ch, []):
        if s <= p <= e:
            dmr_cpg_map[cpg] = r
            break
log("  DMR 内 CpG: %d / %d manifest CpG" % (len(dmr_cpg_map), len(cpg_pos)))

# ---------------------------------------------------------------- 2. 扫描 GoDMC
log("流式扫描 GoDMC cis 文件 (66.6M 行)")
cpg_pairs = defaultdict(int)     # cpg -> 显著 cis-mQTL 对数
t0 = time.time()
with gzip.open(GODMC, "rt") as f:
    for line in f:
        cpg = line.split(",", 1)[0]
        if cpg in cpg_pairs or cpg in cpg_pos:
            cpg_pairs[cpg] += 1
log("  完成: %d 个 CpG 具显著 cis-mQTL (%.1f min)" % (len(cpg_pairs), (time.time() - t0) / 60))

# ---------------------------------------------------------------- 3. DMR 级 + 覆盖度富集
dmr_rows = []
for r in dmrs:
    cpgs = [c for c, d in dmr_cpg_map.items() if d is r]
    n_mqtl = sum(1 for c in cpgs if c in cpg_pairs)
    n_pairs = sum(cpg_pairs.get(c, 0) for c in cpgs)
    dmr_rows.append({
        "region_id": r["region_id"], "chr": r["chr"], "start": r["start"], "end": r["end"],
        "symbol": r.get("symbol", ""), "direction": r["direction"],
        "delta_beta": r.get("delta_beta", ""), "n_cpg": len(cpgs),
        "n_cpg_with_mqtl": n_mqtl, "pct_cpg_mqtl": round(100.0 * n_mqtl / len(cpgs), 1) if cpgs else "",
        "n_mqtl_pairs": n_pairs})

with open("results/GSE282512_godmc_dmr_mqtl.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(dmr_rows[0].keys()))
    w.writeheader()
    w.writerows(dmr_rows)

# 富集: DMR CpG vs 非 DMR HM450 背景
bg_cpgs = [c for c in cpg_pos if c not in dmr_cpg_map]
a = sum(1 for c in dmr_cpg_map if c in cpg_pairs)            # DMR CpG 有 mQTL
b = len(dmr_cpg_map) - a
c_ = sum(1 for c in bg_cpgs if c in cpg_pairs)               # 背景 CpG 有 mQTL
d = len(bg_cpgs) - c_
# Fisher exact 双侧 (超几何 pmf, lgamma 加速)
from math import lgamma, exp
def fisher(a, b, c, d):
    r1, r2 = a + b, c + d
    col1 = a + c
    n = r1 + r2
    if col1 == 0 or col1 == n or r1 == 0 or r2 == 0:
        return 1.0
    def logpmf(x):
        return (lgamma(r1 + 1) - lgamma(x + 1) - lgamma(r1 - x + 1)
                + lgamma(r2 + 1) - lgamma(col1 - x + 1) - lgamma(r2 - col1 + x + 1)
                - (lgamma(n + 1) - lgamma(col1 + 1) - lgamma(n - col1 + 1)))
    lo = max(0, col1 - r2)
    hi = min(col1, r1)
    lp_obs = logpmf(a)
    tot = 0.0
    for x in range(lo, hi + 1):
        lp = logpmf(x)
        if lp <= lp_obs + 1e-9:
            tot += exp(lp)
    return min(1.0, tot)
OR = (a * d) / (b * c_) if b and c_ else float("inf")
p_fisher = fisher(a, b, c_, d)
if p_fisher is None:
    p_fisher = float("nan")
log("富集: DMR %d/%d vs 背景 %d/%d, OR=%.3f, p=%.3g" % (a, len(dmr_cpg_map), c_, len(bg_cpgs), OR, p_fisher))

# hyper vs hypo
hyper = [r for r in dmr_rows if r["direction"] == "hyper"]
hypo = [r for r in dmr_rows if r["direction"] == "hypo"]
mh = sum(r["n_cpg_with_mqtl"] for r in hyper); nh = sum(r["n_cpg"] for r in hyper)
mo = sum(r["n_cpg_with_mqtl"] for r in hypo); no_ = sum(r["n_cpg"] for r in hypo)
OR_hh = (mh * (no_ - mo)) / ((nh - mh) * mo) if (nh - mh) and mo else float("inf")
p_hh = fisher(mh, nh - mh, mo, no_ - mo)

# ---------------------------------------------------------------- 4. MR CpG 结构表
log("MR CpG mQTL 结构 (从 mqtl_full7)")
mr_cpg_rows = []
all_pairs = defaultdict(list)
with open("results/_tmp_smr/mqtl_full7.csv") as f:
    for r in csv.reader(f):
        if len(r) >= 8 and r[0].startswith("cg"):
            loc = r[1].replace('"', "").replace(":SNP", "").split(":")
            try:
                all_pairs[r[0]].append({"chr": loc[0], "pos37": loc[1], "beta": float(r[2]),
                                        "se": float(r[3]), "p": float(r[4]),
                                        "a1": r[5], "a2": r[6], "maf": float(r[7]) if r[7] else None})
            except (ValueError, IndexError):
                pass
for cpg in MR_CPGS:
    pairs = all_pairs.get(cpg, [])
    lead = min(pairs, key=lambda x: x["p"]) if pairs else None
    dmr = dmr_cpg_map.get(cpg)
    mr_cpg_rows.append({
        "cpg": cpg, "region_id": dmr["region_id"] if dmr else "", "symbol": dmr.get("symbol", "") if dmr else "",
        "dmr_direction": dmr["direction"] if dmr else "", "dmr_delta_beta": dmr.get("delta_beta", "") if dmr else "",
        "cpg_chr": cpg_pos[cpg][0] if cpg in cpg_pos else "",
        "cpg_pos38": cpg_pos[cpg][1] if cpg in cpg_pos else "",
        "n_mqtl_pairs": len(pairs),
        "lead_snp_pos37": lead["pos37"] if lead else "",
        "lead_p": lead["p"] if lead else "",
        "lead_beta": lead["beta"] if lead else "",
        "lead_alleles": "%s/%s" % (lead["a1"], lead["a2"]) if lead else "",
        "lead_maf": lead["maf"] if lead else ""})

# ---------------------------------------------------------------- 5. VEP 注释 (160 IV)
log("Ensembl VEP 注释 IV SNP")
iv = {}
with open(WALD) as f:
    for r in csv.DictReader(f):
        iv[r["rsid"]] = r
log("  IV SNP: %d" % len(iv))

vep_res = {}
ids = list(iv.keys())
for i in range(0, len(ids), 100):
    batch = ids[i:i + 100]
    payload = json.dumps({"ids": batch})
    for attempt in range(4):
        try:
            out = subprocess.run(
                ["curl", "-s", "--ssl-no-revoke", "-m", "120", "-X", "POST",
                 "-H", "Content-Type: application/json",
                 "-H", "Accept: application/json",
                 "-d", payload, "https://rest.ensembl.org/vep/human/id"],
                capture_output=True, text=True, timeout=150).stdout
            data = json.loads(out)
            if isinstance(data, list):
                for item in data:
                    vep_res[item.get("input", "")] = item
                break
        except Exception as e:
            log("  VEP 重试 %d: %s" % (attempt + 1, e))
            time.sleep(5)

iv_rows = []
for rs, r in iv.items():
    v = vep_res.get(rs, {})
    msc = v.get("most_severe_consequence", "")
    gene = ""
    if v.get("transcript_consequences"):
        for tc in v["transcript_consequences"]:
            if tc.get("gene_symbol"):
                gene = tc["gene_symbol"]; break
    dist = ""
    if v.get("intergenic_consequences"):
        for ic in v["intergenic_consequences"]:
            if "distance" in ic:
                dist = ic["distance"]; break
    iv_rows.append({"rsid": rs, "cpg": r["cpg"], "region_id": r["region_id"],
                    "pos38": r["pos38"], "consequence": msc, "gene": gene,
                    "distance_to_gene": dist, "beta_m": r["beta_m"], "p_m": r["p_m"]})
with open("results/GSE282512_godmc_iv_vep.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(iv_rows[0].keys()))
    w.writeheader()
    w.writerows(iv_rows)
log("  VEP 注释: %d/%d" % (sum(1 for x in iv_rows if x["consequence"]), len(iv_rows)))

# 后果归类
def cat(c):
    if not c: return "unknown"
    if c in ("missense_variant", "synonymous_variant", "stop_gained", "stop_lost",
             "start_lost", "protein_altering_variant", "inframe_insertion",
             "inframe_deletion", "coding_sequence_variant"):
        return "coding"
    if "splice" in c: return "splice"
    if "UTR" in c: return "UTR"
    if c in ("intron_variant", "intron_variant,non_coding_transcript_variant",
             "non_coding_transcript_exon_variant", "non_coding_transcript_variant"):
        return "intronic/ncRNA"
    if "regulatory" in c or "TF_binding" in c or c == "enhancer_variant": return "regulatory"
    if "upstream" in c or "downstream" in c: return "flanking"
    if "intergenic" in c: return "intergenic"
    return "other"
cat_count = defaultdict(int)
for x in iv_rows:
    cat_count[cat(x["consequence"])] += 1

# ---------------------------------------------------------------- 6. 方向三角验证
log("方向三角验证")
tri_rows = []
with open(IVW) as f:
    ivw = {(r["cpg"], r["pheno"]): r for r in csv.DictReader(f)}
for row in mr_cpg_rows:
    cpg = row["cpg"]
    obs_dir = row["dmr_direction"]  # hyper = PE 中 cfDNA 甲基化更高
    for ph in ("PE_FinnGen", "PE_Tyrmi", "GH_FinnGen"):
        iv = ivw.get((cpg, ph))
        if not iv:
            continue
        ORv = float(iv["OR"]); pv = float(iv["p"])
        # MR: OR>1 = 遗传预测甲基化升高 -> 风险升高
        mr_dir = "risk_up" if ORv > 1 else "risk_down"
        if not obs_dir:
            consist = ""
        else:
            # 观测: hyper(PE 高甲基化) 且 MR OR>1(高甲基化升风险) -> 方向一致(因果解释成立)
            # 观测 hyper 但 MR OR<1 -> 反向 (支持反向因果/混杂, 甲基化为 PE 结果或标志)
            if obs_dir == "hyper":
                consist = "consistent" if ORv > 1 else "opposite"
            else:
                consist = "consistent" if ORv < 1 else "opposite"
        tri_rows.append({"cpg": cpg, "region_id": row["region_id"], "symbol": row["symbol"],
                         "pheno": ph, "obs_cfDNA_direction": obs_dir,
                         "dmr_delta_beta": row["dmr_delta_beta"],
                         "MR_OR": ORv, "MR_p": pv, "MR_direction": mr_dir,
                         "triangulation": consist,
                         "significant": "yes" if pv < 0.05 else "no"})
with open("results/GSE282512_godmc_direction_triangulation.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(tri_rows[0].keys()))
    w.writeheader()
    w.writerows(tri_rows)

# MR CpG 表写盘 (附 VEP 基因计数)
for row in mr_cpg_rows:
    cpg_iv = [x for x in iv_rows if x["cpg"] == row["cpg"]]
    row["n_iv"] = len(cpg_iv)
    row["iv_coding_or_splice"] = sum(1 for x in cpg_iv if cat(x["consequence"]) in ("coding", "splice"))
    row["iv_regulatory"] = sum(1 for x in cpg_iv if cat(x["consequence"]) == "regulatory")
with open("results/GSE282512_godmc_mr_cpg_annotation.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(mr_cpg_rows[0].keys()))
    w.writeheader()
    w.writerows(mr_cpg_rows)

# ---------------------------------------------------------------- 7. 汇总
with open("results/GSE282512_godmc_annotation_summary.txt", "w") as f:
    f.write("== GoDMC mQTL 深度注释 (P1-6) ==\n\n")
    f.write("[A] cis-mQTL 覆盖度 (GoDMC 显著对 p<~1e-5, N=32,000)\n")
    f.write("  DMR CpG (166 DMR):       %d/%d 具显著 cis-mQTL (%.1f%%)\n" % (a, len(dmr_cpg_map), 100.0 * a / len(dmr_cpg_map)))
    f.write("  HM450 背景非 DMR CpG:    %d/%d (%.1f%%)\n" % (c_, len(bg_cpgs), 100.0 * c_ / len(bg_cpgs)))
    f.write("  Fisher OR=%.3f, p=%.3g\n" % (OR, p_fisher))
    f.write("  hyper DMR: %d/%d (%.1f%%) vs hypo DMR: %d/%d (%.1f%%), OR=%.3f, p=%.3g\n\n"
            % (mh, nh, 100.0 * mh / nh if nh else 0, mo, no_, 100.0 * mo / no_ if no_ else 0, OR_hh, p_hh))
    f.write("[B] 7 个 MR CpG cis-mQTL 结构\n")
    for row in mr_cpg_rows:
        f.write("  %-12s %s %-10s %4d 对 | lead p=%.2g beta=%.3f | DMR %s Δβ=%s | IV %d (coding/splice %d, reg %d)\n"
                % (row["cpg"], row["region_id"], row["symbol"], row["n_mqtl_pairs"],
                   float(row["lead_p"]), float(row["lead_beta"]), row["dmr_direction"],
                   row["dmr_delta_beta"][:7], row["n_iv"], row["iv_coding_or_splice"], row["iv_regulatory"]))
    f.write("\n[C] IV SNP 功能后果 (Ensembl VEP, GRCh38, n=%d)\n" % len(iv_rows))
    for k in sorted(cat_count, key=lambda x: -cat_count[x]):
        f.write("  %-20s %d\n" % (k, cat_count[k]))
    f.write("\n[D] cfDNA 观测方向 x MR 因果方向 (显著 MR 结果)\n")
    for t in tri_rows:
        if t["significant"] == "yes":
            f.write("  %-12s %-9s obs=%-5s OR=%.3f p=%.2g -> %s\n"
                    % (t["cpg"], t["pheno"], t["obs_cfDNA_direction"], t["MR_OR"], t["MR_p"], t["triangulation"]))
    f.write("\n判读要点:\n")
    f.write("  - 若显著 MR 的方向与 cfDNA 观测方向相反 (opposite), 支持 DMR 甲基化为 PE 的结果/伴随标志\n")
    f.write("    (反向因果或白细胞构成混杂), 而非 PE 的因果驱动因素\n")
log("完成")
