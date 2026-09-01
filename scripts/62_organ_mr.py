# -*- coding: utf-8 -*-
"""
62_organ_mr.py - 器官孟德尔化: cis-meQTL IV -> 血细胞/器官表型 MR (P1-7)

输入: results/_tmp_organ/organ_outcome_effects.csv (61 号产出)
方法: 等位基因定向 (回文模糊剔除) + IVW (Q 显著切 DL 随机效应) / 单 SNP Wald
输出: results/GSE282512_organ_mr_results.csv
      results/GSE282512_organ_mr_summary.txt
      figures/GSE282512_organ_mr_forest.png
"""
import csv
import math
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

IN = "results/_tmp_organ/organ_outcome_effects.csv"
OUT = "results/GSE282512_organ_mr_results.csv"
SUMMARY = "results/GSE282512_organ_mr_summary.txt"
FIG = "figures/GSE282512_organ_mr_forest.png"

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def harmonize(a1, a2, f1, ea, oa, af):
    """返回 (beta_m 定向, 是否有效) — 暴露等位基因 a1=效应"""
    if not ea or not oa or not a1 or not a2:
        return None, False
    ea, oa, a1, a2 = ea.upper(), oa.upper(), a1.upper(), a2.upper()
    if len(ea) != 1 or len(oa) != 1 or len(a1) != 1 or len(a2) != 1:
        return None, False
    # 正常匹配
    if ea == a1 and oa == a2:
        return 1.0, True
    if ea == a2 and oa == a1:
        return -1.0, True
    # 互补链匹配
    ca1, ca2 = COMP.get(a1), COMP.get(a2)
    if ea == ca1 and oa == ca2:
        return 1.0, True
    if ea == ca2 and oa == ca1:
        return -1.0, True
    # 回文: 用频率消歧
    pal_exp = {a1, a2} == {ea, oa} or {a1, a2} == {COMP.get(ea), COMP.get(oa)}
    if pal_exp and f1 is not None and af not in ("", None):
        try:
            af = float(af)
            if abs(f1 - af) < 0.04 and abs(f1 - (1 - af)) > 0.04:
                return 1.0, True   # 方向一致
            if abs(f1 - (1 - af)) < 0.04 and abs(f1 - af) > 0.04:
                return -1.0, True  # 方向相反
        except ValueError:
            pass
    return None, False


# ---------------------------------------------------------------- 读取
rows = []
with open(IN) as f:
    for r in csv.DictReader(f):
        rows.append(r)
print("载入 IV-结局配对: %d" % len(rows))

# ---------------------------------------------------------------- 追加本地 FinnGen PE/GH (管线验证)
# finngen 窗口文件: chrom pos ref alt rsids ... pval beta sebeta af_alt (GRCh38, alt=effect)
def load_finngen(path, pheno):
    by_pos = {}
    by_rsid = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            try:
                pos = int(r["pos"])
            except (ValueError, KeyError):
                continue
            rec = {"pos": pos, "rsid": (r.get("rsids") or "").split(";")[0],
                   "ea": r.get("alt", ""), "oa": r.get("ref", ""),
                   "beta": r.get("beta", ""), "se": r.get("sebeta", ""),
                   "p": r.get("pval", ""), "af": r.get("af_alt", "")}
            by_pos[(str(r.get("chrom")), pos)] = rec
            if rec["rsid"]:
                by_rsid[rec["rsid"]] = rec
    return by_pos, by_rsid

IV_BY_RSID = {}
for r in rows:
    if r["rsid"]:
        IV_BY_RSID[r["rsid"]] = r
# IV 元数据从 organ_outcome_effects 任取一行 (同 rsid 同暴露)
iv_meta = {}
for r in rows:
    key = (r["cpg"], r["rsid"])
    if key not in iv_meta:
        iv_meta[key] = r

n_add = 0
for pheno, path in [("PE_FinnGen", "results/_tmp_smr/finngen_pe_window.tsv"),
                    ("GH_FinnGen", "results/_tmp_smr/finngen_gh_window.tsv")]:
    by_pos, by_rsid = load_finngen(path, pheno)
    for (cpg, rsid), iv in iv_meta.items():
        rec = by_rsid.get(rsid)
        if rec is None:
            continue
        try:
            bo = float(rec["beta"]); so = float(rec["se"])
        except (ValueError, TypeError):
            continue
        if so is None or so <= 0:
            continue
        rows.append({"family": pheno, "accession": "FinnGen_R10_local", "trait": pheno,
                     "build": "GRCh38", "cpg": cpg, "region_id": iv["region_id"],
                     "rsid": rsid, "pos": rec["pos"], "ea": rec["ea"], "oa": rec["oa"],
                     "af": rec["af"], "beta_o": bo, "se_o": so, "p_o": rec["p"],
                     "a1": iv["a1"], "a2": iv["a2"], "freq_a1": iv["freq_a1"],
                     "beta_m": iv["beta_m"], "se_m": iv["se_m"], "p_m": iv["p_m"]})
        n_add += 1
print("追加 FinnGen PE/GH 配对: %d" % n_add)

# ---------------------------------------------------------------- IVW
tests = defaultdict(list)
for r in rows:
    sign, ok = harmonize(r["a1"], r["a2"], float(r["freq_a1"]) if r["freq_a1"] else None,
                         r["ea"], r["oa"], r["af"])
    if not ok:
        continue
    try:
        bo = float(r["beta_o"]); so = float(r["se_o"])
        bm = sign * float(r["beta_m"]); sm = float(r["se_m"])
    except (ValueError, TypeError):
        continue
    if so is None or so <= 0 or abs(bm) < 1e-9:
        continue
    wald = bo / bm
    se_w = math.sqrt(so ** 2 / bm ** 2 + bo ** 2 * sm ** 2 / bm ** 4)
    F = (bm / sm) ** 2
    tests[(r["family"], r["accession"], r["trait"], r["cpg"], r["region_id"])].append(
        {"rsid": r["rsid"], "wald": wald, "se": se_w, "F": F})

print("可分析 (family, cpg) 组合: %d" % len(tests))


def chi2_sf(x, df):
    """卡方上尾概率: df=1 精确, 否则 Wilson-Hilferty 近似"""
    from math import erf, sqrt
    if x <= 0:
        return 1.0
    if df == 1:
        return math.erfc(math.sqrt(x / 2.0))
    z = ((x / df) ** (1 / 3.0) - (1 - 2.0 / (9 * df))) / math.sqrt(2.0 / (9 * df))
    return 0.5 * (1 - erf(z / sqrt(2.0)))


def ivw(walds, ses):
    """IVW-FE; Q 显著 (p<0.05) 切 DL 随机效应"""
    k = len(walds)
    if k == 1:
        return walds[0], ses[0], None, None, "wald"
    w = [1.0 / (s * s) for s in ses]
    sw = sum(w)
    b = sum(wi * x for wi, x in zip(w, walds)) / sw
    se_fe = math.sqrt(1.0 / sw)
    Q = sum(wi * (x - b) ** 2 for wi, x in zip(w, walds))
    Qp = chi2_sf(Q, k - 1)
    if Qp < 0.05:
        # DL 随机效应
        tau2 = max(0.0, (Q - (k - 1)) / (sw - sum(wi * wi for wi in w) / sw))
        w_re = [1.0 / (s * s + tau2) for s in ses]
        swre = sum(w_re)
        b = sum(wi * x for wi, x in zip(w_re, walds)) / swre
        se = math.sqrt(1.0 / swre)
        return b, se, Q, Qp, "ivw_re"
    return b, se_fe, Q, Qp, "ivw_fe"


def norm_sf(z):
    return 0.5 * math.erfc(abs(z) / math.sqrt(2.0))


res = []
for (fam, acc, trait, cpg, region), ivs in sorted(tests.items()):
    walds = [x["wald"] for x in ivs]
    ses = [x["se"] for x in ivs]
    b, se, Q, Qp, method = ivw(walds, ses)
    z = b / se if se > 0 else 0.0
    p = norm_sf(z) * 2
    meanF = sum(x["F"] for x in ivs) / len(ivs)
    res.append({"family": fam, "accession": acc, "trait": trait, "cpg": cpg,
                "region_id": region, "k_iv": len(ivs), "method": method,
                "beta": b, "se": se, "lo": b - 1.96 * se, "hi": b + 1.96 * se,
                "p": p, "Qp": Qp if Qp is not None else "", "mean_F": meanF})

# 验证结局 (本地 FinnGen) 不参与器官 MR 的 BH
res_val = [r for r in res if r["accession"] == "FinnGen_R10_local"]
res = [r for r in res if r["accession"] != "FinnGen_R10_local"]

# BH
m = len(res)
order = sorted(range(m), key=lambda i: res[i]["p"])
prev = 1.0
for rank in range(m - 1, -1, -1):
    i = order[rank]
    q = min(prev, res[i]["p"] * m / (rank + 1))
    res[i]["p_bh"] = q
    prev = q

with open(OUT, "w", newline="") as f:
    allr = res + res_val
    w = csv.DictWriter(f, fieldnames=list(allr[0].keys()))
    w.writeheader()
    w.writerows(allr)

# ---------------------------------------------------------------- 验证: vs 58 号 MR (本地 FinnGen)
import math as _m
val_expect = {("cg19490609", "PE_FinnGen"): 0.898, ("cg15445000", "PE_FinnGen"): 1.089,
              ("cg15445000", "GH_FinnGen"): 1.048}
val_lines = []
for r in res_val:
    OR = _m.exp(r["beta"])
    exp_v = val_expect.get((r["cpg"], r["family"]))
    tag = ""
    if exp_v:
        tag = " | 58号: OR=%.3f %s" % (exp_v, "一致" if abs(OR - exp_v) / exp_v < 0.10 else "偏差>10%")
    val_lines.append("  %-12s %-10s k=%-3d OR=%.3f [%.3f,%.3f] p=%.3g%s"
                     % (r["cpg"], r["family"], r["k_iv"], OR,
                        _m.exp(r["lo"]), _m.exp(r["hi"]), r["p"], tag))

# ---------------------------------------------------------------- 汇总
sig = [r for r in res if r["p_bh"] < 0.05]
nom = [r for r in res if r["p"] < 0.05]
with open(SUMMARY, "w") as f:
    f.write("== 器官孟德尔化: cis-meQTL MR -> 血细胞/器官表型 (P1-7) ==\n\n")
    f.write("暴露: GoDMC cis-meQTL IV (7 CpG, 与 58 号 MR 同池, LD clump r2<0.1)\n")
    f.write("结局: GWAS Catalog harmonised 全量汇总 (远程 tabix 按位点提取, EBI FTP)\n")
    f.write("  血细胞 12 性状: Sakaue 2021 Nat Genet 图谱 (EUR ~350k, Vuckovic 单核/嗜酸/嗜碱补充)\n")
    f.write("  eGFR: Sakaue 图谱 (EUR 342k); SBP/DBP/PP: Evangelou 2024 Nat Genet (EUR 1,028,980, mmHg)\n")
    f.write("  (肝酶/CRP/尿酸/肌酐: GWAS Catalog 无带 tabix 索引的欧洲大样本全量文件, 未纳入)\n")
    f.write("方法: 等位基因定向 (回文 AF 消歧, 模糊剔除) + IVW (Q 显著切 DL) / 单 SNP Wald; F 统计量\n")
    f.write("效应量: 每 1-SD 甲基化变化的结局变化 (血细胞为 SD 单位, 血压为 mmHg, PE/GH 为 logOR)\n\n")
    f.write("共 %d 项检验 (%d CpG x %d 性状域, 平均可分析 IV %.1f)\n\n"
            % (len(res), len(set(r['cpg'] for r in res)), len(set(r['family'] for r in res)),
               sum(r['k_iv'] for r in res) / max(1, len(res))))
    f.write("---- 管线验证: 本地 FinnGen PE/GH vs 58 号 MR ----\n")
    for ln in val_lines:
        f.write(ln + "\n")
    f.write("\n---- 名义显著 (p<0.05, n=%d, 不含验证结局) ----\n" % len(nom))
    for r in sorted(nom, key=lambda x: x["p"]):
        f.write("  %-12s %-14s %-12s k=%-3d beta=%+.4f [%.4f,%.4f] p=%.3g q=%.3g F=%.0f\n"
                % (r["cpg"], r["family"], r["accession"], r["k_iv"], r["beta"], r["lo"], r["hi"],
                   r["p"], r["p_bh"], r["mean_F"]))
    f.write("\n---- BH 显著 (q<0.05, n=%d) ----\n" % len(sig))
    for r in sorted(sig, key=lambda x: x["p"]):
        f.write("  %-12s %-14s %s k=%d beta=%+.4f p=%.3g q=%.3g\n"
                % (r["cpg"], r["family"], r["accession"], r["k_iv"], r["beta"], r["p"], r["p_bh"]))
    f.write("\n---- 全部结果 (按性状) ----\n")
    for fam in sorted(set(r["family"] for r in res)):
        sub = [r for r in res if r["family"] == fam]
        tr = sub[0]["trait"]; acc = sub[0]["accession"]
        f.write("  [%s] %s (%s)\n" % (fam, tr, acc))
        for r in sorted(sub, key=lambda x: x["p"]):
            star = "*" if r["p"] < 0.05 else " "
            f.write("     %-12s%s k=%-3d beta=%+.4f [%.4f,%.4f] p=%.3g q=%.3g\n"
                    % (r["cpg"], star, r["k_iv"], r["beta"], r["lo"], r["hi"], r["p"], r["p_bh"]))
print("结果: %d 项, 名义显著 %d, BH 显著 %d" % (len(res), len(nom), len(sig)))

# ---------------------------------------------------------------- 森林图 (显著 + 代表性)
show = sorted(nom, key=lambda x: x["p"])[:25]
if len(show) < 10:
    show = sorted(res, key=lambda x: x["p"])[:20]
n = len(show)
fig, ax = plt.subplots(figsize=(9, 0.32 * n + 1.2))
ys = list(range(n))[::-1]
for i, r in enumerate(show):
    y = ys[i]
    ax.plot([r["lo"], r["hi"]], [y, y], color="#c0392b" if r["p"] < 0.05 else "#7f8c8d", lw=1.6)
    ax.plot([r["beta"]], [y], "o", ms=4.5, color="#c0392b" if r["p"] < 0.05 else "#7f8c8d")
ax.axvline(0, color="k", lw=0.8, ls="--")
ax.set_yticks(ys)
ax.set_yticklabels(["%s | %s (k=%d)" % (r["cpg"], r["family"], r["k_iv"]) for r in show], fontsize=7.5)
ax.set_xlabel("Outcome change per 1-SD methylation (outcome units)")
ax.set_title("Organ MR: cis-mQTL instrumented DMR CpG methylation → blood cell & organ traits")
fig.tight_layout()
fig.savefig(FIG, dpi=200)
print("图: %s" % FIG)
