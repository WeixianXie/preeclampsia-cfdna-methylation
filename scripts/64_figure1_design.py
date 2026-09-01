# -*- coding: utf-8 -*-
"""64_figure1_design.py - Figure 1 study design / evidence layers flowchart"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(11.5, 14.5))
ax.set_xlim(0, 100)
ax.set_ylim(0, 128)
ax.axis("off")

C_DISC = "#DCE9F7"   # discovery
C_DMR = "#C4D7EE"
C_QUEST = "#FCE8C8"
C_NEG = "#F6D5D5"    # placenta-negative layers
C_POS = "#D8EAD3"    # leukocyte-supportive layers
C_GEN = "#E7DCEF"
C_CONC = "#B7CDE8"
EDGE = "#4A4A4A"

def box(x, y, w, h, text, fc, fs=8.4, bold=False, ec=EDGE, title=None):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.6",
                       linewidth=1.0, edgecolor=ec, facecolor=fc)
    ax.add_patch(p)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fs, fontweight=("bold" if bold else "normal"), linespacing=1.35)

def arrow(x1, y1, x2, y2, style="-|>", lw=1.4, color=EDGE):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                 mutation_scale=16, linewidth=lw, color=color))

# ---------------- 顶部: 发现队列 ----------------
box(22, 119, 56, 8,
    "Discovery: GSE282512 plasma cfDNA WGBS\n"
    "369 profiles / 279 women -> gestational-age-matched 32 PE vs 32 controls",
    C_DISC, fs=8.6, bold=False)
arrow(50, 119, 50, 114.5)

box(24, 106.5, 52, 7.5,
    "166 candidate DMRs (105 hypermethylated / 61 hypomethylated)\n"
    "robust to collection tube (100% direction consistency) and enriched in early-onset PE",
    C_DMR, fs=8.2)

# 侧支: 分类器
box(3, 106.5, 18, 7.5,
    "Elastic net classifier\nLOOCV AUC 0.932\n(internally validated only)",
    "#FFFFFF", fs=7.8)
arrow(24, 110.2, 21.5, 110.2, style="-")

arrow(50, 106.5, 50, 102.5)

# ---------------- 核心问题 ----------------
box(20, 95, 60, 7,
    "Central question: what is the cellular origin of the\n"
    "PE-associated plasma cfDNA methylation signal?",
    C_QUEST, fs=9.0, bold=True)

# 两个分支箭头
arrow(35, 95, 20, 90)
arrow(65, 95, 80, 90)
ax.text(16.5, 91.8, "Placental\norigin?", fontsize=8.2, ha="center", style="italic")
ax.text(83.5, 91.8, "Maternal\nleukocytes?", fontsize=8.2, ha="center", style="italic")

# ---------------- 五个证据层 ----------------
# 层 1: 表达
box(4, 80, 44, 9.5,
    "Layer 1 - Expression meta-analysis (8 cohorts)\n"
    "6 placental + 2 maternal blood RNA cohorts\n"
    "0/66 DMR genes differentially expressed in placenta (meta p<0.05);\n"
    "direction consistency at chance level",
    C_NEG, fs=7.9)
ax.text(49.5, 84.8, "evidence against\nplacental origin", fontsize=7.2, ha="center",
        color="#943634", style="italic")

box(52, 80, 44, 9.5,
    "Layer 2 - Placental tissue methylome (3 cohorts, 105 samples)\n"
    "1,030 HM450 CpGs in 127 DMRs: mean |z| = 0.35, 57.8th background percentile;\n"
    "permutation p = 1.0; direction agreement 46% (chance); rho = -0.04\n"
    "positive control (top-500 placental CpGs): mean |z| = 3.99",
    C_NEG, fs=7.9)
ax.text(50, 78.2, "the same pipeline detects real placental signal - the DMR set has none",
        fontsize=7.4, ha="center", color="#943634", style="italic")

# 层 3: 调控注释
box(4, 65.5, 44, 10.5,
    "Layer 3 - Regulatory grammar\n"
    "Ensembl regulatory build: enhancer OR 1.47 (FDR 0.037),\n"
    "promoter depletion OR 0.48\n"
    "JASPAR: FOS / AP-1 motif enrichment 15-fold (FDR 0.022)\n"
    "a myeloid inflammatory signature, not trophoblast-response genes",
    C_POS, fs=7.9)
ax.text(49.5, 70.8, "consistent with\nleukocyte origin", fontsize=7.2, ha="center",
        color="#2E6B2E", style="italic")

box(52, 65.5, 44, 10.5,
    "Layer 4 - Longitudinal trajectories and deconvolution\n"
    "GSE154378 cfDNA (134 samples) + GSE37722 leukocytes (84 samples)\n"
    "placental markers rise with gestation (rho +0.54) while blood-cell markers fall;\n"
    "T-cell fraction declines in both cohorts, postpartum reversal\n"
    "PreX subgroup: neutrophil fraction +0.070 already in first trimester\n"
    "at delivery +0.063 despite lower placental fraction -0.126; GDM positive control",
    C_POS, fs=7.9)
ax.text(50, 63.6, "composition shift precedes clinical disease and reverses postpartum",
        fontsize=7.4, ha="center", color="#2E6B2E", style="italic")

# 层 5: 遗传
box(4, 49.5, 92, 12,
    "Layer 5 - Genetic architecture\n"
    "cis-mQTL coverage: DMR CpGs 28.8% vs background 37.7% (OR 0.669, p=2.1e-9; hyper-DMRs 22.2%) - dynamically, not genetically, regulated\n"
    "Two-sample MR (FinnGen R10, Tyrmi 2023): cg19490609 in SLC17A1 OR 0.898 per SD methylation (p=1.0e-4, BH 7.3e-4);\n"
    "cg15445000 OR 1.089 (p=2.2e-4, BH 7.6e-4), replicated in Tyrmi 2023 (OR 1.059)\n"
    "Organ-wide MR (16 traits): 14 BH-significant effects - cg19490609 raises monocyte / lymphocyte / WBC / eosinophil counts,\n"
    "cg24308560 / cg27360567 lower WBC and neutrophils; cg25753631 raises SBP (+0.27 mmHg); pipeline check reproduced FinnGen estimates",
    C_GEN, fs=7.9)
ax.text(50, 47.6, "weak cis-regulation of the DMR set overall; where instruments exist, effects run through blood cell counts",
        fontsize=7.4, ha="center", color="#5B3A75", style="italic")

# 层间箭头
arrow(26, 80, 26, 76.5)
arrow(74, 80, 74, 76.5)
arrow(26, 65.5, 26, 62)
arrow(74, 65.5, 74, 62)
arrow(26, 49.5, 26, 45)
arrow(74, 49.5, 74, 45)

# ---------------- 结论 ----------------
box(8, 33, 84, 11,
    "Conclusion\n"
    "PE-associated plasma cfDNA methylation differences behave as a dynamic fingerprint of\n"
    "maternal leukocyte composition: absent in placental tissue (at 8.7% of positive-control strength),\n"
    "tracking blood-cell trajectories, shifted toward neutrophils before clinical onset, weakly genetically\n"
    "regulated as a set, and causally connected to blood cell counts where instruments exist.\n"
    "cfDNA methylation biomarkers for PE should be validated as readouts of maternal immune state,\n"
    "not as a placental biopsy in a tube.",
    C_CONC, fs=8.6, bold=True)
arrow(50, 33, 50, 28.5)

# ---------------- 底部注释: 外部验证状态 ----------------
box(12, 19, 76, 8.5,
    "Validation status\n"
    "Classifier: internally validated candidate signature (LOOCV AUC 0.932); definitive external\n"
    "validation in a first-trimester cfDNA methylome cohort (EGAS00001007071) under access application.\n"
    "Origin evidence: replicated across cohort, platform and modality (plasma cfDNA, leukocyte layer, tissue, genetics).",
    "#FFFFFF", fs=7.8)

ax.text(50, 14.5,
        "Data: GSE282512 (discovery) | GSE57767, GSE73375, GSE75196 (placenta) | GSE154378 (cfDNA trajectories) | GSE37722 (leukocytes)\n"
        "GoDMC cis-mQTLs | FinnGen R10 | Tyrmi 2023 meta-GWAS | Sakaue 2021 blood cell atlas | Evangelou 2024 blood pressure GWAS",
        fontsize=7.2, ha="center", va="center", color="#555555", linespacing=1.5)

plt.tight_layout()
out = "figures/Figure1_study_design.png"
fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
print("saved", out)
