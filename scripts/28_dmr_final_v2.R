# 28_dmr_final_v2.R -------------------------------------------------------------
# 合并位点级验证 + 敏感性分析 -> 最终 v2 判定表
# 输入:
#   results/GSE282512_dmr_final.csv               166 候选 (含 validation 分级)
#   results/GSE282512_dmr_sensitivity_candidates.csv  敏感性分级 (sens_flag/eo_flag)
# 输出:
#   results/GSE282512_dmr_final_v2.csv
#   results/GSE282512_dmr_final_v2_summary.txt
# 约定: UTF-8 保存, 相对路径, LANG=en_US.UTF-8 运行
# ------------------------------------------------------------------------------
suppressWarnings(suppressMessages(library(data.table)))

fin  <- fread("results/GSE282512_dmr_final.csv")
sens <- fread("results/GSE282512_dmr_sensitivity_candidates.csv")

v2 <- merge(fin, sens[, .(region_id, delta_beta_edta, p_edta, fdr_edta,
                          delta_beta_pax, p_pax, delta_beta_eope, p_eope, fdr_eope,
                          p_inter, fdr_inter, sens_flag, eo_flag)],
            by = "region_id", all.x = TRUE, sort = FALSE)

# 综合证据分级:
#   high_confidence: site_confirmed + tube_robust (FDR 或 nominal)
#   core_tube_robust: tube_robust_fdr + (direction_consistent 或 site_confirmed)
#   tube_robust:      EDTA 内名义显著且方向一致
#   sensitivity_only: 敏感性中方向一致但位点级未复现
v2[, final_call := fcase(
  validation == "site_confirmed" & sens_flag %in% c("tube_robust_fdr", "tube_robust_nominal"),
    "high_confidence",
  sens_flag == "tube_robust_fdr" & validation %in% c("direction_consistent", "site_confirmed"),
    "core_tube_robust",
  sens_flag == "tube_robust_nominal",
    "tube_robust",
  sens_flag == "direction_only",
    "sensitivity_only",
  default = "deprioritized")]

setorder(v2, final_call, fdr_limma)
fwrite(v2, "results/GSE282512_dmr_final_v2.csv")

sink("results/GSE282512_dmr_final_v2_summary.txt", split = TRUE)
cat("===== GSE282512 DMR 最终判定 v2 (位点级验证 + 敏感性合并) =====\n\n")
cat("敏感性分析关键结论:\n")
cat("  1. tube_type 混杂被排除: 166 候选 EDTA-only 方向一致率 100%, 0 个方向翻转;\n")
cat("     EDTA 子集内 63 个 FDR<0.05 复现, 且 EDTA-only 全基因组 FDR<0.05 达 2212 个;\n")
cat("     group x tube 交互检验 0 个区域显著 -> 信号非采血管类型驱动。\n")
cat("  2. PAXgene-only 方向一致率 99.4% (小样本探索性)。\n")
cat("  3. EOPE 子集: 166 候选方向一致率 100%, 27 个 FDR + 135 个名义复现。\n\n")
cat("-- 综合证据分级 --\n")
print(v2[, .N, by = final_call][order(-N)])
cat("\n-- high_confidence / core_tube_robust 区域 --\n")
print(v2[final_call %in% c("high_confidence", "core_tube_robust"),
  .(region_id, chr, start, end, symbol, type, direction, delta_beta,
    delta_beta_edta, fdr_edta, delta_beta_eope, p_eope,
    validation, sens_flag, eo_flag, final_call)])
cat("\n-- 结论 --\n")
cat("区域级信号稳健(不依赖 tube_type、在 EOPE 富集), 但位点级(每 CpG 均值)复现弱,\n")
cat("提示效应集中在覆盖度构成层面; 后续解释与功能注释应围绕 tube-robust 区域,\n")
cat("优先 high_confidence (位点级+敏感性双通过) 核心 DMR。\n")
sink()

cat("DONE\n")
