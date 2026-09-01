# 27_dmr_sensitivity.R ---------------------------------------------------------
# Phase 1 敏感性分析: tube_type 混杂检验 + EOPE 子集
# 输入:
#   results/GSE282512_region_beta.csv.gz    区域x样本 beta 宽矩阵
#   results/GSE282512_subcohort.csv         子队列分组
#   results/GSE282512_region_annot.csv      区域注释
#   results/GSE282512_dmr_candidates.csv    原区域级统计表(含 166 候选)
# 分析:
#   A) EDTA-only  (PE 20 vs Control 25): limma ~ group + batch
#   B) PAXgene-only (PE 12 vs Control 7): 探索性 Wilcoxon + limma ~ group
#   C) EOPE (PE 21 vs GA 匹配对照): limma ~ group + tube + batch
#   D) 全队列交互检验: limma ~ group * tube + batch (group:tube 项)
# 与原 166 候选对比 -> tube_robust / tube_driven 分级
# 输出:
#   results/GSE282512_dmr_sensitivity.csv
#   results/GSE282512_dmr_sensitivity_summary.txt
#   figures/dmr_sensitivity_scatter.png
# 约定: UTF-8 保存, 相对路径, LANG=en_US.UTF-8 运行
# ------------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  library(data.table); library(limma)
}))

MIN_COVER <- 0.80
FDR_CUT   <- 0.05

## 1. 读取数据 ------------------------------------------------------------------
wide <- fread("results/GSE282512_region_beta.csv.gz")
ann  <- fread("results/GSE282512_region_annot.csv")
sub  <- fread("results/GSE282512_subcohort.csv")
cand <- fread("results/GSE282512_dmr_candidates.csv")   # 原分析全区域统计表

gsms <- names(wide)[-1]
mat  <- as.matrix(wide[, -"region_id"])
storage.mode(mat) <- "numeric"
rownames(mat) <- wide$region_id
cat(sprintf("矩阵: %d 区域 x %d 样本\n", nrow(mat), ncol(mat)))

## 2. 样本对齐 (按矩阵列序) ------------------------------------------------------
sub <- sub[match(gsms, gsm)]
stopifnot(nrow(sub) == ncol(mat), !anyNA(sub$gsm))
sub[, tube2 := ifelse(tube_type == "PAXgene DNA", "PAXgene", "EDTA")]
sub[, batch := as.factor(batch)]  # 与脚本 24 保持一致: 因子型批次

## 3. 通用函数 -------------------------------------------------------------------
# 子集矩阵: 限定区域(原分析测过的区域) + 子集内覆盖>=80%, NA 用中位数填充
build_subset <- function(idx) {
  m <- mat[, idx, drop = FALSE]
  keep <- rowMeans(!is.na(m)) >= MIN_COVER
  m <- m[keep, , drop = FALSE]
  na_fill <- function(x) { x[is.na(x)] <- median(x, na.rm = TRUE); x }
  t(apply(m, 1, na_fill))
}
# 安全 limma: 设计矩阵秩亏时自动去掉 batch
safe_fit <- function(m, design_full, coef_name) {
  ok <- qr(design_full)$rank == ncol(design_full)
  design <- if (ok) design_full else design_full[, !grepl("^batch", colnames(design_full)), drop = FALSE]
  fit <- eBayes(lmFit(m, design))
  tt  <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  data.table(region_id = as.integer(rownames(m)),
             p = tt$P.Value, fdr = tt$adj.P.Val)
}
spearman_ok <- function(x, y) {
  ct <- tryCatch(cor.test(x, y, method = "spearman"), error = function(e) NULL)
  if (is.null(ct)) return(c(NA, NA))
  c(ct$estimate, ct$p.value)
}

orig <- cand[, .(region_id, delta_beta_orig = delta_beta,
                 direction_orig = direction, candidate_orig = candidate)]

## 4. 分析 A: EDTA-only ----------------------------------------------------------
idxA <- which(sub$tube2 == "EDTA")
subA  <- sub[idxA]
cat("\n== A) EDTA 子集: PE", sum(subA$group == "PE"), "vs Control",
    sum(subA$group == "Control"), "==\n")
print(table(subA$group, subA$batch))
mA <- build_subset(idxA)
gA <- factor(subA$group, levels = c("Control", "PE"))
dA <- model.matrix(~ gA + batch, data = subA)
resA <- safe_fit(mA, dA, "gAPE")
resA[, delta_beta := rowMeans(mA[, gA == "PE", drop = FALSE]) -
                  rowMeans(mA[, gA == "Control", drop = FALSE])]
setnames(resA, c("p", "fdr"), c("p_edta", "fdr_edta"))
resA[, delta_beta_edta := delta_beta][, delta_beta := NULL]
cat(sprintf("EDTA 分析区域: %d; FDR<0.05: %d\n", nrow(resA), sum(resA$fdr_edta < FDR_CUT)))

## 5. 分析 B: PAXgene-only (探索性) ----------------------------------------------
idxB <- which(sub$tube2 == "PAXgene")
subB <- sub[idxB]
cat("\n== B) PAXgene 子集: PE", sum(subB$group == "PE"), "vs Control",
    sum(subB$group == "Control"), "(探索性) ==\n")
print(table(subB$group, subB$batch))
mB <- build_subset(idxB)
gB <- factor(subB$group, levels = c("Control", "PE"))
dB <- model.matrix(~ gB)
fitB <- eBayes(lmFit(mB, dB))
ttB  <- topTable(fitB, coef = "gBPE", number = Inf, sort.by = "none")
resB <- data.table(region_id = as.integer(rownames(mB)),
                   p_pax = ttB$P.Value, fdr_pax = ttB$adj.P.Val)
resB[, delta_beta_pax := rowMeans(mB[, gB == "PE", drop = FALSE]) -
                    rowMeans(mB[, gB == "Control", drop = FALSE])]
cat(sprintf("PAXgene 分析区域: %d; FDR<0.05: %d\n", nrow(resB), sum(resB$fdr_pax < FDR_CUT)))

## 6. 分析 C: EOPE 子集 ----------------------------------------------------------
pe_eope <- sub[group == "PE" & onset == "EOPE"]
cat("\n== C) EOPE 子集 ==\n")
cat("EOPE PE:", nrow(pe_eope), "; ga 范围:",
    min(pe_eope$ga_weeks), "-", max(pe_eope$ga_weeks), "\n")
ga_e <- pe_eope$ga_weeks
ctrl_all <- sub[group == "Control"]
dmin <- vapply(ctrl_all$ga_weeks, function(g) min(abs(g - ga_e)), numeric(1))
idxC <- which(sub$group == "PE" & sub$onset == "EOPE" |
              (sub$group == "Control" & dmin[match(sub$gsm, ctrl_all$gsm)] <= 3))
subC <- sub[idxC]
cat("GA 匹配对照:", sum(subC$group == "Control"), "\n")
print(table(subC$group, subC$tube2))
mC <- build_subset(idxC)
gC <- factor(subC$group, levels = c("Control", "PE"))
dC <- model.matrix(~ gC + tube2 + batch, data = subC)
resC <- safe_fit(mC, dC, "gCPE")
resC[, delta_beta_eope := rowMeans(mC[, gC == "PE", drop = FALSE]) -
                     rowMeans(mC[, gC == "Control", drop = FALSE])]
setnames(resC, c("p", "fdr"), c("p_eope", "fdr_eope"))
cat(sprintf("EOPE 分析区域: %d; FDR<0.05: %d\n", nrow(resC), sum(resC$fdr_eope < FDR_CUT)))

## 7. 分析 D: 全队列 group x tube 交互检验 ----------------------------------------
mD <- build_subset(seq_len(ncol(mat)))   # 全样本覆盖过滤 + NA 填充
dD <- model.matrix(~ group * tube2 + batch, data = sub[, .(group, tube2, batch)])
fitD <- eBayes(lmFit(mD, dD))
ttD  <- topTable(fitD, coef = "groupPE:tube2PAXgene", number = Inf, sort.by = "none")
resD <- data.table(region_id = as.integer(rownames(mD)),
                   p_inter = ttD$P.Value, fdr_inter = ttD$adj.P.Val)
cat(sprintf("交互项 FDR<0.05 区域: %d\n", sum(resD$fdr_inter < FDR_CUT)))

## 8. 合并与候选分级 --------------------------------------------------------------
sens <- Reduce(function(x, y) merge(x, y, by = "region_id", all = TRUE),
               list(orig, resA, resB, resC, resD))
sens <- merge(sens, ann[, .(region_id, chr, start, end, type, symbol)],
              by = "region_id", all.x = TRUE, sort = FALSE)

candS <- sens[candidate_orig == TRUE]
# tube_robust: EDTA 内方向一致且名义 p<0.05; FDR 复现更强
candS[, sens_flag := fcase(
  !is.na(fdr_edta) & fdr_edta < FDR_CUT & sign(delta_beta_edta) == sign(delta_beta_orig),
    "tube_robust_fdr",
  !is.na(p_edta) & p_edta < 0.05 & sign(delta_beta_edta) == sign(delta_beta_orig),
    "tube_robust_nominal",
  !is.na(p_edta) & sign(delta_beta_edta) == sign(delta_beta_orig),
    "direction_only",
  default = "tube_driven")]
candS[, eo_flag := fcase(
  !is.na(fdr_eope) & fdr_eope < FDR_CUT & sign(delta_beta_eope) == sign(delta_beta_orig),
    "eope_fdr",
  !is.na(p_eope) & p_eope < 0.05 & sign(delta_beta_eope) == sign(delta_beta_orig),
    "eope_nominal",
  default = "not_eope")]

setorder(sens, region_id)
fwrite(sens, "results/GSE282512_dmr_sensitivity.csv")
fwrite(candS, "results/GSE282512_dmr_sensitivity_candidates.csv")

## 9. 汇总报告 ---------------------------------------------------------------------
sink("results/GSE282512_dmr_sensitivity_summary.txt", split = TRUE)
cat("===== GSE282512 Phase 1 敏感性分析汇总 =====\n\n")
cat(sprintf("原分析: %d 候选 DMR (limma, tube+batch 调整)\n\n", sum(cand$candidate)))

cat("-- A) EDTA-only (PE 20 vs Control 25) --\n")
cat(sprintf("测过区域: %d; FDR<0.05: %d\n", nrow(resA), sum(resA$fdr_edta < FDR_CUT)))
sc <- spearman_ok(sens[!is.na(delta_beta_edta), delta_beta_orig],
                  sens[!is.na(delta_beta_edta), delta_beta_edta])
cat(sprintf("全区域 delta_beta Spearman rho = %.3f (p = %.2e)\n", sc[1], sc[2]))
cat(sprintf("166 候选中: EDTA FDR<0.05 且方向一致: %d; 名义 p<0.05 且方向一致: %d; 方向一致: %d; 方向翻转: %d\n",
    sum(candS$sens_flag == "tube_robust_fdr"),
    sum(candS$sens_flag == "tube_robust_nominal"),
    sum(candS$sens_flag %in% c("direction_only")),
    sum(candS$sens_flag == "tube_driven")))
cat(sprintf("候选方向一致率 (EDTA): %.1f%%\n",
    100 * mean(sign(candS$delta_beta_orig) == sign(candS$delta_beta_edta), na.rm = TRUE)))
cat(sprintf("候选 delta_beta 中位收缩比 (EDTA/原): %.2f\n",
    median(candS$delta_beta_edta / candS$delta_beta_orig, na.rm = TRUE)))

cat("\n-- B) PAXgene-only (PE 12 vs Control 7, 探索性) --\n")
cat(sprintf("测过区域: %d; FDR<0.05: %d\n", nrow(resB), sum(resB$fdr_pax < FDR_CUT)))
cat(sprintf("166 候选中 PAXgene 方向一致率: %.1f%%\n",
    100 * mean(sign(candS$delta_beta_orig) == sign(candS$delta_beta_pax), na.rm = TRUE)))

cat("\n-- C) EOPE (PE 21 vs GA 匹配对照) --\n")
cat(sprintf("对照数: %d; 测过区域: %d; FDR<0.05: %d\n",
    sum(subC$group == "Control"), nrow(resC), sum(resC$fdr_eope < FDR_CUT)))
cat(sprintf("166 候选中: EOPE FDR<0.05 且方向一致: %d; 名义且方向一致: %d\n",
    sum(candS$eo_flag == "eope_fdr"), sum(candS$eo_flag == "eope_nominal")))
cat(sprintf("候选方向一致率 (EOPE): %.1f%%\n",
    100 * mean(sign(candS$delta_beta_orig) == sign(candS$delta_beta_eope), na.rm = TRUE)))

cat("\n-- D) 全队列 group x tube 交互检验 --\n")
cat(sprintf("交互项 FDR<0.05 区域: %d (占测过 %.2f%%)\n",
    sum(resD$fdr_inter < FDR_CUT), 100 * mean(resD$fdr_inter < FDR_CUT)))
cat(sprintf("166 候选中交互项 FDR<0.05: %d\n",
    sum(candS$fdr_inter < FDR_CUT, na.rm = TRUE)))

cat("\n-- 分级汇总 (166 候选) --\n")
print(candS[, .N, by = sens_flag][order(-N)])
cat("\n")
print(candS[, .N, by = eo_flag][order(-N)])
cat("\n-- Top 15 tube-robust 候选 --\n")
print(candS[sens_flag %in% c("tube_robust_fdr", "tube_robust_nominal")][
  order(p_edta)][1:15,
  .(region_id, chr, start, end, symbol, direction_orig, delta_beta_orig,
    delta_beta_edta, p_edta, fdr_edta, delta_beta_eope, p_eope, sens_flag)])
sink()

## 10. 图: 原 vs EDTA delta_beta 散点 ---------------------------------------------
png("figures/dmr_sensitivity_scatter.png", width = 1400, height = 1300, res = 150)
d <- sens[!is.na(delta_beta_edta)]
colv <- ifelse(d$candidate_orig,
               ifelse(sign(d$delta_beta_orig) == sign(d$delta_beta_edta), "firebrick", "purple2"),
               "grey80")
plot(d$delta_beta_orig, d$delta_beta_edta, pch = 20, cex = ifelse(d$candidate_orig, 0.8, 0.4),
     col = colv, xlab = expression(Delta*beta~"full cohort (n = 64)"),
     ylab = expression(Delta*beta~"EDTA-only subset"),
     main = "Tube-type sensitivity: full cohort vs EDTA-only")
abline(0, 1, lty = 2, col = "grey40"); abline(h = 0, v = 0, lty = 3, col = "grey70")
legend("topleft", legend = c("candidate DMR - same direction", "candidate DMR - opposite direction", "non-candidate"),
       col = c("firebrick", "purple2", "grey80"), pch = 20, bty = "n", cex = 0.8)
dev.off()

cat("\nDONE\n")
