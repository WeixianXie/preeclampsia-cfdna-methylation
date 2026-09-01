# 54_marker_quasideconv.R — 变化量拟去卷积: 白细胞构成变化跨队列互证 (任务 #50 下半, 接 53)
# 动机:
#  (1) 单一 sign 对齐在 3 种血细胞混合下混淆靶细胞与其余细胞的变化;
#  (2) 4x4 NNLS 绝对组分受 HM27 vs WGBS 平台偏移污染 (观测 β 偏离参考单纯形)。
# 方法: 变化量模型 — β_mix(set,s) = c_set + Σ_k f_k(s)·R[set,k], 时间点作差后 c_set 精确抵消:
#   Δβ[set] = Δf_N·(R_N − R_B) + Δf_T·(R_T − R_B) [+ Δf_Pla·(R_Pla − R_B)]
#   GSE37722 白细胞层: 无胎盘贡献 (Δf_Pla≡0), 4 方程 2 未知量 (超定, LS);
#     其中 Placenta-HIGH 集一行成为内恰性约束 (预测 Δβ≈0 可检验)。
#   GSE154378 cfDNA: 4 方程 3 未知量, Δf_Pla = 胎盘 cfDNA 分数增量。
#   bootstrap (标记集内重抽 bin, 400 次) 给每阶段 Δf 的 95% CI。
# 输入: results/GSE37722_probe_marker_delta.csv / results/GSE154378_bin_marker_delta.csv (53 产出)
#       data/geo_methylation/GSE154378/sun2015_markers_hg38_lift.tsv
# 输出: results/GSE37722_GSE154378_quasideconv.csv (每阶段 Δf + CI + 拟合)
#       figures/GSE37722_GSE154378_cell_fraction.png
suppressPackageStartupMessages(library(data.table))
set.seed(42)

SETS  <- c("Neutrophils", "T-cells", "B-cells", "Placenta-HIGH")
CELLS <- c("Neutrophils", "T-cells", "B-cells", "Placenta")
TISSUE_COLS <- c('Liver','Lungs','Colon','Small intestines','Pancreas',
                 'Adrenal glands','Esophagus','Adipose tissues','Heart','Brain',
                 'T-cells','B-cells','Neutrophils','Placenta')
lift <- fread("data/geo_methylation/GSE154378/sun2015_markers_hg38_lift.tsv")
for (t in TISSUE_COLS) lift[[t]] <- lift[[t]] / 100
lift[, mset := fifelse(type == "I", tissue, "Placenta-HIGH")]

## 设计矩阵构造:
## 对白细胞层: Δf_N + Δf_T + Δf_B = 0  =>  Δβ = Δf_N(R_N − R_B) + Δf_T(R_T − R_B)
## 对 cfDNA:   血内 Δf 之和 = −Δf_Pla  =>  Δβ = Δf_N(R_N−R_B) + Δf_T(R_T−R_B) + Δf_Pla(R_Pla−R_B)
mkX <- function(R, with_placenta) {
  X <- cbind(N = R[, "Neutrophils"] - R[, "B-cells"],
             T = R[, "T-cells"] - R[, "B-cells"])
  if (with_placenta) X <- cbind(X, Placenta = R[, "Placenta"] - R[, "B-cells"])
  X
}

## ============ 1. 载入两队列 bin 级观测 ============
pl <- fread("results/GSE37722_probe_marker_delta.csv")[flank == 5000]
stages <- c("nulligravid","early","middle","delivery","postpartum")
b37 <- pl[, lapply(.SD, mean, na.rm = TRUE),
          .SDcols = paste0("leuk_", stages), by = .(bin, set)]
setnames(b37, paste0("leuk_", stages), stages)
b37 <- b37[set %in% SETS]
R37 <- lift[bin %in% unique(b37$bin) & mset %in% SETS,
            lapply(.SD, mean, na.rm = TRUE), .SDcols = CELLS, by = .(set = mset)]
setkey(R37, set)

bm <- fread("results/GSE154378_bin_marker_delta.csv")
tps <- c("NP_NP","Normal_1stT","Normal_2ndT","Normal_3rdT","Normal_delivery")
bm15 <- bm[set %in% SETS]
R15 <- lift[bin %in% unique(bm15$bin) & mset %in% SETS,
            lapply(.SD, mean, na.rm = TRUE), .SDcols = CELLS, by = .(set = mset)]
setkey(R15, set)

R37m <- as.matrix(R37[SETS][, CELLS, with = FALSE])   # 行按 SETS 顺序
R15m <- as.matrix(R15[SETS][, CELLS, with = FALSE])
X37 <- mkX(R37m, FALSE); X15 <- mkX(R15m, TRUE)
cat("设计矩阵 (GSE37722):\n"); print(round(X37, 3))
cat("\n设计矩阵 (GSE154378):\n"); print(round(X15, 3))

## 标记级参考表 (bin, set, CELLS)
Rlm37 <- lift[bin %in% unique(b37$bin) & mset %in% SETS,
              c("bin", "mset", CELLS), with = FALSE]
setnames(Rlm37, "mset", "set")
Rlm15 <- lift[bin %in% unique(bm15$bin) & mset %in% SETS,
              c("bin", "mset", CELLS), with = FALSE]
setnames(Rlm15, "mset", "set")

## ============ 2. 拟合函数 ============
fit_delta <- function(db, X) {
  fit <- lm.fit(X, db)
  cf <- if (!is.null(fit$coefficients)) fit$coefficients else rep(NA_real_, ncol(X))
  cf[is.na(cf)] <- 0
  pr <- as.vector(X %*% cf)
  list(cf = cf, r2 = 1 - sum((db - pr)^2) / sum((db - mean(db))^2),
       pred = pr)
}
## 每队列逐阶段: db = 集 β(s) − 集 β(baseline)
run_cohort <- function(obs_bin, Rlm, cols, with_pla, boot = TRUE, nboot = 400) {
  obs_set <- obs_bin[, lapply(.SD, mean, na.rm = TRUE), .SDcols = cols, by = set]
  setkey(obs_set, set)
  om <- as.matrix(obs_set[SETS][, cols, with = FALSE])
  Rset <- Rlm[, lapply(.SD, mean, na.rm = TRUE), .SDcols = CELLS, by = set]
  setkey(Rset, set)
  Rm <- as.matrix(Rset[SETS][, CELLS, with = FALSE])
  X <- mkX(Rm, with_pla)
  res <- vector("list", length(cols) - 1)
  names(res) <- cols[-1]
  for (j in seq_along(res)) {
    db <- om[, j + 1] - om[, 1]
    f <- fit_delta(db, X)
    out <- data.table(cell = names(f$cf), delta_f = f$cf,
                      stage = cols[j + 1], r2 = f$r2)
    if (boot) {
      B <- matrix(NA_real_, nboot, ncol(X))
      bins_by_set <- obs_bin[, .(bins = list(bin)), by = set]
      for (i in seq_len(nboot)) {
        samp <- rbindlist(lapply(SETS, function(s) {
          bb <- sample(bins_by_set[set == s, bins][[1]], replace = TRUE)
          data.table(set = s, bin = bb)
        }))
        o <- merge(samp, obs_bin, by = c("set","bin"), all.x = TRUE)
        r <- merge(samp, Rlm, by = c("set","bin"), all.x = TRUE)
        om2 <- as.matrix(o[, lapply(.SD, mean, na.rm = TRUE), .SDcols = cols, by = set][match(SETS, set)][, cols, with = FALSE])
        Rm2 <- as.matrix(r[, lapply(.SD, mean, na.rm = TRUE), .SDcols = CELLS, by = set][match(SETS, set)][, CELLS, with = FALSE])
        X2 <- mkX(Rm2, with_pla)
        db2 <- om2[, j + 1] - om2[, 1]
        B[i, ] <- fit_delta(db2, X2)$cf
      }
      q <- t(apply(B, 2, quantile, c(.025, .975), na.rm = TRUE))
      out[, `:=`(ci_lo = q[, 1], ci_hi = q[, 2])]
    }
    res[[j]] <- out
  }
  rbindlist(res)
}

res37 <- run_cohort(b37, Rlm37, stages, FALSE)
res37[, cohort := "GSE37722_leukocyte"]
res15 <- run_cohort(bm15[, c("bin","set", tps), with = FALSE], Rlm15, tps, TRUE)
res15[, cohort := "GSE154378_cfDNA"]

out <- rbind(res37, res15)
out[, cell := factor(cell, levels = c("N","T","Placenta"),
                     labels = c("Neutrophils","T-cells","Placenta"))]
setorder(out, cohort, cell, stage)
cat("\n== 变化量拟去卷积 Δf (相对基线) ==\n")
print(out)
fwrite(out, "results/GSE37722_GSE154378_quasideconv.csv")

## ============ 3. 图 ============
library(ggplot2)
out[, stage := factor(stage, levels = c(stages[-1], tps[-1]))]
p <- ggplot(out, aes(as.integer(stage), delta_f, color = cell, group = cell)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = cell), alpha = 0.15, color = NA) +
  geom_point(size = 2.2) + geom_line(linewidth = 0.9) +
  facet_wrap(~ cohort, nrow = 2, scales = "free_x") +
  scale_x_continuous(breaks = 1:8,
    labels = c("early","middle","delivery","postpartum","1stT","2ndT","3rdT","delivery")) +
  scale_color_manual(values = c(Neutrophils = "#DD8452", "T-cells" = "#4C72B0", Placenta = "#C44E52")) +
  scale_fill_manual(values = c(Neutrophils = "#DD8452", "T-cells" = "#4C72B0", Placenta = "#C44E52")) +
  labs(x = NULL, y = "Estimated delta fraction (vs baseline)",
       title = "Change-based quasi-deconvolution: leukocyte layer (GSE37722) vs cfDNA (GSE154378)",
       color = NULL) +
  theme_bw(base_size = 11)
png("figures/GSE37722_GSE154378_cell_fraction.png", width = 2400, height = 2000, res = 300)
print(p); dev.off()

cat("\n[完成] quasideconv csv + 图 1\n")
