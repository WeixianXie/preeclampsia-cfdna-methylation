# 55_gse154378_subgroup_quasideconv.R — GSE154378 亚组级变化量去卷积 (任务 #55-58)
# 动机: 54 只做了 Normal 孕期轨迹。本脚本把同一套 Sun 2015 标记 + 变化量 LS 去卷积
#   推广到 PreX (PE/妊娠高血压) / GDM / cHTN 亚组, 检验 "相同胎盘分数窗口 (如 delivery)
#   下 PE 的母体白细胞构成是否与 Normal 偏移" —— 直接检验 "PE cfDNA DMR = 白细胞构成指纹"。
# 设计:
#   - 共同基线 = NP (7 非孕样本), 与 54 完全同构 → Normal 结果可交叉验证
#   - Δβ[set](group,t) = β[set](group,t) − β[set](NP)
#   - Δβ = Δf_N(R_N−R_B) + Δf_T(R_T−R_B) + Δf_Pla(R_Pla−R_B)   (4 方程 3 未知量, LS)
#   - bootstrap: 标记集内重抽 bin (NP + 目标亚组独立重抽, 400 次) → 95% CI
#   - 核心检验: 差值 Δ = Δf_亚组(t) − Δf_Normal(t), bootstrap 差值 CI (400 次)
# 输入: results/GSE154378_bin_mc_long.csv.gz (gsm,bin,m,c)
#       data/geo_methylation/GSE154378/gse154378_samples.tsv
#       data/geo_methylation/GSE154378/sun2015_markers_hg38_lift.tsv
# 输出: results/GSE154378_subgroup_bin_beta.csv       (亚组×时点 bin 级 β, c_min=30)
#       results/GSE154378_subgroup_deltaf.csv          (每亚组每时点 Δf + CI + r2)
#       results/GSE154378_subgroup_deltaf_contrast.csv (vs Normal 差值 + CI)
#       figures/GSE154378_subgroup_deltaf_trajectory.png
#       figures/GSE154378_subgroup_deltaf_forest.png
suppressPackageStartupMessages(library(data.table))
set.seed(42)

SETS  <- c("Neutrophils", "T-cells", "B-cells", "Placenta-HIGH")
CELLS <- c("Neutrophils", "T-cells", "B-cells", "Placenta")
TISSUE_COLS <- c('Liver','Lungs','Colon','Small intestines','Pancreas',
                 'Adrenal glands','Esophagus','Adipose tissues','Heart','Brain',
                 'T-cells','B-cells','Neutrophils','Placenta')
tp_order <- c('1stT'=1, '2ndT'=2, '3rdT'=3, 'delivery'=4)
GROUPS <- c("Normal", "PreX", "GDM", "cHTN")
STAGES <- names(tp_order)

## ---------- 1. 载入与聚合 ----------
sm <- fread('data/geo_methylation/GSE154378/gse154378_samples.tsv')   # gsm,group,subtype,patient,timepoint
mc <- fread('results/GSE154378_bin_mc_long.csv.gz')                    # gsm,bin,m,c
mc[, beta := m / c]
lift <- fread('data/geo_methylation/GSE154378/sun2015_markers_hg38_lift.tsv')
for (t in TISSUE_COLS) lift[[t]] <- lift[[t]] / 100
lift[, mset := fifelse(type == "I", tissue, "Placenta-HIGH")]

sm[, tp2 := fifelse(group == 'NP', 'NP', timepoint)]
sm <- sm[tp2 %in% c('NP', STAGES)]          # 排除 cordB 等
mc_m <- merge(mc, sm[, .(gsm, group, tp2)], by = 'gsm')

## 池化: (group, tp2, bin) → sum(m), sum(c), beta; c_min=30
agg <- mc_m[, .(m = sum(m), c = sum(c)), by = .(group, tp2, bin)][c >= 30]
agg[, beta := m / c]
setkey(agg, bin)

## 只保留 4 个标记集内的 bin, 并附 set
mk <- lift[mset %in% SETS, .(bin, set = mset)]
obs <- merge(agg, mk, by = 'bin')
cat('标记 bin 覆盖 (c>=30):\n')
print(obs[, .(n_bin = uniqueN(bin)), by = .(group, tp2)], digits = 4)

## 宽表: NP 基线 + 每亚组 (bin, set, NP_NP, <g>_1stT, ...)
dw_np <- dcast(obs[group == 'NP', .(bin, set, tp2, beta)], bin + set ~ tp2, value.var = 'beta')
setnames(dw_np, 'NP', 'NP_NP')
dcast_obs <- dw_np
for (g in GROUPS) {
  sub <- obs[group == g, .(bin, set, tp2, beta)]
  dw  <- dcast(sub, bin + set ~ tp2, value.var = 'beta')
  setnames(dw, STAGES, paste0(g, "_", STAGES))
  dcast_obs <- merge(dcast_obs, dw, by = c('bin','set'), all = TRUE)
}
fwrite(dcast_obs, 'results/GSE154378_subgroup_bin_beta.csv')
cat('\n宽表行数:', nrow(dcast_obs), ' (应为 582 左右)\n')

## 参考表 (bin 级, 只取观测到的 bin)
Rlm <- lift[mset %in% SETS, c('bin','mset', CELLS), with = FALSE]
setnames(Rlm, 'mset', 'set')

## ---------- 2. 拟合函数 (沿用 54) ----------
mkX <- function(R, with_placenta) {
  X <- cbind(N = R[, "Neutrophils"] - R[, "B-cells"],
             T = R[, "T-cells"] - R[, "B-cells"])
  if (with_placenta) X <- cbind(X, Placenta = R[, "Placenta"] - R[, "B-cells"])
  X
}
fit_delta <- function(db, X) {
  fit <- lm.fit(X, db)
  cf <- if (!is.null(fit$coefficients)) fit$coefficients else rep(NA_real_, ncol(X))
  cf[is.na(cf)] <- 0
  pr <- as.vector(X %*% cf)
  list(cf = cf, r2 = 1 - sum((db - pr)^2) / sum((db - mean(db))^2))
}

## 单亚组全时点去卷积 (基线 = NP_NP 列)
run_group <- function(dw, Rlm, grp, nboot = 400) {
  cols <- c('NP_NP', paste0(grp, "_", STAGES))
  obs_set <- dw[, lapply(.SD, mean, na.rm = TRUE), .SDcols = cols, by = set]
  setkey(obs_set, set)
  om <- as.matrix(obs_set[SETS][, cols, with = FALSE])
  Rset <- Rlm[, lapply(.SD, mean, na.rm = TRUE), .SDcols = CELLS, by = set]
  setkey(Rset, set)
  Rm <- as.matrix(Rset[SETS][, CELLS, with = FALSE])
  X <- mkX(Rm, TRUE)
  res <- vector('list', length(STAGES)); names(res) <- STAGES
  for (j in seq_along(STAGES)) {
    db <- om[, j + 1] - om[, 1]            # Δβ = β(group,t) − β(NP)
    f <- fit_delta(db, X)
    out <- data.table(cell = names(f$cf), delta_f = f$cf,
                      stage = STAGES[j], r2 = f$r2)
    ## bootstrap: 重抽 NP 与目标亚组的 bin (集内)
    B <- matrix(NA_real_, nboot, ncol(X))
    bins_by_set <- dw[, .(bins = list(bin)), by = set]
    for (i in seq_len(nboot)) {
      samp <- rbindlist(lapply(SETS, function(s) {
        bb <- sample(bins_by_set[set == s, bins][[1]], replace = TRUE)
        data.table(set = s, bin = bb)
      }))
      o2 <- merge(samp, dw, by = c('set','bin'), all.x = TRUE)
      r2t <- merge(samp, Rlm, by = c('set','bin'), all.x = TRUE)
      om2 <- as.matrix(o2[, lapply(.SD, mean, na.rm = TRUE), .SDcols = cols, by = set][match(SETS, set)][, cols, with = FALSE])
      Rm2 <- as.matrix(r2t[, lapply(.SD, mean, na.rm = TRUE), .SDcols = CELLS, by = set][match(SETS, set)][, CELLS, with = FALSE])
      X2 <- mkX(Rm2, TRUE)
      db2 <- om2[, j + 1] - om2[, 1]
      B[i, ] <- fit_delta(db2, X2)$cf
    }
    q <- t(apply(B, 2, quantile, c(.025, .975), na.rm = TRUE))
    out[, `:=`(ci_lo = q[, 1], ci_hi = q[, 2])]
    res[[j]] <- out
  }
  out <- rbindlist(res)
  out[, group := grp]
  out
}

## ---------- 3. 核心检验: 亚组 vs Normal 差值 (bootstrap 差值 CI) ----------
run_contrast <- function(dw, Rlm, grp, ref = 'Normal', nboot = 400) {
  cols_ref <- c('NP_NP', paste0(ref, "_", STAGES))
  cols_tst <- c('NP_NP', paste0(grp, "_", STAGES))
  bins_by_set <- dw[, .(bins = list(bin)), by = set]
  res <- vector('list', length(STAGES)); names(res) <- STAGES
  for (j in seq_along(STAGES)) {
    ## 点估计: 差值 = Δf(grp) − Δf(ref)
    est_deltaf <- function(cc, jj) {
      os <- dw[, lapply(.SD, mean, na.rm = TRUE), .SDcols = cc, by = set]
      setkey(os, set)
      om <- as.matrix(os[SETS][, cc, with = FALSE])
      Rs <- Rlm[, lapply(.SD, mean, na.rm = TRUE), .SDcols = CELLS, by = set]
      setkey(Rs, set)
      Rm <- as.matrix(Rs[SETS][, CELLS, with = FALSE])
      X <- mkX(Rm, TRUE)
      fit_delta(om[, jj + 1] - om[, 1], X)$cf
    }
    est_r <- est_deltaf(cols_ref, j); est_t <- est_deltaf(cols_tst, j)
    diff_pt <- est_t - est_r
    ## bootstrap
    B <- matrix(NA_real_, nboot, 3L)
    for (i in seq_len(nboot)) {
      samp <- rbindlist(lapply(SETS, function(s) {
        bb <- sample(bins_by_set[set == s, bins][[1]], replace = TRUE)
        data.table(set = s, bin = bb)
      }))
      fit1 <- function(cc, jj) {
        o2 <- merge(samp, dw, by = c('set','bin'), all.x = TRUE)
        r2t <- merge(samp, Rlm, by = c('set','bin'), all.x = TRUE)
        om2 <- as.matrix(o2[, lapply(.SD, mean, na.rm = TRUE), .SDcols = cc, by = set][match(SETS, set)][, cc, with = FALSE])
        Rm2 <- as.matrix(r2t[, lapply(.SD, mean, na.rm = TRUE), .SDcols = CELLS, by = set][match(SETS, set)][, CELLS, with = FALSE])
        X2 <- mkX(Rm2, TRUE)
        fit_delta(om2[, jj + 1] - om2[, 1], X2)$cf
      }
      B[i, ] <- fit1(cols_tst, j) - fit1(cols_ref, j)
    }
    q <- t(apply(B, 2, quantile, c(.025, .975), na.rm = TRUE))
    res[[j]] <- data.table(cell = names(est_r), diff = diff_pt,
                           stage = STAGES[j], ci_lo = q[, 1], ci_hi = q[, 2])
  }
  out <- rbindlist(res)
  out[, `:=`(contrast = paste0(grp, "_vs_", ref), group = grp)]
  out
}

## ---------- 4. 执行 ----------
cat('\n== 各亚组 Δf (相对 NP) ==\n')
lst <- lapply(GROUPS, function(g) run_group(dcast_obs, Rlm, g))
res <- rbindlist(lst)
res[, stage_ord := match(stage, STAGES)]
setorder(res, group, cell, stage_ord)
res[, stage_ord := NULL]
print(res, digits = 3)
fwrite(res, 'results/GSE154378_subgroup_deltaf.csv')

cat('\n== vs Normal 差值检验 ==\n')
lstc <- lapply(c('PreX','GDM','cHTN'), function(g) run_contrast(dcast_obs, Rlm, g))
resc <- rbindlist(lstc)
resc[, stage_ord := match(stage, STAGES)]
setorder(resc, contrast, cell, stage_ord)
resc[, stage_ord := NULL]
print(resc, digits = 3)
fwrite(resc, 'results/GSE154378_subgroup_deltaf_contrast.csv')

## ---------- 5. 图 ----------
library(ggplot2)
res[, cell := factor(cell, levels = c('N','T','Placenta'),
                     labels = c('Neutrophils','T-cells','Placenta'))]
res[, stage_f := factor(stage, levels = STAGES)]
res[, group := factor(group, levels = GROUPS)]

p1 <- ggplot(res, aes(stage_f, delta_f, color = cell, group = cell)) +
  geom_hline(yintercept = 0, linetype = 2, color = 'grey60') +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = cell), alpha = 0.15, color = NA) +
  geom_point(size = 2.2) + geom_line(linewidth = 0.9) +
  facet_wrap(~ group, nrow = 2) +
  scale_color_manual(values = c('Neutrophils' = '#DD8452', 'T-cells' = '#4C72B0', 'Placenta' = '#C44E52')) +
  scale_fill_manual(values = c('Neutrophils' = '#DD8452', 'T-cells' = '#4C72B0', 'Placenta' = '#C44E52')) +
  labs(x = NULL, y = 'Estimated delta fraction (vs NP baseline)',
       title = 'GSE154378 subgroup change-based quasi-deconvolution (Sun 2015 markers)') +
  theme_bw(base_size = 11)
png('figures/GSE154378_subgroup_deltaf_trajectory.png', width = 2400, height = 2000, res = 300)
print(p1); dev.off()

## 森林图: delivery 时点 vs Normal 差值
resc2 <- resc[stage == 'delivery']
resc2[, cell := factor(cell, levels = c('N','T','Placenta'),
                       labels = c('Neutrophils','T-cells','Placenta'))]
resc2[, contrast := factor(contrast, levels = c('PreX_vs_Normal','GDM_vs_Normal','cHTN_vs_Normal'))]
p2 <- ggplot(resc2, aes(diff, cell, color = contrast)) +
  geom_vline(xintercept = 0, linetype = 2, color = 'grey60') +
  geom_point(size = 2.6, position = position_dodge(width = 0.6)) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25,
                 position = position_dodge(width = 0.6), linewidth = 0.8) +
  scale_color_manual(values = c('PreX_vs_Normal' = '#C44E52', 'GDM_vs_Normal' = '#55A868', 'cHTN_vs_Normal' = '#8172B2')) +
  labs(x = 'ΔΔf (subgroup − Normal) at delivery', y = NULL,
       title = 'White-blood-cell composition shift vs Normal at delivery (bootstrap 95% CI)') +
  theme_bw(base_size = 11)
png('figures/GSE154378_subgroup_deltaf_forest.png', width = 2400, height = 1400, res = 300)
print(p2); dev.off()

cat('\n[完成] subgroup quasideconv csv + 2 图\n')
