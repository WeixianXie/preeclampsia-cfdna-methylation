# 39_cellscore_deconv.R — 判别标记集评分 + 去卷积（任务 #37，标记 ±250bp 窗口级）
# 动机: 区域注释中位宽 18kb, 区域级代理稀释标记信息 (脚本 38 去卷积病态);
#      改用 39a 从 cov 直接抽取的 333 标记 ±250bp 窗口 m/u -> 窗口级 beta (保真, 可得性 ~84%)
# 分析: 1) 判别标记集评分 (免求逆): 每类细胞特异开放标记的平均 β, 低 = 占比高
#      2) EpiDISH RPC/CP 去卷积 (>=50% 样本覆盖的标记)
#      3) 12 亚型参考同法评分
#      4) 细胞比例/评分 x DMR 甲基化联动

suppressPackageStartupMessages({
  library(data.table); library(EpiDISH); library(ggplot2)
})
set.seed(42)

## ============ 1. 读入标记窗口级 beta ============
sub <- fread("results/GSE282512_subcohort.csv")
mk_files <- list.files("results/_tmp_cellmark", pattern = "\\.tsv$", full.names = TRUE)
long <- rbindlist(lapply(mk_files, function(fp) {
  gsm <- sub("\\.tsv$", "", basename(fp))
  x <- fread(fp, header = FALSE)
  if (!nrow(x)) return(NULL)
  setnames(x, c("cpg","m","u"))
  x[, gsm := gsm][, depth := m + u][, beta := ifelse(depth > 0, m / depth, NA_real_)]
}))
long <- long[!is.na(beta)]
wide <- dcast(long, cpg ~ gsm, value.var = "beta")
mm <- as.matrix(wide[, -1]); rownames(mm) <- wide$cpg
gsm_order <- intersect(names(wide)[-1], sub$gsm)
mm <- mm[, gsm_order, drop = FALSE]
sub <- sub[match(gsm_order, gsm)]
grp <- sub$group
cat(sprintf("标记矩阵: %d 标记 x %d 样本; 平均可测率: %.1f%%\n",
            nrow(mm), ncol(mm), 100 * mean(!is.na(mm))))
cov_rate <- rowMeans(!is.na(mm))

## ============ 2. 判别标记集评分 (7 亚型) ============
data("centDHSbloodDMC.m")
types <- colnames(centDHSbloodDMC.m)
mk <- rownames(mm)
score_rows <- rbindlist(lapply(types, function(t) {
  o <- setdiff(types, t)
  disc <- intersect(mk, rownames(centDHSbloodDMC.m)[centDHSbloodDMC.m[, t] < 0.2 &
    apply(centDHSbloodDMC.m[, o, drop = FALSE], 1, max) > 0.6])
  if (length(disc) < 5) return(NULL)
  sc <- colMeans(mm[disc, , drop = FALSE], na.rm = TRUE)
  d <- data.frame(y = sc, pe = as.integer(grp == "PE"), tube = factor(sub$tube_type),
                  batch = factor(sub$batch))
  cf <- summary(lm(y ~ pe + tube + batch, d))$coefficients
  data.table(cell = t, n_marker = length(disc), avail = mean(cov_rate[disc]),
             score_pe = median(sc[grp == "PE"]), score_ct = median(sc[grp == "Control"]),
             diff = median(sc[grp == "PE"]) - median(sc[grp == "Control"]),
             p_wilcox = suppressWarnings(wilcox.test(sc ~ factor(grp)))$p.value,
             beta_pe = cf["pe","Estimate"], p_adj = cf["pe","Pr(>|t|)"])
}))
score_rows[, p_bh := p.adjust(p_wilcox, "BH")]
score_rows <- score_rows[order(p_wilcox)]

## ============ 3. 12 亚型同法 ============
data("cent12CT450k.m")
types12 <- colnames(cent12CT450k.m)
score12 <- rbindlist(lapply(types12, function(t) {
  o <- setdiff(types12, t)
  disc <- intersect(mk, rownames(cent12CT450k.m)[cent12CT450k.m[, t] < 0.2 &
    apply(cent12CT450k.m[, o, drop = FALSE], 1, max) > 0.6])
  if (length(disc) < 5) return(NULL)
  sc <- colMeans(mm[disc, , drop = FALSE], na.rm = TRUE)
  d <- data.frame(y = sc, pe = as.integer(grp == "PE"), tube = factor(sub$tube_type),
                  batch = factor(sub$batch))
  cf <- summary(lm(y ~ pe + tube + batch, d))$coefficients
  data.table(cell = t, n_marker = length(disc), avail = mean(cov_rate[disc]),
             diff = median(sc[grp=="PE"]) - median(sc[grp=="Control"]),
             p_wilcox = suppressWarnings(wilcox.test(sc ~ factor(grp)))$p.value,
             beta_pe = cf["pe","Estimate"], p_adj = cf["pe","Pr(>|t|)"])
}))
score12[, p_bh := p.adjust(p_wilcox, "BH")]
score12 <- score12[order(p_wilcox)]

## ============ 4. 去卷积 (RPC/CP, >=50% 覆盖标记) ============
keep <- names(cov_rate)[cov_rate >= 0.5]
b <- mm[keep, , drop = FALSE]
ref <- centDHSbloodDMC.m[keep, , drop = FALSE]
med <- apply(b, 1, median, na.rm = TRUE)
idx <- which(is.na(b), arr.ind = TRUE); b[idx] <- med[idx[, "row"]]
cat(sprintf("去卷积标记: %d (覆盖>=50%%)\n", length(keep)))
frac_list <- list()
for (mth in c("RPC","CP")) {
  est <- tryCatch(epidish(b, ref, method = mth)$estF, error = function(e) NULL)
  if (!is.null(est)) frac_list[[mth]] <- data.frame(gsm = gsm_order, group = grp, est,
                                                    check.names = FALSE)
}
conv_tab0 <- if (length(frac_list)) rbindlist(lapply(names(frac_list), function(mn) {
  fr <- frac_list[[mn]]
  rbindlist(lapply(setdiff(names(fr), c("gsm","group")), function(cc) {
    x <- fr[[cc]]
    d <- data.frame(y = x, pe = as.integer(grp == "PE"), tube = factor(sub$tube_type),
                    batch = factor(sub$batch))
    cf <- summary(lm(y ~ pe + tube + batch, d))$coefficients
    data.table(method = mn, cell = cc, med_pe = median(x[grp=="PE"]),
               med_ct = median(x[grp=="Control"]),
               p_wilcox = suppressWarnings(wilcox.test(x ~ factor(grp)))$p.value,
               beta_pe = cf["pe","Estimate"], p_adj = cf["pe","Pr(>|t|)"])
  }))
}))
conv_tab <- if (nrow(conv_tab0)) conv_tab0 else data.table()
if (nrow(conv_tab)) { conv_tab[, p_bh := p.adjust(p_wilcox, "BH"), by = method]
  conv_tab <- conv_tab[order(method, p_wilcox)] }

## ============ 5. DMR 联动 ============
beta_reg <- fread("results/GSE282512_region_beta.csv.gz")
beta_reg <- beta_reg[, region_id := as.character(region_id)]
bmat <- as.matrix(beta_reg[, -1]); rownames(bmat) <- beta_reg$region_id
v2 <- fread("results/GSE282512_dmr_final_v2.csv", select = c("region_id","direction"))
v2[, region_id := as.character(region_id)]
v2m <- v2[region_id %in% rownames(bmat)]
g2 <- sub[match(gsm_order, gsm)] # 对齐样本列
stopifnot(identical(names(beta_reg)[-1][match(gsm_order, names(beta_reg)[-1])], gsm_order))
hyper_mean <- colMeans(bmat[v2m[direction=="hyper"]$region_id,
                            match(gsm_order, colnames(bmat)), drop = FALSE], na.rm = TRUE)
hypo_mean  <- colMeans(bmat[v2m[direction=="hypo"]$region_id,
                            match(gsm_order, colnames(bmat)), drop = FALSE], na.rm = TRUE)

link_rows <- rbindlist(lapply(types, function(t) {
  o <- setdiff(types, t)
  disc <- intersect(mk, rownames(centDHSbloodDMC.m)[centDHSbloodDMC.m[, t] < 0.2 &
    apply(centDHSbloodDMC.m[, o, drop = FALSE], 1, max) > 0.6])
  if (length(disc) < 5) return(NULL)
  sc <- colMeans(mm[disc, , drop = FALSE], na.rm = TRUE)
  rbind(data.table(cell = t, dmr = "hyper", rho = suppressWarnings(
          cor.test(sc, hyper_mean, method="spearman"))$estimate,
        p = suppressWarnings(cor.test(sc, hyper_mean, method="spearman"))$p.value),
        data.table(cell = t, dmr = "hypo", rho = suppressWarnings(
          cor.test(sc, hypo_mean, method="spearman"))$estimate,
        p = suppressWarnings(cor.test(sc, hypo_mean, method="spearman"))$p.value))
}))
link_rows[, p_bh := p.adjust(p, "BH")]

## ============ 6. 输出 ============
fwrite(score_rows, "results/GSE282512_cellscore_7type.csv")
fwrite(score12, "results/GSE282512_cellscore_12ct.csv")
fwrite(conv_tab, "results/GSE282512_deconv_refined.csv")
fwrite(link_rows, "results/GSE282512_cellscore_dmr_link.csv")

sink("results/GSE282512_cellscore_summary.txt")
cat("===== 判别标记集评分 + 去卷积 (标记 ±250bp 窗口级, cov 直接抽取) =====\n\n")
cat(sprintf("样本: PE %d / Control %d; 标记 %d 个, 平均可测率 %.1f%%\n\n",
            sum(grp=="PE"), sum(grp=="Control"), nrow(mm), 100*mean(!is.na(mm))))
cat("评分 = 细胞类型特异开放标记平均 β (越低 = 该类来源占比越高); beta_pe>0 = PE 中占比降低\n\n")
cat("---- 7 亚型标记集评分 ----\n"); print(score_rows)
cat("\n---- 12 亚型标记集评分 ----\n"); print(score12)
cat("\n---- 去卷积 (RPC/CP) ----\n"); if (nrow(conv_tab)) print(conv_tab) else cat("去卷积失败\n")
cat("\n---- 评分 x DMR 甲基化 Spearman ----\n"); print(link_rows)
sink()

lg <- rbindlist(lapply(types, function(t) {
  o <- setdiff(types, t)
  disc <- intersect(mk, rownames(centDHSbloodDMC.m)[centDHSbloodDMC.m[, t] < 0.2 &
    apply(centDHSbloodDMC.m[, o, drop = FALSE], 1, max) > 0.6])
  if (length(disc) < 5) return(NULL)
  sc <- colMeans(mm[disc, , drop = FALSE], na.rm = TRUE)
  data.table(cell = t, score = sc, group = grp)
}))
p1 <- ggplot(lg, aes(group, score, fill = group)) + geom_boxplot(width = 0.6, outlier.size = 0.6) +
  facet_wrap(~cell, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c(Control = "steelblue", PE = "firebrick")) +
  labs(x = NULL, y = "细胞类型特异标记平均 β (低=占比高)") + theme_classic(base_size = 10)
ggsave("figures/cellscore_box.png", p1, width = 8, height = 5, dpi = 300)
cat("完成: results/GSE282512_cellscore_*.csv/txt + figures/cellscore_box.png\n")
