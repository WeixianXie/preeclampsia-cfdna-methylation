# 25_dmr_validate.R -----------------------------------------------------------
# 候选 DMR 位点级验证 (v3: 纯矩阵方案避免 worker 内 data.table 列选择环境问题)
# 输入:
#   results/GSE282512_dmr_candidates.csv / results/GSE282512_subcohort.csv
#   data/geo_methylation/GSE282512_raw/*.cov.gz
# 输出:
#   results/GSE282512_dmr_site_level.csv / results/GSE282512_dmr_final.csv
#   results/GSE282512_dmr_final_summary.txt / figures/dmr_site_top6.png
# 注意: 相对路径 + GBK 编码
# -----------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  library(data.table); library(parallel)
}))

RAW_DIR   <- "data/geo_methylation/GSE282512_raw"
N_WORKER  <- 6
MIN_GROUP_N <- 10
SITE_FDR  <- 0.05

## 1. 候选区域与样本 -----------------------------------------------------------
cand <- fread("results/GSE282512_dmr_candidates.csv")[candidate == TRUE]
subcoh <- fread("results/GSE282512_subcohort.csv")
f_actual <- list.files(RAW_DIR, pattern = "\\.cov\\.gz$")
cov_map  <- data.table(gsm = sub("_DNA.*$", "", f_actual), cov_file = f_actual)
subcoh <- merge(subcoh, cov_map, by = "gsm")
cat(sprintf("候选区域: %d; 样本: %d\n", nrow(cand), nrow(subcoh)))

reg <- cand[, .(region_id, chr, start, end, type, symbol)]
setkey(reg, chr, start, end)

## 2. 每样本流式提取候选区域内位点 (链级重复行按 m/u 求和聚合) -----------------
extract_sites <- function(cov_file) {
  dt <- fread(file.path(RAW_DIR, cov_file), sep = "\t", header = FALSE,
              select = c(1L, 2L, 3L, 5L, 6L),
              colClasses = list(character = 1L, integer = 2L, integer = 3L,
                                integer = 4L, integer = 5L))
  setnames(dt, c("qchr", "qstart", "qend", "m", "u"))
  setkey(dt, qchr, qstart, qend)
  ov <- foverlaps(dt, reg, type = "within", nomatch = NULL)
  # 聚合同一样本内同一位点的多条链记录; site 必须用查询侧坐标 (qstart),
  # 否则与区域侧 start 撞名后会折叠成每区域 1 个位点
  ov[, .(m = sum(as.numeric(m)), u = sum(as.numeric(u))),
     by = .(region_id, site = paste(qchr, qstart, sep = "_"))]
}

cl <- makeCluster(N_WORKER, type = "PSOCK")
clusterEvalQ(cl, suppressWarnings(suppressMessages(library(data.table))))
clusterExport(cl, c("extract_sites", "RAW_DIR", "reg"))
res <- parLapply(cl, subcoh$cov_file, function(f) {
  tryCatch(extract_sites(f), error = function(e) data.table(region_id = integer(),
    site = character(), m = numeric(), u = numeric()))
})
stopCluster(cl)
names(res) <- subcoh$gsm
long <- rbindlist(lapply(seq_along(res), function(i) {
  x <- res[[i]]; if (nrow(x)) x[, gsm := names(res)[i]]
}), use.names = TRUE, fill = TRUE)
long[, depth := m + u][, beta := ifelse(depth > 0, m / depth, NA_real_)]
cat(sprintf("候选区域内位点x样本记录: %s\n", format(nrow(long), big.mark = ",")))

## 3. 位点级统计 (纯矩阵: 主进程转 matrix, worker 只做基础索引) -----------------
gsm_vec <- subcoh$gsm
grp_vec <- subcoh$group
pe_idx  <- which(grp_vec == "PE")
ctr_idx <- which(grp_vec == "Control")

wide <- dcast(long, region_id + site ~ gsm, value.var = "beta", fill = NA_real_)
mat <- as.matrix(wide[, ..gsm_vec])
if (is.character(mat)) mat <- apply(mat, 2, as.numeric)
storage.mode(mat) <- "double"
rid_vec  <- wide$region_id
site_vec <- wide$site
cat(sprintf("位点矩阵: %d 位点 x %d 样本\n", nrow(mat), ncol(mat)))

site_stat <- function(i) {
  pe  <- mat[i, pe_idx, drop = TRUE]; pe <- pe[!is.na(pe)]
  ctr <- mat[i, ctr_idx, drop = TRUE]; ctr <- ctr[!is.na(ctr)]
  if (length(pe) < MIN_GROUP_N || length(ctr) < MIN_GROUP_N)
    return(list(region_id = rid_vec[i], site = site_vec[i],
                n_pe = length(pe), n_ctrl = length(ctr), beta_pe = NA_real_,
                beta_ctrl = NA_real_, delta = NA_real_, p = NA_real_))
  list(region_id = rid_vec[i], site = site_vec[i],
       n_pe = length(pe), n_ctrl = length(ctr),
       beta_pe = mean(pe), beta_ctrl = mean(ctr), delta = mean(pe) - mean(ctr),
       p = suppressWarnings(wilcox.test(pe, ctr)$p.value))
}

cl2 <- makeCluster(N_WORKER, type = "PSOCK")
clusterExport(cl2, c("mat", "rid_vec", "site_vec", "pe_idx", "ctr_idx",
                     "MIN_GROUP_N", "site_stat"))
st <- rbindlist(parLapply(cl2, seq_len(nrow(mat)), site_stat))
stopCluster(cl2)
st[, fdr := p.adjust(p, method = "BH")]
st[, dir_site := ifelse(delta > 0, "hyper", "hypo")]
st[, pos := as.integer(tstrsplit(site, "_", fixed = TRUE)[[2]])]
st <- merge(st, reg[, .(region_id, chr, type, symbol)], by = "region_id",
            all.x = TRUE, sort = FALSE, suffixes = c("", ".reg"))
fwrite(st, "results/GSE282512_dmr_site_level.csv")
cat(sprintf("位点级检验完成: %d 位点 (%d 显著 FDR<%.2f)\n",
            nrow(st), sum(st$fdr < SITE_FDR, na.rm = TRUE), SITE_FDR))

## 4. 区域支持度 ----------------------------------------------------------------
sup <- st[!is.na(p), .(
  n_sites = .N,
  n_sig   = sum(fdr < SITE_FDR, na.rm = TRUE),
  pct_sig = 100 * sum(fdr < SITE_FDR, na.rm = TRUE) / .N,
  dir_agree = {
    s <- sign(delta)
    if (all(is.na(s))) NA_integer_ else as.integer(sum(s > 0, na.rm = TRUE) == .N |
                                                   sum(s < 0, na.rm = TRUE) == .N)
  }
), by = region_id]
final <- merge(cand, sup, by = "region_id", all.x = TRUE)
final[, site_support := fifelse(!is.na(pct_sig) & pct_sig >= 20, "strong",
                          fifelse(!is.na(pct_sig) & pct_sig >= 5, "moderate", "weak"))]
final[, final_dmr := !is.na(pct_sig) & pct_sig >= 20 & dir_agree == 1]
setorder(final, fdr_limma)
fwrite(final, "results/GSE282512_dmr_final.csv")
cat(sprintf("最终 DMR: %d / %d 候选 (位点显著比例>=20%% 且方向一致)\n",
            sum(final$final_dmr), nrow(final)))

sink("results/GSE282512_dmr_final_summary.txt", split = TRUE)
cat("===== GSE282512 位点级验证汇总 =====\n\n")
cat(sprintf("候选区域: %d; 位点总数: %d; 位点显著(FDR<%.2f): %d\n",
            nrow(cand), nrow(st), SITE_FDR, sum(st$fdr < SITE_FDR, na.rm = TRUE)))
cat(sprintf("最终 DMR (位点显著比例>=20%% 且方向一致): %d\n\n", sum(final$final_dmr)))
cat("-- 最终 DMR 类型分布 --\n")
print(final[final_dmr == TRUE, .N, by = type][order(-N)])
cat("\n-- 最终 DMR Top 20 --\n")
print(final[final_dmr == TRUE][order(fdr_limma)][1:min(20, sum(final$final_dmr)),
  .(region_id, chr, start, end, type, symbol, beta_pe, beta_ctrl, delta_beta,
    fdr_limma, n_sites, n_sig, pct_sig)])
sink()

## 5. Top6 区域位点曲线 ---------------------------------------------------------
top6 <- final[final_dmr == TRUE][order(fdr_limma)][1:min(6, sum(final$final_dmr))]$region_id
if (length(top6)) {
  png("figures/dmr_site_top6.png", width = 2400, height = 1600, res = 150)
  par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
  for (rid in top6) {
    dd <- long[region_id == as.integer(rid)]
    dd <- merge(dd, subcoh[, .(gsm, group)], by = "gsm")
    dd[, pos := as.integer(tstrsplit(site, "_", fixed = TRUE)[[2]])]
    sm_pe  <- dd[group == "PE", .(beta = mean(beta, na.rm = TRUE)), by = pos]
    sm_ct  <- dd[group == "Control", .(beta = mean(beta, na.rm = TRUE)), by = pos]
    setorder(sm_pe, pos); setorder(sm_ct, pos)
    rr <- reg[region_id == as.integer(rid)][1]
    plot(NA, xlim = range(c(sm_pe$pos, sm_ct$pos)), ylim = c(0, 1),
         xlab = sprintf("%s:%s-%s", rr$chr, rr$start, rr$end),
         ylab = "methylation beta",
         main = sprintf("region %d (%s)", rid, ifelse(is.na(rr$symbol), "NA", rr$symbol)))
    points(sm_ct$pos, sm_ct$beta, pch = 16, cex = 0.7, col = "steelblue")
    points(sm_pe$pos, sm_pe$beta, pch = 16, cex = 0.7, col = "firebrick")
    lines(sm_ct$pos, sm_ct$beta, col = "steelblue", lwd = 1.5)
    lines(sm_pe$pos, sm_pe$beta, col = "firebrick", lwd = 1.5)
    legend("topright", legend = c("Control", "PE"), col = c("steelblue", "firebrick"),
           pch = 16, bty = "n", cex = 0.8)
  }
  dev.off()
}
cat("DONE\n")
