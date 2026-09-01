# 26_dmr_final_call.R ---------------------------------------------------------
# 最终 DMR 判定 (v2): 基于已有位点级结果重新分级验证结论
# 背景: 位点级 FDR 无一通过 (信号弥散 + 低覆盖), 且位点方向一致率中位仅 ~44%,
#       说明区域级差异主要由覆盖度构成 (cfDNA 片段来源) 而非单 CpG 甲基化驱动。
# 判定分级 (透明报告, 不做二元剔除):
#   site_confirmed      : 位点方向一致率>=70% 且 nominal p<0.05 位点>=3
#   direction_consistent: 位点方向一致率>=70% (方向支持, 未达位点显著)
#   unconfirmed         : 位点级不支持区域级信号
# 输出:
#   results/GSE282512_dmr_final.csv / results/GSE282512_dmr_final_summary.txt
#   figures/dmr_site_top6.png (Top6 区域位点级 PE vs Control 曲线)
# 注意: 相对路径 + UTF-8 编码 (运行时 LANG=en_US.UTF-8)
# -----------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  library(data.table); library(parallel)
}))

RAW_DIR <- "data/geo_methylation/GSE282512_raw"
N_WORKER <- 6

st   <- fread("results/GSE282512_dmr_site_level.csv")
cand <- fread("results/GSE282512_dmr_candidates.csv")[candidate == TRUE]
cat(sprintf("候选区域: %d; 位点级记录: %d\n", nrow(cand), nrow(st)))

## 1. 区域级位点方向一致率与名义显著计数 ---------------------------------------
agg <- st[!is.na(p), .(
  n_tested   = .N,
  n_sig_nom  = sum(p < 0.05),
  min_p      = min(p),
  n_pos      = sum(delta > 0)
), by = region_id]
final <- merge(cand, agg, by = "region_id", all.x = TRUE)
final[, agree_pct := fifelse(!is.na(n_tested) & n_tested > 0,
                             fifelse(direction == "hyper",
                                     100 * n_pos / n_tested,
                                     100 * (n_tested - n_pos) / n_tested),
                             NA_real_)]
final[, validation :=
        fifelse(!is.na(agree_pct) & agree_pct >= 70 & n_sig_nom >= 3, "site_confirmed",
         fifelse(!is.na(agree_pct) & agree_pct >= 70, "direction_consistent",
          "unconfirmed"))]
setorder(final, fdr_limma)
fwrite(final, "results/GSE282512_dmr_final.csv")
cat(sprintf("验证分级: site_confirmed=%d, direction_consistent=%d, unconfirmed=%d\n",
            sum(final$validation == "site_confirmed"),
            sum(final$validation == "direction_consistent"),
            sum(final$validation == "unconfirmed")))

sink("results/GSE282512_dmr_final_summary.txt", split = TRUE)
cat("===== GSE282512 DMR 位点级验证汇总 (v2 分级判定) =====\n\n")
cat("一、区域级发现 (limma, tube_type+batch 协变量调整)\n")
cat(sprintf("  候选 DMR: %d (hyper %d / hypo %d), FDR<0.05\n",
            nrow(final), sum(final$direction == "hyper"), sum(final$direction == "hypo")))
cat("\n二、位点级验证 (Wilcoxon, 每 CpG 均值 beta, 未调整协变量)\n")
cat(sprintf("  可测位点: %d (每组>=10 样本覆盖), 位点 FDR<0.05: 0\n", sum(!is.na(st$p))))
cat(sprintf("  名义 p<0.05 位点: %d; 位点方向一致率中位: %.0f%%\n",
            sum(st$p < 0.05, na.rm = TRUE), median(final$agree_pct, na.rm = TRUE)))
cat("\n三、分级结论\n")
print(final[, .N, by = validation])
cat("\n四、关键方法学发现\n")
cat("  区域级(覆盖度加权)信号未在单 CpG 均值层面复现:\n")
cat("  提示差异主要由区域内覆盖度构成 (cfDNA 片段来源混合/提取方式) 驱动,\n")
cat("  而非逐位点甲基化变化。后续须做 tube_type 分层敏感性分析。\n")
cat("\n五、Top 20 候选 (按区域级 FDR)\n")
print(final[1:min(20, .N),
  .(region_id, chr, start, end, type, symbol, direction, delta_beta,
    fdr_limma, n_tested, agree_pct, n_sig_nom, validation)])
sink()

## 2. Top6 区域位点级曲线 (仅提取 Top6 区域的位点, 重新流式读取) ----------------
top6 <- final[order(fdr_limma)][1:min(6, .N)]
reg6 <- top6[, .(region_id, chr, start, end, type, symbol)]
setkey(reg6, chr, start, end)

extract_one <- function(cov_file) {
  dt <- fread(file.path(RAW_DIR, cov_file), sep = "\t", header = FALSE,
              select = c(1L, 2L, 3L, 5L, 6L),
              colClasses = list(character = 1L, integer = 2L, integer = 3L,
                                integer = 4L, integer = 5L))
  setnames(dt, c("qchr", "qstart", "qend", "m", "u"))
  setkey(dt, qchr, qstart, qend)
  ov <- foverlaps(dt, reg6, type = "within", nomatch = NULL)
  ov[, .(m = sum(as.numeric(m)), u = sum(as.numeric(u))),
     by = .(region_id, site = paste(qchr, qstart, sep = "_"))]
}

subcoh <- fread("results/GSE282512_subcohort.csv")
f_actual <- list.files(RAW_DIR, pattern = "\\.cov\\.gz$")
cov_map  <- data.table(gsm = sub("_DNA.*$", "", f_actual), cov_file = f_actual)
subcoh <- merge(subcoh, cov_map, by = "gsm")

cl <- makeCluster(N_WORKER, type = "PSOCK")
clusterEvalQ(cl, suppressWarnings(suppressMessages(library(data.table))))
clusterExport(cl, c("extract_one", "RAW_DIR", "reg6"))
res <- parLapply(cl, subcoh$cov_file, function(f) {
  tryCatch(extract_one(f), error = function(e) NULL)
})
stopCluster(cl)
names(res) <- subcoh$gsm
long6 <- rbindlist(lapply(seq_along(res), function(i) {
  x <- res[[i]]; if (!is.null(x) && nrow(x)) x[, gsm := names(res)[i]]
}), use.names = TRUE, fill = TRUE)
long6[, beta := ifelse(m + u > 0, m / (m + u), NA_real_)]
long6[, pos := as.integer(tstrsplit(site, "_", fixed = TRUE)[[2]])]
long6 <- merge(long6, subcoh[, .(gsm, group)], by = "gsm")

png("figures/dmr_site_top6.png", width = 2400, height = 1600, res = 150)
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
for (rid in top6$region_id) {
  dd <- long6[region_id == rid]
  sm_pe <- dd[group == "PE", .(beta = mean(beta, na.rm = TRUE)), by = pos]
  sm_ct <- dd[group == "Control", .(beta = mean(beta, na.rm = TRUE)), by = pos]
  setorder(sm_pe, pos); setorder(sm_ct, pos)
  rr <- top6[region_id == rid][1]
  plot(NA, xlim = range(c(sm_pe$pos, sm_ct$pos)), ylim = c(0, 1),
       xlab = sprintf("%s:%s-%s", rr$chr, rr$start, rr$end),
       ylab = "site-level mean beta",
       main = sprintf("%s (%s, %s, agree %.0f%%)",
                       rr$symbol, rr$type, rr$direction,
                       ifelse(is.na(rr$agree_pct), 0, rr$agree_pct)))
  points(sm_ct$pos, sm_ct$beta, pch = 16, cex = 0.7, col = "steelblue")
  points(sm_pe$pos, sm_pe$beta, pch = 16, cex = 0.7, col = "firebrick")
  lines(sm_ct$pos, sm_ct$beta, col = "steelblue", lwd = 1.5)
  lines(sm_pe$pos, sm_pe$beta, col = "firebrick", lwd = 1.5)
  legend("topright", legend = c("Control", "PE"), col = c("steelblue", "firebrick"),
         pch = 16, bty = "n", cex = 0.8)
}
dev.off()
cat("DONE\n")
