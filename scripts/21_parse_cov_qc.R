# 21_parse_cov_qc.R -------------------------------------------------------
# GSE282512 cfDNA WGBS .cov.gz 流式 QC 汇总 (369 文件, 8 核并行)
# cov 格式: chr  start  end  meth%  countM  countU   (无表头)
# 输出: results/GSE282512_qc_summary.csv + results/GSE282512_qc_chrom.csv
# 用法: Rscript 21_parse_cov_qc.R [n_workers]     (默认 6)
# -------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  library(data.table); library(parallel)
}))

args <- commandArgs(trailingOnly = TRUE)
n_workers <- if (length(args) >= 1) as.integer(args[1]) else 6L

proj_dir <- "E:/妊娠期高血压甲基化研究方案/hdp-methylation-project"
raw_dir  <- file.path(proj_dir, "data/geo_methylation/GSE282512_raw")
meta_csv <- file.path(proj_dir, "results/GSE282512_samples_clean.csv")
out_qc   <- file.path(proj_dir, "results/GSE282512_qc_summary.csv")
out_ch   <- file.path(proj_dir, "results/GSE282512_qc_chrom.csv")
prog_log <- file.path(proj_dir, "logs/qc_progress.log")

STD_CHR <- c(paste0("chr", 1:22), "chrX", "chrY")

meta <- fread(meta_csv)
files <- list.files(raw_dir, pattern = "\\.cov\\.gz$", full.names = TRUE)
stopifnot(length(files) > 0)
message(sprintf("待 QC 文件: %d | 并行 worker: %d", length(files), n_workers))

# 按元数据顺序对齐 (文件名含 GSM)
gsm_of <- function(f) sub("^.*(GSM\\d+)_.*\\.cov\\.gz$", "\\1", basename(f))
files <- files[match(meta$gsm, gsm_of(files))]
stopifnot(!any(is.na(files)))

qc_one <- function(f) {
  # 流式读取, 只取必要列: chr / cov 深度 / 甲基化比例
  dt <- tryCatch(
    fread(f, sep = "\t", header = FALSE, showProgress = FALSE,
          colClasses = list(character = 1L, integer = 5L, integer = 6L),
          select = c(1L, 5L, 6L)),
    error = function(e) NULL)
  if (is.null(dt)) return(list(ok = FALSE))
  setnames(dt, c("chr", "m", "u"))
  depth <- dt$m + dt$u
  is_std <- dt$chr %in% STD_CHR
  # 全局指标
  out <- list(
    gsm        = gsm_of(f),
    n_cpg_all  = nrow(dt),
    n_cpg_std  = sum(is_std),
    cov_total  = sum(as.numeric(depth)),
    mean_depth = mean(depth),
    med_depth  = as.numeric(median(depth)),
    n_cov1     = sum(depth >= 1L),
    n_cov5     = sum(depth >= 5L),
    n_cov10    = sum(depth >= 10L),
    mean_beta  = sum(dt$m) / max(sum(depth), 1L)   # 覆盖度加权全局甲基化
  )
  # 常染色体分区指标 (chr1-22, 排除 XY 与随机/替代 contigs)
  std <- dt[is_std & !(chr %in% c("chrX", "chrY"))]
  out$mean_beta_std <- sum(std$m) / max(sum(std$m + std$u), 1L)
  # 分染色体 CpG 计数 (标准染色体)
  ch_tab <- dt[is_std, .(n_cpg = .N, n_cov5 = sum(depth >= 5L)),
               by = chr]
  ch_tab[, gsm := gsm_of(f)]
  list(ok = TRUE, summary = out, chrom = ch_tab)
}

cl <- makeCluster(n_workers)
on.exit(stopCluster(cl), add = TRUE)
clusterExport(cl, c("qc_one", "gsm_of", "STD_CHR", "prog_log"))
clusterEvalQ(cl, suppressMessages(library(data.table)))

t0 <- Sys.time()
res <- parLapplyLB(cl, files, function(f) {
  r <- tryCatch(qc_one(f), error = function(e) list(ok = FALSE))
  cat(file = prog_log, append = TRUE,
      sprintf("%s %s %s\n", format(Sys.time(), "%H:%M:%S"),
              basename(f), if (isTRUE(r$ok)) "OK" else "FAIL"))
  r
})
message(sprintf("QC 完成, 耗时 %.1f 分钟", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

ok_idx <- vapply(res, function(x) isTRUE(x$ok), logical(1))
message(sprintf("成功 %d / 失败 %d", sum(ok_idx), sum(!ok_idx)))
if (any(!ok_idx)) message("失败文件: ",
  paste(basename(files[!ok_idx]), collapse = ", "))

qc <- rbindlist(lapply(res[ok_idx], function(x)
  as.data.table(x$summary)), use.names = TRUE)
qc$gsm <- gsm_of(files[ok_idx])
qc[, sample_id := meta$sample_id[match(gsm, meta$gsm)]]
qc[, cov5x_pct := round(100 * n_cov5 / pmax(n_cpg_all, 1L), 2)]
qc[, cov10x_pct := round(100 * n_cov10 / pmax(n_cpg_all, 1L), 2)]
qc[, mean_depth := round(mean_depth, 1)]

# 与元数据合并, 输出
meta_small <- meta[, .(gsm, sample_id, patient, category, severity,
                       onset, ga_weeks, tube_type, batch)]
qc <- merge(meta_small, qc, by = c("gsm", "sample_id"), all.x = TRUE)
setorder(qc, category, ga_weeks)
fwrite(qc, out_qc)
message(sprintf("输出: %s (%d 行)", out_qc, nrow(qc)))

chrom <- rbindlist(lapply(res[ok_idx], function(x) x$chrom), use.names = TRUE)
fwrite(chrom, out_ch)
message(sprintf("输出: %s (%d 行)", out_ch, nrow(chrom)))

## QC 阈值体检 --------------------------------------------------------------
message("\n===== QC 概览 =====")
print(qc[, .(
  n            = .N,
  cpg_med      = format(round(median(n_cpg_std, na.rm = TRUE), 0), big.mark = ","),
  depth_mean   = round(mean(mean_depth, na.rm = TRUE), 1),
  cov5x_pct    = round(median(cov5x_pct, na.rm = TRUE), 1),
  cov10x_pct   = round(median(cov10x_pct, na.rm = TRUE), 1),
  beta_mean    = round(mean(mean_beta_std, na.rm = TRUE), 4),
  flag_lowcpg  = sum(n_cpg_std < 1e6, na.rm = TRUE),
  flag_cov5    = sum(cov5x_pct < 50, na.rm = TRUE)
)])
message("参考: cfDNA WGBS 建议阈值 n_cpg_std>=1M 且 cov5x_pct>=50; 低于者进入剔除候选。")
