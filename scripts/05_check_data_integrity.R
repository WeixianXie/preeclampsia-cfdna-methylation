# =====================================================================
# 05_check_data_integrity.R
# Phase 0 数据落地核验：遍历全部下载数据，检查完整性并生成核验表
# 运行: E:/R/R-4.6.1/bin/Rscript.exe scripts/05_check_data_integrity.R
# 输出: results/phase0_data_check.csv + 控制台报告
# =====================================================================

base <- "E:/妊娠期高血压甲基化研究方案/hdp-methylation-project"
dirs <- c(
  meth = file.path(base, "data/geo_methylation"),
  tx   = file.path(base, "data/geo_transcriptome"),
  gwas = file.path(base, "data/gwas"),
  mqtl = file.path(base, "data/mqtl")
)

cat("==== Phase 0 数据完整性核验 ====\n")
cat("时间:", format(Sys.time()), "\n\n")

# ---- gzip 完整性检查（读头部 + 尝试解压首行）----
check_gz <- function(path) {
  info <- file.info(path)
  size_mb <- round(info$size / 1048576, 2)
  # 文件可能正被下载进程写入（锁定）；读完立即关闭避免独占锁
  magic <- raw(2)
  open_ok <- tryCatch({
    con <- file(path, "rb")
    magic <- readBin(con, "raw", n = 2)
    close(con)
    TRUE
  }, error = function(e) FALSE)
  if (!open_ok) {
    return(list(size_mb = size_mb, is_gz = NA, stream_ok = NA,
                first_line = "文件被占用(下载中)", mtime = info$mtime))
  }
  is_gz <- (as.integer(magic[1]) == 0x1F) && (as.integer(magic[2]) == 0x8B)
  # 尝试读第一行（只捕获 error；Windows 下编码 warning 属正常）
  first_line <- ""
  stream_ok <- FALSE
  if (is_gz) {
    res <- tryCatch({
      zz <- gzfile(path, "rt"); l <- readLines(zz, n = 1); close(zz)
      first_line <- substr(l, 1, 60); TRUE
    }, error = function(e) FALSE)
    stream_ok <- suppressWarnings(res)
  }
  list(size_mb = size_mb, is_gz = is_gz, stream_ok = stream_ok,
       first_line = first_line, mtime = info$mtime)
}

# ---- series matrix 元数据提取 ----
probe_series <- function(path) {
  res <- tryCatch({
    zz <- gzfile(path, "rt"); lines <- readLines(zz, n = 200); close(zz)
    n_samples <- 0; has_disease <- FALSE; has_ga <- FALSE; platform <- ""
    for (l in lines) {
      if (grepl("^!Sample_geo_accession", l)) {
        n_samples <- length(gregexpr("\\t", l)[[1]])
      }
      if (grepl("disease state|status|diagnosis|condition", l, ignore.case = TRUE) &&
          grepl("^!Sample_characteristics", l)) has_disease <- TRUE
      if (grepl("gestational|pregnancy week|GA", l, ignore.case = TRUE) &&
          grepl("^!Sample_characteristics", l)) has_ga <- TRUE
      if (grepl("^!Series_platform_id|^!Sample_platform_id", l)) platform <- l
    }
    list(n_samples = n_samples, has_disease = has_disease, has_ga = has_ga,
         platform = substr(platform, 1, 50))
  }, error = function(e) list(n_samples = NA, has_disease = NA, has_ga = NA, platform = NA))
  res
}

rows <- list()
add_row <- function(...) rows[[length(rows) + 1]] <<- data.frame(...)

for (dname in names(dirs)) {
  d <- dirs[[dname]]
  if (!dir.exists(d)) next
  files <- list.files(d, pattern = "\\.gz$|\\.tar$", full.names = TRUE)
  for (f in files) {
    fname <- basename(f)
    cat(sprintf("[%s] %s ... ", dname, fname))
    if (grepl("\\.tar$", fname)) {
      sz <- round(file.info(f)$size / 1073741824, 2)
      cat(sprintf("%.2f GB (tar, 详细核验待解压)\n", sz))
      add_row(layer = dname, file = fname, size = paste0(sz, "GB"),
              gz_ok = "tar", readable = "待解压", n_samples = NA,
              disease_field = NA, ga_field = NA)
      next
    }
    chk <- check_gz(f)
    status <- if (is.na(chk$is_gz)) "LOCKED"
              else if (chk$is_gz && isTRUE(chk$stream_ok)) "PASS"
              else "FAIL"
    n <- NA; dis <- NA; ga <- NA
    if (identical(status, "PASS") && grepl("series_matrix", fname)) {
      pr <- tryCatch(probe_series(f),
                     error = function(e) list(n_samples = NA, has_disease = NA,
                                              has_ga = NA, platform = NA))
      n <- pr$n_samples; dis <- ifelse(isTRUE(pr$has_disease), "YES",
                                ifelse(is.na(pr$has_disease), NA, "NO"))
      ga <- ifelse(isTRUE(pr$has_ga), "YES", ifelse(is.na(pr$has_ga), NA, "NO"))
    }
    cat(sprintf("%s (%.1fMB%s)\n", status, chk$size_mb,
                if (identical(status, "PASS")) "" else
                ifelse(identical(status, "LOCKED"), ", 下载中",
                       ", 流损坏")))
    add_row(layer = dname, file = fname, size = paste0(chk$size_mb, "MB"),
            gz_ok = status, readable = chk$stream_ok, n_samples = n,
            disease_field = dis, ga_field = ga)
  }
}

if (length(rows) > 0) {
  df <- do.call(rbind, rows)
  outdir <- file.path(base, "results")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  write.csv(df, file.path(outdir, "phase0_data_check.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  cat("\n==== 核验汇总 ====\n")
  print(df, row.names = FALSE)
  n_pass <- sum(df$gz_ok == "PASS", na.rm = TRUE)
  n_fail <- sum(df$gz_ok == "FAIL", na.rm = TRUE)
  cat(sprintf("\n通过: %d  失败: %d  tar待解压: %d\n", n_pass, n_fail,
              sum(grepl("tar", df$gz_ok))))
  cat("核验表已写入: results/phase0_data_check.csv\n")
} else {
  cat("\n[WARN] 未发现数据文件\n")
}
