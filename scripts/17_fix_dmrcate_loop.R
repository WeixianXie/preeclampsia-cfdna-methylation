# 17_fix_dmrcate_loop.R - 自动循环：解析报错中缺失的包名 -> 安装 -> 重试（最多 15 轮）
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"),
        BioC_mirror = "https://mirrors.ustc.edu.cn/bioc/",
        timeout = 3600,
        Ncpus = max(1, parallel::detectCores() - 1))

log_file <- "E:/妊娠期高血压甲基化研究方案/hdp-methylation-project/logs/fix_dmrcate_loop.log"

for (round in 1:15) {
  cat(sprintf("[ROUND %d] attempt DMRcate\n", round))
  # 用 callr 风格的子进程安装，捕获输出
  out <- tryCatch({
    utils::capture.output(
      BiocManager::install("DMRcate", ask = FALSE, update = FALSE)
    )
  }, error = function(e) c("ERROR:", conditionMessage(e)))

  if (requireNamespace("DMRcate", quietly = TRUE)) {
    cat("[SUCCESS] DMRcate installed at round", round, "\n")
    break
  }

  combined <- paste(out, collapse = "\n")
  # 匹配 "there is no package called 'X'"
  m <- regmatches(combined, regexpr("there is no package called '([^']+)'", combined))
  if (length(m) > 0) {
    misspkg <- sub("there is no package called '([^']+)'", "\\1", m)
    cat(sprintf("[ROUND %d] missing dep: %s -> installing\n", round, misspkg))
    tryCatch(install.packages(misspkg), error = function(e) {
      # CRAN 失败则尝试 Bioconductor
      cat("  CRAN fail, try Bioc\n")
      tryCatch(BiocManager::install(misspkg, ask = FALSE, update = FALSE),
               error = function(e2) cat("  Bioc fail too:", conditionMessage(e2), "\n"))
    })
  } else {
    cat("[ROUND %d] no recognizable missing dep; output tail:\n", round)
    cat(tail(out, 15), sep = "\n")
  }
}

cat("[FINAL VERIFY]\n")
ip <- rownames(installed.packages())
cat("DMRcate", ifelse("DMRcate" %in% ip, "OK", "MISSING"), "\n")
cat("[DONE]\n")
