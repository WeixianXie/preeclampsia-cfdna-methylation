# 15_fix_dmrcate.R - 一次性补齐 DMRcate 全部 CRAN 依赖再安装
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"),
        BioC_mirror = "https://mirrors.ustc.edu.cn/bioc/",
        timeout = 3600,
        Ncpus = max(1, parallel::detectCores() - 1))

# DMRcate 全部 Imports/Depends 类 CRAN 依赖（宁多勿缺）
cran_deps <- c("gtools", "stringr", "plyr", "reshape2", "data.table",
               "latticeExtra", "dichromat", "RColorBrewer", "ggplot2")
for (p in cran_deps) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("[CRAN]", p, "\n")
    tryCatch(install.packages(p), error = function(e) cat("  fail:", conditionMessage(e), "\n"))
  }
}

# DMRcate 的 Bioc 依赖
bioc_deps <- c("bsseq", "DSS", "minfi", "limma", "GenomicRanges", "DelayedArray")
for (p in bioc_deps) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("[BIOC]", p, "\n")
    tryCatch(BiocManager::install(p, ask = FALSE, update = FALSE),
             error = function(e) cat("  fail:", conditionMessage(e), "\n"))
  }
}

cat("[INSTALL] DMRcate final\n")
tryCatch(BiocManager::install("DMRcate", ask = FALSE, update = FALSE),
         error = function(e) cat("  fail:", conditionMessage(e), "\n"))

cat("[VERIFY]\n")
ip <- rownames(installed.packages())
for (p in c(cran_deps, bioc_deps, "DMRcate")) {
  cat(p, ifelse(p %in% ip, "OK", "MISSING"), "\n")
}
cat("[DONE]\n")
