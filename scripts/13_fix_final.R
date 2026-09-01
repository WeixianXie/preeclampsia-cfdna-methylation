# 13_fix_final.R - 安装最后两个缺失包 DMRcate 和 TwoSampleMR
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"),
        BioC_mirror = "https://mirrors.ustc.edu.cn/bioc/",
        timeout = 3600,
        Ncpus = max(1, parallel::detectCores() - 1))

cat("[STEP] dichromat (DMRcate 依赖)\n")
if (!requireNamespace("dichromat", quietly = TRUE)) install.packages("dichromat")

cat("[INSTALL] DMRcate\n")
if (!requireNamespace("DMRcate", quietly = TRUE)) BiocManager::install("DMRcate", ask = FALSE, update = FALSE)

cat("[INSTALL] TwoSampleMR CRAN\n")
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  tryCatch(install.packages("TwoSampleMR"),
           error = function(e) cat("CRAN fail:", conditionMessage(e), "\n"))
}

cat("[INSTALL] TwoSampleMR GitHub fallback\n")
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  tryCatch(remotes::install_github("MRCIEU/TwoSampleMR"),
           error = function(e) cat("GitHub fail:", conditionMessage(e), "\n"))
}

cat("[VERIFY]\n")
ip <- rownames(installed.packages())
for (p in c("dichromat", "DMRcate", "TwoSampleMR")) {
  cat(p, ifelse(p %in% ip, "OK", "MISSING"), "\n")
}
cat("[DONE]\n")
