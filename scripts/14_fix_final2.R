# 14_fix_final2.R - 补装 latticeExtra -> DMRcate -> TwoSampleMR
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"),
        BioC_mirror = "https://mirrors.ustc.edu.cn/bioc/",
        timeout = 3600,
        Ncpus = max(1, parallel::detectCores() - 1))

cat("[STEP] latticeExtra\n")
if (!requireNamespace("latticeExtra", quietly = TRUE)) install.packages("latticeExtra")

cat("[INSTALL] DMRcate\n")
if (!requireNamespace("DMRcate", quietly = TRUE)) BiocManager::install("DMRcate", ask = FALSE, update = FALSE)

cat("[INSTALL] TwoSampleMR (GitHub MRCIEU)\n")
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  tryCatch(remotes::install_github("MRCIEU/TwoSampleMR@main"),
           error = function(e) cat("GitHub fail:", conditionMessage(e), "\n"))
}

cat("[VERIFY]\n")
ip <- rownames(installed.packages())
for (p in c("latticeExtra", "DMRcate", "TwoSampleMR")) {
  cat(p, ifelse(p %in% ip, "OK", "MISSING"), "\n")
}
cat("[DONE]\n")
