# =====================================================================
# 09_reinstall_missing.R
# 补装缺失包（在 Rtools 就绪、PATH 已包含编译器后运行）
# =====================================================================
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"),
        BioC_mirror = "https://mirrors.ustc.edu.cn/bioc/",
        timeout = 3600,
        Ncpus = max(1, parallel::detectCores() - 1))

cat("=== R version:", R.version.string, "===\n")
cat("=== Rtools PATH check ===\n")
cat("make: ", Sys.which("make"), "\n")
cat("gcc:  ", Sys.which("gcc"), "\n")

bioc_missing <- c("GO.db", "KEGGREST", "DMRcate", "GSVA", "MungeSumstats")
cran_missing <- c("shapforxgboost", "TwoSampleMR")

inst <- rownames(installed.packages())

for (p in bioc_missing) {
  if (p %in% inst) { cat("[SKIP]", p, "\n"); next }
  cat("[INSTALL-BIOC]", p, "\n")
  try(BiocManager::install(p, update = FALSE, ask = FALSE, force = TRUE))
}

for (p in cran_missing) {
  if (p %in% inst) { cat("[SKIP]", p, "\n"); next }
  cat("[INSTALL-CRAN]", p, "\n")
  try(install.packages(p))
}

cat("=== FINAL CHECK ===\n")
ip <- rownames(installed.packages())
key <- c(bioc_missing, cran_missing)
for (p in key) cat(p, ifelse(p %in% ip, "OK", "MISSING"), "\n")
cat("R_INSTALL2_EXIT\n")
