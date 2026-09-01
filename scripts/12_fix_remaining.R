# =====================================================================
# 12_fix_remaining.R - 修正剩余包安装
# GO.db(依赖已就绪) / DMRcate / SHAPforxgboost(修正大小写) / TwoSampleMR(CRAN失败改GitHub)
# =====================================================================
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"),
        BioC_mirror = "https://mirrors.ustc.edu.cn/bioc/",
        timeout = 3600,
        Ncpus = max(1, parallel::detectCores() - 1))

cat("=== PATH check ===\n")
cat("gcc:", Sys.which("gcc"), "\n")

ip <- rownames(installed.packages())

# 1) GO.db（KEGGREST 已装，这次应能 lazy load）
if (!("GO.db" %in% ip)) {
  cat("[INSTALL] GO.db\n")
  try(BiocManager::install("GO.db", update = FALSE, ask = FALSE))
}

# 2) DMRcate（依赖 GO.db + DSS）
if (!("DMRcate" %in% ip)) {
  cat("[INSTALL] DMRcate\n")
  try(BiocManager::install("DMRcate", update = FALSE, ask = FALSE))
}

# 3) SHAPforxgboost（正确大小写）
if (!("SHAPforxgboost" %in% ip)) {
  cat("[INSTALL] SHAPforxgboost\n")
  try(install.packages("SHAPforxgboost"))
}

# 4) TwoSampleMR：CRAN 不可用则走 GitHub
if (!("TwoSampleMR" %in% ip)) {
  cat("[INSTALL] TwoSampleMR (CRAN)\n")
  ok <- tryCatch({ install.packages("TwoSampleMR"); TRUE },
                 error = function(e) FALSE)
  ip2 <- rownames(installed.packages())
  if (!("TwoSampleMR" %in% ip2)) {
    cat("[INSTALL] TwoSampleMR (GitHub MRCIEU)\n")
    try(remotes::install_github("MRCIEU/TwoSampleMR", upgrade = "never"))
  }
}

cat("=== FINAL CHECK ===\n")
ip3 <- rownames(installed.packages())
for (p in c("GO.db", "DMRcate", "SHAPforxgboost", "TwoSampleMR", "clusterProfiler", "GSVA", "EpiDISH", "MungeSumstats")) {
  cat(p, ifelse(p %in% ip3, "OK", "MISSING"), "\n")
}
cat("R_FIX_EXIT\n")
