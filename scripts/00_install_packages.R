# =====================================================================
# 00_install_packages.R
# Phase 0 环境搭建：全部 R 包安装（国内镜像源）
# 项目：妊娠期高血压甲基化研究（HDP-Methylation）
# 运行：Rscript 00_install_packages.R
# 说明：逐包检测已装则跳过，支持中断后续装；日志写入 logs/
# =====================================================================

## ---- 0. 镜像配置（清华 CRAN + 中科大 Bioconductor）----
cran_mirror <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"
bioc_mirror <- "https://mirrors.ustc.edu.cn/bioc/"
options(repos = c(CRAN = cran_mirror), BioC_mirror = bioc_mirror,
        timeout = 3600, Ncpus = max(1, parallel::detectCores() - 1))

## ---- 1. 包清单 ----
cran_pkgs <- c(
  "BiocManager", "remotes",
  # 差异分析 / 建模
  "glmnet", "xgboost", "pROC", "rms", "dcurves",
  "shapforxgboost", "caret",
  # 网络分析
  "WGCNA",
  # MR
  "TwoSampleMR", "coloc",
  # 免疫浸润辅助
  "reshape2", "ggpubr", "tidyverse"
)

bioc_pkgs <- c(
  # 数据获取与预处理
  "GEOquery", "minfi", "ChAMP", "limma", "DMRcate",
  "sva", "missMethyl", "preprocessCore",
  # 富集
  "clusterProfiler", "GSVA",
  # 免疫去卷积（甲基化）
  "EpiDISH",
  # GWAS 工具
  "MungeSumstats"
)

## ---- 2. 安装函数 ----
install_if_missing <- function(pkgs, source = c("cran", "bioc")) {
  source <- match.arg(source)
  result <- list()
  for (p in pkgs) {
    if (requireNamespace(p, quietly = TRUE)) {
      cat(sprintf("[SKIP] %s 已安装 (v%s)\n", p,
                  as.character(packageVersion(p))))
      result[[p]] <- "installed_old"
      next
    }
    cat(sprintf("[INSTALL] %s (%s) ...\n", p, source))
    ok <- tryCatch({
      if (source == "cran") {
        install.packages(p, repos = cran_mirror)
      } else {
        BiocManager::install(p, ask = FALSE, update = FALSE)
      }
      requireNamespace(p, quietly = TRUE)
    }, error = function(e) {
      cat(sprintf("[ERROR] %s: %s\n", p, conditionMessage(e))); FALSE
    })
    result[[p]] <- if (ok) "installed_new" else "FAILED"
    cat(sprintf("[%s] %s -> %s\n",
                if (ok) "OK" else "FAIL", p, result[[p]]))
  }
  unlist(result)
}

## ---- 3. 执行安装 ----
t0 <- Sys.time()
cat("==== 开始安装:", format(t0), "====\n")
cat("R version:", R.version.string, "\n\n")

cat("\n---- CRAN 包 ----\n")
res_cran <- install_if_missing(cran_pkgs, "cran")

cat("\n---- Bioconductor 包 ----\n")
res_bioc <- install_if_missing(bioc_pkgs, "bioc")

## ---- 4. MR-PRESSO（GitHub 源，可能失败不阻塞）----
if (!requireNamespace("MRPRESSO", quietly = TRUE)) {
  cat("\n---- MR-PRESSO (GitHub, 可选) ----\n")
  tryCatch(remotes::install_github("rondolab/MR-PRESSO"),
           error = function(e) cat("[WARN] MR-PRESSO 安装失败，可稍后手动安装\n"))
}

## ---- 5. 汇总 ----
all_res <- c(res_cran, res_bioc)
df <- data.frame(
  package = names(all_res),
  status = all_res,
  version = sapply(names(all_res), function(p) {
    tryCatch(as.character(packageVersion(p)), error = function(e) NA)
  })
)
cat("\n==== 安装汇总 ====\n")
print(df, row.names = FALSE)
cat("耗时:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "分钟\n")

failed <- df$package[df$status == "FAILED"]
if (length(failed) > 0) {
  cat("\n[注意] 以下包安装失败，请重跑本脚本续装：\n -",
      paste(failed, collapse = "\n - "), "\n")
  quit(status = 1)
}
cat("\n全部核心包就绪！\n")
