# 35_cross_cohort_expr.R —— 机制层 ① 续: 跨队列 DMR 基因表达 meta 验证
# 胎盘队列: GSE25906 / GSE35574 / GSE44711 / GSE10588 / GSE54618 (+GSE98224 主队列结果并入)
# 母血队列: GSE85307 (早孕期外周血) / GSE48424 (PE 血)
# 输出: results/GSE282512_expr_meta.csv / GSE282512_expr_meta_summary.txt

suppressPackageStartupMessages({library(data.table); library(limma)})
set.seed(1)

## ---- 通用函数 ----
read_series <- function(path) {
  con <- gzfile(path, "rt"); meta <- list(); expr <- c(); in_tab <- FALSE
  repeat {
    ln <- readLines(con, n = 1)
    if (length(ln) == 0) break
    if (grepl("^!series_matrix_table_begin", ln)) { in_tab <- TRUE; next }
    if (grepl("^!series_matrix_table_end", ln)) { in_tab <- FALSE; next }
    if (in_tab) expr <- c(expr, ln) else {
      p <- strsplit(ln, "\t")[[1]]
      if (length(p) > 1) {
        k <- p[1]
        while (k %in% names(meta)) k <- paste0(k, "|", length(grep(paste0("^", sub("\\|", "\\\\|", p[1])), names(meta))) + 1)
        meta[[k]] <- p[-1]
      }
    }
  }
  close(con)
  cl <- function(x) if (is.null(x)) rep(NA_character_, 0) else gsub('^"|"$', "", x)
  gsm <- cl(meta[["!Sample_geo_accession"]])
  title <- cl(meta[["!Sample_title"]])
  # 所有 characteristics 行展成 field -> value 向量
  chf <- list()
  for (k in grep("!Sample_characteristics_ch1", names(meta), value = TRUE)) {
    v <- cl(meta[[k]]); fld <- sub(":.*", "", v[1])
    chf[[fld]] <- sub("^[^:]*:", "", v) |> trimws() |> (\(x) sub('"', "", x))()
  }
  mh <- strsplit(expr[1], "\t")[[1]]
  cols <- gsub('^"|"$', "", mh[-1])
  m <- t(sapply(expr[-1], function(l) {
    p <- strsplit(l, "\t")[[1]]; suppressWarnings(as.numeric(p[-1])) }))
  rownames(m) <- gsub('^"|"$', "", sapply(expr[-1], function(l) strsplit(l, "\t")[[1]][1]))
  colnames(m) <- cols
  list(gsm = gsm, title = title, char = chf, mat = m)
}

read_annot <- function(path) {
  lns <- readLines(gzfile(path))
  i <- grep("^!platform_table_begin", lns); stopifnot(length(i) == 1)
  dt <- fread(paste(lns[(i + 1):length(lns)], collapse = "\n"),
              sep = "\t", header = TRUE, quote = "", fill = TRUE)
  dt[, .(probe = as.character(ID), symbol = as.character(`Gene symbol`))]
}

run_cohort <- function(name, path, annot, group_fun, exclude = NULL, log2fix = TRUE) {
  s <- read_series(path)
  n <- length(s$gsm)
  grp <- group_fun(s, n)
  keep <- !is.na(grp) & grp %in% c("PE", "CT")
  if (!is.null(exclude)) {  # 返回需剔除样本的索引
    keep[exclude(s, n)] <- FALSE
  }
  grp <- grp[keep]; mat <- s$mat[, keep, drop = FALSE]
  if (sum(grp == "PE") < 4 || sum(grp == "CT") < 4) {
    cat(sprintf("[%s] 样本不足 (PE=%d, CT=%d), 跳过\n", name, sum(grp == "PE"), sum(grp == "CT")))
    return(NULL)
  }
  if (log2fix && max(mat, na.rm = TRUE) > 100) mat <- log2(mat + 1)
  mat <- normalizeBetweenArrays(mat, method = "quantile")
  an <- annot[probe %in% rownames(mat)]
  mat <- mat[an$probe, , drop = FALSE]
  me <- rowMeans(mat, na.rm = TRUE)
  an[, me := me]; setorder(an, -me); an <- an[!duplicated(symbol)]
  gmat <- mat[an$probe, , drop = FALSE]; rownames(gmat) <- an$symbol
  gmat <- gmat[rownames(gmat) != "" & !is.na(rownames(gmat)), , drop = FALSE]
  design <- model.matrix(~ factor(grp, levels = c("CT", "PE")))
  fit <- eBayes(lmFit(gmat, design))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  out <- data.table(symbol = rownames(gmat), logFC = tt$logFC,
                    se = abs(tt$logFC / tt$t), t = tt$t,
                    p = tt$P.Value, fdr = tt$adj.P.Val)
  # t=0 保护
  out[!is.finite(se), se := 1e6]
  cat(sprintf("[%s] PE=%d CT=%d, %d 基因\n", name, sum(grp == "PE"), sum(grp == "CT"), nrow(out)))
  out[, cohort := name][]
}

## ---- 注释 ----
an_6102  <- read_annot("data/annot/GPL6102.annot.gz")
an_10558 <- read_annot("data/annot/GPL10558.annot.gz")
an_2986  <- read_annot("data/annot/GPL2986.annot.gz")
an_6480  <- read_annot("data/annot/GPL6480.annot.gz")

GT <- "data/geo_transcriptome"
R <- function(f) file.path(GT, f)

cohs <- list()

# 胎盘队列
cohs[["GSE25906"]] <- run_cohort("GSE25906", R("GSE25906_series_matrix.txt.gz"), an_6102,
  function(s, n) {
    v <- s$char[["classification"]]; g <- rep(NA, n)
    g[grepl("preeclamptic|PE", v, ignore.case = TRUE)] <- "PE"
    g[grepl("control", v, ignore.case = TRUE)] <- "CT"; g })

cohs[["GSE35574"]] <- run_cohort("GSE35574", R("GSE35574_series_matrix.txt.gz"), an_6102,
  function(s, n) {
    v <- s$char[["classification"]]; g <- rep(NA, n)
    g[grepl("^PE$", v)] <- "PE"; g[grepl("^CONTROL$", v)] <- "CT"; g })

cohs[["GSE44711"]] <- run_cohort("GSE44711", R("GSE44711_series_matrix.txt.gz"), an_10558,
  function(s, n) {
    v <- s$char[["condition"]]; g <- rep(NA, n)
    g[grepl("EOPET|PE", v, ignore.case = TRUE)] <- "PE"
    g[grepl("control", v, ignore.case = TRUE)] <- "CT"; g })

cohs[["GSE10588"]] <- run_cohort("GSE10588", R("GSE10588_series_matrix.txt.gz"), an_2986,
  function(s, n) {
    t <- tolower(s$title); g <- rep(NA, n)
    g[grepl("preeclamps", t)] <- "PE"; g[grepl("normal", t)] <- "CT"; g })

cohs[["GSE54618"]] <- run_cohort("GSE54618", R("GSE54618_series_matrix.txt.gz"), an_10558,
  function(s, n) {
    t <- toupper(s$title); g <- rep(NA, n)
    g[grepl("PE", t)] <- "PE"; g[grepl("NORMOTENSIVE", t)] <- "CT"; g })

# 母血队列
# GSE85307 平台为 GPL6244 (Affy HuGene 1.0 ST), 用 Bioconductor 注释包映射
suppressPackageStartupMessages(library(hugene10sttranscriptcluster.db))
probes85 <- rownames(read_series(R("GSE85307-GPL6244_series_matrix.txt.gz"))$mat)
map6244 <- select(hugene10sttranscriptcluster.db, keys = probes85, columns = "SYMBOL")
an_6244 <- as.data.table(map6244)[!is.na(SYMBOL)]
setnames(an_6244, c("probe", "symbol")); an_6244[, symbol := as.character(symbol)]
cohs[["GSE85307_blood"]] <- run_cohort("GSE85307_blood", R("GSE85307-GPL6244_series_matrix.txt.gz"), an_6244,
  function(s, n) {
    v <- s$char[["pregnancy_condition"]]; g <- rep(NA, n)
    g[grepl("preeclampsia", v, ignore.case = TRUE)] <- "PE"
    g[grepl("normal", v, ignore.case = TRUE)] <- "CT"; g })

cohs[["GSE48424_blood"]] <- run_cohort("GSE48424_blood", R("GSE48424_series_matrix.txt.gz"), an_6480,
  function(s, n) {
    v <- s$char[["disease status"]]; g <- rep(NA, n)
    g[grepl("preeclampsia", v, ignore.case = TRUE)] <- "PE"
    g[grepl("healthy", v, ignore.case = TRUE)] <- "CT"; g })

## ---- 并入 GSE98224 主队列 (脚本 34 结果, 含协变量) ----
de98 <- fread("results/GSE982512_expr_gse98224_de_all.csv")
cohs[["GSE98224_placenta"]] <- de98[, .(symbol, logFC, se = abs(logFC / (t)), t, p, fdr,
                                         cohort = "GSE98224_placenta")]

## ---- DMR 基因表 ----
dmr <- fread("results/GSE282512_dmr_final_v2.csv")[candidate == TRUE]
cand <- dmr[!is.na(symbol) & symbol != ""]
genes <- unique(cand$symbol)
cat(sprintf("\nDMR 唯一基因 %d 个, 进入 meta\n", length(genes)))

## ---- 逐队列提取 DMR 基因 ----
placenta_cohorts <- c("GSE98224_placenta", "GSE25906", "GSE35574", "GSE44711", "GSE10588", "GSE54618")
blood_cohorts <- c("GSE85307_blood", "GSE48424_blood")

per_gene <- rbindlist(lapply(names(cohs), function(cn) {
  d <- cohs[[cn]][symbol %in% genes]
  if (nrow(d) == 0) return(NULL)
  d[, .(symbol, cohort, logFC, se, p)]
}), use.names = TRUE)

## ---- 固定效应 meta (胎盘队列 / 全部) ----
meta_fit <- function(d, label) {
  d <- d[is.finite(se) & se > 0]
  if (nrow(d) == 0) return(NULL)
  w <- 1 / d$se^2
  d[, w := w]
  agg <- d[, .(lf = sum(w * logFC) / sum(w), var = 1 / sum(w),
               n_coh = .N, n_pos = sum(logFC > 0), p_min = min(p)),
           by = symbol]
  agg[, z := lf / sqrt(var)]
  agg[, p_meta := 2 * pnorm(abs(z))]
  agg[, label := label][]
}
meta_pl <- meta_fit(per_gene[cohort %in% placenta_cohorts], "placenta_meta")
meta_all <- meta_fit(per_gene, "all_meta")
meta_bd <- meta_fit(per_gene[cohort %in% blood_cohorts], "blood_meta")

## ---- 联回 DMR 表 ----
cand2 <- copy(cand)
for (m in list(meta_pl, meta_bd)) {
  if (is.null(m)) next
  lab <- unique(m$label)
  cand2[, paste0("lf_", lab) := m$lf[match(symbol, m$symbol)]]
  cand2[, paste0("p_", lab) := m$p_meta[match(symbol, m$symbol)]]
  cand2[, paste0("n_", lab) := m$n_coh[match(symbol, m$symbol)]]
  cand2[, paste0("npos_", lab) := m$n_pos[match(symbol, m$symbol)]]
}
# 方向预期 (与脚本 34 相同规则)
cand2[, region_class := fifelse(grepl("promoter", type, ignore.case = TRUE), "promoter",
                         fifelse(grepl("gene|body|exon|intron", type, ignore.case = TRUE), "genebody", "other"))]
cand2[, expr_expected := fcase(
  region_class == "promoter" & direction == "hyper", -1,
  region_class == "promoter" & direction == "hypo", 1,
  region_class == "genebody" & direction == "hypo", 1,
  region_class == "genebody" & direction == "hyper", -1, default = NA_real_)]
cand2[, dir_match_pl := as.integer(expr_expected == sign(lf_placenta_meta))]
fwrite(cand2, "results/GSE282512_expr_meta.csv")
fwrite(per_gene, "results/GSE282512_expr_meta_per_cohort.csv")

## ---- 汇总 ----
sink("results/GSE282512_expr_meta_summary.txt")
cat("===== 跨队列 DMR 基因表达 meta (机制层 ①) =====\n\n")
cat("队列清单:\n")
for (cn in names(cohs)) {
  d <- cohs[[cn]]
  cat(sprintf("  %-20s %d 基因 (DMR 基因可测 %d)\n", cn, nrow(d), sum(d$symbol %in% genes)))
}
cat("\n-- 胎盘 meta --\n")
if (!is.null(meta_pl)) {
  tested <- meta_pl[symbol %in% cand2[!is.na(p_placenta_meta)]$symbol]
  cat(sprintf("DMR 基因有胎盘 meta 估计: %d\n", nrow(tested)))
  cat(sprintf("meta p<0.05: %d; FDR<0.05: %d\n",
              sum(tested$p_meta < 0.05), sum(p.adjust(tested$p_meta, "BH") < 0.05)))
  cat(sprintf("方向一致率 (多队列 logFC 同号比例 > 0.5 且 n>=2): %d/%d\n",
              sum(tested[n_coh >= 2]$npos > tested[n_coh >= 2]$n_coh / 2, na.rm = TRUE) +
              sum(tested[n_coh >= 2]$npos < tested[n_coh >= 2]$n_coh / 2, na.rm = TRUE),
              tested[n_coh >= 2, .N]))
  cat("\nTop 15 (按 meta p):\n")
  print(head(tested[order(p_meta)], 15))
}
cat("\n-- 母血 meta (GSE85307 早孕血 + GSE48424 PE 血) --\n")
if (!is.null(meta_bd)) {
  cat(sprintf("DMR 基因可测 %d; meta p<0.05: %d\n",
              nrow(meta_bd[symbol %in% cand2$symbol]), sum(meta_bd[symbol %in% cand2$symbol]$p_meta < 0.05)))
  print(head(meta_bd[order(p_meta)][1:10]))
}
cat("\n-- SLC17A1 跨队列 --\n")
print(per_gene[symbol == "SLC17A1"])
sink()
cat("DONE\n")
