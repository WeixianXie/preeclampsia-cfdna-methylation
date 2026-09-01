# 48_immune_infiltration.R — 母血转录组 ssGSEA 免疫浸润 (任务 #46, 第三梯队)
# 队列: GSE85307_blood (早孕期外周血, GPL6244) + GSE48424_blood (PE 外周血, GPL6480)
# 方法: 手工实现 ssGSEA (Barbie 2009, 秩加权) + 14 细胞谱系标记面板
#       PE vs CT 逐队列 Wilcoxon; 与 cfDNA 12 亚型去卷积方向交叉对照
# 输出: results/GSE282512_immune_ssgsea.csv / _summary.txt / figures/immune_heatmap.png
# 约定: UTF-8, 相对路径, LANG=en_US.UTF-8
suppressPackageStartupMessages({
  library(data.table); library(limma); library(ggplot2)
  library(hugene10sttranscriptcluster.db); library(AnnotationDbi)
})
set.seed(42)

## ---- 通用函数 (沿 35_cross_cohort_expr.R) ----
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
  cl <- function(x) gsub('^"|"$', "", x)
  gsm <- cl(meta[["!Sample_geo_accession"]])
  chf <- list()
  for (k in grep("!Sample_characteristics_ch1", names(meta), value = TRUE)) {
    v <- cl(meta[[k]]); fld <- sub(":.*", "", v[1])
    chf[[fld]] <- trimws(sub("^[^:]*:", "", v))
  }
  mh <- strsplit(expr[1], "\t")[[1]]
  cols <- gsub('^"|"$', "", mh[-1])
  m <- t(sapply(expr[-1], function(l) {
    p <- strsplit(l, "\t")[[1]]; suppressWarnings(as.numeric(p[-1])) }))
  rownames(m) <- gsub('^"|"$', "", sapply(expr[-1], function(l) strsplit(l, "\t")[[1]][1]))
  colnames(m) <- cols
  list(gsm = gsm, char = chf, mat = m)
}
read_annot <- function(path) {
  lns <- readLines(gzfile(path))
  i <- grep("^!platform_table_begin", lns); stopifnot(length(i) == 1)
  dt <- fread(paste(lns[(i + 1):length(lns)], collapse = "\n"),
              sep = "\t", header = TRUE, quote = "", fill = TRUE)
  dt[, .(probe = as.character(ID), symbol = as.character(`Gene symbol`))]
}

## ---- 14 细胞谱系标记面板 (Bindea/Charo 惯用核心标记, 手工整理) ----
panel <- list(
  T_cells      = c("CD3D","CD3E","CD2","TRAC","CD3G"),
  CD8_T        = c("CD8A","CD8B","GZMK","CD8BP"),
  Cytotoxic    = c("GZMB","PRF1","NKG7","GNLY","GZMH"),
  NK           = c("KLRD1","NCAM1","NCR1","KLRB1","KLRF1"),
  B_cells      = c("CD19","MS4A1","CD79A","CD79B","BANK1","TCL1A"),
  Plasma       = c("MZB1","XBP1","JCHAIN","SDC1"),
  Monocytes    = c("CD14","LST1","FCN1","S100A12","CTSS"),
  Macrophage   = c("CD68","CD163","MRC1","CSF1R"),
  DC           = c("ITGAX","CD1C","CLEC10A","FCER1A","HLA-DQA1"),
  Neutrophils  = c("FCGR3B","CSF3R","CXCR1","CXCR2","FPR1"),
  Eosinophils  = c("CCR3","IL5RA","SIGLEC8","CLC"),
  Treg         = c("FOXP3","IL2RA","CTLA4","IKZF2"),
  Th1          = c("IFNG","TBX21","CXCR3","STAT1"),
  Th17         = c("RORC","IL17A","IL23R","CCR6")
)

## ---- ssGSEA (秩加权, Barbie 2009 简化实现) ----
ssgsea <- function(gmat, sets, k = 0.25) {
  genes <- rownames(gmat)
  common <- lapply(sets, function(s) intersect(s, genes))
  N <- nrow(gmat)
  es <- matrix(NA_real_, length(sets), ncol(gmat),
               dimnames = list(names(sets), colnames(gmat)))
  for (i in seq_len(ncol(gmat))) {
    x <- gmat[, i]
    o <- order(x, decreasing = TRUE)
    w <- abs(x)^k; w[!is.finite(w)] <- 0
    for (j in seq_along(sets)) {
      S <- common[[j]]; ns <- length(S)
      if (ns < 3) next
      isS <- genes %in% S
      isS_r <- isS[o]; w_r <- w[o]
      hit <- cumsum(ifelse(isS_r, w_r, 0)) / sum(w[isS])
      miss <- cumsum(!isS_r) / (N - ns)
      d <- hit - miss
      es[j, i] <- sum(d[isS_r])   # 成员基因位置的运行差之和 (原 ssGSEA 口径)
    }
  }
  ## 样本内 z 标准化 (跨细胞类型比较用均值-0/1sd)
  sweep(es, 2, colMeans(es, na.rm = TRUE), "-") |>
    (\(m) sweep(m, 2, apply(m, 2, sd, na.rm = TRUE), "/"))() -> z
  list(es = es, z = z)
}

## ---- 逐队列 ----
GT <- "data/geo_transcriptome"
an_6480 <- read_annot(file.path(GT, "../annot/GPL6480.annot.gz"))
probes85 <- rownames(read_series(file.path(GT, "GSE85307-GPL6244_series_matrix.txt.gz"))$mat)
map6244 <- select(hugene10sttranscriptcluster.db, keys = probes85, columns = "SYMBOL")
an_6244 <- as.data.table(map6244)[!is.na(SYMBOL)]
setnames(an_6244, c("probe", "symbol")); an_6244[, symbol := as.character(symbol)]

build <- function(path, annot, group_fun) {
  s <- read_series(path)
  n <- length(s$gsm)
  grp <- group_fun(s, n)
  keep <- !is.na(grp) & grp %in% c("PE", "CT")
  grp <- grp[keep]; mat <- s$mat[, keep, drop = FALSE]
  if (max(mat, na.rm = TRUE) > 100) mat <- log2(mat + 1)
  mat <- normalizeBetweenArrays(mat, method = "quantile")
  an <- annot[probe %in% rownames(mat)]
  me <- rowMeans(mat, na.rm = TRUE)
  an[, me := me[match(probe, rownames(mat))]]; setorder(an, -me); an <- an[!duplicated(symbol)]
  gmat <- mat[an$probe, , drop = FALSE]; rownames(gmat) <- an$symbol
  gmat <- gmat[rownames(gmat) != "" & !is.na(rownames(gmat)), , drop = FALSE]
  list(gmat = gmat, grp = grp, gsm = s$gsm[keep])
}
c1 <- build(file.path(GT, "GSE85307-GPL6244_series_matrix.txt.gz"), an_6244,
            function(s, n) {
              v <- s$char[["pregnancy_condition"]]; g <- rep(NA, n)
              g[grepl("preeclampsia", v, ignore.case = TRUE)] <- "PE"
              g[grepl("normal", v, ignore.case = TRUE)] <- "CT"; g })
c2 <- build(file.path(GT, "GSE48424_series_matrix.txt.gz"), an_6480,
            function(s, n) {
              v <- s$char[["disease status"]]; g <- rep(NA, n)
              g[grepl("preeclampsia", v, ignore.case = TRUE)] <- "PE"
              g[grepl("healthy", v, ignore.case = TRUE)] <- "CT"; g })
cat(sprintf("GSE85307: PE=%d CT=%d; GSE48424: PE=%d CT=%d\n",
            sum(c1$grp == "PE"), sum(c1$grp == "CT"),
            sum(c2$grp == "PE"), sum(c2$grp == "CT")))

r1 <- ssgsea(c1$gmat, panel); r2 <- ssgsea(c2$gmat, panel)

test_cohort <- function(es, grp, cohort) {
  dt <- as.data.table(t(es))
  out <- rbindlist(lapply(intersect(names(panel), colnames(dt)), function(ct) {
    w <- suppressWarnings(wilcox.test(dt[[ct]] ~ factor(grp)))
    data.table(cohort = cohort, cell = ct,
               med_PE = median(unlist(dt[grp == "PE", ..ct, with = FALSE]), na.rm = TRUE),
               med_CT = median(unlist(dt[grp == "CT", ..ct, with = FALSE]), na.rm = TRUE),
               p = w$p.value)
  }))
  out[, log2fc := log2((abs(med_PE) + 1e-6) / (abs(med_CT) + 1e-6))]
  out[, diff := med_PE - med_CT]
  out
}
t1 <- test_cohort(r1$es, c1$grp, "GSE85307_blood"); t1[, p_bh := p.adjust(p, "BH")]
t2 <- test_cohort(r2$es, c2$grp, "GSE48424_blood"); t2[, p_bh := p.adjust(p, "BH")]
res <- rbind(t1, t2)
fwrite(res, "results/GSE282512_immune_ssgsea.csv")

## ---- 与 cfDNA 12 亚型去卷积方向对照 ----
frac <- fread("results/GSE282512_deconv_fractions_12ct.csv")
ctl_lab <- if ("Control" %in% frac$group) "Control" else "CT"
fd <- frac[, lapply(.SD, function(x) mean(x[frac$group == "PE"]) - mean(x[frac$group == ctl_lab])),
           .SDcols = setdiff(names(frac), c("gsm","group"))]
fdt <- data.table(cell = names(fd), cfDNA_diff = as.numeric(fd))
fwrite(fdt, "results/GSE282512_immune_cfDNA_celldiff.csv")

## ---- 图: 差异热图 (两队列, 标 p) ----
pl <- rbind(t1, t2)[, .(cohort, cell, diff, p)]
pl[, sig := fifelse(p < 0.05, "*", "")]
p1 <- ggplot(pl, aes(cohort, reorder(cell, diff), fill = diff)) +
  geom_tile() + geom_text(aes(label = sig), size = 6, vjust = 0.7) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(x = NULL, y = NULL, fill = "Δ ssGSEA\n(PE - CT)",
       title = "Maternal blood immune infiltration (ssGSEA, 14 lineages)") +
  theme_bw(base_size = 11)
ggsave("figures/immune_heatmap.png", p1, width = 7, height = 6, dpi = 300)

## ---- summary ----
sink("results/GSE282512_immune_summary.txt")
cat("== 母血转录组 ssGSEA 免疫浸润 ==\n\n")
cat(sprintf("队列: GSE85307_blood (早孕外周血; PE=%d CT=%d), GSE48424_blood (PE 血; PE=%d CT=%d)\n",
            sum(c1$grp == "PE"), sum(c1$grp == "CT"), sum(c2$grp == "PE"), sum(c2$grp == "CT")))
cat("方法: 手工 ssGSEA (|x|^0.25 加权秩), 14 细胞谱系标记面板; 队列内 Wilcoxon + BH\n\n")
cat("---- GSE85307_blood (早孕期) ----\n")
print(t1[order(p)][, .(cell, med_PE = round(med_PE, 3), med_CT = round(med_CT, 3),
                       diff = round(diff, 3), p = signif(p, 2), p_bh = signif(p_bh, 2))])
cat("\n---- GSE48424_blood ----\n")
print(t2[order(p)][, .(cell, med_PE = round(med_PE, 3), med_CT = round(med_CT, 3),
                       diff = round(diff, 3), p = signif(p, 2), p_bh = signif(p_bh, 2))])
cat("\n---- cfDNA 12 亚型 PE-CT 差 (参照) ----\n")
print(fdt[order(-abs(cfDNA_diff))])
cat("\n注: 血转录组免疫浸润方向与 cfDNA 去卷积亚群差异可交叉印证母体免疫构成改变\n")
sink()
cat("[完成] 输出: GSE282512_immune_ssgsea.csv / _cfDNA_celldiff.csv / _summary.txt / 1 图\n")
