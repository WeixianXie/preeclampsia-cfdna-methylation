# 34_gse98224_expr_validate.R —— 机制层 ①：GSE98224 胎盘表达验证 166 候选 DMR 基因
# 输入: data/geo_transcriptome/GSE98224-GPL6244_series_matrix.txt.gz (48 胎盘, log2)
#       data/annot/GPL6244.annot.gz (探针->基因)
#       results/GSE282512_dmr_final_v2.csv (166 候选 DMR)
# 输出: results/GSE282512_expr_gse98224.csv / GSE982512_expr_validation_summary.txt

suppressPackageStartupMessages({library(data.table); library(limma)})
set.seed(1)

## ---- 1. 读 series matrix: 元数据 + 表达矩阵 ----
con <- gzfile("data/geo_transcriptome/GSE98224-GPL6244_series_matrix.txt.gz", "rt")
meta_raw <- list(); expr_lines <- c(); in_tab <- FALSE
repeat {
  ln <- readLines(con, n = 1)
  if (length(ln) == 0) break
  if (grepl("^!series_matrix_table_begin", ln)) { in_tab <- TRUE; next }
  if (grepl("^!series_matrix_table_end", ln)) { in_tab <- FALSE; next }
  if (in_tab) expr_lines <- c(expr_lines, ln) else {
    parts <- strsplit(ln, "\t")[[1]]
    if (length(parts) > 1) {
      # 同名 key (多条 characteristics) 不覆盖: 自动追加编号
      k <- parts[1]
      while (k %in% names(meta_raw)) k <- paste0(k, "|", length(grep(paste0("^", sub("\\|", "\\\\|", parts[1])), names(meta_raw))) + 1)
      meta_raw[[k]] <- parts[-1]
    }
  }
}
close(con)

gsm  <- sub('^"|"$', "", unlist(meta_raw[["!Sample_geo_accession"]]))
get_char <- function(key) {
  v <- unlist(meta_raw[[key]])
  if (is.null(v)) return(rep(NA_character_, length(gsm)))
  sub('^"|"$', "", v)
}
get_num <- function(key, pat) {
  v <- get_char(key); out <- rep(NA_real_, length(v))
  for (i in seq_along(v)) {
    m <- regmatches(v[i], regexpr(pat, v[i]))
    if (length(m) == 1) out[i] <- suppressWarnings(as.numeric(sub("[^0-9.-]", "", m)))
  }
  out
}
pheno <- data.table(
  gsm     = gsm,
  title   = get_char("!Sample_title"),
  dx      = get_char("!Sample_characteristics_ch1") |> (\(x) x[seq_along(gsm)])(),  # 占位, 下面重新解析
  ga      = rep(NA_real_, length(gsm)),
  matage  = rep(NA_real_, length(gsm)),
  sex     = rep(NA_character_, length(gsm)),
  iugr    = rep(NA_character_, length(gsm))
)
# 逐条 characteristics 行解析
char_keys <- grep("!Sample_characteristics_ch1", names(meta_raw), value = TRUE)
for (k in char_keys) {
  v <- get_char(k)
  field <- sub(":.*", "", v[1])
  vv <- sub("^[^:]*:", "", v); vv <- trimws(sub('"', "", vv))
  if (field == "diagnosis")  pheno[, dx := vv]
  if (field == "ga (week)")  pheno[, ga := suppressWarnings(as.numeric(vv))]
  if (field == "maternal age") pheno[, matage := suppressWarnings(as.numeric(vv))]
  if (field == "infant gender") pheno[, sex := vv]
  if (field == "iugr diagnosis") pheno[, iugr := vv]
}
pheno[, pe := as.integer(dx == "PE")]
cat(sprintf("GSE98224: N=%d (PE=%d, non-PE=%d)\n", nrow(pheno), sum(pheno$pe), sum(1 - pheno$pe)))

## ---- 2. 表达矩阵 + 探针->基因 ----
mhead <- strsplit(expr_lines[1], "\t")[[1]]
stopifnot(mhead[1] == '"ID_REF"')
cols <- sub('^"|"$', "", mhead[-1])
mat <- t(sapply(expr_lines[-1], function(l) {
  p <- strsplit(l, "\t")[[1]]
  suppressWarnings(as.numeric(p[-1]))
}))
rownames(mat) <- sapply(expr_lines[-1], function(l) strsplit(l, "\t")[[1]][1])
colnames(mat) <- cols
stopifnot(all(cols %in% pheno$gsm))
mat <- mat[, pheno$gsm]
cat(sprintf("表达矩阵: %d 探针 x %d 样本\n", nrow(mat), ncol(mat)))

# GPL annot 读取辅助: 定位 platform_table_begin 后以文本读入
read_annot <- function(path) {
  lns <- readLines(gzfile(path))
  i <- grep("^!platform_table_begin", lns)
  stopifnot(length(i) == 1)
  dt <- fread(paste(lns[(i + 1):length(lns)], collapse = "\n"),
              sep = "\t", header = TRUE, quote = "", fill = TRUE)
  dt
}
an <- read_annot("data/annot/GPL6244.annot.gz")
an <- an[, .(probe = ID, symbol = `Gene symbol`)]
an <- an[!is.na(symbol) & symbol != ""]
cat(sprintf("GPL6244 注释: %d 探针有基因符号\n", nrow(an)))

# 基因层折叠: 取均值表达最高的探针
keep <- rownames(mat) %in% an$probe
mat <- mat[keep, , drop = FALSE]
mean_expr <- rowMeans(mat, na.rm = TRUE)
an2 <- an[match(rownames(mat), an$probe)]
an2[, mean_expr := mean_expr]
setorder(an2, -mean_expr)
an2 <- an2[!duplicated(symbol)]
gmat <- mat[an2$probe, , drop = FALSE]
rownames(gmat) <- an2$symbol
gmat <- gmat[!is.na(rownames(gmat)) & rownames(gmat) != "", , drop = FALSE]
cat(sprintf("基因层矩阵: %d 基因\n", nrow(gmat)))

## ---- 3. limma: PE vs non-PE ----
library(data.table)
sub_ph <- copy(pheno)
sub_ph[, ga2 := ga]; sub_ph[, age2 := matage]
sub_ph[, sex2 := factor(sex)]
design <- model.matrix(~ pe + ga2 + age2 + sex2, data = sub_ph)
fit <- lmFit(gmat, design); fit <- eBayes(fit)
tt <- topTable(fit, coef = "pe", number = Inf, sort.by = "none")
de <- data.table(symbol = rownames(gmat), logFC = tt$logFC,
                 t = tt$t, p = tt$P.Value, fdr = tt$adj.P.Val)
cat(sprintf("全基因组 DE (FDR<0.05): %d / %d 基因\n",
            sum(de$fdr < 0.05), nrow(de)))
fwrite(de, "results/GSE982512_expr_gse98224_de_all.csv")

## ---- 4. DMR 基因联表 ----
dmr <- fread("results/GSE282512_dmr_final_v2.csv")
dmr[, region_id := as.integer(region_id)]
cand <- dmr[candidate == TRUE]
cand[symbol == "" | is.na(symbol), symbol := NA_character_]
cat(sprintf("候选 DMR 166 个, 有基因符号 %d (唯一基因 %d)\n",
            sum(!is.na(cand$symbol)), uniqueN(cand$symbol)))

res <- copy(cand)
res[, logFC_expr := de$logFC[match(res$symbol, de$symbol)]]
res[, p_expr := de$p[match(res$symbol, de$symbol)]]
res[, fdr_expr := de$fdr[match(res$symbol, de$symbol)]]
# 表达方向预期: 启动子区 hyper->下调, hypo->上调; 基因体 hypo->上调(常见)
res[, region_class := fifelse(grepl("promoter", type, ignore.case = TRUE), "promoter",
                       fifelse(grepl("gene|body|exon|intron", type, ignore.case = TRUE), "genebody", "other"))]
res[, expr_expected := fcase(
  region_class == "promoter" & direction == "hyper", -1,
  region_class == "promoter" & direction == "hypo", 1,
  region_class == "genebody" & direction == "hypo", 1,
  region_class == "genebody" & direction == "hyper", -1,
  default = NA_real_)]
res[, expr_sign := fifelse(is.na(logFC_expr), NA_real_, sign(logFC_expr))]
res[, expr_dir_match := as.integer(expr_expected == expr_sign)]
res[, expr_sig := as.integer(!is.na(p_expr) & p_expr < 0.05)]
res[, expr_fdr_sig := as.integer(!is.na(fdr_expr) & fdr_expr < 0.05)]
fwrite(res, "results/GSE982512_expr_gse98224.csv")

## ---- 5. 富集检验 ----
tested <- res[!is.na(p_expr)]
n_dmr_tested <- nrow(tested)
n_meas <- nrow(tested)  # 有表达测量的候选
de_any <- de[fdr < 0.05]
tab <- table(DMR_DE = tested$expr_fdr_sig == 1,
             GENOME_DE = tested$symbol %in% de_any$symbol)
fisher_all <- if (nrow(tab) > 1 && ncol(tab) > 1) fisher.test(tab) else list(estimate = NA, p.value = NA)
# 方向一致富集: 候选中方向匹配率 vs 全基因组名义 p<0.05 方向平衡
sig_nom <- tested[p_expr < 0.05]
dir_rate <- if (nrow(sig_nom) > 0) mean(sig_nom$expr_dir_match, na.rm = TRUE) else NA
bin_p <- tryCatch(binom.test(sum(sig_nom$expr_dir_match, na.rm = TRUE),
                             sum(!is.na(sig_nom$expr_dir_match)), p = 0.5)$p.value,
           error = function(e) NA)

## ---- 6. SLC17A1 专题 ----
slc <- de[symbol %in% c("SLC17A1", "SLC17A1-AS1")]
dmr_slc <- cand[!is.na(symbol) & grepl("SLC17A1", symbol)]

## ---- 7. 汇总 ----
sink("results/GSE982512_expr_validation_summary.txt")
cat("===== 机制层 ① GSE98224 胎盘表达验证 (48 例: PE vs non-PE) =====\n\n")
cat(sprintf("模型: ~ pe + 孕周 + 母龄 + 胎儿性别; 基因层矩阵 %d 基因\n", nrow(gmat)))
cat(sprintf("全基因组 DE FDR<0.05: %d; 名义 p<0.05: %d\n\n",
            sum(de$fdr < 0.05), sum(de$p < 0.05)))
cat(sprintf("候选 DMR 166 个; 注释到基因符号 %d; 表达可测 %d\n",
            sum(!is.na(cand$symbol)), n_meas))
cat(sprintf("DMR 基因中: 名义 p<0.05 %d (%.1f%%); FDR<0.05 %d\n",
            sum(tested$expr_sig), 100 * mean(tested$expr_sig),
            sum(tested$expr_fdr_sig)))
cat(sprintf("名义显著者方向符合预期 (启动子hyper->下调等): %d/%d (%.0f%%), binomial p=%.3g\n\n",
            sum(sig_nom$expr_dir_match, na.rm = TRUE),
            sum(!is.na(sig_nom$expr_dir_match)),
            100 * dir_rate, bin_p))
cat("-- top 候选 DMR 的表达 (按 p_expr) --\n")
print(tested[order(p_expr)][1:15,
      .(region_id, symbol, type, direction, delta_beta, final_call,
        logFC_expr, p_expr, fdr_expr, expr_dir_match)])
cat("\n-- 分 final_call 的表达可测/方向一致 --\n")
print(res[, .(n = .N, n_meas = sum(!is.na(p_expr)),
              n_sig = sum(expr_sig, na.rm = TRUE),
              n_fdr = sum(expr_fdr_sig, na.rm = TRUE)), by = final_call])
cat("\n-- SLC17A1 --\n")
if (nrow(slc) > 0) print(slc) else cat("SLC17A1 在 GSE98224 不可测\n")
if (nrow(dmr_slc) > 0) {
  cat("对应 DMR:\n"); print(dmr_slc[, .(region_id, chr, start, end, type, direction, delta_beta, final_call)])
}
cat("\n-- 富集检验 --\n")
cat(sprintf("DMR 基因 DE (FDR<0.05) vs 全基因组 DE: Fisher OR=%.2f, p=%.3g\n",
            fisher_all$estimate, fisher_all$p.value))
cat(sprintf("DMR 基因富集于 PE 胎盘 DE: %s\n",
            ifelse(fisher_all$p.value < 0.05, "是", "否")))
cat("\n解读要点:\n")
cat("1. 表达方向符合启动子/基因体常规预期者比例 >50% 即有方向性趋势;\n")
cat("2. 阴性/低复现本身有信息: cfDNA 甲基化变化可能反映胎盘细胞构成,\n")
cat("   而非经典启动子甲基化->转录调控回路。\n")
sink()

cat("DONE\n")
