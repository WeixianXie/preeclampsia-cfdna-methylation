# 23_build_region_matrix.R ---------------------------------------------------
# 构建 hg38 区域注释 + GSE282512 子队列区域聚合甲基化矩阵
# 区域类型:
#   gene_body         TxDb gene body
#   promoter_TSS2kb   TxDb gene TSS ±2kb
#   regbuild_promoter Ensembl Regulatory Build promoter 特征 (hg38)
# 输入: data/annot/regbuild.gff.gz; results/GSE282512_subcohort.csv
#       results/GSE282512_samples_clean.csv; data/geo_methylation/GSE282512_raw/*.cov.gz
# 输出:
#   results/GSE282512_region_annot.csv          区域注释 (region_id/chr/start/end/type/gene_id/symbol/reg_id)
#   results/GSE282512_region_counts_long.csv.gz 长表 (region_id/gsm/m/u)
#   results/GSE282512_region_beta.csv.gz        宽矩阵 β=(m+1)/(m+u+2)
#   results/GSE282512_region_counts.RDS         R 对象 (含 m/u 矩阵)
# 注意: Windows R 无法 fread 中文绝对路径 -> 全部相对路径; 脚本以 GBK 保存
# -----------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  library(data.table)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(parallel)
}))

CHR_MAIN <- paste0("chr", c(1:22, "X", "Y"))
RAW_DIR  <- "data/geo_methylation/GSE282512_raw"
N_WORKER <- 6

## 1. TxDb 基因注释 -> gene body + promoter TSS±2kb ---------------------------
txdb  <- TxDb.Hsapiens.UCSC.hg38.knownGene
gn    <- keepSeqlevels(genes(txdb), CHR_MAIN, pruning.mode = "coarse")
gn    <- sort(gn)

gb <- data.table(
  chr = as.character(seqnames(gn)), start = start(gn), end = end(gn),
  type = "gene_body", gene_id = names(gn))
pr <- promoters(gn, upstream = 2000, downstream = 2000)
prd <- data.table(
  chr = as.character(seqnames(pr)), start = start(pr), end = end(pr),
  type = "promoter_TSS2kb", gene_id = names(pr))

## 2. Ensembl Regulatory Build promoter 特征 ----------------------------------
gff <- fread("data/annot/regbuild.gff.gz", sep = "\t", header = FALSE, fill = TRUE,
             colClasses = rep("character", 9))
gff <- gff[!grepl("^#", V1)]
gff <- gff[V3 == "promoter"]
gff_dt <- data.table(
  chr = paste0("chr", gff$V1),
  start = as.integer(gff$V4), end = as.integer(gff$V5),
  type = "regbuild_promoter", gene_id = NA_character_,
  reg_id = sub("^ID=promoter:", "", sub(";.*", "", gff$V9)))
gff_dt <- gff_dt[chr %in% CHR_MAIN]

## 3. 合并区域表 --------------------------------------------------------------
regions <- rbind(gb, prd, gff_dt, use.names = TRUE, fill = TRUE)
regions[, region_id := .I]
setkey(regions, chr, start, end)
cat(sprintf("区域总数: %d (gene_body %d, promoter_TSS2kb %d, regbuild_promoter %d)\n",
            nrow(regions),
            sum(regions$type == "gene_body"),
            sum(regions$type == "promoter_TSS2kb"),
            sum(regions$type == "regbuild_promoter")))

## gene symbol 注释
sym <- AnnotationDbi::select(org.Hs.eg.db,
                             keys = unique(regions[!is.na(gene_id)]$gene_id),
                             columns = "SYMBOL", keytype = "ENTREZID")
sym <- data.table(gene_id = sym$ENTREZID, symbol = sym$SYMBOL)
sym <- sym[!duplicated(gene_id)]
regions <- merge(regions, sym, by = "gene_id", all.x = TRUE, sort = FALSE)
setkey(regions, chr, start, end)
fwrite(regions, "results/GSE282512_region_annot.csv")
cat("区域注释已写: results/GSE282512_region_annot.csv\n")

## 4. 子队列样本与 cov 文件映射 ----------------------------------------------
# 注意: tar 内实际文件名 GSMxxxx_DNAyyyy.cov.gz 与 series matrix 的
# supplementary 文件名 (GSMxxxx.cov.gz) 不一致 -> 按 GSM 前缀映射实际文件
meta     <- fread("results/GSE282512_samples_clean.csv")
f_actual <- list.files(RAW_DIR, pattern = "\\.cov\\.gz$")
cov_map  <- data.table(gsm = sub("_DNA.*$", "", f_actual), cov_file = f_actual)
sub <- merge(fread("results/GSE282512_subcohort.csv"), cov_map, by = "gsm")
cat(sprintf("子队列样本: %d (PE %d, Control %d)\n", nrow(sub),
            sum(sub$group == "PE"), sum(sub$group == "Control")))

## 5. 单样本区域聚合 ----------------------------------------------------------
agg_one <- function(cov_file) {
  dt <- fread(file.path(RAW_DIR, cov_file), sep = "\t", header = FALSE,
              select = c(1L, 2L, 3L, 5L, 6L),
              colClasses = list(character = 1L, integer = 2L, integer = 3L,
                                integer = 4L, integer = 5L))
  setnames(dt, c("chr", "start", "end", "m", "u"))
  dt <- dt[chr %in% CHR_MAIN]
  if (!nrow(dt)) return(data.table(region_id = integer(), m = numeric(), u = numeric()))
  setkey(dt, chr, start, end)
  ov <- foverlaps(dt, regions, type = "within", nomatch = NULL)
  ov[, .(m = sum(as.numeric(m)), u = sum(as.numeric(u))), by = region_id]
}

cl <- makeCluster(N_WORKER, type = "PSOCK")
clusterEvalQ(cl, { suppressWarnings(suppressMessages(library(data.table))) })
clusterExport(cl, c("agg_one", "RAW_DIR", "CHR_MAIN", "regions"))
res <- parLapply(cl, sub$cov_file, function(f) {
  tryCatch(agg_one(f), error = function(e) data.table(region_id = integer(), m = NA_real_, u = NA_real_))
})
stopCluster(cl)
names(res) <- sub$gsm
ok <- vapply(res, nrow, integer(1)) > 0
cat(sprintf("聚合成功样本: %d / %d\n", sum(ok), length(res)))
if (!all(ok)) cat("失败样本:", paste(sub$gsm[!ok], collapse = ", "), "\n")

## 6. 长表 + 宽 beta 矩阵 ------------------------------------------------------
long <- rbindlist(lapply(seq_along(res), function(i) {
  x <- res[[i]]; x[, gsm := names(res)[i]]
}), use.names = TRUE, fill = TRUE)
long <- long[!is.na(m)]
long[, beta := (m + 1) / (m + u + 2)]
fwrite(long, "results/GSE282512_region_counts_long.csv.gz")
cat(sprintf("长表已写: %s 行\n", format(nrow(long), big.mark = ",")))

wide <- dcast(long, region_id ~ gsm, value.var = "beta")
fwrite(wide, "results/GSE282512_region_beta.csv.gz")
cat("宽矩阵 beta 已写: results/GSE282512_region_beta.csv.gz\n")

m_mat <- dcast(long, region_id ~ gsm, value.var = "m")
u_mat <- dcast(long, region_id ~ gsm, value.var = "u")
saveRDS(list(regions = regions,
             long    = long,
             m_mat   = m_mat,
             u_mat   = u_mat,
             subcohort = sub),
        "results/GSE282512_region_counts.RDS")
cat("RDS 已写: results/GSE282512_region_counts.RDS\n")
cat("DONE\n")
