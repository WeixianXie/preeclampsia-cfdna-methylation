# 32_phase3_lift_tri.R ----------------------------------------------------------
# Phase 3 补充: FinnGen 变异 GRCh38->GRCh37 回映 (Ensembl REST), 与 meQTL GRCh37 精确联表
# 背景: GoDMC 坐标为 GRCh37 (cis 距离检验 + Ensembl 映射双重确认),
#       Phase 2 的 CpG-in-DMR 重叠用 hg38 manifest 不受影响, 但 SNP 联表须同版本
# 输入:
#   results/phase3_lift_targets.csv  DMR±1.5Mb 窗口内 FinnGen p<1e-4 变异 (hg38)
#   results/GSE282512_phase2_dmr_meqtl.csv  DMR 内 meQTL (SNP 坐标 GRCh37)
#   data/gwas/finngen_pe_p1e4.tsv / finngen_gh_p1e4.tsv
# 输出:
#   results/phase3_lift_map.csv          变异 hg38 -> GRCh37 映射表
#   results/phase3_shared_mqtl_gwas.csv  meQTL 与 GWAS 共享变异 (三角证据)
# 约定: UTF-8, 相对路径, LANG=en_US.UTF-8
# ------------------------------------------------------------------------------
suppressWarnings(suppressMessages({ library(data.table); library(jsonlite); library(httr) }))

tg <- fread("results/phase3_lift_targets.csv")
cat(sprintf("待回映变异: %d\n", nrow(tg)))

## 1. Ensembl REST 逐条映射 (限速, 已有回映表则跳过) --------------------------------
if (file.exists("results/phase3_lift_map.csv")) {
  lift <- fread("results/phase3_lift_map.csv")
  lift[, chrom_hg38 := as.character(chrom_hg38)]
  ok <- nrow(lift); fail <- 0
  cat(sprintf("已存在回映表 %d 行, 跳过查询\n", nrow(lift)))
} else {
BASE <- "https://rest.ensembl.org/map/human/GRCh38/%s:%d:%d/GRCh37?content-type=application/json"
fetch_json <- function(url) {
  tryCatch(fromJSON(content(GET(url, timeout(30)), as = "text")),
           error = function(e) NULL)
}
lift <- data.table()
ok <- 0; fail <- 0
for (i in seq_len(nrow(tg))) {
  ch <- tg$chrom[i]; p <- tg$pos[i]
  url <- sprintf(BASE, ch, p, p + 1)
  js <- fetch_json(url)
  if (is.null(js)) { Sys.sleep(0.5); js <- fetch_json(url) }
  if (!is.null(js) && !is.null(js$mappings)) {
    m <- js$mappings$mapped
    lift <- rbind(lift, data.table(chrom_hg38 = ch, pos_hg38 = p,
                                   chrom_hg19 = as.character(m$seq_region_name[1]),
                                   pos_hg19 = as.integer(m$start[1])))
    ok <- ok + 1
  } else fail <- fail + 1
  if (i %% 200 == 0) cat(sprintf("  %d / %d (ok %d, fail %d)\n", i, nrow(tg), ok, fail))
  Sys.sleep(0.07)
}
cat(sprintf("回映完成: ok %d / fail %d\n", ok, fail))
fwrite(lift, "results/phase3_lift_map.csv")
lift[, chrom_hg38 := as.character(chrom_hg38)]
}

## 2. 与 meQTL 联表 (GRCh37 空间) ----------------------------------------------------
me <- fread("results/GSE282512_phase2_dmr_meqtl.csv")
me[, snp_chr37 := sub('"', "", snp_chr)]
me[, snp_chr37 := sub("chr", "", snp_chr37)]
lift[, chrom_hg38 := as.character(chrom_hg38)]
shared <- merge(me, lift, by.x = c("snp_chr37", "snp_pos"),
                by.y = c("chrom_hg19", "pos_hg19"), allow.cartesian = TRUE)
cat(sprintf("meQTL-GWAS 共享变异: %d 对 (变异 %d, DMR %d)\n",
            nrow(shared), uniqueN(shared$pos_hg38), uniqueN(shared$region_id)))

# 加回 GWAS 统计 (PE 与 GH 分别)
pe <- fread("data/gwas/finngen_pe_p1e4.tsv",
            col.names = c("chrom", "pos", "rsids", "nearest_genes", "pval", "beta", "af_alt"))
gh <- fread("data/gwas/finngen_gh_p1e4.tsv",
            col.names = c("chrom", "pos", "rsids", "nearest_genes", "pval", "beta", "af_alt"))
pe[, `:=`(chrom = as.character(chrom), pheno = "PE")]
gh[, `:=`(chrom = as.character(chrom), pheno = "GH")]
gw <- rbind(pe, gh)
setkey(gw, chrom, pos)
shared <- merge(shared, gw[, .(chrom, pos, rsids, pheno, pval, beta)],
                by.x = c("chrom_hg38", "pos_hg38"), by.y = c("chrom", "pos"),
                allow.cartesian = TRUE)
fwrite(shared, "results/phase3_shared_mqtl_gwas.csv")

## 3. 摘要 ---------------------------------------------------------------------------
sink("results/phase3_shared_summary.txt", split = TRUE)
cat("===== Phase 3 meQTL-GWAS 三角证据 (坐标回映版) =====\n\n")
cat(sprintf("回映变异: %d/%d; meQTL-GWAS 共享 SNP-CpG 对: %d (DMR %d, 变异 %d)\n\n",
            nrow(lift), nrow(tg), nrow(shared), uniqueN(shared$region_id),
            uniqueN(paste(shared$chrom_hg38, shared$pos_hg38))))
if (nrow(shared) > 0) {
  cat("-- 共享变异 Top 30 (按 GWAS p) --\n")
  ss <- shared[, .(pheno, pval_gwas = pval, beta_gwas = beta,
                   region_id, dmr_symbol, final_call, cpg,
                   snp = rsids, pos_hg38, pval_mqtl, beta_mqtl)]
  ss <- ss[, .SD[which.min(pval_gwas)], by = .(pheno, region_id, cpg, snp)]
  print(ss[order(pval_gwas)][1:30])
  cat("\n-- 按 DMR 汇总 --\n")
  print(shared[, .(n_shared_pairs = .N, n_shared_var = uniqueN(paste(chrom_hg38, pos_hg38)),
                   min_gwas_p = min(pval)), by = .(region_id, dmr_symbol, final_call)][order(min_gwas_p)])
}
sink()

cat("DONE\n")
