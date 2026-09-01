# 30_phase2_mqtl_overlap.R ------------------------------------------------------
# Phase 2 甲基化-QTL 整合: GoDMC cis-meQTL x DMR 交叉 (内存受控版)
# 数据流: data/mqtl/chunks/chr* (按 SNP 染色体分块, 每块独立处理)
# 输入:
#   data/mqtl/chunks/chrN               cis-meQTL (cpg,snp_chr,snp_pos,pval,beta_a1,se)
#   data/annot/HM450.hg38.manifest.tsv.gz 450K cg -> hg38 坐标
#   results/GSE282512_dmr_final_v2.csv  166 候选 DMR (v2 判定)
#   results/GSE282512_region_annot.csv   全部区域坐标 (背景)
#   results/GSE282512_dmr_candidates.csv 原分析测过的全部区域 (背景)
# 输出:
#   results/GSE282512_phase2_dmr_meqtl.csv  候选 DMR 内全部 meQTL 行
#   results/GSE282512_phase2_summary.txt
# 约定: UTF-8, 相对路径, LANG=en_US.UTF-8
# ------------------------------------------------------------------------------
suppressWarnings(suppressMessages({ library(data.table) }))

COLS <- c("cpg", "snp_chr", "snp_pos", "pval", "beta_a1", "se")

## 注释与小表 ---------------------------------------------------------------------
man <- fread("data/annot/HM450.hg38.manifest.tsv.gz",
             select = c("Probe_ID", "CpG_chrm", "CpG_beg"))
setnames(man, c("cpg", "cchr", "cbeg"))
man <- man[!is.na(cchr) & grepl("^chr[0-9XY]{1,2}$", cchr)]
man[, cchr := sub("chr", "", cchr)]

reg <- fread("results/GSE282512_region_annot.csv",
             select = c("region_id", "chr", "start", "end"))
reg[, chr := sub("chr", "", chr)]
tested <- fread("results/GSE282512_dmr_candidates.csv", select = c("region_id"))
reg <- reg[region_id %in% tested$region_id]
v2 <- fread("results/GSE282512_dmr_final_v2.csv",
            select = c("region_id", "chr", "start", "end", "symbol", "type",
                       "direction", "delta_beta", "final_call"))
v2[, chr := sub("chr", "", chr)]
reg <- reg[!region_id %in% v2$region_id]
cat(sprintf("候选 DMR %d; 背景区域 %d\n", nrow(v2), nrow(reg)))

## 逐染色体块处理 -------------------------------------------------------------------
chunks <- list.files("data/mqtl/chunks", full.names = TRUE)
tot_rows <- 0
bg_counts <- data.table(region_id = integer(), N = integer())
v2_rows <- list()

for (cf in chunks) {
  ch <- sub(".*[/\\\\]", "", cf)          # chr1 .. chrX
  chrn <- sub("chr", "", ch)
  mqc <- tryCatch(fread(cf, header = FALSE, quote = "", col.names = COLS),
                  error = function(e) NULL)
  if (is.null(mqc) || nrow(mqc) == 0) next
  tot_rows <- tot_rows + nrow(mqc)
  # 只留 CpG 在该染色体的行 (cis 时 SNP 与 CpG 同链, 减少映射表 join 量)
  mc <- merge(mqc, man[cchr == chrn, .(cpg, cbeg)], by = "cpg")
  if (nrow(mc) == 0) { rm(mqc); next }
  # 背景: CpG 落入区域
  rc <- reg[chr == chrn]
  if (nrow(rc) > 0) {
    h <- rc[mc, on = .(start <= cbeg, end >= cbeg), nomatch = NULL,
            .(region_id)]
    if (nrow(h) > 0) bg_counts <- rbindlist(list(bg_counts, h[, .N, by = region_id]))
  }
  # 候选 DMR
  vc <- v2[chr == chrn]
  if (nrow(vc) > 0) {
    hv <- vc[mc, on = .(start <= cbeg, end >= cbeg), nomatch = NULL,
             .(region_id, cpg, snp_chr, snp_pos, pval, beta_a1, se)]
    if (nrow(hv) > 0) v2_rows[[length(v2_rows) + 1]] <- hv
  }
  n_row <- nrow(mqc)
  rm(mqc, mc); gc()
  cat(sprintf("  %s: %d 行\n", ch, n_row))
}
cat(sprintf("GoDMC cis-meQTL 总行数: %d\n", tot_rows))

bg_counts <- bg_counts[, .(N = sum(N)), by = region_id]
mq_v2 <- rbindlist(v2_rows)
mq_v2 <- merge(mq_v2, v2[, .(region_id, dmr_symbol = symbol, final_call)], by = "region_id")
setorder(mq_v2, region_id, pval)
cat(sprintf("meQTL CpG 落入候选 DMR: %d 对 (DMR %d); 背景命中区域 %d\n",
            nrow(mq_v2), uniqueN(mq_v2$region_id), nrow(bg_counts)))

fwrite(mq_v2[, .(region_id, dmr_symbol, final_call, cpg,
                 snp_chr = sub('"', "", snp_chr), snp_pos,
                 pval_mqtl = pval, beta_mqtl = beta_a1, se_mqtl = se)],
       "results/GSE282512_phase2_dmr_meqtl.csv")

## 富集检验 -------------------------------------------------------------------------
cand_hit <- uniqueN(mq_v2$region_id)
bg_hit   <- nrow(bg_counts)
tab <- matrix(c(cand_hit, nrow(v2) - cand_hit,
                bg_hit, nrow(reg) - bg_hit), nrow = 2, byrow = TRUE,
              dimnames = list(group = c("candidate", "background"),
                              hit = c("with_meQTL", "without")))
fisher <- fisher.test(tab, alternative = "greater")

# meQTL 密度 (每 kb)
idx_c <- match(bg_counts$region_id, reg$region_id)
dens_bg <- rep(0, nrow(reg))
dens_bg[idx_c] <- bg_counts$N * 1000 / (reg$end[idx_c] - reg$start[idx_c] + 1)
d_c <- mq_v2[, .N, by = region_id]
idx_v <- match(d_c$region_id, v2$region_id)
dens_cand <- rep(0, nrow(v2))
dens_cand[idx_v] <- d_c$N * 1000 / (v2$end[idx_v] - v2$start[idx_v] + 1)
dens_w <- wilcox.test(dens_cand, dens_bg)

## 汇总 -----------------------------------------------------------------------------
sink("results/GSE282512_phase2_summary.txt", split = TRUE)
cat("===== Phase 2 甲基化-QTL 整合 (GoDMC x DMR) =====\n\n")
cat(sprintf("GoDMC cis-meQTL (p<1e-8): %d 对 (22 染色体分块流式处理)\n\n", tot_rows))
cat(sprintf("候选 DMR (%d) 中含 meQTL CpG 的区域: %d (%.1f%%)\n",
            nrow(v2), cand_hit, 100 * cand_hit / nrow(v2)))
cat(sprintf("背景测过区域 (%d) 中含 meQTL CpG: %d (%.1f%%)\n\n",
            nrow(reg), bg_hit, 100 * bg_hit / nrow(reg)))
cat("-- Fisher 精确检验 (候选 vs 背景, 单侧) --\n")
print(tab)
cat(sprintf("OR = %.2f, p = %.3g\n\n",
            (cand_hit * (nrow(reg) - bg_hit)) / ((nrow(v2) - cand_hit) * bg_hit),
            fisher$p.value))
cat(sprintf("meQTL 密度 (个/kb): 候选中位 %.3f vs 背景 %.3f, Wilcoxon p = %.3g\n\n",
            median(dens_cand), median(dens_bg), dens_w$p.value))
cat("-- 含 meQTL 的候选 DMR (按 final_call) --\n")
print(v2[region_id %in% unique(mq_v2$region_id), .N, by = final_call])
cat("\n-- 每个含 meQTL 候选 DMR 的 Top meQTL --\n")
top <- mq_v2[, .SD[which.min(pval)], by = region_id]
top2 <- v2[, .(region_id, chr, start, end, symbol, direction, delta_beta, final_call)]
print(merge(top2, top[, .(region_id, cpg, snp_chr, snp_pos, pval_mqtl = pval,
                          beta_mqtl = beta_a1)], by = "region_id", sort = FALSE))
sink()

cat("DONE\n")
