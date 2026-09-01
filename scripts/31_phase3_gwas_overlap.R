# 31_phase3_gwas_overlap.R ------------------------------------------------------
# Phase 3 GWAS 整合: FinnGen R10 PE GWAS x DMR 定位交叉 + 与 Phase 2 meQTL 三角整合
# 输入:
#   data/gwas/finngen_pe_p1e4.tsv          PE GWAS p<1e-4 变异
#   data/gwas/finngen_gh_p1e4.tsv          妊娠期高血压 GWAS p<1e-4 变异
#   results/GSE282512_dmr_final_v2.csv     166 候选 DMR
#   results/GSE282512_phase2_dmr_meqtl.csv DMR 内 meQTL (Phase 2 输出)
#   results/GSE282512_region_annot.csv     背景区域坐标
#   results/GSE282512_dmr_candidates.csv   背景区域 (测过)
# 分析:
#   1) PE GWAS lead loci (p<5e-8, 500kb 合并)
#   2) 候选 DMR ± 1Mb 内 lead / 亚阈值信号
#   3) 富集检验 (候选 vs 背景靠近 lead 的比例)
#   4) meQTL SNP x GWAS p 三角联表
# 输出:
#   results/GSE282512_phase3_gwas_overlap.csv
#   results/GSE282512_phase3_summary.txt
#   figures/phase3_integration.png
# 约定: UTF-8, 相对路径, LANG=en_US.UTF-8
# ------------------------------------------------------------------------------
suppressWarnings(suppressMessages({ library(data.table) }))

LEAD_P    <- 5e-8
SUB_P     <- 1e-4
FLANK     <- 1e6
MERGE_GAP <- 5e5

## 1. GWAS 数据 ------------------------------------------------------------------
pe <- fread("data/gwas/finngen_pe_p1e4.tsv",
            col.names = c("chrom", "pos", "rsids", "nearest_genes", "pval", "beta", "af_alt"))
gh <- fread("data/gwas/finngen_gh_p1e4.tsv",
            col.names = c("chrom", "pos", "rsids", "nearest_genes", "pval", "beta", "af_alt"))
cat(sprintf("PE GWAS p<1e-4: %d; GH GWAS p<1e-4: %d\n", nrow(pe), nrow(gh)))

pe[, chrom := as.character(chrom)]
gh[, chrom := as.character(chrom)]

# lead loci: p<5e-8, 500kb gap 合并
make_loci <- function(dt) {
  sig <- dt[pval < LEAD_P][order(chrom, pos)]
  if (nrow(sig) == 0) return(data.table())
  sig[, gap := c(Inf, as.numeric(diff(pos)))]
  sig[, newloc := gap > MERGE_GAP | chrom != shift(chrom)]
  sig[1, newloc := TRUE]
  sig[, loc := cumsum(newloc)]
  sig[, .(lead_pos = pos[which.min(pval)], lead_rsids = rsids[which.min(pval)],
          lead_p = min(pval), lead_gene = nearest_genes[which.min(pval)],
          n_sig = .N), by = .(loc, chrom)]
}
loc_pe <- make_loci(pe); loc_gh <- make_loci(gh)
cat(sprintf("PE lead loci: %d; GH lead loci: %d\n", nrow(loc_pe), nrow(loc_gh)))

## 2. DMR 与背景区域 ---------------------------------------------------------------
v2 <- fread("results/GSE282512_dmr_final_v2.csv",
            select = c("region_id", "chr", "start", "end", "symbol", "type",
                       "direction", "delta_beta", "final_call"))
v2[, chr := sub("chr", "", chr)]
v2[, mid := (start + end) %/% 2]
reg <- fread("results/GSE282512_region_annot.csv",
             select = c("region_id", "chr", "start", "end"))
tested <- fread("results/GSE282512_dmr_candidates.csv", select = c("region_id"))
reg <- reg[region_id %in% tested$region_id & !region_id %in% v2$region_id]
reg[, chr := sub("chr", "", chr)]
reg[, mid := (start + end) %/% 2]

## 3. DMR ± 1Mb 内 GWAS 信号 --------------------------------------------------------
near_join <- function(dt, loci) {
  if (nrow(loci) == 0) return(data.table())
  m <- merge(dt[, .(region_id, chr, mid)],
             loci[, .(chrom, lead_pos, lead_rsids, lead_p, lead_gene)],
             by.x = "chr", by.y = "chrom", allow.cartesian = TRUE)
  m[abs(lead_pos - mid) <= FLANK]
}
v2_lead <- near_join(v2, loc_pe)
bg_lead <- near_join(reg, loc_pe)

# 亚阈值: DMR ± 1Mb 内最小 GWAS p (用 p<1e-4 抽取集)
sub_join <- function(dt) {
  m <- merge(dt[, .(region_id, chr, mid)],
             pe[, .(chrom, pos, pval)], by.x = "chr", by.y = "chrom", allow.cartesian = TRUE)
  m[, min_p_sub := min(pval), by = region_id]
  m[, .(region_id, min_p_sub)][!duplicated(region_id)]
}
v2_sub <- sub_join(v2)

## 4. 富集检验 -----------------------------------------------------------------------
cand_near <- uniqueN(v2_lead$region_id)
bg_near   <- uniqueN(bg_lead$region_id)
tab <- matrix(c(cand_near, nrow(v2) - cand_near,
                bg_near, nrow(reg) - bg_near), nrow = 2, byrow = TRUE,
              dimnames = list(group = c("candidate", "background"),
                              near = c("near_lead", "not_near")))
fisher <- fisher.test(tab, alternative = "greater")

## 5. meQTL x GWAS 三角联表 -----------------------------------------------------------
me <- fread("results/GSE282512_phase2_dmr_meqtl.csv")
me[, snp_chr2 := sub("chr", "", snp_chr)]
pe_pos <- pe[, .(min_p_gwas = min(pval), beta_gwas = beta[which.min(pval)]),
             by = .(chrom, pos)]
me_g <- merge(me, pe_pos, by.x = c("snp_chr2", "snp_pos"), by.y = c("chrom", "pos"),
              all.x = TRUE, allow.cartesian = TRUE)
me_g[!is.na(min_p_gwas) & min_p_gwas < SUB_P, has_gwas := TRUE]
me_gw <- me_g[has_gwas == TRUE]
# 每 DMR: 有 GWAS 证据的 meQTL SNP 数与最强信号
me_dmr <- me_g[, .(n_mqtl = .N), by = .(region_id, final_call)]
if (nrow(me_gw) > 0) {
  me_gw_dmr <- me_gw[, .(n_mqtl_gwas = .N, best_gwas_p = min(min_p_gwas)), by = region_id]
} else me_gw_dmr <- data.table(region_id = integer(), n_mqtl_gwas = integer(),
                               best_gwas_p = numeric())

## 6. 汇总输出 -------------------------------------------------------------------------
res <- merge(v2, me_dmr[, .(region_id, n_mqtl)], by = "region_id", all.x = TRUE)
res <- merge(res, me_gw_dmr, by = "region_id", all.x = TRUE)
res <- merge(res, v2_sub, by = "region_id", all.x = TRUE)
res[is.na(n_mqtl), n_mqtl := 0]
res[is.na(n_mqtl_gwas), n_mqtl_gwas := 0]
res[, has_lead := region_id %in% v2_lead$region_id]
res[, evidence := fcase(
  has_lead & n_mqtl_gwas > 0, "gwas_lead_and_mqtl",
  has_lead, "gwas_lead_only",
  n_mqtl_gwas > 0, "mqtl_gwas_subthreshold",
  n_mqtl > 0, "mqtl_only",
  min_p_sub < SUB_P, "gwas_subthreshold_only",
  default = "none")]
setorder(res, evidence, min_p_sub, na.last = TRUE)
fwrite(res, "results/GSE282512_phase3_gwas_overlap.csv")

sink("results/GSE282512_phase3_summary.txt", split = TRUE)
cat("===== Phase 3 GWAS 整合 (FinnGen R10 x DMR) =====\n\n")
cat(sprintf("FinnGen PE (O15_PRE_OR_ECLAMPSIA) p<1e-4 变异: %d\n", nrow(pe)))
cat(sprintf("PE lead loci (p<5e-8, 500kb 合并): %d\n", nrow(loc_pe)))
if (nrow(loc_pe) > 0) print(loc_pe[order(lead_p)])
cat(sprintf("\n妊娠期高血压 (O15_GESTAT_HYPERT) lead loci: %d\n", nrow(loc_gh)))
if (nrow(loc_gh) > 0) print(loc_gh[order(lead_p)])

cat("\n-- DMR x GWAS lead loci (±1Mb) --\n")
cat(sprintf("候选 DMR 靠近 PE lead: %d / %d\n", cand_near, nrow(v2)))
if (nrow(v2_lead) > 0) {
  print(merge(v2_lead[, .(region_id, symbol, final_call, direction,
                           lead_rsids, lead_p, lead_gene,
                           dist_kb = round(abs(lead_pos - mid) / 1000))],
              v2[, .(region_id, delta_beta)], by = "region_id")[order(lead_p)])
}
cat("\n-- 富集检验 (候选 vs 背景, ±1Mb 内有 lead) --\n")
print(tab)
cat(sprintf("OR = %.2f, p = %.3g\n\n",
            (cand_near * (nrow(reg) - bg_near)) / ((nrow(v2) - cand_near) * bg_near),
            fisher$p.value))

cat("-- meQTL x GWAS 三角联表 --\n")
cat(sprintf("DMR 内 meQTL: %d 个 SNP-CpG 对; 其中 SNP 在 FinnGen PE p<1e-4: %d 对 (涉及 DMR %d)\n",
            nrow(me), nrow(me_gw), uniqueN(me_gw$region_id)))
if (nrow(me_gw) > 0) {
  cat("\nTop 15 meQTL-GWAS 双证据 SNP (按 GWAS p):\n")
  print(me_gw[order(min_p_gwas)][1:15,
    .(region_id, dmr_symbol, final_call, cpg, snp = paste0(snp_chr, ":", snp_pos),
      pval_mqtl, beta_mqtl, min_p_gwas, beta_gwas)])
}
cat("\n-- 综合证据分级 (166 候选) --\n")
print(res[, .N, by = evidence])
sink()

## 7. 图: 综合证据分布 ---------------------------------------------------------------
png("figures/phase3_integration.png", width = 1500, height = 900, res = 150)
par(mar = c(8, 5, 3, 1))
tb <- res[, .N, by = evidence][order(N)]
bp <- barplot(tb$N, names.arg = gsub("_", "\n", tb$evidence), col = "steelblue",
              ylab = "候选 DMR 数", main = "Phase 2+3 综合证据分布 (166 候选 DMR)",
              las = 2, cex.names = 0.8)
text(bp, tb$N + 1, tb$N, cex = 0.8)
dev.off()

cat("DONE\n")
