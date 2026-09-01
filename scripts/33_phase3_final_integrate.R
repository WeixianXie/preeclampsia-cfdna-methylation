# 33_phase3_final_integrate.R ---------------------------------------------------
# Phase 3 收尾: 整合 lead loci 重叠 + 富集 + 修正版 (坐标回映) meQTL-GWAS 三角证据
# 输入:
#   results/GSE282512_phase3_gwas_overlap.csv   (含 evidence 旧列, 忽略)
#   results/phase3_shared_mqtl_gwas.csv         修正版三角联表
#   results/GSE282512_dmr_final_v2.csv
#   data/gwas/finngen_pe_p1e4.tsv / finngen_gh_p1e4.tsv
# 输出:
#   results/GSE282512_phase3_final.csv
#   results/GSE282512_phase3_final_summary.txt
#   figures/phase3_integration.png (重绘)
# 约定: UTF-8, 相对路径, LANG=en_US.UTF-8
# ------------------------------------------------------------------------------
suppressWarnings(suppressMessages({ library(data.table) }))

v2 <- fread("results/GSE282512_dmr_final_v2.csv")
p3 <- fread("results/GSE282512_phase3_gwas_overlap.csv")
sh <- fread("results/phase3_shared_mqtl_gwas.csv")

pe <- fread("data/gwas/finngen_pe_p1e4.tsv",
            col.names = c("chrom", "pos", "rsids", "nearest_genes", "pval", "beta", "af_alt"))
pe[, chrom := as.character(chrom)]
sig <- pe[pval < 5e-8][order(chrom, pos)]
sig[, gap := c(Inf, as.numeric(diff(pos)))]
sig[, newloc := gap > 5e5 | chrom != shift(chrom)]
sig[1, newloc := TRUE]; sig[, loc := cumsum(newloc)]
loc_pe <- sig[, .(lead_pos = pos[which.min(pval)], lead_rsids = rsids[which.min(pval)],
                  lead_p = min(pval), lead_gene = nearest_genes[which.min(pval)],
                  n_sig = .N), by = .(loc, chrom)]

# 共享变异按 DMR 汇总 (PE 与 GH 分开)
sh_dmr <- sh[, .(n_shared_var = uniqueN(paste(chrom_hg38, pos_hg38)),
                 min_p_pe = min(pval[pheno == "PE"], na.rm = TRUE) *
                   as.integer(any(pheno == "PE")),
                 min_p_gh = min(pval[pheno == "GH"], na.rm = TRUE) *
                   as.integer(any(pheno == "GH"))),
             by = .(region_id)]
sh_dmr[is.infinite(min_p_pe) | min_p_pe == 0, min_p_pe := NA_real_]
sh_dmr[is.infinite(min_p_gh) | min_p_gh == 0, min_p_gh := NA_real_]

res <- merge(v2, p3[, .(region_id, has_lead, min_p_sub)], by = "region_id", all.x = TRUE)
res <- merge(res, sh_dmr, by = "region_id", all.x = TRUE)
res[is.na(n_shared_var), n_shared_var := 0]
me2 <- fread("results/GSE282512_phase2_dmr_meqtl.csv")
cnt <- me2[, .N, by = region_id]
res[, n_mqtl := cnt$N[match(region_id, cnt$region_id)]]
res[is.na(n_mqtl), n_mqtl := 0]
res[, final_evidence := fcase(
  has_lead == TRUE & n_shared_var > 0, "gwas_lead_and_shared_mqtl",
  has_lead == TRUE, "gwas_lead_only",
  n_shared_var > 0, "shared_mqtl_gwas_subthreshold",
  n_mqtl > 0 & !is.na(min_p_sub) & min_p_sub < 1e-4, "mqtl_and_gwas_regional",
  n_mqtl > 0, "mqtl_only",
  !is.na(min_p_sub) & min_p_sub < 1e-4, "gwas_regional_only",
  default = "none")]
setorder(res, -n_shared_var, min_p_sub, na.last = TRUE)
fwrite(res, "results/GSE282512_phase3_final.csv")

sink("results/GSE282512_phase3_final_summary.txt", split = TRUE)
cat("===== Phase 3 GWAS 整合最终汇总 (FinnGen R10 x DMR + meQTL 三角) =====\n\n")
cat("重要方法学说明: GoDMC 坐标为 GRCh37, 与 FinnGen (GRCh38) 联表前已用\n")
cat("Ensembl REST 将窗口内 FinnGen 变异回映至 GRCh37 (1359/1359 成功);\n")
cat("直接精确联表因此前版本错位得到 0 命中, 回映后恢复统计功效。\n\n")
cat("-- FinnGen PE (O15_PRE_OR_ECLAMPSIA) lead loci (p<5e-8) --\n")
print(loc_pe[, .(chrom, lead_pos, lead_rsids, lead_gene, lead_p, n_sig)])
cat("\nGH (O15_GESTAT_HYPERT) lead loci: 14 个 (详见 phase3_summary)\n")
cat("\n-- DMR x GWAS --\n")
cat("候选 DMR 靠近 PE lead (±1Mb): 0/166; 背景区域 586/72982 (0.8%),\n")
cat("Fisher p=1 -> DMR 位点与 FinnGen PE 全显著基因组位点无直接重叠\n")
cat("DMR ±1Mb 内存在 PE 亚阈值 (p<1e-4) 信号: 50/166\n\n")
cat("-- meQTL-GWAS 三角证据 (修正版) --\n")
cat(sprintf("meQTL-GWAS 共享变异: 235 对, 变异 90 个, 涉及 DMR %d 个:\n", nrow(sh_dmr)))
print(merge(sh_dmr, v2[, .(region_id, chr, start, end, symbol, direction,
                            delta_beta, final_call)], by = "region_id"))
cat("\nTop 共享证据 (按 PE p):\n")
top <- sh[pheno == "PE"][order(pval)][1:min(10, sum(sh$pheno == "PE")),
  .(region_id, dmr_symbol, cpg, snp = rsids, p_gwas = pval,
    beta_gwas = beta, pval_mqtl, beta_mqtl)]
print(top)
cat("\n-- 综合证据分级 (166 候选) --\n")
print(res[, .N, by = final_evidence])
cat("\n-- 解读 --\n")
cat("1. DMR 不富集于 FinnGen PE 全基因组显著位点 (0/166), 亦不富集 meQTL CpG\n")
cat("   (Phase 2: 68.7% vs 75.7%) -> cfDNA 甲基化差异独立于经典遗传易感位点\n")
cat("   与 cis-meQTL 遗传调控, 支持'覆盖度构成/片段来源'生物学解释。\n")
cat("2. 3 个 DMR (chr6 SLC17A1 基因体, chr3:49.9Mb 启动子, chr17:39.45Mb 启动子)\n")
cat("   存在 meQTL-GWAS 共享变异 (GH 信号为主, p~1e-5 量级, 含 PE p=4.3e-5 的\n")
cat("   rs116510165), 是甲基化-基因型-疾病三重交汇的候选位点。\n")
sink()

## 图 -------------------------------------------------------------------------------
png("figures/phase3_integration.png", width = 1500, height = 900, res = 150)
par(mar = c(9, 5, 3, 1))
tb <- res[, .N, by = final_evidence][order(N)]
bp <- barplot(tb$N, names.arg = gsub("_", "\n", tb$final_evidence), col = "steelblue",
              ylab = "候选 DMR 数", las = 2, cex.names = 0.75,
              main = "Phase 2+3 综合证据分布 (166 候选 DMR, meQTL-GWAS 三角修正版)")
text(bp, tb$N + 0.5, tb$N, cex = 0.8)
dev.off()
cat("DONE\n")
