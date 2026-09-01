#!/usr/bin/env Rscript
# 50b_validate_mapping.R — 决定性映射验证: cfDNA 观测 beta vs 参考组织甲基化
suppressMessages(library(data.table))
mk <- fread('data/geo_methylation/GSE154378/sun2015_markers.tsv')
mc <- fread('results/GSE154378_bin_mc_long.csv.gz')
mc[, beta := m / c]
obs <- mc[, .(obs_beta = mean(beta), n_samp = .N), by = bin]
obs <- obs[n_samp >= 30]
d <- merge(obs, mk, by = 'bin')
refs <- c('Neutrophils', 'T-cells', 'B-cells', 'Placenta', 'Liver', 'Brain', 'Colon')
cat('== bins (coverage>=30 samples):', nrow(d), '==\n')
for (tc in refs) {
  ct <- suppressWarnings(cor.test(d[[tc]], d$obs_beta, method = 'spearman'))
  cat(sprintf('%-14s rho=%6.3f  p=%.2e\n', tc, ct$estimate, ct$p.value))
}
for (ty in c('I', 'II')) {
  dd <- d[type == ty]
  if (nrow(dd) > 5) {
    ct <- suppressWarnings(cor.test(dd[['Neutrophils']], dd$obs_beta, method = 'spearman'))
    cat(sprintf('type %s (n=%d): Neutrophils ref rho=%6.3f  p=%.2e\n', ty, nrow(dd), ct$estimate, ct$p.value))
  }
}
# 参考分布描述
cat('\n== 参考甲基化分布 ==\n')
for (tc in c('Neutrophils', 'T-cells', 'B-cells', 'Placenta')) {
  cat(sprintf('%-14s median=%.3f  q25=%.3f  q75=%.3f\n', tc,
              median(mk[[tc]]), quantile(mk[[tc]], .25), quantile(mk[[tc]], .75)))
}
cat('\nType.I bins:', nrow(mk[type=='I']), ' Type.II bins:', nrow(mk[type=='II']), '\n')
