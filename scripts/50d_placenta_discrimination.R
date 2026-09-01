#!/usr/bin/env Rscript
# 50d_placenta_discrimination.R — 决定性映射判别: 胎盘标记 母体cfDNA vs 脐血
# 若 bin->标记 映射正确:
#   placenta HIGH 标记 (胎盘高甲基化): 母体 cfDNA(含胎盘DNA) beta > 脐血(无胎盘) beta
#   placenta LOW  标记 (胎盘低甲基化): 母体 cfDNA beta < 脐血 beta
# 同一患者 delivery(母体) vs cordB(脐血) 配对检验
suppressMessages(library(data.table))
mk <- fread('data/geo_methylation/GSE154378/sun2015_markers.tsv')
sm <- fread('data/geo_methylation/GSE154378/gse154378_samples.tsv')
mc <- fread('results/GSE154378_bin_mc_long.csv.gz')
mc[, beta := m / c]

ref_cols <- c('Liver','Lungs','Colon','Small intestines','Pancreas','Adrenal glands',
              'Esophagus','Adipose tissues','Heart','Brain','T-cells','B-cells','Neutrophils','Placenta')
# 胎盘信息方向 (Type.II 内: 胎盘 vs 其他组织中位数)
mk2 <- mk[type == 'II']
mk2[, med_other := apply(as.matrix(mk2[, ..ref_cols]), 1, function(x) median(x[-14], na.rm=TRUE))]
mk2[, placenta_dir := ifelse(Placenta - med_other > 15, 'HIGH',
                      ifelse(med_other - Placenta > 15, 'LOW', 'FLAT'))]
cat('Type.II placenta_dir: HIGH', sum(mk2$placenta_dir=='HIGH'), ' LOW', sum(mk2$placenta_dir=='LOW'),
    ' FLAT', sum(mk2$placenta_dir=='FLAT'), '\n')

# 配对患者: Normal 且有 delivery + cordB
pats <- sm[group=='Normal' & timepoint %in% c('delivery','cordB'), .(timepoint, gsm), by=patient]
pats <- pats[patient %in% names(which(table(pats$patient) == 2))]
mat <- dcast(pats, patient ~ timepoint, value.var='gsm')
cat('配对患者数:', nrow(mat), '\n')

# 每标记: 组均值 beta (母体 delivery vs 脐血 cordB), 患者配对 delta
res <- lapply(c('HIGH','LOW'), function(dr) {
  bins <- mk2[placenta_dir == dr]$bin
  sub <- mc[bin %in% bins & gsm %in% unlist(mat[, .(delivery, cordB)])]
  sub <- merge(sub, melt(mat, id.vars='patient', variable.name='timepoint', value.name='gsm')[, .(gsm, patient, timepoint)], by='gsm')
  sub <- sub[!is.na(beta)]
  # 患者-标记 级 delta
  dcast_pt <- dcast(sub, patient + bin ~ timepoint, value.var='beta', fun.aggregate=mean)
  dcast_pt <- dcast_pt[!is.na(delivery) & !is.na(cordB)]
  dcast_pt[, delta := delivery - cordB]
  list(dir=dr, n_marker=length(bins), n_pair_pt=nrow(dcast_pt),
       mean_delta=mean(dcast_pt$delta), median_delta=median(dcast_pt$delta),
       wilcox_p=suppressWarnings(wilcox.test(dcast_pt$delta)$p.value),
       t_p=suppressWarnings(t.test(dcast_pt$delta)$p.value),
       frac_pos=mean(dcast_pt$delta > 0))
})
res <- rbindlist(lapply(res, as.data.table))
print(res)

# 组均值对比 (全部标记聚合, 更稳)
cat('\n== 组均值 (全部 HIGH/LOW 标记聚合) ==\n')
for (dr in c('HIGH','LOW')) {
  bins <- mk2[placenta_dir == dr]$bin
  agg <- mc[bin %in% bins & gsm %in% unlist(mat[, .(delivery, cordB)])]
  agg <- merge(agg, melt(mat, id.vars='patient', variable.name='timepoint', value.name='gsm')[, .(gsm, timepoint)], by='gsm')
  g <- agg[, .(m=sum(m), c=sum(c)), by=.(timepoint, bin)][, beta := m/c]
  g2 <- g[, .(m=sum(m), c=sum(c)), by=timepoint][, beta := m/c]
  cat(sprintf('%s: delivery beta=%.4f  cordB beta=%.4f  diff=%.4f\n', dr,
              g2[timepoint=='delivery', beta], g2[timepoint=='cordB', beta],
              g2[timepoint=='delivery', beta] - g2[timepoint=='cordB', beta]))
}
