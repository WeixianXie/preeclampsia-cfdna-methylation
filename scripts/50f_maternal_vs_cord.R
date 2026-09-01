#!/usr/bin/env Rscript
# 50f_maternal_vs_cord.R — delivery vs cordB 配对差异 (数据信号判别)
suppressMessages(library(data.table))
mc <- fread('results/GSE154378_bin_mc_long.csv.gz')
sm <- fread('data/geo_methylation/GSE154378/gse154378_samples.tsv')
mc[, beta := m / c]
m <- merge(mc, sm[, .(gsm, group, patient, timepoint)], by = 'gsm')
sub <- m[group == 'Normal' & timepoint %in% c('delivery', 'cordB')]
dc <- dcast(sub, patient + bin ~ timepoint, value.var = 'beta', fun.aggregate = mean)
dc <- dc[!is.na(delivery) & !is.na(cordB)]
dc[, delta := delivery - cordB]
dc[, n_pat := .N, by = bin]
stats <- dc[n_pat >= 3, .(n_pat = .N, mean_delta = mean(delta), med_delta = median(delta),
                          t_stat = suppressWarnings(t.test(delta)$statistic),
                          p = suppressWarnings(t.test(delta)$p.value)), by = bin]
stats[, padj := p.adjust(p, 'BH')]
cat('bins with >=3 paired patients:', nrow(stats), '\n')
cat('p<0.05:', sum(stats$p < 0.05, na.rm = TRUE), ' p<0.01:', sum(stats$p < 0.01, na.rm = TRUE),
    ' BH<0.05:', sum(stats$padj < 0.05, na.rm = TRUE), '\n')
cat('mean_delta>0:', sum(stats$mean_delta > 0), ' <0:', sum(stats$mean_delta < 0), '\n')
setorder(stats, -abs(mean_delta))
cat('\nTop 12 差异 bin:\n')
print(head(stats[, .(bin, n_pat, mean_delta, p, padj)], 12))

mk <- fread('data/geo_methylation/GSE154378/sun2015_markers.tsv')
sg <- merge(stats, mk[, .(bin, type)], by = 'bin')
cat('\n显著 bin (p<0.05) type 构成:\n'); print(sg[p < 0.05, .(n = .N), by = type])
cat('全部 bin type 构成:\n'); print(sg[, .(n = .N), by = type])
# 胎盘方向富集
refs <- c('Liver','Lungs','Colon','Small intestines','Pancreas','Adrenal glands','Esophagus',
          'Adipose tissues','Heart','Brain','T-cells','B-cells','Neutrophils','Placenta')
mk2 <- mk[type == 'II']
mk2[, med_other := apply(as.matrix(mk2[, ..refs]), 1, function(x) median(x[-14], na.rm = TRUE))]
mk2[, pdir := ifelse(Placenta - med_other > 15, 'HIGH', ifelse(med_other - Placenta > 15, 'LOW', 'FLAT'))]
sg2 <- merge(stats, mk2[, .(bin, pdir)], by = 'bin')
cat('\n显著 bin 胎盘方向构成 (p<0.05):\n'); print(sg2[p < 0.05, .(n = .N), by = pdir])
fwrite(stats, 'results/GSE154378_maternal_vs_cord_bins.csv')
