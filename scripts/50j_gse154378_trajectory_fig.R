#!/usr/bin/env Rscript
# 50j_gse154378_trajectory_fig.R
# GSE154378 DMR 孕期轨迹验证图 (v4 坐标排序标记表)
# 依赖 50_gse154378_trajectory.R 的产物; 独立重算 DMR-bin 轨迹数据
# 运行: LANG=en_US.UTF-8 Rscript scripts/50j_gse154378_trajectory_fig.R
suppressMessages({library(data.table); library(ggplot2)})

mk <- fread('data/geo_methylation/GSE154378/sun2015_markers.tsv')
sm <- fread('data/geo_methylation/GSE154378/gse154378_samples.tsv')
mc <- fread('results/GSE154378_bin_mc_long.csv.gz')
mc[, beta := m / c]

tp_order <- c('1stT'=1, '2ndT'=2, '3rdT'=3, 'delivery'=4)
sm[, tp_num := ifelse(timepoint %in% names(tp_order), tp_order[timepoint], NA)]
sm[group=='NP', tp_num := 0]
tp_lab <- c('0'='NP', '1'='1stT', '2'='2ndT', '3'='3rdT', '4'='delivery')

# ---- DMR-bin 重叠 (v2) ----
ov <- fread('results/GSE154378_dmr_bin_overlap.csv')
ov2 <- ov[, .(bin=unique(bin)), by=region_id]
ov2 <- merge(ov2, unique(fread('results/GSE282512_dmr_final_v2.csv')[, .(region_id, direction, symbol)]), by='region_id')

# Normal 组每 bin 每时间点加权 beta (排除 cordB, 因 tp_num 未定义)
mc_m <- merge(mc, sm[group %in% c('Normal','NP') & !is.na(tp_num), .(gsm, group, tp_num)], by='gsm')
agg_dmr <- mc_m[bin %in% unique(ov$bin), .(m=sum(m), c=sum(c)), by=.(tp_num, bin)]
agg_dmr[, beta := m/c]
agg_dmr <- merge(agg_dmr, unique(ov[, .(bin, direction)]), by='bin', all.x=TRUE)  # direction 来源
agg_dmr[, tp_lab := factor(tp_lab[as.character(tp_num)], levels=tp_lab)]
setorder(agg_dmr, bin, tp_num)

# ---- Type.I 组织标记轨迹 (Neutrophils vs 其他) ----
mc_m2 <- mc_m[group=='Normal']
agg_allb <- mc_m2[, .(m=sum(m), c=sum(c)), by=.(tp_num, bin)]
agg_allb[, beta := m/c]
rho_all <- agg_allb[order(bin, tp_num), {
  if (.N >= 3 && length(unique(beta)) >= 2) list(rho=suppressWarnings(cor(tp_num, beta, method='spearman')), c_min=min(c))
  else list(rho=NA_real_, c_min=min(c))
}, by=bin]
rho_all <- rho_all[!is.na(rho) & c_min >= 10]
mkI <- mk[type=='I']
mkI_rho <- merge(mkI[, .(bin, tissue)], rho_all, by='bin')
# 组织聚合轨迹: Neutrophils vs 全部 Type.I
neu_bins <- mkI_rho[tissue=='Neutrophils', bin]
agg_neu <- agg_allb[bin %in% neu_bins, .(m=sum(m), c=sum(c)), by=tp_num][, beta:=m/c][, grp:='Neutrophil markers']
agg_other <- agg_allb[!bin %in% neu_bins & bin %in% mkI_rho$bin, .(m=sum(m), c=sum(c)), by=tp_num][, beta:=m/c][, grp:='Other Type.I markers']
agg_tis <- rbind(agg_neu, agg_other)
cat('agg_tis rows:', nrow(agg_tis), '\n'); print(head(agg_tis))
agg_tis[, tp_lab := factor(tp_lab[as.character(tp_num)], levels=tp_lab)]

# ---- Figure 1: DMR-bin 轨迹 ----
# DMR bin 注释标签 (ov 已含 symbol; direction 已在 agg_dmr 中)
bin_lab <- ov[, .(sym = paste0(unique(symbol[symbol!='']), collapse='/')), by=bin]
bin_lab[sym=='', sym := 'no gene']
agg_dmr <- merge(agg_dmr, bin_lab, by='bin', all.x=TRUE)

p1 <- ggplot(agg_dmr[!is.na(tp_lab)], aes(x=tp_lab, y=beta, group=bin, color=direction)) +
  geom_line(linewidth=0.9) + geom_point(size=2.2) +
  facet_wrap(~bin, ncol=5, scales='free_y') +
  scale_color_manual(values=c('hyper'='#C0392B', 'hypo'='#27AE60'), name='DMR direction') +
  labs(x=NULL, y='Group-mean methylated fraction (cfDNA)',
       title='GSE154378 Normal pregnancy: DMR-overlapping marker bins across gestation') +
  theme_bw(base_size=10) +
  theme(axis.text.x=element_text(angle=45, hjust=1), strip.background=element_rect(fill='grey90'))
png('figures/GSE154378_dmr_trajectory.png', width=10, height=3.4, units='in', res=300)
print(p1); dev.off()

# ---- Figure 2: Neutrophil vs other Type.I ----
p2 <- ggplot(agg_tis, aes(x=tp_lab, y=beta, group=grp, color=grp)) +
  geom_line(linewidth=1) + geom_point(size=2.5) +
  scale_color_manual(values=c('Neutrophil markers'='#C0392B', 'Other Type.I markers'='grey40')) +
  labs(x=NULL, y='Group-mean methylated fraction (cfDNA)',
       title='GSE154378 Normal pregnancy: Type.I tissue markers across gestation (v4)') +
  theme_bw(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom')
png('figures/GSE154378_typeI_trajectory.png', width=6.5, height=4.2, units='in', res=300)
print(p2); dev.off()

# ---- 汇总表: 每 bin 轨迹 ----
fwrite(agg_dmr, 'results/GSE154378_dmr_bin_trajectory_data.csv')
fwrite(rho_all, 'results/GSE154378_allbin_rho.csv')
cat('figures saved: figures/GSE154378_dmr_trajectory.png, figures/GSE154378_typeI_trajectory.png\n')
cat('data saved: results/GSE154378_dmr_bin_trajectory_data.csv\n')
print(agg_dmr[!is.na(tp_lab), .(bin, direction, sym, tp_lab, beta)])
