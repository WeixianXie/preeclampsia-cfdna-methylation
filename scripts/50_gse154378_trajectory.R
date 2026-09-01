#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# 50_gse154378_trajectory.R
# GSE154378 (UCLA Del Vecchio) DMR 孕期轨迹验证 + bin→标记映射验证
# 依赖 49_gse154378_aggregate.py 与 49b_liftover_dmrs.py 的产物
# 运行: LANG=en_US.UTF-8 Rscript scripts/50_gse154378_trajectory.R

suppressMessages({library(data.table); library(ggplot2)})
# 运行前需 cd 到项目根 (中文绝对路径在 UTF-8 环境下 setwd 失败)
options(datatable.fread.datatable=TRUE)

# ---------- 载入数据 ----------
mk <- fread('data/geo_methylation/GSE154378/sun2015_markers.tsv')      # bin,chr,start,end,type,tissue,placenta_informative + 14 tissues
sm <- fread('data/geo_methylation/GSE154378/gse154378_samples.tsv')    # gsm,group,subtype,patient,timepoint
mc <- fread('results/GSE154378_bin_mc_long.csv.gz')                    # gsm,bin,m,c
mc[, beta := m / c]

# 时间点顺序 (母体; cordB=脐血排除于母体轨迹)
tp_order <- c('1stT'=1, '2ndT'=2, '3rdT'=3, 'delivery'=4)
sm[, tp_num := ifelse(timepoint %in% names(tp_order), tp_order[timepoint], NA)]
sm[group=='NP', tp_num := 0]

# 样本-时间点 汇总 (用于组均值轨迹)
sm_sub <- sm[group %in% c('Normal','NP')]
mc_m <- merge(mc, sm_sub[, .(gsm, group, tp_num)], by='gsm')

# ---------- Part 1: bin→标记映射验证 (Type.II 胎盘标记孕期轨迹) ----------
mkII <- mk[type=='II' & placenta_informative=='Y']
mkII[, placenta_dir := ifelse(Placenta - median(c(Liver,Lungs,Colon,`Small intestines`,Pancreas,`Adrenal glands`,Esophagus,`Adipose tissues`,Heart,Brain,`T-cells`,`B-cells`,Neutrophils)) > 15, 'HIGH',
                       ifelse(median(c(Liver,Lungs,Colon,`Small intestines`,Pancreas,`Adrenal glands`,Esophagus,`Adipose tissues`,Heart,Brain,`T-cells`,`B-cells`,Neutrophils)) - Placenta > 15, 'LOW', 'FLAT')), by=bin]

agg1 <- mc_m[mkII, on='bin'][!is.na(tp_num) & group %in% c('Normal','NP')]
# 组均值: 每 (tp_num, bin, dir) 合并
agg1 <- agg1[, .(m=sum(m), c=sum(c)), by=.(tp_num, bin, placenta_dir)]
agg1[, beta := m/c]
setorder(agg1, placenta_dir, bin, tp_num)

# 每标记在 母体时间点 (NP=0,1stT=1,2ndT=2,3rdT=3,delivery=4) 的 Spearman rho
rho1 <- agg1[, {
  x <- tp_num
  if (length(unique(x)) >= 3 && length(unique(beta)) >= 2) {
    list(rho = suppressWarnings(cor(x, beta, method='spearman')), n_tp = length(unique(x)))
  } else list(rho=NA_real_, n_tp=length(unique(x)))
}, by=.(bin, placenta_dir)]
rho1 <- rho1[!is.na(rho)]
rho1[is.na(placenta_dir) | placenta_dir=='', placenta_dir := 'FLAT']

stat_plac <- rho1[, .(n_bin=.N, mean_rho=mean(rho), n_pos=sum(rho>0), n_neg=sum(rho<0),
                      sign_p=binom.test(sum(rho>0), .N, p=0.5)$p.value), by=placenta_dir]
cat('\n=== Part 1: Type.II 胎盘标记孕期轨迹验证 ===\n')
print(stat_plac)

# 各时间点组均值 (加权, 全部胎盘信息标记)
agg1_all <- agg1[, .(m=sum(m), c=sum(c)), by=.(tp_num, placenta_dir)][, beta:=m/c]
setorder(agg1_all, placenta_dir, tp_num)
cat('\n胎盘信息标记 组均值beta (NP=0,1stT=1,...,delivery=4):\n')
print(dcast(agg1_all, placenta_dir ~ tp_num, value.var='beta'))

# ---------- Part 2: Type.I 组织标记 cfDNA 均值 (判别两种行序解释) ----------
mkI <- mk[type=='I']
# 全部 134 样本在该组织标记的加权平均 beta
mc_all <- merge(mc, sm[, .(gsm)], by='gsm')
aggI <- mc_all[mkI, on='bin'][, .(m=sum(m), c=sum(c)), by=.(bin, tissue)]
aggI[, beta := m/c]
tissue_means <- aggI[, .(n_bin=.N, mean_beta=mean(beta), median_beta=median(beta),
                         mean_methyl_pct=mean(beta)*100), by=tissue]
setorder(tissue_means, -mean_beta)
cat('\n=== Part 2: Type.I 组织标记的 cfDNA 平均甲基化 (判别行序) ===\n')
print(tissue_means)

# ---------- Part 3: DMR→标记重叠 ----------
dmr <- fread('results/GSE282512_dmr_hg19.csv')[status=='OK']
dmr[, c('qchr','qstart','qend') := .(chr_hg19, as.numeric(start_hg19), as.numeric(end_hg19))]
mk2 <- mk[, .(bin, chr, start, end, type, tissue, placenta_informative)]
mk2[, chr := ifelse(startsWith(chr, 'chr'), chr, paste0('chr', chr))]  # v4 表已带前缀, 仅补缺失的
setkey(mk2, chr, start, end)
ov <- foverlaps(dmr[, .(region_id, qchr, qstart, qend)], mk2,
                by.x=c('qchr','qstart','qend'), by.y=c('chr','start','end'),
                type='any', nomatch=NULL)
ov <- merge(ov, dmr[, .(region_id, chr_hg38, start_hg38, end_hg38)], by='region_id')
# 关联 DMR 方向
dmr_dir <- unique(fread('results/GSE282512_dmr_final_v2.csv')[, .(region_id, direction, delta_beta, symbol)])
ov <- merge(ov, dmr_dir, by='region_id')
cat('\n=== Part 3: DMR(hg19) × Sun标记 重叠 ===\n')
cat('重叠标记数:', nrow(ov), ' 涉及DMR数:', length(unique(ov$region_id)), '\n')
if (nrow(ov) > 0) {
  print(ov[, .(region_id, chr_hg19=qchr, start_hg19=qstart, end_hg19=qend, bin, type, tissue, direction, delta_beta)])
  fwrite(ov, 'results/GSE154378_dmr_bin_overlap.csv')
}

# ---------- Part 4: DMR-重叠标记的孕期轨迹 (位点级) ----------
if (nrow(ov) > 0) {
  bins_dmr <- unique(ov$bin)
  agg_dmr <- mc_m[bin %in% bins_dmr & !is.na(tp_num) & group=='Normal',
                  .(m=sum(m), c=sum(c)), by=.(tp_num, bin)]
  agg_dmr[, beta := m/c]
  agg_dmr <- merge(agg_dmr, unique(ov[, .(bin, direction)]), by='bin')
  # 每 DMR-bin 的 rho
  rho_dmr <- agg_dmr[, {
    x <- tp_num
    list(rho=suppressWarnings(cor(x, beta, method='spearman')), n_tp=length(unique(x)))
  }, by=.(bin, direction)]
  cat('\n=== Part 4: DMR-bin 孕期轨迹 (Normal, 1stT→delivery) ===\n')
  print(rho_dmr)
  fwrite(rho_dmr, 'results/GSE154378_dmr_bin_trajectory.csv')
}

# ---------- Part 5: 标记级 WBC 组成轨迹 (全部标记) ----------
# 每标记: Normal 母体轨迹 rho(beta, tp)
agg_allb <- mc_m[!is.na(tp_num) & group=='Normal', .(m=sum(m), c=sum(c)), by=.(tp_num, bin)]
agg_allb[, beta := m/c]
rho_all <- agg_allb[order(bin, tp_num), {
  if (.N >= 3 && length(unique(beta)) >= 2) {
    list(rho = suppressWarnings(cor(tp_num, beta, method='spearman')), c_min = min(c))
  } else list(rho=NA_real_, c_min=min(c))
}, by=bin]
rho_all <- rho_all[!is.na(rho) & c_min >= 10]
mk_merge <- merge(mk, rho_all, by='bin')
# 与参考组织谱的相关: Δrho vs (tissue_ref - median_ref)
ref_cols <- c('Liver','Lungs','Colon','Small intestines','Pancreas','Adrenal glands',
              'Esophagus','Adipose tissues','Heart','Brain','T-cells','B-cells','Neutrophils','Placenta')
for (tc in ref_cols) {
  mk_merge[[paste0('d_', tc)]] <- mk_merge[[tc]] - apply(mk_merge[, ..ref_cols], 1, median)
}
cat('\n=== Part 5: 标记级轨迹 rho 与参考组织谱的 Spearman 相关 ===\n')
res_ref <- data.table(tissue=character(), rho_cor=numeric(), p=numeric())
for (tc in ref_cols) {
  xx <- mk_merge[!is.na(rho)]
  ct <- suppressWarnings(cor.test(xx[[paste0('d_', tc)]], xx$rho, method='spearman'))
  res_ref <- rbind(res_ref, data.table(tissue=tc, rho_cor=ct$estimate, p=ct$p.value))
}
setorder(res_ref, -rho_cor)
print(res_ref)
fwrite(res_ref, 'results/GSE154378_marker_trajectory_tissue_cor.csv')

# 轨迹上升/下降标记 在组织分类上的富集 (Type.I)
mkI_rho <- merge(mkI, rho_all, by='bin')
mkI_rho[is.na(tissue), tissue := 'FLAT']
tissue_traj <- mkI_rho[, .(n_bin=.N, mean_rho=mean(rho), n_up=sum(rho>0.3), n_down=sum(rho< -0.3),
                           sign_p=suppressWarnings(binom.test(sum(rho>0), .N)$p.value)), by=tissue]
setorder(tissue_traj, -mean_rho)
cat('\nType.I 组织标记的孕期轨迹 (Normal):\n')
print(tissue_traj)
fwrite(tissue_traj, 'results/GSE154378_typeI_tissue_trajectory.csv')

# ---------- Part 6: DMR 相关标记的 Δβ (delivery−1stT) 患者配对 ----------
if (nrow(ov) > 0) {
  # 患者配对: Normal 且 同时有 1stT 与 delivery
  pats <- sm[group=='Normal', .(has1=any(timepoint=='1stT'), hasD=any(timepoint=='delivery')), by=patient]
  pats <- pats[has1 & hasD]
  mc_p <- merge(mc, sm[group=='Normal' & patient %in% pats$patient & timepoint %in% c('1stT','delivery'),
                       .(gsm, patient, timepoint)], by='gsm')
  dcast_p <- dcast(mc_p[bin %in% bins_dmr], patient + bin ~ timepoint, value.var='beta', fun.aggregate=mean)
  dcast_p <- dcast_p[!is.na(`1stT`) & !is.na(delivery)]
  dcast_p[, delta := delivery - `1stT`]
  dcast_p <- merge(dcast_p, unique(ov[, .(bin, direction)]), by='bin')
  cat('\n=== Part 6: DMR-bin 患者配对 Δβ (delivery−1stT) ===\n')
  print(dcast_p[, .(n=length(delta), median_delta=median(delta), mean_delta=mean(delta),
                    n_pos=sum(delta>0), wilcox_p=suppressWarnings(wilcox.test(delta)$p.value)), by=.(bin, direction)])
  fwrite(dcast_p, 'results/GSE154378_dmr_bin_paired_delta.csv')
}

# ---------- 输出汇总 ----------
sink('results/GSE154378_trajectory_summary.txt')
cat('GSE154378 DMR 孕期轨迹验证 汇总\n')
cat('日期:', format(Sys.time(), '%Y-%m-%d %H:%M'), '\n\n')
cat('样本: 134 (Normal 44 / PreX 40 / GDM 33 / cHTN 10 / NP 7)\n')
cat('标记:', nrow(mk), '(Type.I', nrow(mk[type=='I']), '+ Type.II', nrow(mk[type=='II']),
    '), bin', min(mk$bin), '..', max(mk$bin), '(v4 终版表, bin=坐标排序排名, 三重验证 rho=0.89/0.82/0.81)\n\n')
cat('[Part 1] Type.II 胎盘信息标记轨迹验证\n'); print(stat_plac)
cat('\n[Part 2] Type.I 组织标记 cfDNA 均值\n'); print(tissue_means)
cat('\n[Part 3] DMR×标记重叠\n'); cat('重叠标记数:', nrow(ov), ' 涉及DMR:', length(unique(ov$region_id)), '\n')
if (nrow(ov)>0) print(ov[, .(region_id, qchr, qstart, qend, bin, type, tissue, direction)])
cat('\n[Part 5] 轨迹 rho × 参考组织谱相关\n'); print(res_ref)
cat('\n[Part 5b] Type.I 组织轨迹\n'); print(tissue_traj)
sink()

cat('\n=== 全部完成. 输出见 results/GSE154378_trajectory_summary.txt ===\n')
