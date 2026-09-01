# -*- coding: utf-8 -*-
"""50n: 阳性对照 — 胎盘特异性标记在母体 cfDNA 中的孕期轨迹 (v4 映射)"""
import gzip
from collections import defaultdict
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu

ROOT = r"E:/妊娠期高血压甲基化研究方案/hdp-methylation-project"
log = lambda *a: print(*a, flush=True)

# v4 markers
TISSUES = None
mk = {}
with open(ROOT + '/data/geo_methylation/GSE154378/sun2015_markers.tsv', encoding='utf-8') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    TISSUES = hdr[7:]
    TI = {t: i for i, t in enumerate(TISSUES)}
    for line in f:
        p = line.rstrip('\n').split('\t')
        b = int(p[0]); refs = [float(x) if x else np.nan for x in p[7:]]
        mk[b] = {'chrom': p[1], 'start': int(p[2]), 'type': p[4], 'tissue': p[5], 'refs': refs}
log(f'markers: {len(mk)}')

# 胎盘特异: Placenta ref 高, 其他 13 组织(含血细胞)低
pla_spec = {}
for b, m in mk.items():
    v = m['refs']
    other = [v[TI[t]] for t in TISSUES if t != 'Placenta']
    if not np.isnan(v[TI['Placenta']]) and v[TI['Placenta']] >= 60 and np.nanmean(other) <= 20:
        pla_spec[b] = v[TI['Placenta']]
log(f'placenta-specific markers (Pla>=60, others<=20): {len(pla_spec)}')

# 样本
smap = {}
with open(ROOT + '/data/geo_methylation/GSE154378/gse154378_samples.tsv', encoding='utf-8') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    gi, gri, ti, pi = hdr.index('gsm'), hdr.index('group'), hdr.index('timepoint'), hdr.index('patient')
    for line in f:
        p = line.rstrip('\n').split('\t')
        smap[p[gi]] = (p[gri], p[ti], p[pi])

TP_ORDER = ['1stT', '2ndT', '3rdT', 'delivery', 'cordB']
obs = defaultdict(lambda: defaultdict(lambda: [0, 0]))  # gsm -> bin -> [m,c]
with gzip.open(ROOT + '/results/GSE154378_bin_mc_long.csv.gz', 'rt') as fh:
    next(fh)
    for line in fh:
        g, b, m, c = line.split(',')
        b, m, c = int(b), float(m), int(c)
        d = obs[g][b]; d[0] += m; d[1] += c

# 每时间点 Normal 组的胎盘标记聚合 beta + 胎盘 cfDNA 分数估计
log('\n== Normal 妊娠: 胎盘特异标记甲基化随孕期 ==')
res = {}
for tp in TP_ORDER:
    ms, cs = 0, 0
    for g, (gr, t, p) in smap.items():
        if gr == 'Normal' and t == tp:
            for b in pla_spec:
                d = obs[g].get(b)
                if d and d[1] > 0:
                    ms += d[0]; cs += d[1]
    res[tp] = ms / cs if cs else np.nan
    log(f'{tp:9s}: pla-marker beta = {res[tp]:.4f}  (≈ placental fraction {res[tp]/0.8:.3f})')

# NP 对照
ms = cs = 0
for g, (gr, t, p) in smap.items():
    if gr == 'NP':
        for b in pla_spec:
            d = obs[g].get(b)
            if d and d[1] > 0:
                ms += d[0]; cs += d[1]
log(f'NP       : pla-marker beta = {ms/cs if cs else float("nan"):.4f}')

# 患者配对: 1stT -> delivery 变化 (paired)
log('\n== 患者配对 1stT -> delivery (Normal) ==')
deltas = []
for g, (gr, t, p) in smap.items():
    if gr == 'Normal' and t == '1stT':
        # 找同患者 delivery
        gd = [gg for gg, (grr, tt, pp) in smap.items() if grr == 'Normal' and pp == p and tt == 'delivery']
        if not gd: continue
        gd = gd[0]
        def pla_beta(gsm):
            m_, c_ = 0, 0
            for b in pla_spec:
                d = obs[gsm].get(b)
                if d and d[1] > 0:
                    m_ += d[0]; c_ += d[1]
            return m_ / c_ if c_ else np.nan
        a, bb = pla_beta(g), pla_beta(gd)
        if not (np.isnan(a) or np.isnan(bb)):
            deltas.append(bb - a)
log(f'n pairs = {len(deltas)}, mean delta = {np.mean(deltas):+.4f}, '
    f'positive: {sum(d>0 for d in deltas)}/{len(deltas)}')

# 血细胞标记对照 (不应随孕期大变)
blood_spec = {}
for b, m in mk.items():
    v = m['refs']
    blood = [v[TI[t]] for t in ['T-cells', 'B-cells', 'Neutrophils']]
    other = [v[TI[t]] for t in TISSUES if t not in ('T-cells', 'B-cells', 'Neutrophils', 'Placenta')]
    if all(not np.isnan(x) for x in blood) and np.mean(blood) >= 70 and np.nanmean(other) <= 40:
        blood_spec[b] = np.mean(blood)
log(f'\nblood-specific markers: {len(blood_spec)}')
log('== Normal 妊娠: 血细胞标记随孕期 ==')
for tp in ['1stT', '2ndT', '3rdT', 'delivery']:
    ms, cs = 0, 0
    for g, (gr, t, p) in smap.items():
        if gr == 'Normal' and t == tp:
            for b in blood_spec:
                d = obs[g].get(b)
                if d and d[1] > 0:
                    ms += d[0]; cs += d[1]
    log(f'{tp:9s}: blood-marker beta = {ms/cs if cs else float("nan"):.4f}')
