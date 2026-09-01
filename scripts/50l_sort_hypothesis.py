# -*- coding: utf-8 -*-
"""50l: 检验 GSE154378 bin 编号 = 基因组坐标排序假设"""
import re, gzip, zipfile
from collections import defaultdict
from xml.etree import ElementTree as ET
import numpy as np
from scipy.stats import spearmanr

ROOT = r"E:/妊娠期高血压甲基化研究方案/hdp-methylation-project"
log = lambda *a: print(*a, flush=True)

NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
z = zipfile.ZipFile(ROOT + '/data/geo_methylation/GSE154378/sun2015_sd01.xlsx')
root = ET.fromstring(z.read('xl/sharedStrings.xml'))
ss = [''.join(x.text or '' for x in si.iter(NS + 't')) for si in root.findall(NS + 'si')]

TISSUES = ['Liver','Lungs','Colon','Small intestines','Pancreas','Adrenal glands',
           'Esophagus','Adipose tissues','Heart','Brain','T-cells','B-cells',
           'Neutrophils','Placenta']

def sheet(fn):
    r = ET.fromstring(z.read(fn)); rows = []
    for row in r.iter(NS + 'row'):
        cells = []
        for c in row.findall(NS + 'c'):
            v = c.find(NS + 'v'); val = v.text if v is not None else ''
            if c.get('t') == 's' and val: val = ss[int(val)]
            cells.append(val)
        rows.append(cells)
    return rows

pat = re.compile(r'chr([0-9XY]+):(\d+)-(\d+)')
def fnum(x):
    try: return float(x)
    except: return np.nan

t1, t2 = sheet('xl/worksheets/sheet1.xml'), sheet('xl/worksheets/sheet2.xml')
# Type.I: loc in col idx1, tissues start at col idx2 (header: ?, locus, tissues...)
# Type.II: loc in col idx0, tissues start col idx1
I, II = [], []
for r in t1[1:]:
    m = pat.match(r[1] or '')
    if not m: continue
    vals = [fnum(v) for v in r[2:16]]
    I.append({'chrom': m.group(1), 'start': int(m.group(2)), 'end': int(m.group(3)),
              'type': 'I', 'tissue': (r[0] or '').strip(), 'refs': vals})
for r in t2[1:]:
    m = pat.match(r[0] or '')
    if not m: continue
    vals = [fnum(v) for v in r[1:15]]
    II.append({'chrom': m.group(1), 'start': int(m.group(2)), 'end': int(m.group(3)),
               'type': 'II', 'tissue': '', 'refs': vals})
log(f'parsed Type.I={len(I)}, Type.II={len(II)}, total={len(I)+len(II)}')

# ---- 检查 sheet 本身是否已按坐标排序 ----
for lab, lst in [('Type.I sheet order', I), ('Type.II sheet order', II)]:
    key = [(c['chrom'], c['start']) for c in lst]
    log(f'{lab}: is_lex_sorted={key == sorted(key, key=lambda t:(t[0], t[1]))}')

# ---- 候选排序 ----
def cnum(c): return (0, int(c)) if c.isdigit() else (1, c)
ALL = I + II
orderings = {
    'row_I_first': ALL[:],
    'row_II_first': II + I,
    'sorted_numeric': sorted(ALL, key=lambda m: (cnum(m['chrom']), m['start'])),
    'sorted_lex': sorted(ALL, key=lambda m: (m['chrom'], m['start'])),
}
# 单类型排序假设: 只有 Type.II? 不可能(5820=1013+4807), 跳过

# ---- obs154: NP 组 (非孕成人血) 与 全部样本 ----
smap = {}
with open(ROOT + '/data/geo_methylation/GSE154378/gse154378_samples.tsv', encoding='utf-8') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    gi, gri, ti = hdr.index('gsm'), hdr.index('group'), hdr.index('timepoint')
    for line in f:
        p = line.rstrip('\n').split('\t')
        smap[p[gi]] = (p[gri], p[ti])
groups = defaultdict(list)
for g, (gr, tp) in smap.items(): groups[gr].append(g)
log('groups: ' + ', '.join(f'{k}={len(v)}' for k, v in sorted(groups.items())))

obs_all = defaultdict(lambda: [0, 0]); obs_np = defaultdict(lambda: [0, 0])
obs_del = defaultdict(lambda: [0, 0])
with gzip.open(ROOT + '/results/GSE154378_bin_mc_long.csv.gz', 'rt') as f:
    next(f)
    for line in f:
        g, b, m, c = line.split(',')
        b, m, c = int(b), float(m), int(c)
        gr, tp = smap.get(g, ('', ''))
        for tgt, cond in ((obs_all, True),
                          (obs_np, gr == 'NP'),
                          (obs_del, tp == 'delivery' and gr == 'Normal')):
            if cond:
                d = tgt[b]; d[0] += m; d[1] += c
beta_all = {b: v[0] / v[1] for b, v in obs_all.items() if v[1] >= 20}
beta_np = {b: v[0] / v[1] for b, v in obs_np.items() if v[1] >= 20}
beta_del = {b: v[0] / v[1] for b, v in obs_del.items() if v[1] >= 20}
log(f'beta_all={len(beta_all)}, beta_np={len(beta_np)}, beta_del={len(beta_del)}')

TI = {t: i for i, t in enumerate(TISSUES)}
def ref_blood(m): 
    v = m['refs']; return np.nanmean([v[TI['T-cells']], v[TI['B-cells']], v[TI['Neutrophils']]])
def ref_plac(m): return m['refs'][TI['Placenta']]
def ref_mean(m): return np.nanmean(m['refs'])

log('\n== 排序假设检验 (obs vs blood-mean reference) ==')
best = []
for name, ol in orderings.items():
    for obs, oname in ((beta_np, 'NP'), (beta_all, 'ALL'), (beta_del, 'DEL')):
        xs, ys = [], []
        for i, m in enumerate(ol, 1):
            if i in obs and not np.isnan(ref_blood(m)):
                xs.append(ref_blood(m)); ys.append(obs[i])
        if len(xs) < 50: continue
        r = spearmanr(xs, ys)
        best.append((r.statistic, name, oname, len(xs)))
        log(f'{name:16s} vs {oname}: rho={r.statistic:+.3f} p={r.pvalue:.3g} n={len(xs)}')

log('\n== 最佳: ' + str(sorted(best, reverse=True)[:3]))
