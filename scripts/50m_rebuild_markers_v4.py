# -*- coding: utf-8 -*-
"""50m: 重建 v4 标记表（bin = 基因组坐标排序行号）+ 决定性验证"""
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

def fnum(x):
    try: return float(x)
    except: return np.nan

pat = re.compile(r'chr([0-9XY]+):(\d+)-(\d+)')
t1, t2 = sheet('xl/worksheets/sheet1.xml'), sheet('xl/worksheets/sheet2.xml')
ALL = []
for r in t1[1:]:
    m = pat.match(r[1] or '')
    if not m: continue
    ALL.append({'chrom': m.group(1), 'start': int(m.group(2)), 'end': int(m.group(3)),
                'type': 'I', 'tissue': (r[0] or '').strip(),
                'refs': [fnum(v) for v in r[2:16]]})
for r in t2[1:]:
    m = pat.match(r[0] or '')
    if not m: continue
    ALL.append({'chrom': m.group(1), 'start': int(m.group(2)), 'end': int(m.group(3)),
                'type': 'II', 'tissue': '',
                'refs': [fnum(v) for v in r[1:15]]})
log(f'total markers parsed: {len(ALL)}')

def cnum(c): return (0, int(c)) if c.isdigit() else (1, ord(c[0]))
ALL.sort(key=lambda m: (cnum(m['chrom'])[0], cnum(m['chrom'])[1], m['start']))
for i, m in enumerate(ALL, 1):
    m['bin'] = i

# ---- 写 v4 标记表（列结构与旧表一致，下游 R 脚本无感切换）----
out = ROOT + '/data/geo_methylation/GSE154378/sun2015_markers.tsv'
with open(out, 'w', encoding='utf-8', newline='') as f:
    f.write('bin\tchr\tstart\tend\ttype\ttissue\tplacenta_informative\t' + '\t'.join(TISSUES) + '\n')
    for m in ALL:
        pla = 'Y' if (m['refs'][TISSUES.index('Placenta')] >= 50 and
                      np.nanmean([m['refs'][TISSUES.index(t)] for t in
                                  ['Liver','Lungs','Colon','Small intestines','Pancreas',
                                   'Adrenal glands','Esophagus','Adipose tissues','Heart','Brain',
                                   'T-cells','B-cells','Neutrophils']]) < 50) else 'N'
        f.write(f"{m['bin']}\tchr{m['chrom']}\t{m['start']}\t{m['end']}\t{m['type']}\t"
                f"{m['tissue']}\t{pla}\t" +
                '\t'.join('' if np.isnan(v) else f'{v:.4g}' for v in m['refs']) + '\n')
log(f'v4 markers written: {out} ({len(ALL)} rows)')

# ---- 决定性验证 ----
smap = {}
with open(ROOT + '/data/geo_methylation/GSE154378/gse154378_samples.tsv', encoding='utf-8') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    gi, gri, ti = hdr.index('gsm'), hdr.index('group'), hdr.index('timepoint')
    for line in f:
        p = line.rstrip('\n').split('\t')
        smap[p[gi]] = (p[gri], p[ti])

obs = {g: defaultdict(lambda: [0, 0]) for g in ('NP', 'cordB', 'delN', 'all')}
with gzip.open(ROOT + '/results/GSE154378_bin_mc_long.csv.gz', 'rt') as fh:
    next(fh)
    for line in fh:
        g, b, m, c = line.split(',')
        b, m, c = int(b), float(m), int(c)
        gr, tp = smap.get(g, ('', ''))
        keys = ['all']
        if gr == 'NP': keys.append('NP')
        if tp == 'cordB': keys.append('cordB')
        if gr == 'Normal' and tp == 'delivery': keys.append('delN')
        for k in keys:
            d = obs[k][b]; d[0] += m; d[1] += c
beta = {k: {b: v[0] / v[1] for b, v in o.items() if v[1] >= 20} for k, o in obs.items()}

TI = {t: i for i, t in enumerate(TISSUES)}
def ref(m, name): return m['refs'][TI[name]]
def ref_blood(m):
    return np.nanmean([ref(m, 'T-cells'), ref(m, 'B-cells'), ref(m, 'Neutrophils')])

log('\n== v4 映射验证: obs vs 参考 (Spearman) ==')
for grp in ('NP', 'cordB', 'delN', 'all'):
    for refn, fn in [('blood_mean', ref_blood), ('Placenta', lambda m: ref(m, 'Placenta')),
                     ('Brain', lambda m: ref(m, 'Brain')), ('Liver', lambda m: ref(m, 'Liver'))]:
        xs, ys = [], []
        for m in ALL:
            b = m['bin']
            if b in beta[grp] and not np.isnan(fn(m)):
                xs.append(fn(m)); ys.append(beta[grp][b])
        r = spearmanr(xs, ys)
        log(f'{grp:6s} vs {refn:12s}: rho={r.statistic:+.3f} p={r.pvalue:.3g} n={len(xs)}')

# 脐血 vs 母血 delivery: 胎盘参考差值方向
log('\n== cordB - delivery (胎盘方向) ==')
db = []
for m in ALL:
    b = m['bin']
    if b in beta['cordB'] and b in beta['delN']:
        d = beta['cordB'][b] - beta['delN'][b]
        hi = ref(m, 'Placenta') >= 50
        db.append((d, hi))
hi = [d for d, h in db if h]; lo = [d for d, h in db if not h]
log(f'placenta-HI markers (n={len(hi)}): cordB-del mean={np.mean(hi):+.4f}')
log(f'placenta-LO markers (n={len(lo)}): cordB-del mean={np.mean(lo):+.4f}')
from scipy.stats import mannwhitneyu
u = mannwhitneyu(hi, lo, alternative='greater')
log(f'Mann-Whitney HI>LO: p={u.pvalue:.3g}')
