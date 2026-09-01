#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
50l_test_sorted.py — 测试 bin 按基因组坐标排序的假设
假设: 5820 个可解析 Sun 标记按 (染色体, start) 排序后 bin=排名
变体: 数值染色体序 vs 字典序(chr1,chr10,chr2..)、含 Type.I+II 合并
判据: GSE154378 obs beta vs Sun 血细胞参考 Spearman
"""
import gzip, re, zipfile
import numpy as np
from collections import defaultdict
from xml.etree import ElementTree as ET
from scipy.stats import spearmanr

ROOT = r"E:/妊娠期高血压甲基化研究方案/hdp-methylation-project"
NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
TISSUES = ['Liver','Lungs','Colon','Small intestines','Pancreas','Adrenal glands',
           'Esophagus','Adipose tissues','Heart','Brain','T-cells','B-cells',
           'Neutrophils','Placenta']

z = zipfile.ZipFile(ROOT + '/data/geo_methylation/GSE154378/sun2015_sd01.xlsx')
root = ET.fromstring(z.read('xl/sharedStrings.xml'))
ss = [''.join(x.text or '' for x in si.iter(NS+'t')) for si in root.findall(NS+'si')]
def sheet(fn):
    r = ET.fromstring(z.read(fn)); rows = []
    for row in r.iter(NS+'row'):
        cells = []
        for c in row.findall(NS+'c'):
            v = c.find(NS+'v'); val = v.text if v is not None else ''
            if c.get('t') == 's' and val: val = ss[int(val)]
            cells.append(val)
        rows.append(cells)
    return rows
t1, t2 = sheet('xl/worksheets/sheet1.xml'), sheet('xl/worksheets/sheet2.xml')

pat = re.compile(r'chr([0-9XY]+):(\d+)-(\d+)')
markers = []  # (chrom, start, end, type, tissue, refs[14])
for r in t1[1:]:
    m = pat.match(r[1] or '')
    if not m: continue
    vals = []
    for j in range(2, 2+len(TISSUES)):
        try: vals.append(float(r[j]))
        except (ValueError, IndexError): vals.append(np.nan)
    markers.append((m.group(1), int(m.group(2)), int(m.group(3)), 'I', r[0], vals))
for r in t2[1:]:
    m = pat.match(r[0] or '')
    if not m: continue
    vals = []
    for j in range(1, 1+len(TISSUES)):
        try: vals.append(float(r[j]))
        except (ValueError, IndexError): vals.append(np.nan)
    markers.append((m.group(1), int(m.group(2)), int(m.group(3)), 'II', '', vals))
print(f'markers: {len(markers)}')

def ckey_num(c):
    try: return (0, int(c), 0)
    except ValueError: return (1, 0, ord(c))

hyps = {
    'sorted_num_start': sorted(markers, key=lambda m: ckey_num(m[0])[:2] + (m[1],)),
    'sorted_num_end':   sorted(markers, key=lambda m: ckey_num(m[0])[:2] + (m[2],)),
    'sorted_lex_start': sorted(markers, key=lambda m: (m[0], m[1])),
}
# Type.II-only 已排序检查
ii = [m for m in markers if m[3]=='II']
print('Type.II sheet sorted by (chr,lex)? :', ii == sorted(ii, key=lambda m:(m[0], m[1])))
ii_num = [m for m in sorted(ii, key=lambda m: ckey_num(m[0])[:2]+(m[1],))]
print('Type.II sheet sorted by (chr,num)? :', ii == ii_num)

# ---- obs154 ----
obs = defaultdict(lambda: [0, 0])
np_gs = set()
with open(ROOT + '/data/geo_methylation/GSE154378/gse154378_samples.tsv', encoding='utf-8') as f:
    next(f)
    for line in f:
        p = line.rstrip('\n').split('\t')
        if p[1] == 'NP': np_gs.add(p[0])
obs_np = defaultdict(lambda: [0, 0])
with gzip.open(ROOT + '/results/GSE154378_bin_mc_long.csv.gz', 'rt') as f:
    next(f)
    for line in f:
        g, b, m, c = line.split(',')
        b, c = int(b), int(c)
        if c < 20: continue
        d = obs[b]; d[0] += int(m); d[1] += c
        if g in np_gs:
            d2 = obs_np[b]; d2[0] += int(m); d2[1] += c
beta = {b: v[0]/v[1] for b, v in obs.items() if v[1] >= 100}
beta_np = {b: v[0]/v[1] for b, v in obs_np.items() if v[1] >= 20}
print(f'bins usable ALL: {len(beta)}, NP: {len(beta_np)}')

REFS = {'Neutrophils': 12, 'T-cells': 10, 'B-cells': 11, 'Placenta': 13}
for name, ml in hyps.items():
    print(f'\n== {name} ==')
    for lab, bb in [('ALL', beta), ('NP', beta_np)]:
        bs = [b for b in bb if 1 <= b <= len(ml)]
        for tc, ti in REFS.items():
            x = [ml[b-1][5][ti] for b in bs]
            y = [bb[b] for b in bs]
            ok = [(a, b2) for a, b2 in zip(x, y) if not np.isnan(a)]
            rho, p = spearmanr([a for a, _ in ok], [b2 for _, b2 in ok])
            print(f'  {lab} vs {tc:12s}: rho={rho:+.3f} p={p:.2e} n={len(ok)}')
