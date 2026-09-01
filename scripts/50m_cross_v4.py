#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""50m_cross_v4.py — v4 映射下的交叉队列确认: obs154 vs GSE282512 (60-bin oracle)"""
import gzip, csv, re, zipfile
import numpy as np
from collections import defaultdict
from scipy.stats import spearmanr
from xml.etree import ElementTree as ET

ROOT = r"E:/妊娠期高血压甲基化研究方案/hdp-methylation-project"
NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'

obs = {}
for row in csv.DictReader(open(ROOT + '/results/GSE282512_60bin_obs.csv')):
    try:
        oc = float(row['obs_Control']); ope = float(row['obs_PE'])
    except ValueError:
        continue
    obs[int(row['bin'])] = (oc, ope)

obs154 = defaultdict(lambda: [0, 0])
with gzip.open(ROOT + '/results/GSE154378_bin_mc_long.csv.gz', 'rt') as f:
    next(f)
    for line in f:
        g, b, m, c = line.split(',')
        b, c = int(b), int(c)
        if c < 20: continue
        d = obs154[b]; d[0] += int(m); d[1] += c
obs154 = {b: v[0]/v[1] for b, v in obs154.items() if v[1] >= 100}

mk = {}
with open(ROOT + '/data/geo_methylation/GSE154378/sun2015_markers.tsv', encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        mk[(row['chr'], int(row['start']), int(row['end']))] = int(row['bin'])

# v2 坐标表重建 (Type.I rows 1..1014 -> bin; Type.II rows 1..4806 -> bin+1014)
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
v2 = {}
for i, r in enumerate(t1[1:], start=1):
    m = pat.match(r[1] or '')
    if m and i <= 1014: v2[i] = ('chr'+m.group(1), int(m.group(2)), int(m.group(3)))
for i, r in enumerate(t2[1:], start=1):
    m = pat.match(r[0] or '')
    if m and i <= 4806: v2[i+1014] = ('chr'+m.group(1), int(m.group(2)), int(m.group(3)))

x_ctrl, x_pe, y154 = [], [], []
for b2, (oc, ope) in obs.items():
    if b2 in v2 and v2[b2] in mk:
        b4 = mk[v2[b2]]
        if b4 in obs154:
            x_ctrl.append(oc); x_pe.append(ope); y154.append(obs154[b4])
print('matched regions:', len(x_ctrl))
for lab, xs in [('GSE282512_Control', x_ctrl), ('GSE282512_PE', x_pe)]:
    r = spearmanr(xs, y154)
    print(f'obs154(v4) vs {lab}: rho={r.statistic:+.3f} p={r.pvalue:.2e} n={len(xs)}')
