#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
50k_test_v3.py — 验证 v3 bin->Sun标记映射假设
v3a: bin 1..1013 = Type.I rows 1..1013; bin 1014..5820 = Type.II rows 1..4807
v3b: bin 1..4807 = Type.II rows 1..4807; bin 4808..5820 = Type.I rows 1..1013
判据: GSE154378 deconv 观测 beta vs Sun 血细胞参考 (Neutrophils/T-cells/B-cells) Spearman
若映射正确, NP(未孕) 纯成人血 cfDNA 应显著正相关
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

def parse_xlsx():
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
    return sheet('xl/worksheets/sheet1.xml'), sheet('xl/worksheets/sheet2.xml')

pat = re.compile(r'chr([0-9XY]+):(\d+)-(\d+)')
t1, t2 = parse_xlsx()

def rows_I():
    out = []
    for r in t1[1:]:
        m = pat.match(r[1] or '')
        if not m: continue
        vals = []
        for j in range(2, 2+len(TISSUES)):
            try: vals.append(float(r[j]))
            except (ValueError, IndexError): vals.append(np.nan)
        out.append((m.group(1), int(m.group(2)), int(m.group(3)), 'I', vals))
    return out

def rows_II():
    out = []
    for r in t2[1:]:
        m = pat.match(r[0] or '')
        if not m: continue
        vals = []
        for j in range(1, 1+len(TISSUES)):
            try: vals.append(float(r[j]))
            except (ValueError, IndexError): vals.append(np.nan)
        out.append((m.group(1), int(m.group(2)), int(m.group(3)), 'II', vals))
    return out

I, II = rows_I(), rows_II()
print(f'parsed Type.I: {len(I)}, Type.II: {len(II)}')

# 候选映射: bin -> marker entry
maps = {
    'v3a (Type.I first)': {i+1: I[i] for i in range(len(I))} | {len(I)+1+j: II[j] for j in range(len(II))},
    'v3b (Type.II first)': {j+1: II[j] for j in range(len(II))} | {len(II)+1+i: I[i] for i in range(len(I))},
}

# ---- 载入 GSE154378 观测 ----
obs_all = defaultdict(lambda: [0, 0])   # 全样本
obs_np  = defaultdict(lambda: [0, 0])   # NP 组 (未孕对照, 纯成人血)
# NP 样本 gsm 列表: 从样本表读
np_gs = set()
with open(ROOT + '/data/geo_methylation/GSE154378/gse154378_samples.tsv', encoding='utf-8') as f:
    next(f)
    for line in f:
        p = line.rstrip('\n').split('\t')
        if p[1] == 'NP': np_gs.add(p[0])

with gzip.open(ROOT + '/results/GSE154378_bin_mc_long.csv.gz', 'rt', encoding='utf-8') as f:
    next(f)
    for line in f:
        gsm, b, m, c = line.split(',')
        b, m, c = int(b), int(m), int(c)
        if c < 5: continue
        obs_all[b][0] += m; obs_all[b][1] += c
        if gsm in np_gs:
            obs_np[b][0] += m; obs_np[b][1] += c

beta_all = {b: v[0]/v[1] for b, v in obs_all.items() if v[1] >= 100}
beta_np  = {b: v[0]/v[1] for b, v in obs_np.items() if v[1] >= 20}
print(f'usable bins (all, c>=100): {len(beta_all)}; NP (c>=20): {len(beta_np)}; NP samples: {len(np_gs)}')

REFS = {'Neutrophils': 12, 'T-cells': 10, 'B-cells': 11, 'Placenta': 13}
for name, mp in maps.items():
    print(f'\n== {name} ==')
    for subset, beta in [('ALL', beta_all), ('NP', beta_np)]:
        bs = [b for b in beta if b in mp]
        for tname, ti in REFS.items():
            x = [mp[b][4][ti] for b in bs]
            y = [beta[b] for b in bs]
            ok = [(a, b2) for a, b2 in zip(x, y) if not np.isnan(a)]
            if len(ok) < 20: continue
            rho, p = spearmanr([a for a, _ in ok], [b2 for _, b2 in ok])
            print(f'  {subset} vs {tname:12s}: rho={rho:+.3f} p={p:.2e} n={len(ok)}')
