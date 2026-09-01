#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
49d_rebuild_markers_v4.py — 重建 sun2015_markers.tsv (v4, 最终正确映射)
v4: bin = 5820 个坐标可解析 Sun 标记按 (染色体数值序, start) 排序的排名 (1..5820)
验证: 50l_test_sorted.py sorted_num_start -> obs154 vs Neu rho=+0.886 (n=5792)
"""
import re, zipfile, csv
import numpy as np
from xml.etree import ElementTree as ET

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
markers = []
for r in t1[1:]:
    m = pat.match(r[1] or '')
    if not m: continue
    vals = []
    for j in range(2, 2+len(TISSUES)):
        try: vals.append(float(r[j]))
        except (ValueError, IndexError): vals.append('')
    markers.append([m.group(1), int(m.group(2)), int(m.group(3)), 'I', r[0],
                    (r[-1] if len(r) > 2+len(TISSUES) else 'N') or 'N'] + vals)
for r in t2[1:]:
    m = pat.match(r[0] or '')
    if not m: continue
    vals = []
    for j in range(1, 1+len(TISSUES)):
        try: vals.append(float(r[j]))
        except (ValueError, IndexError): vals.append('')
    markers.append([m.group(1), int(m.group(2)), int(m.group(3)), 'II', '',
                    (r[-1] if len(r) > 1+len(TISSUES) else 'N') or 'N'] + vals)

def ckey(c):
    try: return (0, int(c), 0)
    except ValueError: return (1, 0, ord(c))

markers.sort(key=lambda m: ckey(m[0]) + (m[1],))
print(f'markers sorted: {len(markers)}, bin 1..{len(markers)}')

hdr = ['bin','chr','start','end','type','tissue','placenta_informative'] + TISSUES
with open(ROOT + '/data/geo_methylation/GSE154378/sun2015_markers.tsv', 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(hdr)
    for i, m in enumerate(markers, start=1):
        w.writerow([i, 'chr'+m[0], m[1], m[2], m[3], m[4], m[5]] + m[6:])
print('written: data/geo_methylation/GSE154378/sun2015_markers.tsv (v4)')

# 快速自检: bin 分布 + Type.I 数
nI = sum(1 for m in markers if m[3]=='I')
print(f'Type.I: {nI}, Type.II: {len(markers)-nI}')
# 首尾 bin
print('bin 1 :', markers[0][:5])
print('bin 5820:', markers[-1][:5])
