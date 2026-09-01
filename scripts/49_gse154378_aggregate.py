#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
49_gse154378_aggregate.py (v2)
GSE154378 (UCLA Del Vecchio) cfDNA WGBS deconv 数据聚合
1) 解析 Sun 2015 PNAS 补充表 (pnas.1508736112.sd01.xlsx) → bin→hg19 坐标映射
   行序假设（唯一可行拼接）：观察 bins 526..5820
     = Type.I rows 1..489（→bin 526..1014，组织标记）
     + Type.II rows 1..4806（→bin 1015..5820，泛组织标记）
   Type.I rows 490..1014 与 Type.II rows 4807..4914 不在观察范围 → 表内保留但不标记
2) 聚合 134 个 *_deconv.txt.gz → per-(sample,bin) 甲基化计数 (m=甲基化CpG, c=总CpG)
3) 输出：sun2015_markers.tsv（bin坐标+组织+14组织参考甲基化值）、
        gse154378_samples.tsv（含 subtype）、GSE154378_bin_mc.npz、long CSV
"""
import gzip, os, re, glob, csv
import numpy as np
from xml.etree import ElementTree as ET
import zipfile

ROOT = r"E:/妊娠期高血压甲基化研究方案/hdp-methylation-project"
DATA = os.path.join(ROOT, "data/geo_methylation/GSE154378")
RES  = os.path.join(ROOT, "results")
os.makedirs(RES, exist_ok=True)

NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
TISSUES = ['Liver','Lungs','Colon','Small intestines','Pancreas','Adrenal glands',
           'Esophagus','Adipose tissues','Heart','Brain','T-cells','B-cells',
           'Neutrophils','Placenta']

def parse_xlsx(path):
    z = zipfile.ZipFile(path)
    root = ET.fromstring(z.read('xl/sharedStrings.xml'))
    ss = []
    for si in root.findall(NS+'si'):
        ss.append(''.join(x.text or '' for x in si.iter(NS+'t')))
    def read_sheet(fn):
        r = ET.fromstring(z.read(fn))
        rows = []
        for row in r.iter(NS+'row'):
            cells = []
            for c in row.findall(NS+'c'):
                v = c.find(NS+'v')
                val = v.text if v is not None else ''
                if c.get('t') == 's' and val:
                    val = ss[int(val)]
                cells.append(val)
            rows.append(cells)
        return rows
    return read_sheet('xl/worksheets/sheet1.xml'), read_sheet('xl/worksheets/sheet2.xml')

def parse_loc(s):
    m = re.match(r'chr([0-9XY]+):(\d+)-(\d+)', s)
    if not m: return None
    return m.group(1), int(m.group(2)), int(m.group(3))

t1, t2 = parse_xlsx(os.path.join(DATA, 'sun2015_sd01.xlsx'))
print(f'Type.I rows={len(t1)} (header+{len(t1)-1}), Type.II rows={len(t2)} (header+{len(t2)-1})')

# observed bins: 526..5820
# Type.I rows 1..489 -> bin row+525 (526..1014)
# Type.II rows 1..4806 -> bin row+1014 (1015..5820)
markers = []
for i, r in enumerate(t1[1:], start=1):
    bin_id = i + 525
    if not (526 <= bin_id <= 1014): continue      # Type.I observed bins only
    loc = parse_loc(r[1])
    if not loc: continue
    vals = []
    for j in range(2, 2+len(TISSUES)):
        try: vals.append(float(r[j]))
        except (ValueError, IndexError): vals.append('')
    markers.append([bin_id, loc[0], loc[1], loc[2], 'I', r[0], r[-1]] + vals)
for i, r in enumerate(t2[1:], start=1):
    bin_id = i + 1014
    if not (1015 <= bin_id <= 5820): continue     # only observed range
    loc = parse_loc(r[0])
    if not loc: continue
    vals = []
    for j in range(1, 1+len(TISSUES)):
        try: vals.append(float(r[j]))
        except (ValueError, IndexError): vals.append('')
    markers.append([bin_id, loc[0], loc[1], loc[2], 'II', '', r[-1]] + vals)

markers.sort(key=lambda x: x[0])
print(f'observed-range markers: {len(markers)} (expect 5295), bin {markers[0][0]}..{markers[-1][0]}')

hdr = ['bin','chr','start','end','type','tissue','placenta_informative'] + TISSUES
with open(os.path.join(DATA, 'sun2015_markers.tsv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(hdr)
    for m in markers:
        w.writerow(m)

# ---------- 2. aggregate deconv ----------
files = sorted(glob.glob(os.path.join(DATA, 'GSM*_deconv.txt.gz')))
print(f'deconv files: {len(files)}')
TPS = {'1stT','2ndT','3rdT','delivery','cordB'}

samples = []
for f in files:
    base = os.path.basename(f)[:-len('_deconv.txt.gz')]  # GSM4669310_Normal_1_2ndT
    parts = base.split('_')
    gsm, group = parts[0], parts[1]
    subtype = ''
    if group == 'PreX' and len(parts) >= 5:
        subtype = parts[2]; patient = parts[3]; timepoint = parts[4]
    elif len(parts) >= 4:
        patient = parts[2]; timepoint = parts[3]
    elif len(parts) == 3:
        patient = parts[2]; timepoint = ''
    else:
        patient = parts[2] if len(parts) > 2 else ''; timepoint = ''
    if timepoint not in TPS: timepoint = ''
    samples.append((gsm, group, subtype, patient, timepoint, f))

n = len(samples)
BINMAX = 5820
M = np.zeros((n, BINMAX+1), dtype=np.int64)
C = np.zeros((n, BINMAX+1), dtype=np.int64)
for si, (gsm, group, subtype, patient, tp, f) in enumerate(samples):
    bins = []; ms = []; cs = []
    with gzip.open(f, 'rt', errors='replace') as fh:
        for line in fh:
            sp = line.split()
            if len(sp) != 2: continue
            b = int(sp[0]); s = sp[1]
            bins.append(b); ms.append(s.count('1')); cs.append(len(s))
    b = np.array(bins, dtype=np.int64)
    M[si] = np.bincount(b, weights=np.array(ms, dtype=np.float64), minlength=BINMAX+1).astype(np.int64)
    C[si] = np.bincount(b, weights=np.array(cs, dtype=np.float64), minlength=BINMAX+1).astype(np.int64)
    if (si+1) % 25 == 0 or si == n-1:
        print(f'  [{si+1}/{n}] {gsm} {group}{"/"+subtype if subtype else ""} P{patient} {tp}: {len(b)} reads')

with open(os.path.join(DATA, 'gse154378_samples.tsv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['gsm','group','subtype','patient','timepoint','file'])
    for s in samples:
        w.writerow([s[0], s[1], s[2], s[3], s[4], os.path.basename(s[5])])

np.savez_compressed(os.path.join(RES, 'GSE154378_bin_mc.npz'),
                    M=M, C=C, sample_names=np.array([s[0] for s in samples]))

with gzip.open(os.path.join(RES, 'GSE154378_bin_mc_long.csv.gz'), 'wt', encoding='utf-8') as f:
    f.write('gsm,bin,m,c\n')
    for si in range(n):
        gsm = samples[si][0]
        idx = np.nonzero(C[si])[0]
        for b in idx:
            f.write(f'{gsm},{b},{M[si,b]},{C[si,b]}\n')
print('done.')
