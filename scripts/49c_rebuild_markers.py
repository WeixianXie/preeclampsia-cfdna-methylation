#!/usr/bin/env python3
# 49c_rebuild_markers.py — 重建 sun2015_markers.tsv (全量 bin 1..5820)
# 修正映射: bin 1..1014 = Type.I rows 1..1014 (bin=row); bin 1015..5820 = Type.II rows 1..4806 (bin=row+1014)
# 旧版仅 526..5820 (遗漏 bin 1..525 = Type.I rows 1..525)
import gzip, os, re, csv, zipfile
from xml.etree import ElementTree as ET

ROOT = r"E:/妊娠期高血压甲基化研究方案/hdp-methylation-project"
DATA = os.path.join(ROOT, "data/geo_methylation/GSE154378")
NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
TISSUES = ['Liver','Lungs','Colon','Small intestines','Pancreas','Adrenal glands',
           'Esophagus','Adipose tissues','Heart','Brain','T-cells','B-cells','Neutrophils','Placenta']

def parse_xlsx(path):
    z = zipfile.ZipFile(path)
    root = ET.fromstring(z.read('xl/sharedStrings.xml'))
    ss = [''.join(x.text or '' for x in si.iter(NS+'t')) for si in root.findall(NS+'si')]
    def read_sheet(fn):
        r = ET.fromstring(z.read(fn))
        rows = []
        for row in r.iter(NS+'row'):
            cells = []
            for c in row.findall(NS+'c'):
                v = c.find(NS+'v'); val = v.text if v is not None else ''
                if c.get('t')=='s' and val: val = ss[int(val)]
                cells.append(val)
            rows.append(cells)
        return rows
    return read_sheet('xl/worksheets/sheet1.xml'), read_sheet('xl/worksheets/sheet2.xml')

def parse_loc(s):
    m = re.match(r'chr([0-9XY]+):(\d+)-(\d+)', s)
    if not m: return None
    return m.group(1), int(m.group(2)), int(m.group(3))

t1, t2 = parse_xlsx(os.path.join(DATA, 'sun2015_sd01.xlsx'))
print(f'Type.I rows={len(t1)-1}, Type.II rows={len(t2)-1}')

markers = []
# Type.I rows 1..1014 -> bin 1..1014
for i, r in enumerate(t1[1:], start=1):
    bin_id = i
    if not (1 <= bin_id <= 1014): continue
    loc = parse_loc(r[1])
    if not loc: continue
    vals = []
    for j in range(2, 2+len(TISSUES)):
        try: vals.append(float(r[j]))
        except (ValueError, IndexError): vals.append('')
    markers.append([bin_id, loc[0], loc[1], loc[2], 'I', r[0], r[-1]] + vals)
# Type.II rows 1..4806 -> bin 1015..5820
for i, r in enumerate(t2[1:], start=1):
    bin_id = i + 1014
    if not (1015 <= bin_id <= 5820): continue
    loc = parse_loc(r[0])
    if not loc: continue
    vals = []
    for j in range(1, 1+len(TISSUES)):
        try: vals.append(float(r[j]))
        except (ValueError, IndexError): vals.append('')
    markers.append([bin_id, loc[0], loc[1], loc[2], 'II', '', r[-1]] + vals)

markers.sort(key=lambda x: x[0])
print(f'markers: {len(markers)} (expect 5820), bin {markers[0][0]}..{markers[-1][0]}')

hdr = ['bin','chr','start','end','type','tissue','placenta_informative'] + TISSUES
with open(os.path.join(DATA, 'sun2015_markers.tsv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(hdr)
    for m in markers:
        w.writerow(m)
print('saved sun2015_markers.tsv')

# 统计 Type.I 的组织构成
from collections import Counter
c = Counter(m[5] for m in markers if m[4]=='I')
print('Type.I tissue composition:', dict(c))
