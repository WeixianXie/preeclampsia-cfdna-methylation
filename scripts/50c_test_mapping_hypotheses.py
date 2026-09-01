#!/usr/bin/env python3
# 50c_test_mapping_hypotheses.py
# 系统测试 bin->标记 拼接假设: 观测 cfDNA beta 与 14 组织参考的秩相关
# 假设空间: bin 526..5820 是 Sun 全表(顺序未知)的连续子段
#   全表 = Type.I(1014) + Type.II(4914)  (I先)  或 Type.II + Type.I (II先)
#   且 bin 526 可对齐全表任意起点 s0 (连续子段)
import gzip, os, re, csv, zipfile
import numpy as np
from xml.etree import ElementTree as ET

ROOT = r"E:/妊娠期高血压甲基化研究方案/hdp-methylation-project"
DATA = os.path.join(ROOT, "data/geo_methylation/GSE154378")
os.chdir(ROOT)
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
# 组织参考矩阵: typeI_rows(1014) x 14, typeII_rows(4914) x 14
def get_ref(rows, col0, use_plac_info_last=False):
    refs = []
    for r in rows[1:]:
        vals = []
        ok = True
        for j in range(col0, col0+14):
            try: vals.append(float(r[j]))
            except: ok = False; vals.append(np.nan)
        refs.append(vals)
    return np.array(refs)
refI = get_ref(t1, 2)     # Type.I: 14 tissues start col2
refII = get_ref(t2, 1)    # Type.II: 14 tissues start col1
print('ref shapes: I', refI.shape, 'II', refII.shape)

# 观测 beta (134 样本均值, 覆盖>=30)
mc = []
with gzip.open('results/GSE154378_bin_mc_long.csv.gz', 'rt') as f:
    next(f)
    for line in f:
        g, b, m, c = line.strip().split(',')
        mc.append((int(b), float(m), float(c)))
import collections
acc = collections.defaultdict(lambda: [0,0,0])  # bin -> [n, sum_m, sum_c]
for b, m, c in mc:
    acc[b][0] += 1; acc[b][1] += m; acc[b][2] += c
bins_obs = sorted(b for b, (n, sm, sc) in acc.items() if n >= 30 and sc > 0)
obs = np.array([acc[b][1]/acc[b][2] for b in bins_obs])
print(f'observed bins: {len(bins_obs)}, range {bins_obs[0]}..{bins_obs[-1]}')

from scipy.stats import spearmanr

def test_hypo(order, s0, label):
    """order='I_II' or 'II_I'; bin b -> row (b - s0 + 1) in concatenated table"""
    n_obs = len(bins_obs)
    ref_all = np.vstack([refI, refII]) if order=='I_II' else np.vstack([refII, refI])
    nI, nII = len(refI), len(refII)
    # bin b -> concat row index (1-based)
    rows = []
    ok_all = True
    for b in bins_obs:
        ri = b - s0 + 1
        if ri < 1 or ri > len(ref_all):
            ok_all = False; break
        rows.append(ri-1)
    if not ok_all: return None
    rr = ref_all[rows]
    out = {}
    for k, tc in enumerate(TISSUES):
        col = rr[:, k]
        m = ~np.isnan(col)
        if m.sum() > 50 and np.nanstd(col) > 1e-9:
            rho, p = spearmanr(col[m], obs[m])
            out[tc] = (rho, p, m.sum())
    return out

# 候选假设
cands = []
for order, nfirst, nsecond in [('I_II', 1014, 4914), ('II_I', 4914, 1014)]:
    for s0 in range(1, 600):
        # bin 526 -> row (526-s0+1) >= 1; bin 5820 -> row (5820-s0+1) <= 5928
        if 526-s0+1 >= 1 and 5820-s0+1 <= 5928:
            cands.append((order, s0))
print(f'candidate (order, s0): {len(cands)}')

best = []
for order, s0 in cands:
    res = test_hypo(order, s0, '')
    if res is None: continue
    # 目标: 血细胞组织 (Neutrophils, T-cells, B-cells) 的 rho 是否显著正
    neut = res.get('Neutrophils'); tcel = res.get('T-cells'); bcel = res.get('B-cells')
    if neut and tcel and bcel:
        score = (neut[0] + tcel[0] + bcel[0]) / 3
        best.append((score, order, s0, neut, tcel, bcel))
best.sort(reverse=True)
print('\n=== Top 10 假设 (按血细胞 rho 均值) ===')
for score, order, s0, n, t, b in best[:10]:
    print(f'score={score:+.3f} order={order} s0={s0} | Neu rho={n[0]:+.3f}(p={n[1]:.1e}) T rho={t[0]:+.3f} B rho={b[0]:+.3f}')

# 也检查最差/当前假设
print('\n=== 当前假设 (I_II, s0=1: Type.I rows 526-1014 + Type.II rows 1-4806) ===')
cur = test_hypo('I_II', 1, 'cur')
for tc in ['Neutrophils','T-cells','B-cells','Placenta']:
    r, p, nn = cur[tc]
    print(f'{tc}: rho={r:+.3f} p={p:.2e} n={nn}')
