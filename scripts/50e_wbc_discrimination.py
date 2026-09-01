#!/usr/bin/env python3
# 50e_wbc_discrimination.py — 血细胞判别 (决定性映射验证)
# 若映射正确: NP(纯成人血) cfDNA 在血细胞特异标记(Z>0)上高甲基化 > 实体组织特异标记(Z>0)
import gzip, os, re, zipfile, numpy as np
from xml.etree import ElementTree as ET
from collections import defaultdict
from scipy.stats import spearmanr

ROOT = r"E:/妊娠期高血压甲基化研究方案/hdp-methylation-project"
DATA = os.path.join(ROOT, "data/geo_methylation/GSE154378")
NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'

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

t1, _ = parse_xlsx(os.path.join(DATA, 'sun2015_sd01.xlsx'))
# Type.I: bin=row, columns: 0 Tissue, 1 Genomic location, 2-15 tissues, 16 Mean, 17 SD, 18 Z, 19 Placenta-info
ti = {}
for i, r in enumerate(t1[1:], start=1):
    if i > 1014: break
    m = re.match(r'chr([0-9XY]+):(\d+)-(\d+)', r[1])
    if not m: continue
    try: z = float(r[18])
    except: z = np.nan
    ti[i] = dict(tissue=r[0], z=z)

# NP 观测
np_gsm = set()
with open(os.path.join(DATA, 'gse154378_samples.tsv'), encoding='utf-8') as f:
    next(f)
    for line in f:
        p = line.strip().split('\t')
        if p[1] == 'NP': np_gsm.add(p[0])
acc = defaultdict(lambda: [0, 0, 0])
with gzip.open(os.path.join(ROOT, 'results/GSE154378_bin_mc_long.csv.gz'), 'rt') as f:
    next(f)
    for line in f:
        g, b, m, c = line.strip().split(',')
        b = int(b)
        if g in np_gsm:
            acc[b][0] += 1; acc[b][1] += float(m); acc[b][2] += float(c)
obs = {b: (s[1]/s[2] if s[2] > 0 else np.nan) for b, s in acc.items() if s[0] >= 4}

BLOOD = {'T-cells', 'B-cells', 'Neutrophils'}
SOLID = {'Liver', 'Lungs', 'Colon', 'Small intestines', 'Pancreas', 'Adrenal glands',
         'Esophagus', 'Adipose tissues', 'Heart', 'Brain'}
rows = []
for b, info in ti.items():
    if b not in obs: continue
    rows.append((b, info['tissue'], info['z'], obs[b]))
rows = [r for r in rows if not np.isnan(r[2]) and r[1] in BLOOD | SOLID]
print(f'Type.I bins with NP obs + Z: {len(rows)}')

# 1) 血细胞 vs 实体 组织特异标记 (Z>0: 高甲基化标记) 的 NP obs beta
for zcond, zlab in [(lambda z: z > 0, 'Z>0 (target HIGH)'), (lambda z: z < 0, 'Z<0 (target LOW)')]:
    bl = [r[3] for r in rows if r[1] in BLOOD and zcond(r[2])]
    sd = [r[3] for r in rows if r[1] in SOLID and zcond(r[2])]
    import statistics
    print(f'\n== {zlab} ==')
    print(f'blood markers n={len(bl)} mean={statistics.mean(bl):.3f} median={statistics.median(bl):.3f}')
    print(f'solid markers n={len(sd)} mean={statistics.mean(sd):.3f} median={statistics.median(sd):.3f}')
    if bl and sd:
        from scipy.stats import mannwhitneyu
        u, p = mannwhitneyu(bl, sd, alternative='greater')
        print(f'Mann-Whitney (blood>solid) p={p:.2e}')

# 2) 组织级均值
import collections, statistics
by_tissue = collections.defaultdict(list)
for b, tc, z, ob in rows:
    by_tissue[(tc, z > 0)].append(ob)
print('\n== 组织级 NP obs beta ==')
for key in sorted(by_tissue):
    tc, zpos = key
    v = by_tissue[key]
    print(f'{tc:16s} Z>0: n={len(v):3d} mean={statistics.mean(v):.3f} | Z<0: n={len(by_tissue[(tc, not zpos)]):3d} mean={statistics.mean(by_tissue[(tc, not zpos)]):.3f}' if (tc, not zpos) in by_tissue else f'{tc:16s} Z>0: n={len(v):3d} mean={statistics.mean(v):.3f}')

# 3) 综合: 全 Type.I bin, obs vs Z score 相关
zs = [r[2] for r in rows]; obs_v = [r[3] for r in rows]
rho, p = spearmanr(zs, obs_v)
print(f'\n== obs vs Z score (全 Type.I) ==\nrho={rho:+.3f} p={p:.2e} n={len(rows)}')
