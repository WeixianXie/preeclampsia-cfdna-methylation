#!/usr/bin/env python3
# 50g_cross_cohort_validation.py
# 交叉队列验证: 抽样 sun 标记 bin, lift hg19->hg38, 在 GSE282512 (hg38 WGBS) 中聚合,
# 检验 GSE282512 cfDNA 观测 vs sun 血细胞参考 (若映射正确应正相关)
# 同时与 GSE154378 观测比较
import gzip, os, re, csv, json, time, glob, random, urllib.request
import numpy as np
from collections import defaultdict
from scipy.stats import spearmanr

ROOT = r"E:/妊娠期高血压甲基化研究方案/hdp-methylation-project"
DATA = os.path.join(ROOT, "data/geo_methylation")
os.chdir(ROOT)
random.seed(42)
BASE = 'https://grch37.rest.ensembl.org'

def fetch(url, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read().decode())
        except Exception:
            time.sleep(4 * (k + 1))
    return None

def liftover(chrom, start, end, frm='GRCh37', to='GRCh38'):
    c = chrom[3:] if chrom.startswith('chr') else chrom
    j = fetch(f'{BASE}/map/human/{frm}/{c}:{start}..{end}/{to}?content-type=application/json')
    if not j or 'mappings' not in j:
        return None
    ms = []
    for m in j['mappings']:
        if 'mapped' in m:
            mm = m['mapped']
            ch = mm['seq_region_name']
            ms.append((('chr'+ch) if not ch.startswith('chr') else ch, mm['start'], mm['end']))
    return ms

# 1) 载入标记
mk = {}
with open(os.path.join(DATA, 'GSE154378/sun2015_markers.tsv'), encoding='utf-8') as f:
    rd = csv.DictReader(f, delimiter='\t')
    for r in rd:
        mk[int(r['bin'])] = r
typeI_bins = [b for b, r in mk.items() if r['type'] == 'I']
typeII_bins = [b for b, r in mk.items() if r['type'] == 'II']
sample = random.sample(typeI_bins, 20) + random.sample(typeII_bins, 40)
print(f'sampled bins: {len(sample)}')

# 2) lift
liftover_map = {}
for b in sample:
    r = mk[b]
    res = liftover(r['chr'], int(r['start']), int(r['end']))
    liftover_map[b] = res
    time.sleep(0.4)
ok = {b: v for b, v in liftover_map.items() if v}
print(f'liftover OK: {len(ok)}/{len(sample)}')

# 3) 聚合 GSE282512 (Control 与 PE)
annot = {}
with gzip.open(os.path.join(DATA, 'GSE282512_sample_annot.csv.gz'), 'rt', encoding='utf-8') as f:
    rd = csv.DictReader(f)
    for r in rd:
        annot[r['Sample_id']] = r['Category']
cov_files = sorted(glob.glob(os.path.join(DATA, 'GSE282512_raw/*.cov.gz')))
# 匹配 cov -> sample_id -> category
cat_of_cov = {}
for fp in cov_files:
    base = os.path.basename(fp)[:-7]  # GSM8644973_DNA031134
    sid = base.split('_', 1)[1]
    if sid in annot:
        cat_of_cov[fp] = annot[sid]
ctrl = [fp for fp, c in cat_of_cov.items() if c == 'Control']
pe = [fp for fp, c in cat_of_cov.items() if c == 'PE']
print(f'cov files: Control {len(ctrl)} PE {len(pe)}')

# 每区域: 收集 [start_hg38, end_hg38] (可能多 block, 取最大)
regions = {}
for b, maps in ok.items():
    best = max(maps, key=lambda x: x[2] - x[1])
    regions[b] = (best[0], best[1], best[2])

def aggregate(files, regions, label):
    acc = {b: [0, 0] for b in regions}  # [m, c]
    for fp in files:
        base = os.path.basename(fp)
        with gzip.open(fp, 'rt', errors='replace') as fh:
            for line in fh:
                p = line.split('\t')
                if len(p) < 5: continue
                ch, st, en = p[0], int(p[1]), int(p[2])
                meth = float(p[3]); cm = int(p[4]); cu = int(p[5])
                if cm + cu < 3: continue
                for b, (rch, rs, re_) in regions.items():
                    if ch == rch and st < re_ and rs < en:
                        acc[b][0] += int(round(meth / 100 * (cm + cu)))
                        acc[b][1] += cm + cu
    return {b: (v[0]/v[1] if v[1] > 0 else np.nan) for b, v in acc.items()}

obs_ctrl = aggregate(ctrl, regions, 'Control')
obs_pe = aggregate(pe, regions, 'PE')
print('aggregation done')

# 4) 相关检验: obs vs sun 血细胞参考
refs = ['Neutrophils', 'T-cells', 'B-cells']
print('\n== GSE282512 Control obs vs 参考 ==')
for tc in refs:
    xs = [mk[b][tc] for b in regions]; ys = [obs_ctrl[b] for b in regions]
    okk = [(float(x), y) for x, y in zip(xs, ys) if x and not np.isnan(y)]
    if len(okk) > 10:
        rho, p = spearmanr([x[0] for x in okk], [x[1] for x in okk])
        print(f'{tc}: rho={rho:+.3f} p={p:.2e} n={len(okk)}')

# 5) GSE282512 Control vs PE 的相关 (数据质量检查: 若 PE vs Control 全无差异则数据可能无信号)
xs = [obs_ctrl[b] for b in regions]; ys = [obs_pe[b] for b in regions]
okk = [(x, y) for x, y in zip(xs, ys) if not np.isnan(x) and not np.isnan(y)]
if len(okk) > 10:
    rho, p = spearmanr([x[0] for x in okk], [x[1] for x in okk])
    print(f'\nControl vs PE obs: rho={rho:+.3f} p={p:.2e} n={len(okk)}')

# 6) GSE282512 Control obs vs GSE154378 obs (同 60 bin)
# 载入 GSE154378 观测
obs154 = {}
with gzip.open(os.path.join(ROOT, 'results/GSE154378_bin_mc_long.csv.gz'), 'rt') as f:
    next(f)
    acc = defaultdict(lambda: [0, 0, 0])
    for line in f:
        g, b, m, c = line.strip().split(',')
        acc[int(b)][0] += 1; acc[int(b)][1] += float(m); acc[int(b)][2] += float(c)
for b, v in acc.items():
    if v[0] >= 30 and v[2] > 0:
        obs154[b] = v[1] / v[2]
xs = [obs154[b] for b in regions if b in obs154]
ys = [obs_ctrl[b] for b in regions if b in obs154]
okk = [(x, y) for x, y in zip(xs, ys) if not np.isnan(y)]
if len(okk) > 10:
    rho, p = spearmanr([x[0] for x in okk], [x[1] for x in okk])
    print(f'GSE154378(all) vs GSE282512(Control): rho={rho:+.3f} p={p:.2e} n={len(okk)}')
