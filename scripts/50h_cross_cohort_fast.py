#!/usr/bin/env python3
# 50h_cross_cohort_fast.py
# 交叉队列映射验证(快速版): 抽样 60 个 sun 标记 bin(同 50g, seed 42),
# lift hg19->hg38, 多进程聚合 GSE282512 全部 369 cov 文件(按染色体索引),
# 检验 GSE282512 Control/PE 观测 vs sun 血细胞参考的相关性
# 若映射正确: Control(正常妊娠母体血) 应与 Neutrophils/T-cells/B-cells 参考正相关
import gzip, os, csv, json, time, glob, random, sys
import numpy as np
from collections import defaultdict
from multiprocessing import Pool
from scipy.stats import spearmanr

ROOT = r"E:/妊娠期高血压甲基化研究方案/hdp-methylation-project"
DATA = os.path.join(ROOT, "data/geo_methylation")
os.chdir(ROOT)
random.seed(42)
BASE = 'https://grch37.rest.ensembl.org'
log = lambda *a: print(*a, flush=True)

def fetch(url, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, timeout=45) as r:
                return json.loads(r.read().decode())
        except Exception:
            time.sleep(3 * (k + 1))
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

# 1) 载入标记(重建版 bin 1..5820)
mk = {}
with open(os.path.join(DATA, 'GSE154378/sun2015_markers.tsv'), encoding='utf-8') as f:
    rd = csv.DictReader(f, delimiter='\t')
    for r in rd:
        mk[int(r['bin'])] = r
typeI_bins = [b for b, r in mk.items() if r['type'] == 'I']
typeII_bins = [b for b, r in mk.items() if r['type'] == 'II']
sample = random.sample(typeI_bins, 20) + random.sample(typeII_bins, 40)
log(f'sampled bins: {len(sample)}')

# 2) lift
liftover_map = {}
for b in sample:
    r = mk[b]
    res = liftover(r['chr'], int(r['start']), int(r['end']))
    liftover_map[b] = res
    time.sleep(0.3)
ok = {b: v for b, v in liftover_map.items() if v}
log(f'liftover OK: {len(ok)}/{len(sample)}')
if len(ok) < 30:
    log('FATAL: too few lifted regions'); sys.exit(1)

# 每 bin 取最大 block
regions = {}
for b, maps in ok.items():
    best = max(maps, key=lambda x: x[2] - x[1])
    regions[b] = (best[0], best[1], best[2])

# 按染色体索引: chr -> [(bin, start, end)]
regions_by_chr = defaultdict(list)
for b, (ch, s, e) in regions.items():
    regions_by_chr[ch].append((b, s, e))
log(f'regions by chr: {len(regions_by_chr)} chrs')

# 3) 样本分组
annot = {}
with gzip.open(os.path.join(DATA, 'GSE282512_sample_annot.csv.gz'), 'rt', encoding='utf-8') as f:
    rd = csv.DictReader(f)
    for r in rd:
        annot[r['Sample_id']] = r['Category']
cov_files = sorted(glob.glob(os.path.join(DATA, 'GSE282512_raw/*.cov.gz')))
cat_of_cov = {}
for fp in cov_files:
    base = os.path.basename(fp)[:-7]
    sid = base.split('_', 1)[1]
    if sid in annot:
        cat_of_cov[fp] = annot[sid]
ctrl = sorted(fp for fp, c in cat_of_cov.items() if c == 'Control')
pe = sorted(fp for fp, c in cat_of_cov.items() if c == 'PE')
log(f'cov files: Control {len(ctrl)} PE {len(pe)}  (total {len(cov_files)})')

# 4) 单文件聚合: 返回 {bin: [m, c]}
def aggregate_one(fp):
    acc = defaultdict(lambda: [0, 0])
    with gzip.open(fp, 'rt', errors='replace') as fh:
        for line in fh:
            p = line.split('\t')
            if len(p) < 6:
                continue
            ch = p[0]
            rl = regions_by_chr.get(ch)
            if not rl:
                continue
            try:
                st = int(p[1]); en = int(p[2])
                meth = float(p[3]); cm = int(p[4]); cu = int(p[5])
            except ValueError:
                continue
            d = cm + cu
            if d < 3:
                continue
            for b, rs, re_ in rl:
                if st < re_ and rs < en:
                    acc[b][0] += int(round(meth / 100.0 * d))
                    acc[b][1] += d
    return dict(acc)

t0 = time.time()
with Pool(8) as pool:
    results = pool.map(aggregate_one, ctrl + pe)
log(f'aggregation done in {time.time()-t0:.0f}s')

ctrl_n, pe_n = len(ctrl), len(pe)
sum_ctrl = defaultdict(lambda: [0, 0])
sum_pe = defaultdict(lambda: [0, 0])
for i, r in enumerate(results):
    tgt = sum_ctrl if i < ctrl_n else sum_pe
    for b, (m, c) in r.items():
        tgt[b][0] += m; tgt[b][1] += c

obs_ctrl = {b: (v[0]/v[1] if v[1] > 0 else np.nan) for b, v in sum_ctrl.items()}
obs_pe = {b: (v[0]/v[1] if v[1] > 0 else np.nan) for b, v in sum_pe.items()}

out_rows = []
def test(label, xfun, yfun, bins):
    xs = [xfun(b) for b in bins]; ys = [yfun(b) for b in bins]
    okk = [(float(x), y) for x, y in zip(xs, ys) if x is not None and x != '' and not np.isnan(y)]
    if len(okk) > 10:
        rho, p = spearmanr([x[0] for x in okk], [x[1] for x in okk])
        log(f'{label}: rho={rho:+.3f} p={p:.2e} n={len(okk)}')
        out_rows.append([label, f'{rho:.4f}', f'{p:.2e}', len(okk)])
    else:
        log(f'{label}: n too small ({len(okk)})')
        out_rows.append([label, 'NA', 'NA', len(okk)])

bins_all = sorted(regions.keys())
refs = ['Neutrophils', 'T-cells', 'B-cells']
log('\n== GSE282512 Control obs vs 参考 ==')
for tc in refs:
    test(f'Control_vs_{tc}', lambda b: mk[b][tc], lambda b: obs_ctrl[b], bins_all)
log('\n== GSE282512 PE obs vs 参考 ==')
for tc in refs:
    test(f'PE_vs_{tc}', lambda b: mk[b][tc], lambda b: obs_pe[b], bins_all)

# Control vs PE obs(数据质量)
log('\n== 数据质量 ==')
test('Control_vs_PE_obs', lambda b: obs_ctrl[b], lambda b: obs_pe[b], bins_all)

# GSE154378 obs(全样本加权) vs GSE282512 Control
obs154 = {}
with gzip.open(os.path.join(ROOT, 'results/GSE154378_bin_mc_long.csv.gz'), 'rt') as f:
    next(f)
    acc = defaultdict(lambda: [0, 0, 0])
    for line in f:
        parts = line.strip().split(',')
        if len(parts) < 4: continue
        b = int(parts[1])
        acc[b][0] += 1; acc[b][1] += float(parts[2]); acc[b][2] += float(parts[3])
for b, v in acc.items():
    if v[0] >= 30 and v[2] > 0:
        obs154[b] = v[1] / v[2]
b154 = [b for b in bins_all if b in obs154]
log(f'\n== GSE154378(all) vs GSE282512(Control) == n_bins={len(b154)}')
if len(b154) > 10:
    xs = [obs154[b] for b in b154]; ys = [obs_ctrl[b] for b in b154]
    okk = [(x, y) for x, y in zip(xs, ys) if not np.isnan(y)]
    rho, p = spearmanr([x[0] for x in okk], [x[1] for x in okk])
    log(f'GSE154378_vs_GSE282512_Control: rho={rho:+.3f} p={p:.2e} n={len(okk)}')
    out_rows.append(['GSE154378_vs_GSE282512_Control', f'{rho:.4f}', f'{p:.2e}', len(okk)])

with open(os.path.join(ROOT, 'results/GSE282512_marker_mapping_validation.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['test', 'rho', 'p', 'n'])
    w.writerows(out_rows)
log('\nsaved results/GSE282512_marker_mapping_validation.csv')

# 存每 bin 观测供后续检查
with open(os.path.join(ROOT, 'results/GSE282512_60bin_obs.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['bin', 'chr_hg38', 'start_hg38', 'end_hg38', 'type',
                'Neutrophils', 'T-cells', 'B-cells', 'Placenta',
                'obs_Control', 'obs_PE'])
    for b in bins_all:
        r = mk[b]; ch, s, e = regions[b]
        w.writerow([b, ch, s, e, r['type'],
                    r['Neutrophils'], r['T-cells'], r['B-cells'], r['Placenta'],
                    f"{obs_ctrl[b]:.4f}" if not np.isnan(obs_ctrl[b]) else 'NA',
                    f"{obs_pe[b]:.4f}" if not np.isnan(obs_pe[b]) else 'NA'])
log('saved results/GSE282512_60bin_obs.csv')
