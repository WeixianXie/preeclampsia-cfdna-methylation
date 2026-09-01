#!/usr/bin/env python3
# 50i_cross_cohort_final.py
# 交叉队列映射验证(最终版): 抽样 60 个 sun 标记 bin(seed 42),
# curl liftover hg19->hg38(修复 urllib SSL 卡死), 多进程聚合 GSE282512 全部 369 cov,
# 检验 GSE282512 Control/PE 观测 vs sun 血细胞参考相关
# 若 bin->标记映射正确: Control(正常妊娠母体血) 应与 Neutrophils/T-cells/B-cells 参考正相关
import gzip, os, csv, json, time, glob, random, sys, subprocess
import numpy as np
from collections import defaultdict
from multiprocessing import Pool
# scipy 延迟导入(仅主进程用; worker 加载 OpenBLAS 会耗尽内存)

ROOT = r"E:/妊娠期高血压甲基化研究方案/hdp-methylation-project"
DATA = os.path.join(ROOT, "data/geo_methylation")
CACHE = os.path.join(ROOT, 'results/GSE282512_60bin_liftover.json')
log = lambda *a: print(*a, flush=True)

# ---------- 全局(供 Pool 子进程使用) ----------
# 注意: Windows spawn 下 worker 重新导入模块, 全局变量不可靠
#       -> regions_by_chr 必须作为参数传入 aggregate_one

def aggregate_one(args):
    fp, rbc = args
    acc = defaultdict(lambda: [0, 0])
    with gzip.open(fp, 'rt', errors='replace') as fh:
        for line in fh:
            p = line.split('\t')
            if len(p) < 6:
                continue
            rl = rbc.get(p[0])
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

def liftover_curl(chrom, start, end, frm='GRCh37', to='GRCh38'):
    c = chrom[3:] if chrom.startswith('chr') else chrom
    url = f'https://grch37.rest.ensembl.org/map/human/{frm}/{c}:{start}..{end}/{to}?content-type=application/json'
    for k in range(3):
        try:
            r = subprocess.run(['curl', '-s', '--ssl-no-revoke', '-m', '30',
                                '-H', 'User-Agent: Mozilla/5.0', url],
                               capture_output=True, text=True, timeout=40)
            if r.returncode == 0 and r.stdout.strip():
                j = json.loads(r.stdout)
                ms = []
                for m in j.get('mappings', []):
                    if 'mapped' in m:
                        mm = m['mapped']
                        ch2 = mm['seq_region_name']
                        ms.append((('chr'+ch2) if not ch2.startswith('chr') else ch2,
                                   mm['start'], mm['end']))
                return ms
        except Exception:
            pass
        time.sleep(2 * (k + 1))
    return None

def main():
    from scipy.stats import spearmanr  # 仅主进程导入
    random.seed(42)

    # 1) 载入标记(v2 重建表 bin 1..5820)
    mk = {}
    with open(os.path.join(DATA, 'GSE154378/sun2015_markers.tsv'), encoding='utf-8') as f:
        rd = csv.DictReader(f, delimiter='\t')
        for r in rd:
            mk[int(r['bin'])] = r
    typeI_bins = [b for b, r in mk.items() if r['type'] == 'I']
    typeII_bins = [b for b, r in mk.items() if r['type'] == 'II']
    sample = random.sample(typeI_bins, 20) + random.sample(typeII_bins, 40)
    log(f'sampled bins: {len(sample)}')

    # 2) lift (curl, 带缓存)
    liftover_map = {}
    if os.path.exists(CACHE):
        with open(CACHE, encoding='utf-8') as f:
            liftover_map = {int(k): v for k, v in json.load(f).items()}
    todo = [b for b in sample if b not in liftover_map]
    log(f'liftover cache: {len(sample)-len(todo)}/{len(sample)} done, {len(todo)} todo')
    for b in todo:
        r = mk[b]
        res = liftover_curl(r['chr'], int(r['start']), int(r['end']))
        liftover_map[b] = res
        time.sleep(0.25)
        if res:
            with open(CACHE, 'w', encoding='utf-8') as f:
                json.dump({str(k): v for k, v in liftover_map.items()}, f)
    ok = {b: v for b, v in liftover_map.items() if v}
    log(f'liftover OK: {len(ok)}/{len(sample)}')
    if len(ok) < 30:
        log('FATAL: too few lifted regions'); sys.exit(1)

    regions = {}
    for b, maps in ok.items():
        best = max(maps, key=lambda x: x[2] - x[1])
        regions[b] = (best[0], best[1], best[2])
    regions_by_chr = defaultdict(list)
    for b, (ch, s, e) in regions.items():
        regions_by_chr[ch].append((b, s, e))
    log(f'regions by chr: {len(regions_by_chr)} chrs')
    rbc = dict(regions_by_chr)  # 传给 worker(避免 spawn 全局丢失)

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

    # 4) 多进程聚合 (结果缓存, 避免重复 14 分钟聚合)
    SUMCACHE = os.path.join(ROOT, 'results/GSE282512_60bin_sums.json')
    sum_ctrl = defaultdict(lambda: [0, 0])
    sum_pe = defaultdict(lambda: [0, 0])
    if os.path.exists(SUMCACHE):
        with open(SUMCACHE, encoding='utf-8') as f:
            cached = json.load(f)
        for b, v in cached['ctrl'].items():
            sum_ctrl[int(b)] = v
        for b, v in cached['pe'].items():
            sum_pe[int(b)] = v
        log(f'aggregation loaded from cache: {len(sum_ctrl)} ctrl bins, {len(sum_pe)} pe bins')
    else:
        t0 = time.time()
        with Pool(4) as pool:
            results = pool.map(aggregate_one, [(fp, rbc) for fp in (ctrl + pe)])
        log(f'aggregation done in {time.time()-t0:.0f}s')
        ctrl_n = len(ctrl)
        for i, r in enumerate(results):
            tgt = sum_ctrl if i < ctrl_n else sum_pe
            for b, (m, c) in r.items():
                tgt[b][0] += m; tgt[b][1] += c
        with open(SUMCACHE, 'w', encoding='utf-8') as f:
            json.dump({'ctrl': {str(k): v for k, v in sum_ctrl.items()},
                       'pe': {str(k): v for k, v in sum_pe.items()}}, f)
        log('aggregation cached')

    obs_ctrl = {b: (v[0]/v[1] if v[1] > 0 else np.nan) for b, v in sum_ctrl.items()}
    obs_pe = {b: (v[0]/v[1] if v[1] > 0 else np.nan) for b, v in sum_pe.items()}

    out_rows = []
    def test(label, xfun, yfun, bins):
        xs = [xfun(b) for b in bins]; ys = [yfun(b) for b in bins]
        okk = [(float(x), y) for x, y in zip(xs, ys)
               if x is not None and x != '' and not np.isnan(y)]
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
        test(f'Control_vs_{tc}', lambda b: mk[b][tc], lambda b: obs_ctrl.get(b, np.nan), bins_all)
    log('\n== GSE282512 PE obs vs 参考 ==')
    for tc in refs:
        test(f'PE_vs_{tc}', lambda b: mk[b][tc], lambda b: obs_pe.get(b, np.nan), bins_all)

    log('\n== 数据质量 ==')
    test('Control_vs_PE_obs', lambda b: obs_ctrl.get(b, np.nan), lambda b: obs_pe.get(b, np.nan), bins_all)

    # GSE154378 obs vs GSE282512 Control
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
        xs = [obs154[b] for b in b154]; ys = [obs_ctrl.get(b, np.nan) for b in b154]
        okk = [(x, y) for x, y in zip(xs, ys) if not np.isnan(y)]
        rho, p = spearmanr([x[0] for x in okk], [x[1] for x in okk])
        log(f'GSE154378_vs_GSE282512_Control: rho={rho:+.3f} p={p:.2e} n={len(okk)}')
        out_rows.append(['GSE154378_vs_GSE282512_Control', f'{rho:.4f}', f'{p:.2e}', len(okk)])

    with open(os.path.join(ROOT, 'results/GSE282512_marker_mapping_validation.csv'), 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['test', 'rho', 'p', 'n'])
        w.writerows(out_rows)
    log('\nsaved results/GSE282512_marker_mapping_validation.csv')

    with open(os.path.join(ROOT, 'results/GSE282512_60bin_obs.csv'), 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['bin', 'chr_hg38', 'start_hg38', 'end_hg38', 'type',
                    'Neutrophils', 'T-cells', 'B-cells', 'Placenta',
                    'obs_Control', 'obs_PE'])
        for b in bins_all:
            r = mk[b]; ch, s, e = regions[b]
            oc = obs_ctrl.get(b, np.nan); op = obs_pe.get(b, np.nan)
            w.writerow([b, ch, s, e, r['type'],
                        r['Neutrophils'], r['T-cells'], r['B-cells'], r['Placenta'],
                        f"{oc:.4f}" if not np.isnan(oc) else 'NA',
                        f"{op:.4f}" if not np.isnan(op) else 'NA'])
    log('saved results/GSE282512_60bin_obs.csv')

if __name__ == '__main__':
    os.chdir(ROOT)
    main()
