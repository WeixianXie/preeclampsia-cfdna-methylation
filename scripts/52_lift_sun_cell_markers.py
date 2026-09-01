# -*- coding: utf-8 -*-
# 52_lift_sun_cell_markers.py — Sun 2015 细胞/胎盘标记 hg19->hg38 回映 (任务 #50 下半: GSE37722 互证准备)
# 背景: GSE37722 = HM27 母体白细胞 (Zhou manifest 为 hg38); Sun v4 标记表为 hg19
#       必须把标记区域回映到 hg38 才能与 HM27 探针重叠
# 选择: Type.I cell markers (Neutrophils/T-cells/B-cells) + Type.II placenta HIGH
#       = 白细胞构成假说互证所需的全部标记
# 实现: grch37.rest.ensembl.org GET /map/human/GRCh37/{region}/GRCh38
#       (POST 不支持该端点); 并发 6 线程, 每线程起 curl --ssl-no-revoke, 失败重试 3 次
# 输出: data/geo_methylation/GSE154378/sun2015_markers_hg38_lift.tsv
import csv, json, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = r'E:/妊娠期高血压甲基化研究方案/hdp-methylation-project'
SRC = ROOT + '/data/geo_methylation/GSE154378/sun2015_markers.tsv'
OUT = ROOT + '/data/geo_methylation/GSE154378/sun2015_markers_hg38_lift.tsv'
TISSUE_COLS = ['Liver', 'Lungs', 'Colon', 'Small intestines', 'Pancreas',
               'Adrenal glands', 'Esophagus', 'Adipose tissues', 'Heart', 'Brain',
               'T-cells', 'B-cells', 'Neutrophils', 'Placenta']


def lift_region(chrom, start, end):
    """hg19 -> hg38, 返回 (chr38, start38, end38) 或 None"""
    c = chrom[3:] if chrom.startswith('chr') else chrom
    url = ('https://grch37.rest.ensembl.org/map/human/GRCh37/'
           f'{c}:{start}..{end}/GRCh38?content-type=application/json')
    for k in range(3):
        try:
            r = subprocess.run(['curl', '-s', '--ssl-no-revoke', '-m', '30',
                                '-H', 'User-Agent: Mozilla/5.0', url],
                               capture_output=True, text=True, timeout=40)
            if r.returncode == 0 and r.stdout.strip():
                j = json.loads(r.stdout)
                ms = [m['mapped'] for m in j.get('mappings', []) if 'mapped' in m]
                if ms:
                    mm = ms[0]
                    ch2 = ('chr' + mm['seq_region_name']) if not mm['seq_region_name'].startswith('chr') else mm['seq_region_name']
                    return ch2, mm['start'], mm['end']
                return None
        except Exception:
            pass
        time.sleep(1.5 * (k + 1))
    return None


def main():
    import statistics
    OTHER = [t for t in TISSUE_COLS if t != 'Placenta']
    rows = []
    with open(SRC, encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            cell_ok = r['type'] == 'I' and r['tissue'] in ('Neutrophils', 'T-cells', 'B-cells')
            pla_dir = ''
            if r['type'] == 'II' and r['placenta_informative'] == 'Y':
                med = statistics.median(float(r[t]) for t in OTHER)
                diff = float(r['Placenta']) - med
                pla_dir = 'HIGH' if diff > 15 else ('LOW' if -diff > 15 else 'FLAT')
            pla_ok = pla_dir == 'HIGH'
            if cell_ok or pla_ok:
                r['pla_dir'] = pla_dir
                rows.append(r)
    print(f'待回映标记: {len(rows)} '
          f"(Neutrophils {sum(1 for r in rows if r['tissue']=='Neutrophils')} / "
          f"T-cells {sum(1 for r in rows if r['tissue']=='T-cells')} / "
          f"B-cells {sum(1 for r in rows if r['tissue']=='B-cells')} / "
          f"Placenta-HIGH {sum(1 for r in rows if r['placenta_informative']=='HIGH')})")

    def job(r):
        res = lift_region(r['chr'], int(r['start']), int(r['end']))
        return r, res

    out_rows, failed = [], []
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = [ex.submit(job, r) for r in rows]
        for i, fut in enumerate(as_completed(futs), 1):
            r, res = fut.result()
            if res is None:
                failed.append(r['bin'])
                continue
            ch38, s38, e38 = res
            rr = dict(bin=r['bin'], chr=r['chr'], start=r['start'], end=r['end'],
                      type=r['type'], tissue=r['tissue'], pla_dir=r['pla_dir'],
                      chr38=ch38, start38=s38, end38=e38)
            # 参考谱数值 (供下游 sign 对齐)
            for t in TISSUE_COLS:
                rr[t] = r[t]
            out_rows.append(rr)
            if i % 100 == 0:
                print(f'  进度 {i}/{len(rows)}', flush=True)

    # 失败重试一轮 (串行)
    for b in list(failed):
        r = next(x for x in rows if x['bin'] == b)
        res = lift_region(r['chr'], int(r['start']), int(r['end']))
        if res:
            failed.remove(b)
            ch38, s38, e38 = res
            rr = dict(bin=r['bin'], chr=r['chr'], start=r['start'], end=r['end'],
                      type=r['type'], tissue=r['tissue'], pla_dir=r['pla_dir'],
                      chr38=ch38, start38=s38, end38=e38)
            for t in TISSUE_COLS:
                rr[t] = r[t]
            out_rows.append(rr)

    out_rows.sort(key=lambda x: int(x['bin']))
    cols = ['bin', 'chr', 'start', 'end', 'type', 'tissue', 'pla_dir',
            'chr38', 'start38', 'end38'] + TISSUE_COLS
    with open(OUT, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(out_rows)
    print(f'完成: {len(out_rows)} 成功 / {len(failed)} 失败 -> {OUT}')
    if failed:
        print('失败 bin:', failed[:20])
        sys.exit(1)


if __name__ == '__main__':
    main()
