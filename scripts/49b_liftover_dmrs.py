#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
49b_liftover_dmrs.py  (v2)
将 GSE282512 的 166 候选 DMR (hg38) 通过 Ensembl REST 映射到 hg19 (GRCh37)
GET /map/human/GRCh38/{chr}:{start}..{end}/GRCh37  （区域级映射，每 DMR 一次调用）

v2 修复：
- 主站点 rest.ensembl.org 的 map 端点当前返回 500 → 改用 grch37.rest.ensembl.org（GRCh37 专用 REST，200）
- 响应字段取 m['mapped']['seq_region_name'/'start'/'end']（v1 误取 m['seq_region_name'] 导致 8 个"OK"也为空）
输出: results/GSE282512_dmr_hg19.csv
"""
import csv, json, os, sys, time, urllib.request, urllib.error

ROOT = r"E:/妊娠期高血压甲基化研究方案/hdp-methylation-project"
os.chdir(ROOT)
BASE = 'https://grch37.rest.ensembl.org'

def fetch(url, tries=5):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            if e.code == 429:
                time.sleep(8 * (k + 1)); continue
            if e.code == 400:
                return None  # 区域本身无法映射
            time.sleep(4 * (k + 1))
        except Exception as e:
            time.sleep(5 * (k + 1))
    return None

rows = list(csv.DictReader(open('results/GSE282512_dmr_final_v2.csv', encoding='utf-8')))
cand = [r for r in rows if str(r.get('candidate', '')).strip().upper() == 'TRUE']
print(f'candidates: {len(cand)}')

out = []
fail = []
for i, r in enumerate(cand):
    chrom, start, end = r['chr'], int(r['start']), int(r['end'])
    chrom_no = chrom[3:] if chrom.startswith('chr') else chrom  # grch37 端点不接受 chr 前缀
    url = f'{BASE}/map/human/GRCh38/{chrom_no}:{start}..{end}/GRCh37?content-type=application/json'
    j = fetch(url)
    if j is None or 'mappings' not in j:
        fail.append((r['region_id'], chrom, start, end))
        out.append([r['region_id'], chrom, start, end, 'FAIL', '', '', ''])
        print(f'  [{i+1}/{len(cand)}] FAIL {r["region_id"]}')
        continue
    maps = j['mappings']
    n = 0
    for m in maps:
        if 'mapped' not in m:
            continue
        mm = m['mapped']
        ch = mm.get('seq_region_name', '')
        if ch and not ch.startswith('chr'):
            ch = 'chr' + ch
        out.append([r['region_id'], chrom, start, end, 'OK', ch,
                    mm.get('start', ''), mm.get('end', '')])
        n += 1
    if n == 0:
        fail.append((r['region_id'], chrom, start, end))
        out.append([r['region_id'], chrom, start, end, 'FAIL', '', '', ''])
        print(f'  [{i+1}/{len(cand)}] FAIL(no mapped) {r["region_id"]}')
        continue
    if (i + 1) % 20 == 0:
        print(f'  [{i+1}/{len(cand)}] ok, accumulated {len(out)} segments')
    time.sleep(0.5)

with open('results/GSE282512_dmr_hg19.csv', 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['region_id','chr_hg38','start_hg38','end_hg38','status','chr_hg19','start_hg19','end_hg19'])
    w.writerows(out)

nok = sum(1 for o in out if o[4] == 'OK')
print(f'done: {len(out)} segments ({nok} OK, {len(fail)} failed DMRs)')
for f in fail[:10]:
    print('  FAIL:', f)
