# 58a_ld_compute.py — 1000G phase3 EUR 局部 LD 表 + GRCh37->GRCh38 位置回映
# 输入: data/ld_1000g/{chr3,chr6,chr17}.vcf.gz + panel.txt + results/_tmp_smr/mqtl_full7.csv
# 输出: results/_tmp_smr/ld_1000g_eur.csv.gz      (chr,pos1,pos2,r2,r; r2>=0.01)
#       results/_tmp_smr/ld_1000g_eur_positions.csv (chr,pos,ref,alt — 1000G 匹配位点)
#       results/_tmp_smr/lift_positions.csv       (chr,pos,pos38 — Ensembl REST 回映)
# 逻辑: 只保留 7 个 DMR CpG 的 cis meQTL (p<1e-5) 位点; 块状 Pearson 相关
# 注意: VCF 与 GoDMC 均为 GRCh37, 直接按位置匹配; 多等位位点跳过
import gzip
import json
import subprocess
import sys
import time
import numpy as np

ROOT = "."
VCF = {"3": "data/ld_1000g/chr3.vcf.gz",
       "6": "data/ld_1000g/chr6.vcf.gz",
       "17": "data/ld_1000g/chr17.vcf.gz"}
PANEL = "data/ld_1000g/panel.txt"
MQTL = "results/_tmp_smr/mqtl_full7.csv"
OUT = "results/_tmp_smr/ld_1000g_eur.csv.gz"
OUT_POS = "results/_tmp_smr/ld_1000g_eur_positions.csv"
OUT_LIFT = "results/_tmp_smr/lift_positions.csv"
P_THR = 1e-5
R2_MIN = 0.01


def curl_json(url):
    """curl 子进程取 JSON (httr/requests 在本环境不稳), 4 次重试"""
    for k in range(4):
        try:
            r = subprocess.run(['curl', '-s', '--ssl-no-revoke', '-m', '60',
                                '-H', 'User-Agent: Mozilla/5.0', url],
                               capture_output=True, text=True, timeout=70)
            if r.returncode == 0 and r.stdout.strip().startswith('{'):
                return json.loads(r.stdout)
        except Exception:
            pass
        time.sleep(2 * (k + 1))
    return None


def lift_positions(chr_, poss):
    """GRCh37 -> GRCh38 分段回映 (每段 <=900kb, Ensembl /map 大区段易失败), 返回 [(chr, pos, pos38)]"""
    pl = sorted(poss)
    out = []
    # 按 <=900kb 分段
    segs = []
    cur = [pl[0]]
    for p in pl[1:]:
        if p - cur[0] <= 900_000:
            cur.append(p)
        else:
            segs.append(cur); cur = [p]
    segs.append(cur)
    for si, seg in enumerate(segs):
        pmin, pmax = seg[0] - 100, seg[-1] + 100
        url = ('https://rest.ensembl.org/map/human/GRCh37/'
               f'{chr_}:{pmin}:{pmax}/GRCh38?content-type=application/json')
        js = curl_json(url)
        if not js or 'mappings' not in js:
            print(f"lift chr{chr_} seg{si} ({pmin}-{pmax}): FAILED")
            continue
        segmap = js['mappings']
        for p in seg:
            for m in segmap:
                o, mp = m['original'], m['mapped']
                if o['seq_region_name'] == chr_ and o['start'] <= p <= o['end']:
                    pos38 = (mp['start'] + (p - o['start']) if mp['strand'] == 1
                             else mp['end'] - (p - o['start']))
                    out.append((chr_, p, pos38))
                    break
        print(f"lift chr{chr_} seg{si} ({pmin}-{pmax}): {sum(1 for x in out if x[0]==chr_)} cumulative mapped")
    print(f"lift chr{chr_}: {len(out)}/{len(pl)} positions mapped")
    return out


def load_eur():
    eur = []
    with open(PANEL) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        si = hdr.index("sample"); pi = hdr.index("super_pop")
        for line in f:
            a = line.rstrip("\n").split("\t")
            if a[pi] == "EUR":
                eur.append(a[si])
    return eur


def load_positions():
    """返回 {chr: set(positions)} — cis 且 p<P_THR 的 meQTL 位点"""
    pos = {}
    with open(MQTL) as f:
        for line in f:
            a = line.rstrip("\n").split(",")
            if len(a) < 10:
                continue
            cpg, snp, beta, se, p, a1, a2, f1, cis, clumped = a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9]
            try:
                p = float(p)
            except ValueError:
                continue
            if cis != "TRUE" or p >= P_THR:
                continue
            # snp 如 "chr6:25983010:SNP"
            parts = snp.replace('"', "").split(":")
            if len(parts) < 2:
                continue
            c, ps = parts[0].replace("chr", ""), parts[1]
            try:
                ppos = int(ps)
            except ValueError:
                continue
            pos.setdefault(c, set()).add(ppos)
    return pos


def main():
    eur = load_eur()
    print(f"EUR samples: {len(eur)}")
    posmap = load_positions()
    for c, s in posmap.items():
        print(f"chr{c}: {len(s)} target positions")

    out_rows = []
    pos_rows = []
    import os
    for c, vcf in VCF.items():
        if c not in posmap:
            continue
        if not os.path.exists(vcf):
            print(f"chr{c}: VCF 尚未下载 ({vcf}), 跳过")
            continue
        targets = posmap[c]
        sites = []      # (pos, dosage array)
        with gzip.open(vcf, "rt") as f:
            for line in f:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    cols = line.rstrip("\n").split("\t")
                    eur_idx = [i for i, s in enumerate(cols) if s in set(eur)]
                    print(f"chr{c}: {len(eur_idx)} EUR columns found")
                    continue
                a = line.rstrip("\n").split("\t")
                # CHROM POS ID REF ALT QUAL FILTER INFO FORMAT ...
                try:
                    ppos = int(a[1])
                except ValueError:
                    continue
                if ppos not in targets:
                    continue
                alt = a[4]
                if "," in alt or alt == ".":
                    continue  # 多等位/缺失
                gts = [a[i] for i in eur_idx]
                if len(gts) != len(eur_idx):
                    print(f"chr{c} pos {ppos}: column mismatch, skip")
                    continue
                dos = np.empty(len(gts), dtype=np.float32)
                ok = True
                for k, g in enumerate(gts):
                    gt = g.split(":")[0]
                    if gt in ("0|0", "0/0"):
                        dos[k] = 0.0
                    elif gt in ("0|1", "1|0", "0/1", "1/0"):
                        dos[k] = 1.0
                    elif gt in ("1|1", "1/1"):
                        dos[k] = 2.0
                    elif gt in (".|.", "./.", "."):
                        ok = False
                        break
                    else:
                        ok = False
                        break
                if not ok:
                    continue
                sites.append((ppos, dos))
                pos_rows.append((c, ppos, a[3], alt))
        print(f"chr{c}: matched {len(sites)}/{len(targets)} positions in 1000G")
        if len(sites) < 2:
            continue
        n = len(sites)
        M = np.stack([s[1] for s in sites])             # n x 503
        # 块状相关
        B = 800
        for i0 in range(0, n, B):
            i1 = min(i0 + B, n)
            X = M[i0:i1]
            Xc = X - X.mean(axis=1, keepdims=True)
            sx = np.sqrt((Xc ** 2).sum(axis=1))
            for j0 in range(0, n, B):
                j1 = min(j0 + B, n)
                if j0 < i0:
                    continue
                Y = M[j0:j1]
                Yc = Y - Y.mean(axis=1, keepdims=True)
                sy = np.sqrt((Yc ** 2).sum(axis=1))
                den = np.outer(sx, sy)
                den[den == 0] = np.nan
                R = (Xc @ Yc.T) / den
                r2 = R ** 2
                ii, jj = np.where(r2 >= R2_MIN)
                for a_, b_ in zip(ii, jj):
                    gi, gj = i0 + a_, j0 + b_
                    if gi == gj:
                        continue
                    if sites[gi][0] > sites[gj][0]:
                        gi, gj = gj, gi
                    out_rows.append((c, sites[gi][0], sites[gj][0], float(r2[a_, b_]), float(R[a_, b_])))
        print(f"chr{c}: cumulative pairs {len(out_rows)}")

    # 去重 (对称块会产生重复对)
    seen = {}
    for c, p1, p2, r2v, rv in out_rows:
        seen[(c, p1, p2)] = (r2v, rv)
    with gzip.open(OUT, "wt", newline="") as f:
        f.write("chr,pos1,pos2,r2,r\n")
        for (c, p1, p2), (r2v, rv) in sorted(seen.items()):
            f.write(f"{c},{p1},{p2},{r2v:.4f},{rv:.4f}\n")
    with open(OUT_POS, "w", newline="") as f:
        f.write("chr,pos,ref,alt\n")
        for c, p, rf, al in sorted(pos_rows, key=lambda x: (x[0], int(x[1]))):
            f.write(f"{c},{p},{rf},{al}\n")
    print(f"DONE: {len(seen)} LD pairs -> {OUT}; {len(pos_rows)} matched positions -> {OUT_POS}")

    # ---- GRCh37 -> GRCh38 回映 (Ensembl REST, 供 58 号 R 脚本读取) ----
    lift_out = []
    for c, ps in posmap.items():
        lift_out.extend(lift_positions(c, ps))
    with open(OUT_LIFT, "w", newline="") as f:
        f.write("chr,pos,pos38\n")
        for c, p, p38 in lift_out:
            f.write(f"{c},{p},{p38}\n")
    print(f"LIFT DONE: {len(lift_out)} -> {OUT_LIFT}")


if __name__ == "__main__":
    main()
