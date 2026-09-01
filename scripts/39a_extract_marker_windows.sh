# 39a_extract_marker_windows.sh — 每样本提取 333 标记 ±250bp 窗口的 m/u 聚合
# 用法: bash scripts/39a_extract_marker_windows.sh
# 输入: results/_tmp_mwin.tsv (chr ws we cpg, 已按 chr,ws 排序)
#        data/geo_methylation/GSE282512_raw/GSM*.cov.gz (染色体内位置有序)
# 输出: results/_tmp_cellmark/GSMxxx.tsv (cpg \t m \t u)

set -u
RAW=data/geo_methylation/GSE282512_raw
WIN=results/_tmp_mwin.tsv
OUTDIR=results/_tmp_cellmark
mkdir -p "$OUTDIR"

for f in "$RAW"/GSM*.cov.gz; do
  gsm=$(basename "$f" | sed 's/_DNA.*$//')
  out="$OUTDIR/${gsm}.tsv"
  [ -s "$out" ] && continue
  gzip -dc "$f" | awk -F'\t' -v win="$WIN" -v outf="$out" '
    BEGIN {
      h = 0; nc = 0
      while ((getline l < win) > 0) {
        if (h++ == 0) continue
        split(l, a, "\t")
        nc++
        wch[nc] = a[1]; wst[nc] = a[2] + 0; wen[nc] = a[3] + 0; wid[nc] = a[4]
      }
      cur = ""
    }
    {
      chr = $1; pos = $2 + 0
      if (chr != cur) { cur = chr; ptr = 1 }
      while (ptr <= nc && wch[ptr] == chr && wst[ptr] <= pos) ptr++
      # 回看起点在 (pos-501, pos] 的窗口 (窗口最长 501bp)
      for (j = ptr - 1; j >= 1 && wch[j] == chr && wst[j] > pos - 501; j--) {
        if (pos >= wst[j] && pos <= wen[j]) { m[wid[j]] += $5; u[wid[j]] += $6 }
      }
    }
    END { for (k in m) print k "\t" m[k] "\t" u[k] > outf }
  '
  echo "done $gsm $( [ -s "$out" ] && wc -l < "$out" || echo 0 )"
done
echo ALL_DONE
