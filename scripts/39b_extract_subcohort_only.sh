# 39b_extract_subcohort_only.sh — 仅子队列 64 样本抽取标记窗口 (跳过已存在)
set -u
RAW=data/geo_methylation/GSE282512_raw
WIN=results/_tmp_mwin.tsv
OUTDIR=results/_tmp_cellmark
mkdir -p "$OUTDIR"

while read gsm; do
  f=$(ls "$RAW"/${gsm}_DNA*.cov.gz 2>/dev/null | head -1)
  [ -z "$f" ] && { echo "MISS $gsm"; continue; }
  out="$OUTDIR/${gsm}.tsv"
  [ -s "$out" ] && { echo "skip $gsm"; continue; }
  gzip -dc "$f" | awk -F'\t' -v win="$WIN" -v outf="$out" '
    BEGIN {
      h = 0; nc = 0
      while ((getline l < win) > 0) {
        if (h++ == 0) continue
        split(l, a, "\t")
        nc++
        wch[nc] = a[1]; wst[nc] = a[2] + 0; wen[nc] = a[3] + 0; wid[nc] = a[4]
        if (!(a[1] in fst)) fst[a[1]] = nc   # 每条染色体窗口块首索引
      }
      cur = ""
    }
    {
      chr = $1; pos = $2 + 0
      if (chr != cur) { cur = chr; ptr = (chr in fst) ? fst[chr] : nc + 1 }
      while (ptr <= nc && wch[ptr] == chr && wst[ptr] <= pos) ptr++
      for (j = ptr - 1; j >= 1 && wch[j] == chr && wst[j] > pos - 501; j--) {
        if (pos >= wst[j] && pos <= wen[j]) { m[wid[j]] += $5; u[wid[j]] += $6 }
      }
    }
    END { for (k in m) print k "\t" m[k] "\t" u[k] > outf }
  '
  echo "done $gsm $( [ -s "$out" ] && wc -l < "$out" || echo 0 )"
done < results/_tmp_subgsms.txt
echo ALL_DONE
