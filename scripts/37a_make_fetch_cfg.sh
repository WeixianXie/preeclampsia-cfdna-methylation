#!/bin/bash
# 37a: 生成 Ensembl REST 序列抓取配置并用 curl --parallel 批量下载
cd "$(dirname "$0")/.." || exit 1
mkdir -p results/_tmp_seq
awk -F'\t' 'NR>1 {
  c=$2; sub(/^chr/,"",c);
  printf "url = \"https://rest.ensembl.org/sequence/region/human/%s:%s:%s:1\"\n", c, $3, $4;
  printf "output = \"results/_tmp_seq/cand_%s.txt\"\n", $1;
}' results/_tmp_cand_regions.tsv > results/_tmp_fetch.cfg

awk -F'\t' 'NR>1 {
  c=$2; sub(/^chr/,"",c);
  printf "url = \"https://rest.ensembl.org/sequence/region/human/%s:%s:%s:1\"\n", c, $3, $4;
  printf "output = \"results/_tmp_seq/bg_%s.txt\"\n", $1;
}' results/_tmp_bg_regions.tsv >> results/_tmp_fetch.cfg

echo "cfg lines: $(wc -l < results/_tmp_fetch.cfg)"

curl -s --ssl-no-revoke --max-time 60 -H "Accept: text/plain" \
  --config results/_tmp_fetch.cfg --parallel --parallel-max 3 \
  --retry 2 --retry-delay 1
echo "batch1 done"

# 第二遍: 顺序补抓失败/非序列输出 (429 或空文件)
bad=0
while IFS=$'\t' read -r rid chr st en tag; do
  f="results/_tmp_seq/${tag}_${rid}.txt"
  if [ ! -s "$f" ] || ! head -c 1 "$f" | grep -q '[ACGT]'; then
    c="$chr"; c="${c#chr}"
    for try in 1 2 3; do
      curl -s --ssl-no-revoke --max-time 30 -H "Accept: text/plain" \
        "https://rest.ensembl.org/sequence/region/human/${c}:${st}:${en}:1" > "$f"
      if head -c 1 "$f" | grep -q '[ACGT]'; then break; fi
      sleep 1
    done
    bad=$((bad+1))
  fi
done < <(awk -F'\t' 'NR>1 {tag="cand"; print $1"\t"$2"\t"$3"\t"$4"\t"tag}' results/_tmp_cand_regions.tsv;
         awk -F'\t' 'NR>1 {tag="bg";   print $1"\t"$2"\t"$3"\t"$4"\t"tag}' results/_tmp_bg_regions.tsv)
echo "retried: $bad"
n_ok=$(ls results/_tmp_seq | wc -l)
echo "files: $n_ok"
