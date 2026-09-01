#!/bin/bash
# 58e_fix_bad_chunks.sh — 自动修复坏块循环: 下载 -> 重组 -> 复扫 -> 直到 ALL OK
# 用法: bash scripts/58e_fix_bad_chunks.sh <chr编号> <坏块逗号列表>
set -o pipefail
CHR="$1"; BADLIST="$2"
CS=8388608
URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
OUT="data/ld_1000g/chr${CHR}.vcf.gz"
PARTS="${OUT}.parts"
PY="C:/Users/xieweixian/.workbuddy/binaries/python/versions/3.13.12/python.exe"

LEN=$(curl -s --ssl-no-revoke -I "$URL" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}' | tail -1)
NCHUNK=$(( (LEN + CS - 1) / CS ))
echo "[fix] chr${CHR} LEN=$LEN NCHUNK=$NCHUNK bad=[$BADLIST]"

dl_one() {
  local i=$1
  local s=$((i * CS))
  local e=$((s + CS - 1))
  [ $e -ge $LEN ] && e=$((LEN - 1))
  local exp=$((e - s + 1))
  local pf="${PARTS}/part_$(printf '%05d' $i)"
  for t in 1 2 3 4 5 6; do
    curl -s --ssl-no-revoke -m 600 -r ${s}-${e} -o "$pf" "$URL"
    local sz
    sz=$(stat -c%s "$pf" 2>/dev/null || echo 0)
    if [ "$sz" -eq "$exp" ]; then echo "[ok] part_$i (try $t)"; return 0; fi
    echo "[retry] part_$i: got $sz / $exp (try $t)"; sleep 3
  done
  echo "[FAIL] part_$i"; return 1
}

for ROUND in 1 2 3; do
  if [ -z "$BADLIST" ]; then echo "[done] 无坏块"; break; fi
  echo "=== ROUND $ROUND: 修复 [$BADLIST] ==="
  FAIL=0
  for i in ${BADLIST//,/ }; do
    while [ "$(jobs -rp | wc -l)" -ge 4 ]; do sleep 1; done
    dl_one "$i" || FAIL=1 &
  done
  wait
  [ "$FAIL" = "1" ] && echo "[warn] 本轮有分块失败"
  cat "$PARTS"/part_* > "$OUT"
  ACT=$(stat -c%s "$OUT")
  echo "[asm] $ACT / $LEN"
  if [ "$ACT" -ne "$LEN" ]; then echo "[err] 尺寸不符, 终止"; exit 1; fi
  SCAN=$("$PY" scripts/58d_scan_all_bad_chunks.py "$OUT" $CS | tail -1)
  echo "[scan] $SCAN"
  if [[ "$SCAN" == ALL*OK* ]]; then echo "[SUCCESS] chr${CHR} 修复完成"; exit 0; fi
  BADLIST="${SCAN#BAD_CHUNKS }"
done
echo "[err] 3 轮后仍有坏块: $BADLIST"
exit 1
