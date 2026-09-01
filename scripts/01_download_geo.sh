#!/bin/bash
# =====================================================================
# 01_download_geo.sh
# Phase 0 数据获取：GEO series matrix 批量下载（14 个数据集）
# 用法: bash 01_download_geo.sh [代理地址(可选)]
#   例: bash 01_download_geo.sh http://127.0.0.1:7890
# 说明: 断点续传 + gzip 完整性校验；GSE282512 若无 matrix 会转查 supplementary
# =====================================================================

PROXY="${1:-}"
if [ -n "$PROXY" ]; then
  export http_proxy="$PROXY" https_proxy="$PROXY"
  echo "[INFO] 使用代理: $PROXY"
fi

BASE="/e/妊娠期高血压甲基化研究方案/hdp-methylation-project/data"
METH_DIR="$BASE/geo_methylation"
TX_DIR="$BASE/geo_transcriptome"
mkdir -p "$METH_DIR" "$TX_DIR"
LOG="$BASE/../logs/download_geo.log"
mkdir -p "$(dirname "$LOG")"

# gse|目录|用途
METH_SETS="GSE282512|cfDNA发现层
GSE37722|全血重现层27K
GSE73375|胎盘机制层
GSE75196|胎盘机制层
GSE57767|胎盘机制层"

TX_SETS="GSE86200|早孕外周血验证
GSE85307|早孕外周血验证
GSE48424|免疫浸润主数据
GSE10588|胎盘DEG合并
GSE25906|胎盘DEG合并
GSE44711|胎盘DEG合并
GSE35574|胎盘DEG合并
GSE98224|胎盘DEG合并
GSE54618|胎盘DEG合并"

geo_url() {
  local gse="$1"
  local num=${gse#GSE}
  local prefix=$(echo "$num" | sed 's/[0-9]\{3\}$/nnn/')
  echo "https://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}/${gse}/matrix/${gse}_series_matrix.txt.gz"
}

download_one() {
  local gse="$1" dir="$2" label="$3"
  local url=$(geo_url "$gse")
  local out="$dir/${gse}_series_matrix.txt.gz"
  if [ -f "$out" ] && gzip -t "$out" 2>/dev/null; then
    echo "[SKIP] $gse ($label) 已存在且完整"
    return 0
  fi
  echo "[GET ] $gse ($label) ..."
  curl -sL --connect-timeout 20 --retry 4 --retry-delay 3 -C - \
       -o "$out" "$url" 2>>"$LOG"
  if [ -f "$out" ] && gzip -t "$out" 2>/dev/null; then
    local sz=$(du -h "$out" | cut -f1)
    echo "[OK  ] $gse 完成 ($sz)"
    return 0
  else
    echo "[FAIL] $gse 下载失败或文件损坏，清理后待重试"
    rm -f "$out"
    # 探测是否存在（区分 404 与网络失败）
    local code=$(curl -s -o /dev/null -w "%{http_code}" --connect-timeout 15 -I "$url")
    if [ "$code" = "404" ]; then
      echo "      -> $gse 无 series matrix（测序型数据集），需从 supplementary 获取:"
      curl -s --connect-timeout 15 "https://ftp.ncbi.nlm.nih.gov/geo/series/${gse#GSE}" >/dev/null 2>&1
      echo "         https://ftp.ncbi.nlm.nih.gov/geo/series/$(echo ${gse#GSE} | sed 's/[0-9]\{3\}$/nnn/')/${gse}/suppl/"
    fi
    return 1
  fi
}

echo "==== GEO 批量下载 $(date '+%F %T') ===="
fail_list=""

echo "--- 甲基化数据集 ---"
echo "$METH_SETS" | while IFS='|' read -r gse label; do
  [ -z "$gse" ] && continue
  download_one "$gse" "$METH_DIR" "$label" || echo "$gse" >> "$BASE/../logs/failed.txt"
done

echo "--- 转录组数据集 ---"
echo "$TX_SETS" | while IFS='|' read -r gse label; do
  [ -z "$gse" ] && continue
  download_one "$gse" "$TX_DIR" "$label" || echo "$gse" >> "$BASE/../logs/failed.txt"
done

echo ""
echo "==== 下载汇总 ===="
echo "甲基化: $(ls "$METH_DIR"/*.gz 2>/dev/null | wc -l) / 5 个文件"
echo "转录组: $(ls "$TX_DIR"/*.gz 2>/dev/null | wc -l) / 9 个文件"
if [ -s "$BASE/../logs/failed.txt" ]; then
  echo "[注意] 失败数据集: $(sort -u "$BASE/../logs/failed.txt" | tr '\n' ' ')"
  echo "       重新运行本脚本即可续传（已完成的自动跳过）"
else
  echo "全部成功！"
fi
