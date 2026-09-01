#!/bin/bash
# =====================================================================
# 02_download_gwas.sh
# Phase 0 数据获取：FinnGen R10 GWAS 汇总统计 + GoDMC mQTL 说明
# 用法: bash 02_download_gwas.sh [代理地址(可选)]
#   例: bash 02_download_gwas.sh http://127.0.0.1:7890
# 注意: FinnGen 文件在 Google Storage，国内需代理；GoDMC 在 EBI FTP
# =====================================================================

PROXY="${1:-}"
if [ -n "$PROXY" ]; then
  export http_proxy="$PROXY" https_proxy="$PROXY"
  echo "[INFO] 使用代理: $PROXY"
fi

BASE="/e/妊娠期高血压甲基化研究方案/hdp-methylation-project/data"
GWAS_DIR="$BASE/gwas"
MQTL_DIR="$BASE/mqtl"
mkdir -p "$GWAS_DIR" "$MQTL_DIR"
LOG="$BASE/../logs/download_gwas.log"

# FinnGen R10 端点（O15 = 妊娠相关）
# 主结局：O15_PRE_OR_ECLAMPSIA（子痫前期，~7,377 例）
# 次结局：O15_GESTAT_HYPERT（妊娠期高血压，~16,417 例）
# 备选：O15_PREEC_OR_FETGRO（PE+胎儿生长受限）、O15_EXIST_HYPERT_COMPLIC
FINNGEN_URL_TPL="https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_%s.gz"

download_finngen() {
  local endpoint="$1"
  local url=$(printf "$FINNGEN_URL_TPL" "$endpoint")
  local out="$GWAS_DIR/finngen_R10_${endpoint}.gz"
  if [ -f "$out" ] && gzip -t "$out" 2>/dev/null; then
    echo "[SKIP] $endpoint 已存在且完整"
    return 0
  fi
  echo "[GET ] FinnGen R10 $endpoint ..."
  # 首先验证 URL（FinnGen 版本目录名可能微调，失败时打印候选入口）
  local code=$(curl -s -o /dev/null -w "%{http_code}" --connect-timeout 15 -I "$url")
  if [ "$code" = "404" ]; then
    echo "      URL 变更，请到 https://r10.finngen.fi/ 页面确认汇总统计下载链接"
    return 1
  fi
  curl -sL --connect-timeout 20 --retry 4 --retry-delay 3 -C - -o "$out" "$url" 2>>"$LOG"
  if [ -f "$out" ] && gzip -t "$out" 2>/dev/null; then
    echo "[OK  ] $endpoint 完成 ($(du -h "$out" | cut -f1))"
  else
    echo "[FAIL] $endpoint 下载失败"; rm -f "$out"; return 1
  fi
}

download_finngen "O15_PRE_OR_ECLAMPSIA"
download_finngen "O15_GESTAT_HYPERT"

# ---------------------------------------------------------------------
# GoDMC mQTL（Phase 10 MR 用，此处仅做入口核验与说明）
# 全量 cis-mQTL BESD 约 15-20 GB，Phase 0 不全量下载；
# Phase 10 时按候选 CpG 清单提取对应染色体的 BESD 分区即可。
# 入口: https://www.godmc.org.uk/ -> Downloads
# 文件结构（EBI FTP）:
#   ftp.ebi.ac.uk/pub/databases/mQTL/  (cis 文件按染色体分卷 BESD)
# 使用方式: SMR 软件 + --besd-file 指定，或 TwoSampleMR 格式化
# ---------------------------------------------------------------------
echo ""
echo "[INFO] GoDMC mQTL 全量文件较大（15-20GB），Phase 0 不下载。"
echo "       入口: https://www.godmc.org.uk/  (Downloads -> cis BESD by chr)"
echo "       Phase 10 按候选 CpG 所在染色体选择性下载到: $MQTL_DIR"
echo ""
echo "==== GWAS 下载汇总 ===="
echo "FinnGen: $(ls "$GWAS_DIR"/*.gz 2>/dev/null | wc -l) / 2 个端点文件"
