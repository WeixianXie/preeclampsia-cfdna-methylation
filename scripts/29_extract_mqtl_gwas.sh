#!/bin/bash
# 29_extract_mqtl_gwas.sh -------------------------------------------------------
# 流式抽取大文件 (bash + awk, 避免全部载入内存):
#   1) GoDMC assoc_meta_all.csv.gz (5.9GB) -> cis-meQTL, pval<1e-8
#      输出列: cpg,snp_chr,snp_pos,pval,beta_a1,se
#   2) FinnGen R10 PE / 妊娠期高血压 GWAS -> p<1e-4 的变异
#      输出列: chrom,pos,rsids,nearest_genes,pval,beta,af_alt
# 约定: 相对路径, 在项目根运行
# ------------------------------------------------------------------------------
set -o pipefail
cd "$(dirname "$0")/.."

echo "[godmc] start $(date)"
zcat data/mqtl/assoc_meta_all.csv.gz | awk -F',' '
NR>1 {
  if ($11 == "TRUE" && $5+0 < 1e-8) {
    split($2, a, ":")
    c = $1; gsub(/"/, "", c)
    printf "%s,%s,%s,%s,%s,%s\n", c, a[1], a[2], $5, $3, $4
  }
}' | gzip -1 > data/mqtl/godmc_cis.csv.gz
echo "[godmc] done: $(zcat data/mqtl/godmc_cis.csv.gz | wc -l) rows, $(date)"

echo "[finngen PE] start $(date)"
zcat data/gwas/finngen_R10_O15_PRE_OR_ECLAMPSIA.gz | awk -F'\t' 'NR>1 && $7+0 < 1e-4 { printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $5, $6, $7, $8, $10 }' > data/gwas/finngen_pe_p1e4.tsv
echo "[finngen PE] done: $(wc -l < data/gwas/finngen_pe_p1e4.tsv) rows, $(date)"

echo "[finngen GH] start $(date)"
zcat data/gwas/finngen_R10_O15_GESTAT_HYPERT.gz | awk -F'\t' 'NR>1 && $7+0 < 1e-4 { printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $5, $6, $7, $8, $10 }' > data/gwas/finngen_gh_p1e4.tsv
echo "[finngen GH] done: $(wc -l < data/gwas/finngen_gh_p1e4.tsv) rows, $(date)"

echo "ALL DONE"
