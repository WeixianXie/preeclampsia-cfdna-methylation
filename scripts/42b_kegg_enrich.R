## ============================================================
## 42b. KEGG 通路富集（REST 批量接口，不依赖 clusterProfiler）
## 输入: dmr_final_v2.csv / region_annot.csv + KEGG REST link/list
## 输出: results/GSE282512_dmr_enrichment_kegg.csv, 汇总追加
## ============================================================
suppressPackageStartupMessages(library(data.table))
setwd("E:/妊娠期高血压甲基化研究方案/hdp-methylation-project")
curl <- "curl -s --ssl-no-revoke"
if (!file.exists("results/_tmp_kegg_link.tsv"))
  system(paste(curl, "-o results/_tmp_kegg_link.tsv https://rest.kegg.jp/link/hsa/pathway"))
if (!file.exists("results/_tmp_kegg_list.tsv"))
  system(paste(curl, "-o results/_tmp_kegg_list.tsv https://rest.kegg.jp/list/pathway/hsa"))

lk <- fread("results/_tmp_kegg_link.tsv", header = FALSE,
            col.names = c("pathway", "gene"), sep = "\t")
nm <- fread("results/_tmp_kegg_list.tsv", header = FALSE,
            col.names = c("pathway", "name"), sep = "\t")
lk[, gene := sub("^hsa:", "", gene)]
lk[, pathway := sub("^path:", "", pathway)]
lk <- lk[nzchar(gene)]
nm[, pathway := sub("^path:", "", pathway)]
cat(sprintf("KEGG 映射: %d 通路, %d 对（%d 基因）\n",
            uniqueN(lk$pathway), nrow(lk), uniqueN(lk$gene)))

## 与 42 相同的基因映射逻辑（direct + nearest）
v2 <- fread("results/GSE282512_dmr_final_v2.csv"); cand <- v2[candidate == TRUE]
ann <- fread("results/GSE282512_region_annot.csv")
ann[, region_id := as.character(region_id)]
genes <- ann[!is.na(symbol) & symbol != "",
             .(chr = chr[1], gstart = min(start), gend = max(end)), by = symbol]
cand[, symbol_clean := fifelse(is.na(symbol) | symbol == "", "", as.character(symbol))]
direct <- cand[symbol_clean != "", .(region_id = as.character(region_id),
                                     direction, symbol = symbol_clean, how = "direct")]
miss <- cand[symbol_clean == "", .(region_id = as.character(region_id),
                                   direction, chr, start, end)]
g2 <- copy(genes); setnames(g2, "chr", "gchr")
j <- merge(miss, g2, by.x = "chr", by.y = "gchr", all.x = TRUE, allow.cartesian = TRUE)
j[, mid := (start + end) / 2]
j[, d := fifelse(mid < gstart, gstart - mid, fifelse(mid > gend, mid - gend, 0))]
j <- setorder(j, region_id, d)[!duplicated(region_id)]
gene_tbl <- unique(rbind(direct, j[, .(region_id, direction, symbol, how = "nearest")]),
                   by = c("region_id", "symbol"))

## symbol -> Entrez（org.Hs.eg.db）
suppressPackageStartupMessages({library(org.Hs.eg.db); library(AnnotationDbi)})
x2e <- as.data.table(select(org.Hs.eg.db, keys = genes$symbol,
                            columns = "ENTREZID", keytype = "SYMBOL"))
setnames(x2e, c("SYMBOL", "ENTREZID"), c("symbol", "gene_id"))
x2e <- x2e[!is.na(gene_id)]
univ_eg <- unique(x2e$gene_id)
cand_eg <- unique(x2e[symbol %in% gene_tbl$symbol, .(symbol, gene_id)])

## 通路集合（限定 universe 内 10-500 基因）
sets <- lk[, .(genes = list(unique(intersect(gene, univ_eg)))), by = pathway]
sets <- sets[sapply(genes, length) >= 10 & sapply(genes, length) <= 500]
sets[, name := nm$name[match(pathway, nm$pathway)]]
cat(sprintf("限定通路数: %d（10-500 基因）\n", nrow(sets)))

hyper_kegg <- function(cg, label) {
  N <- length(univ_eg); k <- length(cg)
  if (k < 5) return(data.table())
  rbindlist(lapply(seq_len(nrow(sets)), function(i) {
    gi <- unlist(sets$genes[[i]]); n <- length(gi); x <- length(intersect(gi, cg))
    if (x < 2) return(NULL)
    data.table(set = label, pathway = sets$pathway[i], name = sets$name[i],
               n_term = n, n_hit = x, n_cand = k, n_univ = N,
               p = phyper(x - 1, n, N - n, k, lower.tail = FALSE),
               OR = (x * (N - n)) / ((k - x) * n + 1e-9))
  }))
}
kegg <- rbindlist(list(
  hyper_kegg(cand_eg$gene_id, "all"),
  hyper_kegg(cand_eg[symbol %in% gene_tbl[direction == "hyper", symbol], gene_id], "hyper"),
  hyper_kegg(cand_eg[symbol %in% gene_tbl[direction == "hypo", symbol], gene_id], "hypo")))
if (nrow(kegg)) kegg[, p_bh := p.adjust(p, "BH"), by = set]
setorder(kegg, set, p)
fwrite(kegg, "results/GSE282512_dmr_enrichment_kegg.csv")

cat("\n## KEGG 富集追加（42b）\n", file = "results/GSE282512_dmr_enrichment_summary.txt", append = TRUE)
cat(sprintf("KEGG 通路显著项 (BH<0.05): %d；Top 名义条目:\n",
            sum(kegg$p_bh < 0.05)),
    file = "results/GSE282512_dmr_enrichment_summary.txt", append = TRUE)
topk <- kegg[p < 0.05][order(p)][1:8]
print(topk[, .(set, pathway, name, n_hit, n_term, OR = round(OR, 2), p = signif(p, 3))])
cat(capture.output(print(topk[, .(set, pathway, name, n_hit, n_term,
                                  OR = round(OR, 2), p = signif(p, 3))])),
    file = "results/GSE282512_dmr_enrichment_summary.txt", append = TRUE, sep = "\n")
cat("完成\n")
