## ============================================================
## 42. Phase 4 功能富集：候选 DMR 基因 GO/KEGG（WGBS 区域级）
## 输入: dmr_final_v2.csv (candidate==TRUE), region_annot.csv
## 设计: 基因集合 = 直接符号基因 + 无符号区域最近基因(flag)
##       背景   = region_annot 全部注释基因（limma 检验宇宙）
##       检验   = 超几何 + BH; GO 项限定 universe 内 10-500 基因
## 输出: results/GSE282512_dmr_enrichment_*.csv, summary.txt
##       figures/enrichment_bubble.png
## ============================================================
suppressPackageStartupMessages({
  library(data.table); library(org.Hs.eg.db); library(GO.db)
})
set.seed(1)

## ---------- 1. 数据 ----------
v2  <- fread("results/GSE282512_dmr_final_v2.csv")
cand <- v2[candidate == TRUE]
ann  <- fread("results/GSE282512_region_annot.csv")
ann[, region_id := as.character(region_id)]

## 基因坐标表（去重基因：chr/start/end 取最小区间）
genes <- ann[!is.na(symbol) & symbol != "",
             .(chr = chr[1], gstart = min(start), gend = max(end)),
             by = symbol]
cat(sprintf("背景注释基因: %d（检验区域 %d）\n", nrow(genes), nrow(ann)))

## ---------- 2. DMR -> 基因（直接符号优先，缺符号取最近基因） ----------
cand[, symbol_clean := fifelse(is.na(symbol) | symbol == "", "", as.character(symbol))]
direct <- cand[symbol_clean != "", .(region_id = as.character(region_id),
                                     direction, symbol = symbol_clean, how = "direct")]
miss <- cand[symbol_clean == "", .(region_id = as.character(region_id),
                                   direction, chr, start, end)]
nearest <- data.table()
if (nrow(miss)) {
  g2 <- copy(genes); setnames(g2, "chr", "gchr")
  j <- merge(miss, g2, by.x = "chr", by.y = "gchr", all.x = TRUE, allow.cartesian = TRUE)
  j[, mid := (start + end) / 2]
  j[, d := fifelse(mid < gstart, gstart - mid, fifelse(mid > gend, mid - gend, 0))]
  j <- setorder(j, region_id, d)[!duplicated(region_id)]
  nearest <- j[, .(region_id, direction, symbol, how = "nearest")]
}
gene_tbl <- unique(rbind(direct, nearest), by = c("region_id", "symbol"))
cat(sprintf("候选 DMR 基因映射: %d 区域 -> %d 基因（direct %d / nearest %d）\n",
            nrow(gene_tbl), uniqueN(gene_tbl$symbol),
            sum(gene_tbl$how == "direct"), sum(gene_tbl$how == "nearest")))

## ---------- 3. symbol -> Entrez ----------
x2e <- as.data.table(select(org.Hs.eg.db, keys = genes$symbol,
                            columns = "ENTREZID", keytype = "SYMBOL"))
setnames(x2e, c("SYMBOL", "ENTREZID"), c("symbol", "gene_id"))
x2e <- x2e[!is.na(gene_id)]
univ_eg <- unique(x2e$gene_id)
cand_eg <- unique(x2e[symbol %in% gene_tbl$symbol, .(symbol, gene_id)])
cat(sprintf("Entrez 映射: 背景 %d / 候选 %d\n", length(univ_eg), nrow(cand_eg)))

## ---------- 4. GO 注释 ----------
## 逐本体的 GO term -> universe 基因（GO2ALLEGS 列表接口，含祖先项）
suppressPackageStartupMessages(library(AnnotationDbi))
g2e <- as.list(org.Hs.egGO2ALLEGS)
onto <- as.data.table(select(GO.db, keys = names(g2e),
                             columns = c("ONTOLOGY", "TERM"), keytype = "GOID"))
go_sets <- list()
for (on in c("BP", "CC", "MF")) {
  ids <- onto[ONTOLOGY == on, GOID]
  gl <- g2e[ids]
  sizes <- vapply(gl, function(g) length(unique(intersect(g, univ_eg))), integer(1))
  keep <- sizes >= 10 & sizes <= 500
  st <- data.table(go_id = ids[keep], term = onto[match(ids[keep], GOID), TERM],
                   genes = lapply(gl[keep], function(g) unique(intersect(g, univ_eg))))
  go_sets[[on]] <- st
  cat(sprintf("GO %s: %d 项（10-500 基因）\n", on, nrow(st)))
}

hyper_enrich <- function(cand_genes, sets, label) {
  N <- length(univ_eg); k <- length(cand_genes)
  if (k < 5) return(data.table())
  rbindlist(lapply(sets, function(st) {
    rbindlist(lapply(seq_len(nrow(st)), function(i) {
      gi <- unlist(st$genes[[i]])
      n <- length(gi); x <- length(intersect(gi, cand_genes))
      if (x < 2) return(NULL)
      p <- phyper(x - 1, n, N - n, k, lower.tail = FALSE)
      data.table(set = label, go_id = st$go_id[i], term = st$term[i],
                 n_term = n, n_hit = x, n_cand = k, n_univ = N, p = p,
                 OR = (x * (N - n)) / ((k - x) * n + 1e-9))
    }))
  }))
}

go_res <- rbindlist(list(
  hyper_enrich(cand_eg$gene_id, go_sets, "all"),
  hyper_enrich(cand_eg[symbol %in% gene_tbl[direction == "hyper", symbol], gene_id], go_sets, "hyper"),
  hyper_enrich(cand_eg[symbol %in% gene_tbl[direction == "hypo",  symbol], gene_id], go_sets, "hypo")
))
if (nrow(go_res)) go_res[, p_bh := p.adjust(p, "BH"), by = set]
setorder(go_res, set, p)
fwrite(go_res, "results/GSE282512_dmr_enrichment_go.csv")

## ---------- 5. KEGG（clusterProfiler 可用则用；否则 REST 拉取通路） ----------
kegg_res <- data.table()
if (requireNamespace("clusterProfiler", quietly = TRUE)) {
  suppressPackageStartupMessages(library(clusterProfiler))
  for (lab in c("all", "hyper", "hypo")) {
    gs <- if (lab == "all") cand_eg$gene_id else
      cand_eg[symbol %in% gene_tbl[direction == lab, symbol], gene_id]
    if (length(gs) < 5) next
    ek <- tryCatch(enrichKEGG(gene = gs, universe = univ_eg, organism = "hsa",
                              pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500),
                   error = function(e) NULL)
    if (!is.null(ek) && nrow(as.data.frame(ek))) {
      d <- as.data.table(ek); d$set <- lab
      kegg_res <- rbindlist(list(kegg_res, d), fill = TRUE)
    }
  }
}
if (nrow(kegg_res)) fwrite(kegg_res, "results/GSE282512_dmr_enrichment_kegg.csv")

## ---------- 6. 汇总与图 ----------
sig <- function(d) if (nrow(d)) d[p_bh < 0.05 | p.adjust < 0.05] else d
top_go <- go_res[p < 0.01][order(p)][1:12]
sink("results/GSE282512_dmr_enrichment_summary.txt")
cat("DMR 基因功能富集（Phase 4）\n==========================\n")
cat(sprintf("候选 DMR: %d（hyper %d / hypo %d）\n", nrow(cand),
            sum(cand$direction == "hyper"), sum(cand$direction == "hypo")))
cat(sprintf("映射基因: %d（direct %d, nearest %d）; 背景基因: %d\n",
            nrow(cand_eg), sum(gene_tbl$how == "direct"),
            sum(gene_tbl$how == "nearest"), length(univ_eg)))
cat(sprintf("GO 项数: BP %d / CC %d / MF %d\n",
            nrow(go_sets$BP), nrow(go_sets$CC), nrow(go_sets$MF)))
cat(sprintf("GO 显著项 (BH<0.05): all %d / hyper %d / hypo %d\n",
            sum(go_res[set == "all"]$p_bh < 0.05),
            sum(go_res[set == "hyper"]$p_bh < 0.05),
            sum(go_res[set == "hypo"]$p_bh < 0.05)))
if (nrow(kegg_res)) cat(sprintf("KEGG 显著项 (BH<0.05): %d\n", sum(kegg_res$p.adjust < 0.05)))
cat("\nTop GO 条目 (p<0.01):\n")
if (nrow(top_go)) print(top_go[, .(set, go_id, term, n_hit, n_term, OR = round(OR, 2), p)]) else cat("（无 p<0.01 条目）\n")
sink()

if (nrow(top_go)) {
  png("figures/enrichment_bubble.png", 2400, 1800, res = 260)
  td <- top_go
  td$lab <- paste0(substr(td$term, 1, 44), " [", td$set, "]")
  par(mar = c(4, 22, 3, 1))
  plot(-log10(td$p), seq_len(nrow(td)), xlim = c(0, max(-log10(td$p)) * 1.15),
       yaxt = "n", xlab = "-log10(p)", ylab = "", pch = 19,
       cex = 0.6 + 2.2 * td$n_hit / max(td$n_hit),
       col = ifelse(td$set == "hyper", "#c0392b", ifelse(td$set == "hypo", "#2471a3", "#5d6d7e")),
       main = "GO enrichment of DMR genes (hyper/hypo/all)")
  axis(2, seq_len(nrow(td)), td$lab, las = 1, cex.axis = 0.72)
  dev.off()
}
cat("完成\n")
