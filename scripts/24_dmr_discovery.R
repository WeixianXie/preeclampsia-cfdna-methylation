# 24_dmr_discovery.R ----------------------------------------------------------
# 区域级 DMR 发现: GSE282512 子队列 PE vs Control
# 输入:
#   results/GSE282512_region_beta.csv.gz   区域x样本 beta 宽矩阵
#   results/GSE282512_subcohort.csv        子队列分组 (group/tube_type/batch)
#   results/GSE282512_region_annot.csv     区域注释
# 方法:
#   1) 区域过滤: 在 >=80% 样本有覆盖的区域, 其余 NA 用中位数填充
#   2) limma 主分析: beta ~ group + tube_type + batch (eBayes)
#   3) Wilcoxon 稳健对照: 每区域 Mann-Whitney U + BH 校正
#   4) 候选: FDR<0.05 且 |delta_beta|>=0.10
# 输出:
#   results/GSE282512_dmr_candidates.csv   全区域统计表(含候选标志)
#   results/GSE282512_dmr_candidates_top.txt 候选汇总
#   figures/dmr_volcano.png  火山图
#   figures/dmr_manhattan.png 曼哈顿图(按染色体 -log10 FDR)
# 注意: 相对路径 + GBK 编码 (Windows R 约定)
# -----------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  library(data.table); library(limma); library(parallel)
}))

MIN_COVER  <- 0.80   # 区域在多少比例样本有覆盖才保留
FDR_CUT    <- 0.05
DELTA_CUT  <- 0.10
N_WORKER   <- 6

## 1. 读取数据 ----------------------------------------------------------------
wide <- fread("results/GSE282512_region_beta.csv.gz")
ann  <- fread("results/GSE282512_region_annot.csv")
sub  <- fread("results/GSE282512_subcohort.csv")

rids <- wide$region_id
gsms <- names(wide)[-1]
mat  <- as.matrix(wide[, -1, with = FALSE])
rownames(mat) <- rids
mat[mat == ""] <- NA
mat <- apply(mat, 2, as.numeric)
rownames(mat) <- rids
cat(sprintf("矩阵: %d 区域 x %d 样本\n", nrow(mat), ncol(mat)))

## 2. 样本分组信息 -------------------------------------------------------------
sub <- sub[gsm %in% gsms]
setkey(sub, gsm)
sub[, tube_type2 := ifelse(tube_type == "PAXgene DNA", "PAXgene", "EDTA")]
sub[, batch := as.factor(batch)]
stopifnot(nrow(sub) == ncol(mat))

## 3. 区域过滤 (覆盖度) --------------------------------------------------------
cov_ratio <- rowMeans(!is.na(mat))
keep <- cov_ratio >= MIN_COVER
cat(sprintf("覆盖>=%.0f%%样本的区域: %d / %d\n", 100 * MIN_COVER, sum(keep), length(keep)))
mat <- mat[keep, , drop = FALSE]

## NA 填充 (区域中位数)
na_fill <- function(x) { x[is.na(x)] <- median(x, na.rm = TRUE); x }
mat <- t(apply(mat, 1, na_fill))

## 4. limma 主分析 ------------------------------------------------------------
group <- factor(sub$group, levels = c("Control", "PE"))
design <- model.matrix(~ group + tube_type2 + batch, data = sub)
fit <- lmFit(mat, design)
fit <- eBayes(fit)
tt <- topTable(fit, coef = "groupPE", number = Inf, sort.by = "none")
tt <- data.table(region_id = as.integer(rownames(mat)),
                 p_limma = tt$P.Value, fdr_limma = tt$adj.P.Val)

## 5. Wilcoxon 稳健对照 --------------------------------------------------------
idx_pe   <- which(group == "PE")
idx_ctrl <- which(group == "Control")
wilcox_one <- function(i) {
  x <- mat[i, idx_pe]; y <- mat[i, idx_ctrl]
  tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)
}
cl <- makeCluster(N_WORKER, type = "PSOCK")
clusterEvalQ(cl, suppressWarnings(suppressMessages(library(stats))))
clusterExport(cl, c("mat", "idx_pe", "idx_ctrl"))
p_w <- parLapply(cl, seq_len(nrow(mat)), wilcox_one)
stopCluster(cl)
tt[, p_wilcox := unlist(p_w)]
tt[, fdr_wilcox := p.adjust(p_wilcox, method = "BH")]

## 6. 汇总与候选 ---------------------------------------------------------------
ann_k <- ann[region_id %in% tt$region_id]
res <- merge(tt, ann_k, by = "region_id", sort = FALSE)
beta_pe   <- rowMeans(mat[, idx_pe, drop = FALSE])
beta_ctrl <- rowMeans(mat[, idx_ctrl, drop = FALSE])
res[, beta_pe := beta_pe][, beta_ctrl := beta_ctrl][, delta_beta := beta_pe - beta_ctrl]
res[, direction := ifelse(delta_beta > 0, "hyper", "hypo")]
res[, candidate := (fdr_limma < FDR_CUT & abs(delta_beta) >= DELTA_CUT)]
setorder(res, fdr_limma)

fwrite(res, "results/GSE282512_dmr_candidates.csv")
cat(sprintf("候选 DMR 区域: %d (hyper %d, hypo %d)\n",
            sum(res$candidate),
            sum(res$candidate & res$direction == "hyper"),
            sum(res$candidate & res$direction == "hypo")))

## 7. 候选汇总文本 -------------------------------------------------------------
sink("results/GSE282512_dmr_candidates_top.txt", split = TRUE)
cat("===== GSE282512 区域级 DMR 发现 (PE vs Control) =====\n\n")
cat(sprintf("矩阵: %d 区域 x %d 样本 (PE %d / Control %d)\n", nrow(mat), ncol(mat),
            length(idx_pe), length(idx_ctrl)))
cat(sprintf("区域过滤: 覆盖>=%.0f%% 保留 %d\n\n", 100 * MIN_COVER, sum(keep)))
cat(sprintf("候选阈值: FDR<%.2f 且 |delta_beta|>=%.2f\n", FDR_CUT, DELTA_CUT))
cat(sprintf("候选数: %d (hyper %d / hypo %d)\n",
            sum(res$candidate), sum(res$candidate & res$direction == "hyper"),
            sum(res$candidate & res$direction == "hypo")))
cat("\n-- 候选区域类型分布 --\n")
print(res[candidate == TRUE, .N, by = type][order(-N)])
cat("\n-- Top 15 候选 --\n")
print(res[candidate == TRUE][order(fdr_limma)][1:min(15, sum(res$candidate)),
  .(region_id, chr, start, end, type, symbol, beta_pe, beta_ctrl, delta_beta,
    p_limma, fdr_limma, fdr_wilcox, direction)])
cat("\n-- Wilcoxon 一致率 (候选内 FDR_wilcox<0.05) --\n")
cat(sprintf("%d / %d\n", sum(res$candidate & res$fdr_wilcox < FDR_CUT),
            sum(res$candidate)))
sink()

## 8. 图 ----------------------------------------------------------------------
png("figures/dmr_volcano.png", width = 1600, height = 1200, res = 150)
plot(res$delta_beta, -log10(res$fdr_limma), pch = 20, cex = 0.5,
     col = ifelse(res$candidate, ifelse(res$direction == "hyper", "firebrick", "steelblue"), "grey70"),
     xlab = "delta beta (PE - Control)", ylab = "-log10 FDR (limma)",
     main = "GSE282512 region-level DMR volcano (PE vs Control)")
abline(v = c(-DELTA_CUT, DELTA_CUT), lty = 2, col = "grey40")
abline(h = -log10(FDR_CUT), lty = 2, col = "grey40")
legend("topright", legend = c("candidate hyper", "candidate hypo"),
       col = c("firebrick", "steelblue"), pch = 20, bty = "n")
dev.off()

## 曼哈顿图: 染色体位置
ord <- order(factor(res$chr, levels = paste0("chr", c(1:22, "X", "Y"))), res$start)
res_o <- res[ord]
chr_pos <- split(seq_len(nrow(res_o)), res_o$chr)
cum_start <- c(0, cumsum(vapply(chr_pos, length, integer(1)))[-length(chr_pos)])
pos <- unlist(lapply(seq_along(chr_pos), function(i) cum_start[i] + seq_along(chr_pos[[i]])),
              use.names = FALSE)
png("figures/dmr_manhattan.png", width = 2400, height = 1200, res = 150)
cols <- rep(c("grey40", "grey70"), ceiling(length(chr_pos) / 2))[seq_along(chr_pos)]
colv <- unlist(lapply(seq_along(chr_pos), function(i) rep(cols[i], length(chr_pos[[i]]))))
plot(pos, -log10(res_o$fdr_limma), pch = 20, cex = 0.5, col = colv,
     xaxt = "n", xlab = "Chromosome", ylab = "-log10 FDR",
     main = "GSE282512 region-level DMR Manhattan")
mid <- cum_start + vapply(chr_pos, function(x) length(x) / 2, numeric(1))
axis(1, at = mid, labels = names(chr_pos), cex.axis = 0.6, las = 2)
abline(h = -log10(FDR_CUT), lty = 2, col = "red")
dev.off()

cat("DONE\n")
