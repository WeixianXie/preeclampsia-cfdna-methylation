# 47_wgcna.R — cfDNA 区域甲基化 WGCNA 共模块网络 (任务 #45, 第三梯队)
# 输入: results/GSE282512_region_beta.csv.gz (99,495 区域 x 64 亚队列样本, hg38)
#       results/GSE282512_subcohort.csv (group), results/GSE282512_deconv_fractions.csv (7 型)
#       results/GSE282512_dmr_final_v2.csv (166 候选 DMR)
# 设计: top 5000 高变异区域 -> WGCNA 模块 -> 模块-性状 (PE / 白细胞比例) 关联
#       + 166 DMR 在模块中的富集
# 输出: results/GSE282512_wgcna_module_trait.csv / _dmr_module.csv / _summary.txt
#       figures/wgcna_module_trait.png
# 约定: UTF-8, 相对路径, LANG=en_US.UTF-8
suppressPackageStartupMessages({
  library(data.table); library(WGCNA); library(ggplot2)
})
options(stringsAsFactors = FALSE)
disableWGCNAThreads()   # PSOCK worker 闭包在本机不稳 (见项目记忆), 单线程 5000x64 规模可接受
set.seed(42)

## ============ 1. 数据 ============
beta <- fread("results/GSE282512_region_beta.csv.gz")
beta[, region_id := as.character(region_id)]
mat <- as.matrix(beta[, -1]); rownames(mat) <- beta$region_id
sub <- fread("results/GSE282512_subcohort.csv")
sub <- sub[match(names(beta)[-1], gsm)]
trait <- data.frame(PE = as.integer(sub$group == "PE"), row.names = sub$gsm)
## 白细胞比例 (12 亚型敏感性参考)
frac <- fread("results/GSE282512_deconv_fractions_12ct.csv")
frac_cols <- setdiff(names(frac), c("gsm","group"))
frac_mat <- as.matrix(frac[, frac_cols, with = FALSE]); rownames(frac_mat) <- frac$gsm
frac_mat <- frac_mat[rownames(trait), , drop = FALSE]
frac_mat <- frac_mat[, apply(frac_mat, 2, function(x) sd(x, na.rm = TRUE) > 1e-8), drop = FALSE]
cat(sprintf("矩阵 %d x %d; PE %d\n", nrow(mat), ncol(mat), sum(trait$PE)))

## ============ 2. top 5000 高变异区域 (缺失 <20% 再行中位数插补) ============
nafrac <- apply(mat, 1, function(r) mean(is.na(r)))
v <- apply(mat, 1, var, na.rm = TRUE); v[is.na(v)] <- 0
cand <- names(which(nafrac <= 0.2))
top <- cand[order(v[cand], decreasing = TRUE)][1:5000]
X <- mat[top, , drop = FALSE]
X <- t(X)  # 样本 x 区域
## 行(区域)中位数插补残余 NA
if (any(is.na(X))) {
  med <- apply(X, 2, median, na.rm = TRUE)
  idx <- which(is.na(X), arr.ind = TRUE)
  X[idx] <- med[idx[, "col"]]
}
good <- goodSamplesGenes(X, verbose = 0)
X <- X[good$goodSamples, good$goodGenes]
cat(sprintf("WGCNA 输入: %d 样本 x %d 区域 (缺失区域剔除后)\n", nrow(X), ncol(X)))

## ============ 3. 网络 ============
sft <- pickSoftThreshold(X, powerVector = 1:12, networkType = "signed", verbose = 0)
pw <- sft$powerEstimate; if (is.na(pw)) pw <- 6
r2 <- sft$fitIndices$SFT.R.sq[match(pw, sft$fitIndices$Power)]
cat(sprintf("soft power = %d (R^2=%.2f)\n", pw, r2))
net <- blockwiseModules(X, power = pw, networkType = "signed", TOMType = "signed",
                        minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        maxBlockSize = 6000, verbose = 0)
colors <- labels2colors(net$colors)
cat(sprintf("模块: %d 个 (大小: %s)\n", length(unique(colors)),
            paste(sort(table(colors), decreasing = TRUE)[1:8], collapse = ", ")))

## ============ 4. 模块-性状关联 ============
MEs <- orderMEs(moduleEigengenes(X, colors)$eigengenes)
allTrait <- cbind(trait, frac_mat)
cors <- suppressWarnings(cor(MEs, allTrait, use = "pairwise.complete.obs"))
ps <- matrix(NA_real_, nrow(cors), ncol(cors), dimnames = dimnames(cors))
for (i in seq_len(nrow(cors))) for (j in seq_len(ncol(cors))) {
  ct <- suppressWarnings(tryCatch(cor.test(MEs[, i], allTrait[, j]), error = function(e) NULL))
  if (!is.null(ct)) ps[i, j] <- ct$p.value
}
mt <- data.table(module = rep(rownames(cors), ncol(cors)),
                 trait = rep(colnames(cors), each = nrow(cors)),
                 r = as.vector(cors), p = as.vector(ps))
setorder(mt, p)
fwrite(mt, "results/GSE282512_wgcna_module_trait.csv")

## ============ 5. DMR 模块归属与富集 ============
dmr <- fread("results/GSE282512_dmr_final_v2.csv")
dmr <- dmr[!is.na(region_id)][, region_id := as.character(region_id)]
dmr_in <- dmr[region_id %in% colnames(X)]
mod_of <- data.table(region_id = colnames(X), module = colors)
dmr_mod <- merge(dmr_in[, .(region_id, direction, final_call)], mod_of, by = "region_id")
fwrite(dmr_mod, "results/GSE282512_wgcna_dmr_module.csv")
## 富集: 各模块内 DMR 比例 vs 全体
tab <- mod_of[, .(n = .N, n_dmr = sum(region_id %in% dmr_in$region_id)), by = module]
tot_dmr <- nrow(dmr_in); tot_reg <- nrow(mod_of)
tab[, p_fisher := mapply(function(a, n)
  fisher.test(matrix(c(a, tot_dmr - a, n - a, (tot_reg - n) - (tot_dmr - a)), 2, 2),
              alternative = "greater")$p.value, n_dmr, n)]
tab[, p_fisher := fifelse(is.finite(p_fisher), p_fisher, 1)]
setorder(tab, -n_dmr)
fwrite(tab, "results/GSE282512_wgcna_dmr_enrichment.csv")

## 图: 模块-性状热图 (r 与 p)
pl <- mt[!grepl("grey", module)]
p1 <- ggplot(pl, aes(trait, module, fill = r)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f\n%.1g", r, p)), size = 2.6) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(x = NULL, y = NULL, fill = "r", title = "WGCNA module-trait relationships") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("figures/wgcna_module_trait.png", p1, width = 8, height = 6, dpi = 300)

## ============ 6. summary ============
sink("results/GSE282512_wgcna_summary.txt")
cat("== WGCNA: cfDNA 区域甲基化共模块网络 ==\n\n")
cat(sprintf("输入: top %d 高变异区域 x %d 样本; soft power %d (signed)\n", ncol(X), nrow(X), pw))
cat(sprintf("模块数: %d (含 grey); 模块大小 top: %s\n\n", length(unique(colors)),
            paste(names(sort(table(colors), decreasing = TRUE)[1:8]),
                  sort(table(colors), decreasing = TRUE)[1:8], sep = "=", collapse = ", ")))
cat("---- 模块-性状 (前 12 行) ----\n")
print(head(mt, 12))
cat("\n---- DMR 模块富集 (前 6) ----\n")
print(head(tab[, .(module, n, n_dmr, p_fisher = signif(p_fisher, 2))], 6))
cat(sprintf("\nDMR 归属: %d/%d 候选 DMR 进入 WGCNA 集合\n", nrow(dmr_mod), nrow(dmr)))
cat(sprintf("hyper/hypo 分布: %s\n", paste(names(table(dmr_mod$direction)), table(dmr_mod$direction), sep = "=", collapse = ", ")))
cat("\n注: 模块-白细胞比例关联可指示哪些共甲基化模块对应细胞构成; DMR 富集模块指示 PE 信号的组织来源\n")
sink()
cat("[完成] WGCNA 输出: 3 csv + summary + 1 图\n")
