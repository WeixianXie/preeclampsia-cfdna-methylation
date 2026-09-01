# 38_leukocyte_deconv.R — 白细胞亚群 cfDNA 去卷积（任务 #37 / 方案 A）
# 思路: EpiDISH 血液参考标记 CpG (HM450) → 映射到本项目区域注释 → 取区域级 beta 作代理
#      (cfDNA WGBS 位点级覆盖稀疏, 直接位点匹配仅 ~5%/样本; 区域级矩阵完整)
# 主参考: centDHSbloodDMC (333 tsDHS-DMC × 7 亚型, 跨平台设计)
# 敏感性: cent12CT450k (600 CpG × 12 亚型) + CBS 方法
# 关键坐标事实: cov 文件位置 = Zhou manifest CpG_beg + 1 (区域级不受影响, 此处仅用于记录)

suppressPackageStartupMessages({
  library(data.table); library(EpiDISH); library(ggplot2)
})
set.seed(42)

## ============ 1. 数据准备 ============
ra <- fread("results/GSE282512_region_annot.csv", select = c("region_id","chr","start","end"))
beta <- fread("results/GSE282512_region_beta.csv.gz")
beta_mat <- as.matrix(beta[, -1]); rownames(beta_mat) <- beta$region_id
sub <- fread("results/GSE282512_subcohort.csv")
stopifnot(all(names(beta)[-1] %in% sub$gsm))
gsm_order <- names(beta)[-1]
sub <- sub[match(gsm_order, gsm)]
grp <- sub$group
cat(sprintf("区域矩阵: %d x %d; PE %d / CT %d\n", nrow(beta_mat), ncol(beta_mat),
            sum(grp=="PE"), sum(grp=="Control")))

man <- fread("data/annot/HM450.hg38.manifest.tsv.gz", select = c("Probe_ID","CpG_chrm","CpG_beg"))

## 标记 CpG -> 区域 (多点区域取最短, 减轻稀释)
map_markers <- function(cg_ids) {
  mm <- man[Probe_ID %in% cg_ids]
  mm[, cend := CpG_beg]
  ra2 <- copy(ra); ra2[, qend := end]; ra2[, width := end - start]
  setkey(ra2, chr, start, qend)
  fo <- foverlaps(mm[, .(CpG_chrm, CpG_beg, cend, Probe_ID)], ra2,
                  by.x = c("CpG_chrm","CpG_beg","cend"),
                  by.y = c("chr","start","qend"), nomatch = NULL)
  fo <- fo[fo[, .I[width == min(width)], by = Probe_ID]$V1]
  fo <- fo[order(Probe_ID, width, region_id)]
  fo <- fo[!duplicated(fo$Probe_ID)]  # 并列最短区域只取一个, 避免标记重复
  ## %in% 会隐式 int->char 强转, 但矩阵 [整数, ] 按位置索引 -> 必须显式转字符
  fo[, region_id := as.character(region_id)]
}

## 矩阵行名是字符型 region_id; 所有用于索引的 region_id 统一转字符
beta <- beta[, region_id := as.character(region_id)]

## ============ 2. 主参考 centDHSbloodDMC (7 亚型) ============
data("centDHSbloodDMC.m")
map7 <- map_markers(rownames(centDHSbloodDMC.m))
cat(sprintf("[主参考] 标记→区域: %d/%d, 落入矩阵: %d\n",
            uniqueN(map7$Probe_ID), nrow(centDHSbloodDMC.m),
            sum(map7$region_id %in% rownames(beta_mat))))
map7 <- map7[region_id %in% rownames(beta_mat)]
ref7 <- centDHSbloodDMC.m[map7$Probe_ID, , drop = FALSE]
b7  <- beta_mat[map7$region_id, , drop = FALSE]
rownames(b7) <- map7$Probe_ID
rownames(ref7) <- map7$Probe_ID

## 区域无覆盖产生的 NA (13.7%) 用标记中位数插补 (RPC/CBS 统一口径)
impute_row_median <- function(m) {
  if (!any(is.na(m))) return(m)
  med <- apply(m, 1, median, na.rm = TRUE)
  idx <- which(is.na(m), arr.ind = TRUE)
  m[idx] <- med[idx[, "row"]]
  m
}
b7 <- impute_row_median(b7)

rpc7 <- epidish(b7, ref7, method = "RPC")
cbs7 <- tryCatch(epidish(b7, ref7, method = "CBS"), error = function(e) NULL)
frac7 <- data.frame(gsm = gsm_order, group = grp, rpc7$estF, check.names = FALSE)
frac7_cbs <- if (!is.null(cbs7)) data.frame(gsm = gsm_order, group = grp, cbs7$estF, check.names = FALSE) else NULL
cat("[主参考 RPC] 各亚型均值:\n"); print(round(colMeans(rpc7$estF), 3))

## ============ 3. 敏感性: cent12CT450k (12 亚型) ============
data("cent12CT450k.m")
map12 <- map_markers(rownames(cent12CT450k.m))
map12 <- map12[region_id %in% rownames(beta_mat)]
cat(sprintf("[12亚型] 标记可用: %d/%d\n", uniqueN(map12$Probe_ID), nrow(cent12CT450k.m)))
ref12 <- cent12CT450k.m[map12$Probe_ID, , drop = FALSE]
b12 <- beta_mat[map12$region_id, , drop = FALSE]
b12 <- impute_row_median(b12)
rownames(b12) <- map12$Probe_ID; rownames(ref12) <- map12$Probe_ID
rpc12 <- epidish(b12, ref12, method = "RPC")
frac12 <- data.frame(gsm = gsm_order, group = grp, rpc12$estF, check.names = FALSE)

## ============ 4. PE vs Control 检验 ============
mk_tests <- function(fr, tag) {
  ct_cols <- setdiff(names(fr), c("gsm","group"))
  rbindlist(lapply(ct_cols, function(cc) {
    x <- fr[[cc]]; pe <- x[grp == "PE"]; ctr <- x[grp == "Control"]
    wt <- suppressWarnings(wilcox.test(pe, ctr))
    ## 调整 tube_type + batch 的线性模型
    d <- data.frame(y = x, pe = as.integer(grp == "PE"),
                    tube = factor(sub$tube_type), batch = factor(sub$batch))
    lm1 <- summary(lm(y ~ pe + tube + batch, d))
    p_adj <- coef(lm1)["pe","Pr(>|t|)"]; b_adj <- coef(lm1)["pe","Estimate"]
    data.table(ref = tag, cell = cc,
               med_pe = median(pe), med_ct = median(ctr),
               diff = median(pe) - median(ctr),
               p_wilcox = wt$p.value, beta_adj = b_adj, p_adj = p_adj)
  }))
}
tests <- rbindlist(list(mk_tests(frac7, "DHS7_RPC"),
                        if (!is.null(frac7_cbs)) mk_tests(frac7_cbs, "DHS7_CBS") else NULL,
                        mk_tests(frac12, "12CT_RPC")))
tests[, p_bh := p.adjust(p_wilcox, "BH"), by = ref]
tests <- tests[order(ref, p_wilcox)]

## ============ 5. 细胞比例与 DMR 甲基化联动 ============
v2 <- fread("results/GSE282512_dmr_final_v2.csv", select = c("region_id","direction"))
v2[, region_id := as.character(region_id)]
rid <- intersect(v2$region_id, rownames(beta_mat))
v2m <- v2[region_id %in% rid]
hyper_mean <- colMeans(beta_mat[v2m[direction=="hyper"]$region_id, , drop = FALSE], na.rm = TRUE)
hypo_mean  <- colMeans(beta_mat[v2m[direction=="hypo"]$region_id, , drop = FALSE], na.rm = TRUE)
dmr_link <- rbindlist(lapply(intersect(names(frac7), c("B","NK","CD4T","CD8T","Mono","Neutro","Eosino")), function(cc) {
  rbind(
    data.table(cell = cc, dmr = "hyper_DMR", rho =
      suppressWarnings(cor.test(frac7[[cc]], hyper_mean, method="spearman"))$estimate,
      p = suppressWarnings(cor.test(frac7[[cc]], hyper_mean, method="spearman"))$p.value),
    data.table(cell = cc, dmr = "hypo_DMR", rho =
      suppressWarnings(cor.test(frac7[[cc]], hypo_mean, method="spearman"))$estimate,
      p = suppressWarnings(cor.test(frac7[[cc]], hypo_mean, method="spearman"))$p.value))
}))
dmr_link[, p_bh := p.adjust(p, "BH")]

## ============ 6. 输出 ============
fwrite(frac7, "results/GSE282512_deconv_fractions.csv")
fwrite(frac12, "results/GSE282512_deconv_fractions_12ct.csv")
fwrite(tests, "results/GSE282512_deconv_tests.csv")
fwrite(dmr_link, "results/GSE282512_deconv_dmr_link.csv")

sink("results/GSE282512_deconv_summary.txt")
cat("===== 白细胞亚群 cfDNA 去卷积 (EpiDISH, 区域级 beta 投影) =====\n\n")
cat(sprintf("样本: PE %d / Control %d; 方法: RPC + CBS 敏感性\n\n", sum(grp=="PE"), sum(grp=="Control")))
cat("---- PE vs Control 各亚型 (主参考 7 亚型 RPC) ----\n")
print(tests[ref == "DHS7_RPC"])
cat("\n---- 敏感性: CBS 方法 ----\n")
if (!is.null(frac7_cbs)) print(tests[ref == "DHS7_CBS"]) else cat("CBS 失败, 已跳过\n")
cat("\n---- 敏感性: 12 亚型参考 ----\n")
print(tests[ref == "12CT_RPC"])
cat("\n---- 细胞比例 x DMR 甲基化 Spearman ----\n")
print(dmr_link)
sink()

## 图: 7 亚型组间箱线 + Neutro vs hyper DMR 散点
lg <- melt(as.data.table(frac7), id.vars = c("gsm","group"),
           variable.name = "cell", value.name = "frac")
p1 <- ggplot(lg, aes(group, frac, fill = group)) +
  geom_boxplot(outlier.size = 0.6, width = 0.6) +
  facet_wrap(~cell, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c(Control = "steelblue", PE = "firebrick")) +
  labs(x = NULL, y = "估计细胞比例") + theme_classic(base_size = 10)
pd <- data.frame(neutro = frac7$Neutro, hyper = hyper_mean, group = grp)
p2 <- ggplot(pd, aes(neutro, hyper, color = group)) +
  geom_point(size = 1.8) + scale_color_manual(values = c(Control="steelblue", PE="firebrick")) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.5) +
  labs(x = "Neutrophil 比例 (RPC)", y = "hyper DMR 平均 β", color = NULL) +
  theme_classic(base_size = 10)
ggsave("figures/deconv_celltypes.png", p1, width = 8, height = 5, dpi = 300)
ggsave("figures/deconv_neutro_dmr.png", p2, width = 5, height = 4, dpi = 300)
cat("\n完成: results/GSE282512_deconv_*.csv/txt + figures/deconv_*.png\n")
