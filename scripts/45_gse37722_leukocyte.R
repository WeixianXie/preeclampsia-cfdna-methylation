# 45_gse37722_leukocyte.R — 白细胞层重现验证 (任务 #43, 第二梯队)
# 思路: GSE37722 = 妊娠各期母体白细胞 HM27 甲基化 (nulligravid/early/middle/delivery/postpartum, 84 样本)
#       若 PE cfDNA DMR 定位于"正常妊娠中随孕期变化的母体白细胞位点" => 支持白细胞构成机制
# 数据: data/geo_methylation/GSE37722_series_matrix.txt.gz (GPL8490, 27,103 探针)
#       data/annot/HM27.hg38.manifest.tsv.gz (Zhou, hg38 坐标, 与 DMR(hg38) 直接重叠)
# 输出: results/GSE282512_gse37722_validation.csv / _summary.txt
#       figures/gse37722_stage_trend.png
# 约定: UTF-8, 相对路径, LANG=en_US.UTF-8
suppressPackageStartupMessages({
  library(data.table); library(limma); library(ggplot2)
})
set.seed(42)

## ============ 1. 解析 series matrix (27k 行, 手工解析避免 strsplit 行尾丢列) ============
ln <- readLines(gzfile("data/geo_methylation/GSE37722_series_matrix.txt.gz"))
mat_b <- grep("^!series_matrix_table_begin", ln); mat_e <- grep("^!series_matrix_table_end", ln)
tab <- ln[(mat_b + 1):(mat_e - 1)]
hdr <- gsub('"', "", strsplit(tab[1], "\t", fixed = TRUE)[[1]])   # ID_REF + 84 GSM
gsm <- hdr[-1]
cells <- tab[-1]
parts <- strsplit(gsub('"', "", cells), "\t", fixed = TRUE)
probe_id <- vapply(parts, `[`, "", 1L)
vals <- do.call(rbind, lapply(parts, function(x) {
  v <- as.numeric(x[-1]); if (length(v) < length(gsm)) v <- c(v, rep(NA, length(gsm) - length(v))); v
}))
rownames(vals) <- probe_id; colnames(vals) <- gsm
cat(sprintf("GSE37722 矩阵: %d 探针 x %d 样本; 缺失 %0.1f%%\n",
            nrow(vals), ncol(vals), 100 * mean(is.na(vals))))

## 元数据: 孕期状态 + 个体
meta_lines <- grep("^!Sample_characteristics_ch1", ln, value = TRUE)
get_field <- function(lines, key) {
  for (l in lines) {
    f <- gsub('"', "", strsplit(l, "\t", fixed = TRUE)[[1]])[-1]
    hit <- grepl(key, f, fixed = TRUE)
    if (any(hit)) return(sub(paste0(key, ": *"), "", f[hit]))
  }
  rep(NA_character_, length(gsm))
}
stage <- get_field(meta_lines, "pregnancy status")
stage <- trimws(sub(" *$", "", stage))
indiv <- get_field(meta_lines, "individual identifier")
title <- gsub('"', "", strsplit(grep("^!Sample_title", ln, value = TRUE), "\t", fixed = TRUE)[[1]])[-1]
pheno <- data.table(gsm = gsm, stage = stage, indiv = ifelse(indiv == "", NA, indiv), title = title)
print(pheno[, .N, by = stage])
stopifnot(!any(is.na(pheno$stage)))

## ============ 2. DMR 探针定位 (HM27 hg38 manifest × 166 DMR) ============
man <- fread("data/annot/HM27.hg38.manifest.tsv.gz",
             select = c("Probe_ID","CpG_chrm","CpG_beg"))
man <- man[!is.na(CpG_beg) & CpG_chrm %in% paste0("chr", c(1:22, "X", "Y"))]
man[, cend := CpG_beg]
setkey(man, CpG_chrm, CpG_beg, cend)

dmr <- fread("results/GSE282512_dmr_final_v2.csv",
             select = c("region_id","chr","start","end","delta_beta","direction","symbol"))
dmr <- dmr[!is.na(region_id)]
dmr[, region_id := as.character(region_id)]
dmr[, qend := end]
setkey(dmr, chr, start, qend)
ov <- foverlaps(man, dmr, by.x = c("CpG_chrm","CpG_beg","cend"),
                by.y = c("chr","start","qend"), nomatch = NULL)
ov[, region_id := as.character(region_id)]
ov <- ov[Probe_ID %in% rownames(vals)]
cat(sprintf("DMR 探针定位: 166 DMR -> %d 探针 (%d DMR 至少 1 探针)\n",
            uniqueN(ov$Probe_ID), uniqueN(ov$region_id)))
## 敏感性: ±2kb flank
dmr2 <- copy(dmr); dmr2[, `:=`(start = start - 2000, qend = end + 2000)]
setkey(dmr2, chr, start, qend)
ov2 <- foverlaps(man, dmr2, by.x = c("CpG_chrm","CpG_beg","cend"),
                 by.y = c("chr","start","qend"), nomatch = NULL)
ov2[, region_id := as.character(region_id)]
ov2 <- ov2[Probe_ID %in% rownames(vals)]

## ============ 3. 孕期趋势全基因组检验 (limma, 重复测量 duplicateCorrelation) ============
sc <- c(nulligravid = 0, early = 1, middle = 2, delivery = 3, postpartum = 4)
pheno[, trend := sc[stage]]
keep <- rowSums(is.na(vals)) == 0
vg <- vals[keep, , drop = FALSE]
pg <- pheno  # 样本顺序 = 矩阵列顺序
design <- model.matrix(~ trend, data = pg)
dupcor <- tryCatch(duplicateCorrelation(vg, design, block = pg$indiv), error = function(e) NULL)
fit <- if (!is.null(dupcor) && is.finite(dupcor$consensus.correlation)) {
  eBayes(dreamless_fit <- lmFit(vg, design, block = pg$indiv, correlation = dupcor$consensus.correlation))
} else eBayes(lmFit(vg, design))
tr <- summary(decideTests(fit))  # not used, avoid print
tt <- topTable(fit, coef = "trend", number = Inf, sort.by = "none")
tt$Probe_ID <- rownames(vg)
trend_res <- data.table(Probe_ID = tt$Probe_ID, t_trend = tt$t, p_trend = tt$P.Value, fdr_trend = tt$adj.P.Val)
cat(sprintf("全基因组孕期趋势: %d 探针, FDR<0.05: %d\n", nrow(trend_res), sum(trend_res$fdr_trend < 0.05)))

## ============ 4. DMR 探针的孕期趋势 ============
dmr_trend <- trend_res[Probe_ID %in% ov$Probe_ID]
dmr_trend <- merge(dmr_trend, ov[, .(Probe_ID, region_id, delta_beta, direction, symbol)],
                   by = "Probe_ID")
cat(sprintf("DMR 探针趋势: n=%d, t>0: %d, t<0: %d; 方向(DMR hyper 的 t 均值 / hypo 的 t 均值):\n",
            nrow(dmr_trend), sum(dmr_trend$t_trend > 0), sum(dmr_trend$t_trend < 0)))
print(dmr_trend[, .(n = .N, mean_t = mean(t_trend), p_stouffer = 2 * pnorm(-abs(sum(t_trend)/sqrt(.N))),
                    p_sign = binom.test(sum(t_trend > 0), .N)$p.value), by = direction])

## ============ 5. 富集检验: DMR 探针是否富集于孕期动态探针 ============
trend_res[, rank_pct := frank(-abs(t_trend)) / .N]
topq <- trend_res[rank_pct <= 0.05, Probe_ID]
flag <- trend_res$Probe_ID %in% ov$Probe_ID
topflag <- trend_res$Probe_ID %in% intersect(trend_res[rank_pct <= 0.05, Probe_ID], ov$Probe_ID)
ctab <- matrix(c(sum(topflag), sum(flag) - sum(topflag),
                 sum(trend_res$rank_pct <= 0.05) - sum(topflag),
                 nrow(trend_res) - sum(trend_res$rank_pct <= 0.05) - (sum(flag) - sum(topflag))),
               2, 2, byrow = TRUE)
fish <- fisher.test(ctab, alternative = "greater")
cat(sprintf("富集: DMR 探针 %d, top5%% 动态探针中 DMR %d; OR=%.2f p=%.3g\n",
            sum(flag), sum(topflag), unname(fish$estimate), fish$p.value))

## ============ 6. 白细胞层效应 vs cfDNA 效应方向一致性 ============
## delivery - nulligravid Δβ (白细胞层) vs DMR delta_beta (cfDNA 层)
dv <- which(pg$stage == "delivery"); ng <- which(pg$stage == "nulligravid")
leuk_delta <- rowMeans(vg[, dv, drop = FALSE]) - rowMeans(vg[, ng, drop = FALSE])
ld <- data.table(Probe_ID = rownames(vg), leuk_delta = leuk_delta)
cmp <- merge(dmr_trend[, .(Probe_ID, region_id, delta_beta, direction, symbol)], ld, by = "Probe_ID")
rho <- suppressWarnings(cor.test(cmp$delta_beta, cmp$leuk_delta, method = "spearman"))
cat(sprintf("cfDNA Δβ vs 白细胞 delivery-nulligravid Δβ: rho=%.3f p=%.3g (n=%d)\n",
            rho$estimate, rho$p.value, nrow(cmp)))

## ============ 7. 输出 ============
out <- merge(dmr_trend, ld, by = "Probe_ID", all.x = TRUE)
setorder(out, region_id, Probe_ID)
fwrite(out, "results/GSE282512_gse37722_validation.csv")
fwrite(cmp, "results/GSE282512_gse37722_cfDNA_vs_leuk.csv")

## 图: DMR 探针平均甲基化随孕期 + 散点
dmr_mat <- vg[intersect(ov$Probe_ID, rownames(vg)), , drop = FALSE]
pd <- data.table(stage = rep(pg$stage, each = nrow(dmr_mat)), beta = as.vector(dmr_mat))
pd[, stage := factor(stage, levels = c("nulligravid","early","middle","delivery","postpartum"))]
mu <- pd[, .(m = mean(beta, na.rm = TRUE), se = sd(beta, na.rm = TRUE)/sqrt(.N)), by = stage]
p1 <- ggplot(mu, aes(stage, m)) +
  geom_col(fill = "#4C72B0", width = 0.65) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.2) +
  labs(x = NULL, y = "Mean methylation (DMR probes)",
       title = "Maternal leukocyte methylation at PE cfDNA DMR loci across pregnancy") +
  theme_bw(base_size = 11)
p2 <- ggplot(cmp, aes(delta_beta, leuk_delta)) +
  geom_point(aes(color = direction), size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30") +
  scale_color_manual(values = c(hyper = "#C44E52", hypo = "#55A868")) +
  labs(x = "cfDNA DMR delta-beta (PE - Control)", y = "Leukocyte delta-beta (delivery - nulligravid)",
       title = sprintf("Effect concordance (Spearman rho = %.2f, p = %.3g)",
                       rho$estimate, rho$p.value)) +
  theme_bw(base_size = 11)
ggsave("figures/gse37722_stage_trend.png", p1, width = 6, height = 4.5, dpi = 300)
ggsave("figures/gse37722_concordance.png", p2, width = 6, height = 4.5, dpi = 300)

## summary
sink("results/GSE282512_gse37722_summary.txt")
cat("== GSE37722 白细胞层重现验证 ==\n\n")
cat(sprintf("数据: 84 母体白细胞样本 (nulligravid %d / early %d / middle %d / delivery %d / postpartum %d), HM27 %d 探针\n",
            pheno[stage=="nulligravid",.N], pheno[stage=="early",.N], pheno[stage=="middle",.N],
            pheno[stage=="delivery",.N], pheno[stage=="postpartum",.N], nrow(vg)))
cat(sprintf("定位: 166 候选 DMR 内 HM27 探针 %d 个 (覆盖 %d DMR); ±2kb 敏感性 %d 探针\n\n",
            uniqueN(ov$Probe_ID), uniqueN(ov$region_id), uniqueN(ov2$Probe_ID)))
cat(sprintf("全基因组孕期趋势 (limma ~ trend, 重复测量相关 %s): FDR<0.05 探针 %d\n",
            if (!is.null(dupcor)) sprintf("%.2f", dupcor$consensus.correlation) else "不可估",
            sum(trend_res$fdr_trend < 0.05)))
cat(sprintf("DMR 探针趋势汇总:\n"))
print(dmr_trend[, .(n = .N, mean_t = round(mean(t_trend), 3),
                    p_stouffer = signif(2 * pnorm(-abs(sum(t_trend)/sqrt(.N))), 3),
                    p_sign = signif(binom.test(sum(t_trend > 0), .N)$p.value, 3)), by = direction])
cat(sprintf("\n富集 (DMR 探针 vs 全基因组 top5%% 动态): OR=%.2f, p=%.3g\n",
            unname(fish$estimate), fish$p.value))
cat(sprintf("cfDNA DMR Δβ vs 白细胞 (delivery-nulligravid) Δβ: rho=%.3f, p=%.3g, n=%d\n",
            rho$estimate, rho$p.value, nrow(cmp)))
sink()

cat("\n[完成] 输出: GSE282512_gse37722_validation.csv / _summary.txt, 2 图\n")
