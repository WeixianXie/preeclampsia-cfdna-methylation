# 53_gse37722_concordance.R — GSE37722 白细胞层 × GSE154378 cfDNA 细胞标记方向互证 (任务 #50 下半)
# 假说: PE cfDNA DMR = 母体白细胞构成变化的表观指纹。
#   若 GSE154378 cfDNA 中细胞标记轨迹由白细胞构成驱动,
#   则同一批 Sun 标记在 GSE37722 真实白细胞层 (混合 DNA) 中应呈同方向孕期变化。
# 设计:
#   [A] GSE37722 (HM27, 84 样本, nulligravid/early/middle/delivery/postpartum):
#       Sun 细胞标记区域 (hg38 回映) ±flank 内探针, 每探针 Δβ = delivery - nulligravid
#   [B] GSE154378 (cfDNA WGBS, read 级): 每 bin Δβ = delivery(Normal) - NP
#   [C] sign 对齐: contrast = ref[靶细胞] - mean(ref[其余两种血细胞)] (Sun 百分比)
#       sign 后正值 = 靶细胞比例上升; Placenta-HIGH 作 cfDNA 特异性负对照 (白细胞层应平坦)
#   [D] set 级方向一致性 + bin 级 Spearman + 轨迹图
# 输出: results/GSE37722_probe_marker_*.csv, results/GSE154378_bin_marker_delta.csv,
#       results/GSE37722_GSE154378_concordance_summary.txt,
#       figures/GSE37722_GSE154378_marker_trajectory.png, figures/GSE37722_GSE154378_delta_scatter.png
suppressPackageStartupMessages({
  library(data.table)
})
set.seed(42)

FLANKS <- c(2000, 5000, 10000)
FL_PRI  <- 5000    # 主分析层
BC      <- c("Neutrophils", "T-cells", "B-cells")

## ============ 1. GSE37722 矩阵 (复用 45 的解析逻辑) ============
ln <- readLines(gzfile("data/geo_methylation/GSE37722_series_matrix.txt.gz"))
mat_b <- grep("^!series_matrix_table_begin", ln); mat_e <- grep("^!series_matrix_table_end", ln)
tab <- ln[(mat_b + 1):(mat_e - 1)]
hdr <- gsub('"', "", strsplit(tab[1], "\t", fixed = TRUE)[[1]])
gsm <- hdr[-1]
cells <- tab[-1]
parts <- strsplit(gsub('"', "", cells), "\t", fixed = TRUE)
probe_id <- vapply(parts, `[`, "", 1L)
vals <- do.call(rbind, lapply(parts, function(x) {
  v <- as.numeric(x[-1]); if (length(v) < length(gsm)) v <- c(v, rep(NA, length(gsm) - length(v))); v
}))
rownames(vals) <- probe_id; colnames(vals) <- gsm
meta_lines <- grep("^!Sample_characteristics_ch1", ln, value = TRUE)
get_field <- function(lines, key) {
  for (l in lines) {
    f <- gsub('"', "", strsplit(l, "\t", fixed = TRUE)[[1]])[-1]
    hit <- grepl(key, f, fixed = TRUE)
    if (any(hit)) return(sub(paste0(key, ": *"), "", f[hit]))
  }
  rep(NA_character_, length(gsm))
}
stage <- trimws(get_field(meta_lines, "pregnancy status"))
pheno <- data.table(gsm = gsm, stage = stage)
STAGES <- c("nulligravid", "early", "middle", "delivery", "postpartum")
print(pheno[, .N, by = stage])
cat(sprintf("GSE37722: %d 探针 x %d 样本\n", nrow(vals), ncol(vals)))
stopifnot(!any(is.na(pheno$stage)))

dv <- which(pheno$stage == "delivery"); ng <- which(pheno$stage == "nulligravid")
## 每探针 stage 均值 + Δβ
smat <- sapply(STAGES, function(s) {
  idx <- which(pheno$stage == s); rowMeans(vals[, idx, drop = FALSE], na.rm = TRUE)
})
dprobe <- data.table(Probe_ID = rownames(vals), smat)
setnames(dprobe, STAGES, paste0("leuk_", STAGES))
dprobe[, leuk_delta := leuk_delivery - leuk_nulligravid]
dprobe <- dprobe[!is.na(leuk_delta) & !is.na(leuk_early)]

## ============ 2. 标记×探针映射 (hg38, flank 分层) ============
mk <- fread("data/geo_methylation/GSE154378/sun2015_markers_hg38_lift.tsv")
mk[, set := ifelse(type == "I", tissue, "Placenta-HIGH")]
## sign: 血细胞对比 (靶细胞 vs 其余两种血细胞, Sun %)
mk[, contrast := fifelse(tissue == "Neutrophils", Neutrophils - (`T-cells` + `B-cells`) / 2,
                   fifelse(tissue == "T-cells",   `T-cells` - (Neutrophils + `B-cells`) / 2,
                   fifelse(tissue == "B-cells",   `B-cells` - (Neutrophils + `T-cells`) / 2,
                           Placenta - (Neutrophils + `T-cells` + `B-cells`) / 3)))]
mk[, sign := sign(contrast)]
mk[, strict := abs(contrast) >= 15]
cat(sprintf("标记: %d (Neutrophils %d / T-cells %d / B-cells %d / Placenta-HIGH %d); contrast 中位数 %s\n",
            nrow(mk), sum(mk$set=="Neutrophils"), sum(mk$set=="T-cells"), sum(mk$set=="B-cells"),
            sum(mk$set=="Placenta-HIGH"),
            paste(sprintf("%s=%.0f", sort(unique(mk$set)),
                  mk[, median(contrast), by=set][, setorder(.SD, set)][, V1]), collapse=" ")))

man <- fread("data/annot/HM27.hg38.manifest.tsv.gz", select = c("Probe_ID","CpG_chrm","CpG_beg"))
man <- man[!is.na(CpG_beg) & CpG_chrm %in% paste0("chr", c(1:22, "X", "Y")) & Probe_ID %in% rownames(vals)]
man[, cend := CpG_beg]

map_list <- list()
for (fl in FLANKS) {
  m <- copy(mk); m[, `:=`(s2 = start38 - fl, e2 = end38 + fl)]
  setkey(m, chr38, s2, e2)
  ov <- foverlaps(man, m, by.x = c("CpG_chrm","CpG_beg","cend"),
                  by.y = c("chr38","s2","e2"), nomatch = NULL)
  ov <- ov[, .(Probe_ID, bin, set, type, tissue, sign, contrast, strict, flank = fl)]
  map_list[[as.character(fl)]] <- ov
  cat(sprintf("flank +/-%d: 探针 %d, bin %d (%s)\n", fl, uniqueN(ov$Probe_ID), uniqueN(ov$bin),
              paste(sprintf("%s=%d", sort(unique(ov$set)),
                    ov[, uniqueN(Probe_ID), by=set][, setorder(.SD, set)][, V1]), collapse=" ")))
}
map_all <- rbindlist(map_list)

## ============ 3. GSE37722 set 级 sign 对齐 Δβ ============
pl <- merge(map_all, dprobe, by = "Probe_ID")
pl[, leuk_delta_sa := leuk_delta * sign]
## set 级: 同一探针可能落入多个同 set bin -> 去重 (探针 x set x flank)
pl_set <- unique(pl, by = c("Probe_ID", "set", "flank"))
res37722 <- pl_set[, .(
  n_probe = .N,
  mean_delta = mean(leuk_delta),
  mean_delta_sa = mean(leuk_delta_sa),
  p_wilcox = tryCatch(wilcox.test(leuk_delta_sa, mu = 0)$p.value, error = function(e) NA),
  p_sign = tryCatch(binom.test(sum(leuk_delta_sa > 0), .N)$p.value, error = function(e) NA)
), by = .(set, flank)]
setorder(res37722, flank, set)
cat("\n== GSE37722 白细胞层 set 级 (sign 对齐 Δβ = delivery - nulligravid) ==\n")
print(res37722[flank == FL_PRI])

## 5 阶段轨迹 (set 级, 原始 β)
traj37722 <- pl_set[, lapply(.SD, mean, na.rm = TRUE),
                    .SDcols = paste0("leuk_", STAGES), by = .(set, flank)]
traj37722 <- melt(traj37722, id.vars = c("set","flank"), variable.name = "stage", value.name = "beta")
traj37722[, stage := sub("leuk_", "", stage)]

## ============ 4. GSE154378 cfDNA bin 级 Δβ ============
sm <- fread("data/geo_methylation/GSE154378/gse154378_samples.tsv")
mc <- fread("results/GSE154378_bin_mc_long.csv.gz")
mc <- merge(mc, sm[, .(gsm, group, timepoint)], by = "gsm")
mc[, tp2 := fifelse(group == "NP", "NP", timepoint)]
TPS <- c("1stT","2ndT","3rdT","delivery")
C_MIN <- 30
agg <- mc[group %in% c("Normal","NP") & tp2 %in% c(TPS, "NP"),
          .(m = sum(m), c = sum(c)), by = .(group, tp2, bin)]
agg[, beta := m / c]
agg <- dcast(agg, bin ~ interaction(group, tp2, sep = "_"), value.var = "beta")
cc <- dcast(mc[group %in% c("Normal","NP") & tp2 %in% c(TPS, "NP"),
              .(c = sum(c)), by = .(group, tp2, bin)],
            bin ~ interaction(group, tp2, sep = "_"), value.var = "c")
cnames <- setdiff(names(cc), "bin")
setnames(cc, cnames, paste0("c_", cnames))
agg <- agg[cc, on = "bin"]
keep <- !is.na(agg$c_Normal_delivery) & agg$c_Normal_delivery >= C_MIN &
        !is.na(agg$c_NP_NP) & agg$c_NP_NP >= C_MIN
agg <- agg[keep]
agg[, cfdna_delta := Normal_delivery - NP_NP]
cat(sprintf("\nGSE154378 bin: c_min=%d 过滤后 %d bin (标记内 %d)\n",
            C_MIN, nrow(agg), nrow(agg[bin %in% mk$bin])))

binmr <- merge(agg[, .(bin, cfdna_delta, Normal_1stT, Normal_2ndT, Normal_3rdT,
                       Normal_delivery, NP_NP, c_Normal_delivery, c_NP_NP)],
               mk[, .(bin, set, sign, contrast, strict)], by = "bin")
binmr[, cfdna_delta_sa := cfdna_delta * sign]
res154 <- binmr[, .(
  n_bin = .N,
  mean_delta = mean(cfdna_delta),
  mean_delta_sa = mean(cfdna_delta_sa),
  p_wilcox = tryCatch(wilcox.test(cfdna_delta_sa, mu = 0)$p.value, error = function(e) NA),
  p_sign = tryCatch(binom.test(sum(cfdna_delta_sa > 0), .N)$p.value, error = function(e) NA)
), by = set]
setorder(res154, set)
cat("\n== GSE154378 cfDNA set 级 (sign 对齐 Δβ = delivery - NP) ==\n")
print(res154)

## ============ 5. 跨队列方向一致性 ============
cmp <- merge(res37722[flank == FL_PRI, .(set, n_probe, leuk_sa = mean_delta_sa, leuk_p = p_wilcox)],
             res154[, .(set, n_bin, cfdna_sa = mean_delta_sa, cfdna_p = p_wilcox)], by = "set")
cmp[, direction_agree := sign(leuk_sa) == sign(cfdna_sa)]
cat("\n== 方向互证 (主层 +/-%d, GSE37722 Δβ=delivery-nulligravid, GSE154378 Δβ=delivery-NP) ==\n", FL_PRI)
print(cmp)

## bin 级 Spearman: GSE37722 每 bin 探针均值 Δβ vs GSE154378 bin Δβ
bl <- pl[flank == FL_PRI, .(leuk_delta_bin = mean(leuk_delta)), by = .(bin, set)]
bcmp <- merge(bl, binmr[, .(bin, cfdna_delta, set)], by = c("bin","set"))
sp <- suppressWarnings(cor.test(bcmp$leuk_delta_bin, bcmp$cfdna_delta, method = "spearman"))
cat(sprintf("\nbin 级一致性 (±%d): n=%d bin, Spearman rho=%.3f p=%.3g\n",
            FL_PRI, nrow(bcmp), sp$estimate, sp$p.value))
## 分 set 也算
sp_set <- bcmp[, .(rho = tryCatch(suppressWarnings(cor.test(leuk_delta_bin, cfdna_delta,
              method="spearman"))$estimate, error=function(e) NA),
              p = tryCatch(suppressWarnings(cor.test(leuk_delta_bin, cfdna_delta,
              method="spearman"))$p.value, error=function(e) NA), n = .N), by = set]
print(sp_set)

## ============ 6. 输出 ============
fwrite(pl[flank == FL_PRI], "results/GSE37722_probe_marker_map.csv")
fwrite(pl_set, "results/GSE37722_probe_marker_delta.csv")
fwrite(binmr, "results/GSE154378_bin_marker_delta.csv")
fwrite(bcmp, "results/GSE37722_GSE154378_bin_concordance.csv")

## ============ 7. 图 ============
library(ggplot2)
## 图 A: 轨迹 (facet = set, color = cohort)
t1 <- traj37722[flank == FL_PRI][, .(set, stage, beta, cohort = "GSE37722 leukocyte")]
t1[, x := match(stage, STAGES)]
t1[, stage2 := stage]
t2 <- binmr[, .(beta = mean(NP_NP), set, stage = "NP", cohort = "GSE154378 cfDNA")]
t2 <- rbind(t2,
  rbindlist(lapply(c("1stT","2ndT","3rdT","delivery"), function(tp)
    binmr[, .(beta = mean(get(paste0("Normal_", tp))), set, stage = tp, cohort = "GSE154378 cfDNA")])))
t2[, x := match(stage, c("NP","1stT","2ndT","3rdT","delivery"))]
t2[, stage2 := stage]
tr <- rbind(t1, t2)
tr[, set := factor(set, levels = c("Neutrophils","T-cells","B-cells","Placenta-HIGH"))]
pA <- ggplot(tr, aes(x, beta, color = cohort, group = cohort)) +
  geom_point(size = 1.8) + geom_line(linewidth = 0.8) +
  facet_grid(cohort ~ set, scales = "free_y") +
  scale_x_continuous(breaks = 1:5,
    labels = c("1","2","3","4","5")) +
  labs(x = "Time point (GSE37722: nulligravid→postpartum; GSE154378: NP→delivery)",
       y = "Mean beta at marker regions",
       title = "Cell-marker trajectories: leukocyte layer (GSE37722) vs cfDNA (GSE154378)",
       color = "Cohort") +
  theme_bw(base_size = 10)
png("figures/GSE37722_GSE154378_marker_trajectory.png", width = 2600, height = 1500, res = 300)
print(pA); dev.off()

## 图 B: bin 级散点
bcmp[, set := factor(set, levels = c("Neutrophils","T-cells","B-cells","Placenta-HIGH"))]
pB <- ggplot(bcmp, aes(leuk_delta_bin, cfdna_delta)) +
  geom_point(aes(color = set), size = 2.2, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  scale_color_manual(values = c(Neutrophils = "#DD8452", "T-cells" = "#4C72B0",
                                "B-cells" = "#55A868", "Placenta-HIGH" = "#C44E52")) +
  labs(x = "GSE37722 leukocyte delta-beta (delivery - nulligravid)",
       y = "GSE154378 cfDNA delta-beta (delivery - NP)",
       title = sprintf("Marker-level concordance (Spearman rho = %.2f, p = %.2g, n = %d)",
                       sp$estimate, sp$p.value, nrow(bcmp)), color = NULL) +
  theme_bw(base_size = 11)
png("figures/GSE37722_GSE154378_delta_scatter.png", width = 2000, height = 1600, res = 300)
print(pB); dev.off()

## ============ 8. summary ============
sink("results/GSE37722_GSE154378_concordance_summary.txt")
cat("== GSE37722 白细胞层 × GSE154378 cfDNA 细胞标记方向互证 ==\n\n")
cat(sprintf("标记: Sun v4 细胞/胎盘标记 hg38 回映 %d 个 (Neutrophils %d / T-cells %d / B-cells %d / Placenta-HIGH %d)\n",
            nrow(mk), sum(mk$set=="Neutrophils"), sum(mk$set=="T-cells"),
            sum(mk$set=="B-cells"), sum(mk$set=="Placenta-HIGH")))
cat("sign 对齐: contrast = ref[靶细胞] - mean(其余两种血细胞); 正值 = 靶细胞比例上升\n\n")
cat(sprintf("[A] GSE37722 白细胞层 (HM27, %d 探针矩阵):\n", nrow(vals)))
print(pheno[, .N, by = stage])
cat(sprintf("    flank 主层 +/-%dbp\n\n", FL_PRI))
print(res37722[flank == FL_PRI])
cat("\n    敏感性 (flank +/-%d 与 +/-%d):\n", FLANKS[1], FLANKS[3])
print(res37722[flank != FL_PRI])
cat("\n[B] GSE154378 cfDNA (bin 级, c_min=%d):\n\n", C_MIN)
print(res154)
cat("\n[C] 方向互证汇总:\n\n")
print(cmp)
cat(sprintf("\n[D] bin 级 Spearman: rho=%.3f, p=%.3g, n=%d bin\n", sp$estimate, sp$p.value, nrow(bcmp)))
print(sp_set)
cat("\n结论要点:\n")
cat("- 若 Neutrophils 白细胞层 sign 对齐 Δβ>0 且 T-cells/B-cells<0, 与 cfDNA 同向 => 跨队列跨模态支持白细胞构成机制\n")
cat("- Placenta-HIGH 在白细胞层应平坦 (无胎盘贡献), 而在 cfDNA 上升 => 证明胎盘信号为 cfDNA 特异\n")
sink()
cat("\n[完成] 输出 3 csv + summary + 2 图\n")
