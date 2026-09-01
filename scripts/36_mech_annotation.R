# 36_mech_annotation.R —— 机制层 ③: DMR 调控元件富集 + 胎盘组织甲基化复现/来源检验
# A) Ensembl regbuild (enhancer/promoter/open chromatin/CTCF/TFBS) 富集: 候选166 vs 背景72982
# B) GSE57767/GSE73375/GSE75196 (HM450 胎盘 PE vs 对照) 组织级复现 (awk 流式抽取)
# C) 胎盘来源检验: cfDNA Δβ vs (胎盘β - 对照cfDNAβ) 相关性
# 输出: results/GSE282512_regbuild_enrichment.csv/.txt, placenta_tissue.csv, placenta_origin.csv

suppressPackageStartupMessages({library(data.table)})
set.seed(1)

## ============ A. regbuild 调控元件富集 ============
rb <- fread("data/annot/regbuild.gff.gz", header = FALSE, sep = "\t",
            select = c(1, 3, 4, 5), col.names = c("chr", "ftype", "fstart", "fend"))
rb[, chr := as.character(chr)]
rb <- rb[chr %in% as.character(1:22) | chr %in% c("X", "Y")]
cat(sprintf("regbuild 特征: %d 个\n", nrow(rb)))

reg <- fread("results/GSE282512_region_annot.csv",
             select = c("region_id", "chr", "start", "end"))
reg[, chr := sub("chr", "", chr)]
tested <- fread("results/GSE282512_dmr_candidates.csv", select = c("region_id"))
reg <- reg[region_id %in% tested$region_id]
v2 <- fread("results/GSE282512_dmr_final_v2.csv",
            select = c("region_id", "chr", "start", "end", "symbol", "type",
                       "direction", "delta_beta", "final_call"))
v2[, chr := sub("chr", "", chr)]
bg <- reg[!region_id %in% v2$region_id]
cat(sprintf("候选 %d, 背景 %d\n", nrow(v2), nrow(bg)))

ovl_hit <- function(regions, feats) {
  qs <- regions[, .(chr, qstart = start, qend = end)]
  sub <- feats[, .(chr, fstart, fend)]
  setkey(sub, chr, fstart, fend)
  fo <- foverlaps(qs, sub, by.x = c("chr", "qstart", "qend"),
                  by.y = c("chr", "fstart", "fend"),
                  nomatch = NULL, which = TRUE)
  length(unique(fo$xid))
}
enrich <- rbindlist(lapply(unique(rb$ftype), function(ft) {
  cat(sprintf("  overlap: %s (%d feats)...\n", ft, sum(rb$ftype == ft)))
  feats <- rb[ftype == ft]
  data.table(feature = ft,
             cand_n = ovl_hit(v2, feats), cand_tot = nrow(v2),
             bg_n = ovl_hit(bg, feats), bg_tot = nrow(bg))
}))
enrich[, prop_cand := cand_n / cand_tot][, prop_bg := bg_n / bg_tot]
enrich[, or := (cand_n / (cand_tot - cand_n)) / (bg_n / (bg_tot - bg_n))]
enrich[, p_fisher := mapply(function(a, b, c, d)
  fisher.test(matrix(c(a, b - a, c, d - c), nrow = 2))$p.value,
  cand_n, cand_tot, bg_n, bg_tot)]
enrich[, p_bh := p.adjust(p_fisher, "BH")]
fwrite(enrich, "results/GSE282512_regbuild_enrichment.csv")
cat("Part A 完成\n")

## ============ B. 胎盘组织甲基化复现 ============
## 先确定 CpG 列表, 再 awk 流式抽取 (避免整表入内存)
man <- fread("data/annot/HM450.hg38.manifest.tsv.gz",
             select = c("Probe_ID", "CpG_chrm", "CpG_beg"))
setnames(man, c("cpg", "cchr", "cbeg"))
man <- man[!is.na(cchr) & grepl("^chr[0-9XY]{1,2}$", cchr)]
man[, cchr := sub("chr", "", cchr)]
dmrreg <- v2[, .(chr, qstart = start, qend = end, region_id)]
setkey(dmrreg, chr, qstart, qend)   # foverlaps 要求 y 前 3 列 = 键列
man[, cend := cbeg]  # foverlaps 不允许 by.x 重复列, 复制作端点
fo <- foverlaps(man[, .(cchr, cbeg, cend, cpg)], dmrreg,
                by.x = c("cchr", "cbeg", "cend"),
                by.y = c("chr", "qstart", "qend"), nomatch = NULL)
dmr_cpg <- unique(fo[, .(cpg, region_id)])
cat(sprintf("DMR 区域内 HM450 CpG: %d (覆盖 %d 个候选)\n",
            nrow(dmr_cpg), uniqueN(dmr_cpg$region_id)))

## 位点级 cfDNA 表 (仅可测位点) -> cpg 联表
sl <- fread("results/GSE282512_dmr_site_level.csv")
sl <- sl[!is.na(beta_pe) & !is.na(beta_ctrl) & n_pe > 0 & n_ctrl > 0]
sl2 <- copy(sl)
sl2[, site_chr := sub("^chr", "", sub("_.*", "", site))][, site_pos := as.integer(sub(".*_", "", site))]
cfdna <- man[sl2, on = .(cchr = site_chr, cbeg = site_pos), nomatch = NULL]
cat(sprintf("位点级 cfDNA 表与 manifest 匹配: %d CpG\n", nrow(cfdna)))

want <- union(dmr_cpg$cpg, cfdna$cpg)
writeLines(want, "results/_tmp_cpgs.txt")
cat(sprintf("需从胎盘矩阵抽取 CpG: %d\n", length(want)))

## awk 抽取写成独立 bash 脚本执行 (R system() 走 cmd.exe, 引号语义与 bash 不兼容)
filt <- 'BEGIN{while((getline l < "results/_tmp_cpgs.txt")>0) want[l]=1} /^!/{print; next} /^"?ID_REF/{print; next} {id=$1; gsub(/"/,"",id); if (id in want) print}'
pl_files <- c(GSE57767 = "data/geo_methylation/GSE57767_series_matrix.txt.gz",
              GSE73375 = "data/geo_methylation/GSE73375_series_matrix.txt.gz",
              GSE75196 = "data/geo_methylation/GSE75196_series_matrix.txt.gz")
cmds <- sprintf("gzip -dc %s | awk -F'\\t' '%s' > results/_tmp_%s.txt",
                pl_files, filt, names(pl_files))
writeLines(c("#!/bin/bash", cmds), "results/_tmp_extract.sh")
system("bash results/_tmp_extract.sh")

read_beta_filtered <- function(out, group_pat) {
  lns <- readLines(out)
  meta <- list()
  for (ln in lns[grepl("^!", lns) & !grepl("^!series_matrix_table", lns)]) {
    p <- strsplit(ln, "\t")[[1]]
    if (length(p) > 1) {
      k <- p[1]
      while (k %in% names(meta)) k <- paste0(k, "|", length(grep(paste0("^", sub("\\|", "\\\\|", p[1])), names(meta))) + 1)
      meta[[k]] <- p[-1]
    }
  }
  cl <- function(x) gsub('^"|"$', "", x)
  chf <- list()
  for (k in grep("!Sample_characteristics_ch1", names(meta), value = TRUE)) {
    v <- cl(meta[[k]]); fld <- sub(":.*", "", v[1])
    chf[[fld]] <- trimws(sub("^[^:]*:", "", v))
  }
  gf <- names(chf)[sapply(chf, function(v) any(grepl(group_pat[1], v, ignore.case = TRUE)) |
                                      any(grepl(group_pat[2], v, ignore.case = TRUE)))]
  v <- chf[[gf[1]]]
  g <- rep(NA_character_, length(v))
  g[grepl(group_pat[1], v, ignore.case = TRUE)] <- "PE"
  g[grepl(group_pat[2], v, ignore.case = TRUE)] <- "CT"
  tab <- lns[!grepl("^!", lns)]
  mh <- strsplit(tab[1], "\t")[[1]]
  ncol_m <- length(mh) - 1
  m <- t(vapply(tab[-1], function(l) {
    p <- strsplit(l, "\t", fixed = TRUE)[[1]]
    if (length(p) < length(mh)) p <- c(p, rep("", length(mh) - length(p)))
    suppressWarnings(as.numeric(p[-1]))
  }, numeric(ncol_m)))
  rownames(m) <- gsub('^"|"$', "", sapply(tab[-1],
    function(l) strsplit(l, "\t", fixed = TRUE)[[1]][1]))
  keep <- !is.na(g) & g %in% c("PE", "CT") & seq_along(g) <= ncol_m
  list(mat = m[, keep, drop = FALSE], grp = g[keep])
}

pl <- list()
pl[["GSE57767"]] <- read_beta_filtered("results/_tmp_GSE57767.txt",
  c("preeclamp", "normal"))
pl[["GSE73375"]] <- read_beta_filtered("results/_tmp_GSE73375.txt",
  c("preeclamp", "normotensive"))
pl[["GSE75196"]] <- read_beta_filtered("results/_tmp_GSE75196.txt",
  c("preeclamp", "normal|healthy"))
for (nm in names(pl)) {
  cat(sprintf("[%s] PE=%d CT=%d, %d CpG\n", nm, sum(pl[[nm]]$grp == "PE"),
              sum(pl[[nm]]$grp == "CT"), nrow(pl[[nm]]$mat)))
}

## 组织 Δβ (每队列) + 池化
tissue_list <- list(); pl_all_beta <- list()
for (nm in names(pl)) {
  m <- pl[[nm]]$mat; g <- pl[[nm]]$grp
  cgs <- intersect(rownames(m), cfdna$cpg)
  mm <- m[cgs, , drop = FALSE]
  bpe <- rowMeans(mm[, g == "PE", drop = FALSE], na.rm = TRUE)
  bct <- rowMeans(mm[, g == "CT", drop = FALSE], na.rm = TRUE)
  ball <- rowMeans(mm, na.rm = TRUE)
  tissue_list[[nm]] <- data.table(cpg = cgs, beta_pe_t = bpe, beta_ct_t = bct,
                                  delta_t = bpe - bct, cohort = nm)
  pl_all_beta[[nm]] <- data.table(cpg = cgs, beta_all = ball)
}
tis <- rbindlist(tissue_list)
tis[, z := scale(delta_t), by = cohort]
tis_pool <- tis[, .(delta_t_pool = mean(z), n_coh = .N), by = cpg]
pl_beta <- Reduce(function(a, b) merge(a, b, by = "cpg", all = TRUE), pl_all_beta)
pl_beta[, beta_placenta := rowMeans(.SD, na.rm = TRUE), .SDcols = patterns("^beta_all")]
pl_beta <- pl_beta[, .(cpg, beta_placenta)]

## ============ C. 相关性检验 ============
cf <- cfdna[, .(cpg, beta_pe, beta_ctrl, delta = beta_pe - beta_ctrl)]
dat <- merge(cf, tis_pool, by = "cpg")
dat <- merge(dat, pl_beta, by = "cpg")
cat(sprintf("\n联表后可检验 CpG: %d\n", nrow(dat)))

cor_org <- cor.test(dat$delta, dat$beta_placenta - dat$beta_ctrl,
                    method = "spearman", exact = FALSE)
perm <- sapply(1:10000, function(i) {
  cor(dat$delta, sample(dat$beta_placenta - dat$beta_ctrl), method = "spearman",
      use = "complete.obs")
})
perm_p <- mean(abs(perm) >= abs(cor_org$estimate), na.rm = TRUE)

sink("results/GSE282512_regbuild_enrichment.txt")
cat("===== 机制层 ③-A: DMR 调控元件富集 (Ensembl regbuild) =====\n\n")
print(enrich[, .(feature, cand_n, cand_tot, prop_cand, bg_n, bg_tot, prop_bg,
                 or = round(or, 2), p_fisher, p_bh)])
cat("\n注: regbuild 的 TF_binding_site 不区分具体转录因子, TF 特异性富集需 JASPAR motif 扫描。\n\n")
cat("===== 机制层 ③-B: 胎盘组织甲基化复现 (3 个 HM450 队列) =====\n\n")
for (nm in names(pl)) cat(sprintf("  %s: PE=%d CT=%d\n", nm,
                                  sum(pl[[nm]]$grp == "PE"), sum(pl[[nm]]$grp == "CT")))
cat(sprintf("DMR 内 HM450 CpG %d; 位点级可测 %d\n", nrow(dmr_cpg), nrow(dat)))
cor_tis <- cor.test(dat$delta, dat$delta_t_pool, method = "spearman", exact = FALSE)
sgn <- sign(dat$delta) == sign(dat$delta_t_pool)
binom_p <- binom.test(sum(sgn, na.rm = TRUE), sum(!is.na(sgn)), 0.5)$p.value
cat(sprintf("\n组织 Δβ vs cfDNA Δβ: Spearman rho=%.3f (p=%.3g); 符号一致 %d/%d (%.0f%%), binom p=%.3g\n",
            cor_tis$estimate, cor_tis$p.value, sum(sgn, na.rm = TRUE),
            sum(!is.na(sgn)), 100 * mean(sgn, na.rm = TRUE), binom_p))
for (nm in names(tissue_list)) {
  d2 <- merge(cf, tissue_list[[nm]][, .(cpg, delta_t)], by = "cpg")
  ct <- suppressWarnings(cor.test(d2$delta, d2$delta_t, method = "spearman", exact = FALSE))
  cat(sprintf("  %s: rho=%.3f (p=%.3g), n=%d\n", nm, ct$estimate, ct$p.value, nrow(d2)))
}
cat("\n===== 机制层 ③-C: 胎盘来源检验 =====\n\n")
cat(sprintf("cfDNA Δβ vs (胎盘β - 对照cfDNAβ): Spearman rho=%.3f (p=%.3g),\n",
            cor_org$estimate, cor_org$p.value))
cat(sprintf("置换检验 p=%.4f (n=%d CpG)\n", perm_p, nrow(dat)))
cat("正相关 = PE cfDNA 甲基化向胎盘谱偏移 -> 支持胎盘来源成分变化;\n")
cat("无相关/负相关 = 不支持胎盘成分变化, 需考虑母体造血细胞来源。\n")
sink()

fwrite(dat, "results/GSE282512_placenta_origin.csv")
fwrite(tis, "results/GSE282512_placenta_tissue.csv")

# 清理临时文件
file.remove("results/_tmp_cpgs.txt", "results/_tmp_extract.sh",
            "results/_tmp_GSE57767.txt", "results/_tmp_GSE73375.txt",
            "results/_tmp_GSE75196.txt")
cat("DONE\n")
