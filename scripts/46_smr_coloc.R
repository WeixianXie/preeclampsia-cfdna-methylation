# 46_smr_coloc.R — 3 共享 DMR 的 SMR + coloc 因果定位 (任务 #44, 第二梯队)
# 对象: region 106528 (chr3:49.9Mb), 12268 (chr6 SLC17A1 gene_body), 84604 (chr17:39.45Mb)
# 暴露: GoDML cis-meQTL (GRCh37, beta/se 为 SD 尺度, 无等位基因列)
# 结局: FinnGen R10 PE (7965/211852) 与 GH (9535/211957), GRCh38
# 方法:
#   1) Ensembl REST 区段映射 GRCh37->GRCh38 (分片, 含反向片段)
#   2) SMR: b_zy = b_gy/b_gx (Wald), T = b_zy^2 / se(b_zy)^2 ~ chi2(1)
#      [T 对单侧等位基因翻转不变 (b_zy^2 与 se 均为平方项), 仅方向符号有歧义]
#   3) coloc.abf p 值模式 (免等位基因方向): 甲基化 N=32000 (GoDML meta), GWAS N=PE 219817 / GH 221492
# 限制: GoDML 精简文件无等位基因列、无 LD 参考面板 -> 正式 HEIDI 不可行,
#       以 coloc PP4 (共享因果变异后验) 替代连锁-多效性判别证据
# 约定: UTF-8, 相对路径, LANG=en_US.UTF-8
suppressPackageStartupMessages({
  library(data.table); library(jsonlite); library(httr); library(coloc)
})

N_MQTL <- 32000
GWAS_N <- c(PE = 7965 + 211852, GH = 9535 + 211957)
gwas_s <- c(PE = 7965 / (7965 + 211852), GH = 9535 / (9535 + 211957))
loci <- list(
  list(region_id = "106528", chr = "3",  cpgs = c("cg18129748","cg24308560","cg27360567"),
       symbol = "", direction = "hyper", dmr_delta = 0.1019),
  list(region_id = "12268",  chr = "6",  cpgs = c("cg00387872","cg19490609","cg25753631"),
       symbol = "SLC17A1", direction = "hypo", dmr_delta = -0.1186),
  list(region_id = "84604",  chr = "17", cpgs = "cg15445000",
       symbol = "", direction = "hypo", dmr_delta = -0.1392)
)
gwas_files <- list(PE = "results/_tmp_smr/finngen_pe_window.tsv",
                   GH = "results/_tmp_smr/finngen_gh_window.tsv")

## ============ 1. 工具函数 ============
map_region <- function(chr, pmin, pmax) {
  url <- sprintf("https://rest.ensembl.org/map/human/GRCh37/%s:%d:%d/GRCh38?content-type=application/json",
                 chr, pmin, pmax)
  js <- tryCatch(fromJSON(content(GET(url, timeout(60)), as = "text")),
                 error = function(e) NULL)
  if (is.null(js) || is.null(js$mappings)) return(NULL)
  m <- js$mappings
  data.table(chr37 = as.character(m$original$seq_region_name),
             s37 = as.integer(m$original$start), e37 = as.integer(m$original$end),
             s38 = as.integer(m$mapped$start),   e38 = as.integer(m$mapped$end),
             strand = as.integer(m$mapped$strand))
}
lift_snps <- function(dt, seg) {
  seg <- copy(seg); setkey(seg, chr37, s37, e37)
  q <- dt[, .(chr = as.character(chr), pos = as.integer(pos))]
  q[, `:=`(s = pos, e = pos)]
  ov <- foverlaps(q, seg, by.x = c("chr","s","e"),
                  by.y = c("chr37","s37","e37"), nomatch = NULL)
  ov[, pos38 := fifelse(strand == 1, pos - s37 + s38, e38 - (pos - s37))]
  unique(ov[, .(pos, pos38)])
}

## ============ 2. 主循环 ============
smr_rows <- list(); coloc_rows <- list(); lead_rows <- list()
for (L in loci) {
  cat(sprintf("\n===== region %s chr%s (%s, %s) =====\n",
              L$region_id, L$chr, L$symbol, L$direction))
  ## godmc 行形如 cg02518338,"chr17,47929884,...  — 引号未闭合, 直接按逗号切 6 列
  mq <- rbindlist(lapply(L$cpgs, function(cg)
    fread(sprintf("results/_tmp_smr/mqtl_%s.csv", cg), header = FALSE, quote = "",
          col.names = c("cpg","chr","pos","p_m","beta_m","se_m"))))
  mq[, chr := sub('"', "", sub("chr", "", chr))]
  mq[, pos := as.integer(pos)]
  ## godmc cis 文件仅收录显著 meQTL 对 (p<~1e-5); p 极小/下溢 (含 denormal ~5e-324) 一律封底
  mq[, p_m := fifelse(p_m < 1e-250, 1e-250, p_m)]
  cat(sprintf("meQTL: %d 对 (%d CpG)\n", nrow(mq), uniqueN(mq$cpg)))

  seg <- map_region(L$chr, min(mq$pos) - 500, max(mq$pos) + 500)
  if (is.null(seg)) { cat("区段映射失败, 跳过\n"); next }
  lf <- lift_snps(mq, seg)
  mq <- merge(mq, lf, by = "pos", all.x = FALSE)
  cat(sprintf("映射到 GRCh38: %d meQTL 行 / %d 位点\n", nrow(mq), uniqueN(mq$pos)))

  for (ph in names(gwas_files)) {
    gw <- fread(gwas_files[[ph]])
    gw[, chrom := as.character(chrom)]
    gw <- gw[chrom == L$chr]
    setnames(gw, c("pos","rsids","pval","beta","sebeta","af_alt"),
             c("pos38","rsid","p_gw","beta_gw","se_gw","af"))
    j <- merge(mq, gw[, .(pos38, rsid, p_gw, beta_gw, se_gw, af)], by = "pos38")
    j <- j[!is.na(af) & af > 0.005 & af < 0.995 & !is.na(se_gw) & !is.na(se_m)]
    j <- j[j[, .I[which.min(p_gw)], by = pos38]$V1]   # FinnGen 多等位同位点多行 -> 取 GWAS p 最小行
    cat(sprintf("  [%s] 共享可分析 SNP: %d\n", ph, nrow(j)))

    for (cg in L$cpgs) {
      sub <- j[cpg == cg]
      if (nrow(sub) < 10) { cat(sprintf("  [%s] %s SNP 不足, 跳过\n", ph, cg)); next }
      sub[, z_m := beta_m / se_m]
      lead <- sub[order(p_m, -abs(z_m))][1]

      ## ---- SMR ----
      b_gx <- lead$beta_m; s_gx <- lead$se_m
      b_gy <- lead$beta_gw; s_gy <- lead$se_gw
      b_zy <- b_gy / b_gx
      se_zy <- sqrt(s_gy^2 / b_gx^2 + b_gy^2 * s_gx^2 / b_gx^4)
      T_smr <- (b_zy / se_zy)^2
      p_smr <- pchisq(T_smr, 1, lower.tail = FALSE)
      smr_rows[[length(smr_rows) + 1]] <- data.table(
        region_id = L$region_id, symbol = L$symbol, direction = L$direction,
        pheno = ph, cpg = cg, n_snp = nrow(sub),
        lead_rsid = lead$rsid, lead_af = lead$af,
        p_mqtl = lead$p_m, beta_m = b_gx, se_m = s_gx,
        p_gw = lead$p_gw, beta_gw = b_gy, se_gw = s_gy,
        b_zy = b_zy, se_zy = se_zy, T_smr = T_smr, p_smr = p_smr)
      lead_rows[[length(lead_rows) + 1]] <- data.table(
        region_id = L$region_id, pheno = ph, cpg = cg,
        rsid = lead$rsid, pos38 = lead$pos38, af = lead$af,
        beta_m = b_gx, se_m = s_gx, p_m = lead$p_m,
        beta_gw = b_gy, se_gw = s_gy, p_gw = lead$p_gw)

      ## ---- coloc (p 值模式: pvalues + MAF + N + s) ----
      snp_id <- sprintf("%s:%d", L$chr, sub$pos38)
      d1 <- list(pvalues = sub$p_m, MAF = pmin(sub$af, 1 - sub$af), N = N_MQTL,
                 type = "quant", snp = snp_id)
      d2 <- list(pvalues = sub$p_gw, MAF = pmin(sub$af, 1 - sub$af), N = GWAS_N[ph],
                 s = as.numeric(gwas_s[ph]), type = "cc", snp = snp_id)
      cres <- tryCatch(suppressMessages(coloc.abf(d1, d2)), error = function(e) NULL)
      if (!is.null(cres)) {
        s <- cres$summary
        coloc_rows[[length(coloc_rows) + 1]] <- data.table(
          region_id = L$region_id, symbol = L$symbol, pheno = ph, cpg = cg,
          nsnps = s["nsnps"], PP0 = s["PP.H0.abf"], PP1 = s["PP.H1.abf"],
          PP2 = s["PP.H2.abf"], PP3 = s["PP.H3.abf"], PP4 = s["PP.H4.abf"])
        cat(sprintf("  [%s] %s: PP4=%.3f (n=%d)\n", ph, cg, s["PP.H4.abf"], s["nsnps"]))
      }
    }
  }
}

## ============ 3. 汇总输出 ============
smr <- rbindlist(smr_rows)
smr[, p_smr_bh := p.adjust(p_smr, "BH")]
setorder(smr, region_id, pheno, cpg)
fwrite(smr, "results/GSE282512_smr_results.csv")
cloc <- rbindlist(coloc_rows)
fwrite(cloc, "results/GSE282512_coloc_results.csv")
fwrite(rbindlist(lead_rows), "results/GSE282512_smr_lead_snps.csv")

sink("results/GSE282512_smr_coloc_summary.txt")
cat("== SMR + coloc: 3 个 meQTL-GWAS 共享 DMR 因果定位 ==\n\n")
cat("方法: SMR (Wald 比值, chi2 检验; T 对等位基因方向不敏感) + coloc.abf (p 值模式)\n")
cat(sprintf("暴露: GoDML cis-meQTL (N=32000); 结局: FinnGen R10 PE (7965+211852) / GH (9535+211957)\n\n"))
cat("---- SMR 结果 ----\n")
print(smr[, .(region = region_id, symbol, pheno, cpg, lead_rsid,
              p_mqtl = signif(p_mqtl, 2), p_gw = signif(p_gw, 2),
              b_zy = round(b_zy, 3), T_smr = round(T_smr, 2),
              p_smr = signif(p_smr, 2), p_bh = signif(p_smr_bh, 2))])
cat("\n---- coloc PP4 ----\n")
print(cloc[, .(region = region_id, pheno, cpg, nsnps,
               PP3 = round(PP3, 3), PP4 = round(PP4, 3))])
cat("\n说明:\n")
cat("- GoDML 精简 cis 文件不含等位基因列 -> b_zy 方向符号不可判读 (SMR 检验统计量不受影响)\n")
cat("- 无 LD 参考面板, 正式 HEIDI 不可行; 以 coloc PP4 判定共享因果变异 (PP4>0.8 强, 0.5-0.8 提示)\n")
cat("- GoDML cis 文件仅含显著 meQTL 对 -> coloc 的 SNP 集为 meQTL 显著集 (H0/H1 分量受限), PP3/PP4 相对比较仍有效\n")
cat("- 注: coloc 的 H4(共享) vs H3(两个不同因果变异, 均有信号) 由 PP3/PP4 比值区分\n")
sink()

cat("\n[完成] 输出: GSE282512_smr_results.csv / _coloc_results.csv / _smr_lead_snps.csv / _smr_coloc_summary.txt\n")
