# 58_mr_sensitivity.R — MR 补强: 多 IV IVW + 敏感性套件 + Tyrmi 2023 meta-GWAS 替代结局
# 对象: 3 个共享 DMR 的 7 个 CpG (region 106528 chr3, 12268 SLC17A1 chr6, 84604 chr17)
# 暴露: GoDMC cis-meQTL (GRCh37, allele1=effect, N~32000)
# 结局: FinnGen R10 PE / GH (GRCh38) + Tyrmi 2023 JAMA Cardiol PE meta-GWAS (GCST90269903, GRCh37, 16743/280081)
# LD: 1000G phase3 EUR 本地相关表 (58a_ld_compute.py 产出; Ensembl REST /ld 已退役)
# 方法:
#   1) 等位基因定向 (GoDMC allele1 vs GWAS effect allele), 模糊回文剔除
#   2) 每 CpG: p_m<1e-5 cis-meQTL 为 IV 池; 1000G EUR LD 贪婪 clumping r2<0.1 (LD 缺失且距 lead<50kb 保守剔除)
#   3) 多 IV MR: IVW-FE/RE + Wald ratio + MR-Egger 截距 + 加权中位数 (bootstrap) + leave-one-out + F 统计量
#   4) HEIDI 近似: top IV 的 LD 伙伴 (r2 0.05-0.9) Wald 比值同质性 chi2 (单共享因果变异零假设)
#   5) 输出: CSV + 森林图 + 汇总 txt (BH 校正)
# 约定: UTF-8, 相对路径, LANG=en_US.UTF-8
suppressPackageStartupMessages({
  library(data.table)
})

N_MQTL <- 32000
GWAS_N  <- c(PE_FinnGen = 7965 + 211852, GH_FinnGen = 9535 + 211957, PE_Tyrmi = 16743 + 280081)
loci <- list(
  list(region_id = "106528", chr = "3",  cpgs = c("cg18129748","cg24308560","cg27360567"), symbol = "", direction = "hyper"),
  list(region_id = "12268",  chr = "6",  cpgs = c("cg00387872","cg19490609","cg25753631"), symbol = "SLC17A1", direction = "hypo"),
  list(region_id = "84604",  chr = "17", cpgs = "cg15445000", symbol = "", direction = "hypo")
)
PVAL_IV <- 1e-5      # IV 池阈值
R2_CLUMP <- 0.1      # LD clumping 阈值
R2_HEIDI_LO <- 0.05; R2_HEIDI_HI <- 0.9
NEAR <- 50000        # LD 缺失时的保守距离剔除

## ============ 1. 工具函数 ============
## (GRCh37->GRCh38 回映由 58a_ld_compute.py 通过 Ensembl REST 完成, 本脚本读取 lift_positions.csv)

## 等位基因定向: meQTL (allele1=effect, freq1) vs GWAS (effect allele=ea, other=oa, freq=ef)
harmonize <- function(a1, a2, f1, ea, oa, ef, pal_drop_maf = 0.42) {
  a1 <- toupper(a1); a2 <- toupper(a2); ea <- toupper(ea); oa <- toupper(oa)
  n <- length(a1); sign_v <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    if (is.na(a1[i]) || is.na(ea[i])) next
    if (a1[i] == ea[i] && a2[i] == oa[i]) { sign_v[i] <- 1;  next }
    if (a1[i] == oa[i] && a2[i] == ea[i]) { sign_v[i] <- -1; next }
    pal <- setequal(c(a1[i], a2[i]), c("A","T")) || setequal(c(a1[i], a2[i]), c("C","G"))
    if (pal) {
      maf1 <- min(f1[i], 1 - f1[i], na.rm = TRUE); mafg <- min(ef[i], 1 - ef[i], na.rm = TRUE)
      if (maf1 > pal_drop_maf || mafg > pal_drop_maf) next
      if (abs(f1[i] - ef[i]) < abs(1 - f1[i] - ef[i])) sign_v[i] <- 1 else sign_v[i] <- -1
    }
  }
  sign_v
}

## 贪婪 LD clumping: 按 p_m 升序; 与已选 lead r2>R2_CLUMP 剔除; 不在 1000G 且 |pos差|<NEAR 保守剔除
## 性能: 每个 lead 只做一次 LD 表查询 (伙伴一次拿全), 不逐对扫描
## 注意: 参数名用 chr_ 避免与 ldt 列 chr 同名 (data.table 作用域: ldt[chr == chr] 恒 TRUE)
clump_ivs <- function(iv, ldt, chr_, in1000g) {
  iv <- copy(iv); setorder(iv, p_m)
  iv[, drop := FALSE]
  iv[, known := pos %in% in1000g]
  sel <- integer(0)
  for (i in seq_len(nrow(iv))) {
    if (iv$drop[i]) next
    sel <- c(sel, i)
    pi <- iv$pos[i]
    part <- ldt[chr == chr_ & (pos1 == pi | pos2 == pi)]
    if (nrow(part)) {
      pp <- fifelse(part$pos1 == pi, part$pos2, part$pos1)
      iv$drop[iv$pos %in% pp[part$r2 > R2_CLUMP]] <- TRUE
    }
    near_pos <- iv$pos[!iv$drop & iv$pos != pi & abs(iv$pos - pi) < NEAR]
    iv$drop[iv$pos %in% setdiff(near_pos, iv$pos[iv$known])] <- TRUE   # 1000G 无该位点且距离近 -> 保守剔除
  }
  iv[sel][, `:=`(drop = NULL, known = NULL)]
}

## MR 估计量
mr_estimates <- function(b, se, bx) {
  w <- 1 / se^2
  k <- length(b)
  if (k == 1) {
    return(data.table(method = "Wald_ratio", beta = b, se = se, p = 2 * pnorm(-abs(b / se)),
                      Q = NA_real_, Qp = NA_real_, egger_int = NA_real_, egger_p = NA_real_,
                      wmed_beta = NA_real_, wmed_p = NA_real_))
  }
  b_fe <- sum(w * b) / sum(w); se_fe <- sqrt(1 / sum(w))
  Q <- sum(w * (b - b_fe)^2); Qp <- pchisq(Q, k - 1, lower.tail = FALSE)
  tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
  w_re <- 1 / (se^2 + tau2); b_re <- sum(w_re * b) / sum(w_re); se_re <- sqrt(1 / sum(w_re))
  beta_use <- ifelse(Qp < 0.05, b_re, b_fe); se_use <- ifelse(Qp < 0.05, se_re, se_fe)
  egger_int <- NA_real_; egger_p <- NA_real_
  if (k >= 4 && sd(bx) > 0) {
    fit <- lm(b ~ bx, weights = w)
    cf <- summary(fit)$coefficients
    egger_int <- cf[1, 1]; egger_p <- cf[1, 4]
  }
  wmed_beta <- NA_real_; wmed_p <- NA_real_
  if (k >= 3) {
    o <- order(b); ws <- w[o]; bs <- b[o]
    wmed_beta <- bs[which.min(abs(cumsum(ws) / sum(ws) - 0.5))]
    set.seed(20260831)
    boot <- replicate(2000, {
      idx <- sample.int(k, k, replace = TRUE)
      oo <- order(b[idx]); ww <- w[idx][oo]; bb <- b[idx][oo]
      bb[which.min(abs(cumsum(ww) / sum(ww) - 0.5))]
    })
    wmed_se <- sd(boot, na.rm = TRUE)
    wmed_p <- 2 * pnorm(-abs(wmed_beta / wmed_se))
  }
  data.table(method = "IVW", beta = beta_use, se = se_use, p = 2 * pnorm(-abs(beta_use / se_use)),
             Q = Q, Qp = Qp, egger_int = egger_int, egger_p = egger_p,
             wmed_beta = wmed_beta, wmed_p = wmed_p)
}

## ============ 2. 载入数据 ============
cat("[load] meQTL / GWAS / LD\n")
mq <- fread("results/_tmp_smr/mqtl_full7.csv", header = FALSE,
            col.names = c("cpg","snp","beta_m0","se_m0","p_m","allele1","allele2","freq_a1","cistrans","clumped"))
mq[, `:=`(chr = sub("^chr([0-9XY]+):.*", "\\1", snp),
          pos = as.integer(sub("^chr[0-9XY]+:(\\d+):.*", "\\1", snp)))]
mq <- mq[cistrans == "TRUE" & p_m < PVAL_IV]
cat(sprintf("  meQTL (cis, p<%.0e): %d 行; clumped=TRUE: %d\n", PVAL_IV, nrow(mq), sum(mq$clumped == "TRUE")))

ldt <- fread("results/_tmp_smr/ld_1000g_eur.csv.gz")
ldpos <- fread("results/_tmp_smr/ld_1000g_eur_positions.csv")
liftmap <- fread("results/_tmp_smr/lift_positions.csv")
ldt[, chr := as.character(chr)]
ldpos[, chr := as.character(chr)]
liftmap[, chr := as.character(chr)]
cat(sprintf("  LD 表: %d 对; 1000G 匹配位点: %d; 回映位点: %d\n", nrow(ldt), nrow(ldpos), nrow(liftmap)))

fg_pe <- fread("results/_tmp_smr/finngen_pe_window.tsv")
fg_gh <- fread("results/_tmp_smr/finngen_gh_window.tsv")
setnames(fg_pe, c("pos","rsids","pval","beta","sebeta","af_alt"), c("pos38","rsid","p_gw","beta_gw","se_gw","af"))
setnames(fg_gh, c("pos","rsids","pval","beta","sebeta","af_alt"), c("pos38","rsid","p_gw","beta_gw","se_gw","af"))
fg_pe <- fg_pe[, .(chr = as.character(chrom), pos38, rsid, p_gw, beta_gw, se_gw, af, ref, alt)]
fg_gh <- fg_gh[, .(chr = as.character(chrom), pos38, rsid, p_gw, beta_gw, se_gw, af, ref, alt)]

ty <- fread("data/gwas/tyrmi2023_PE_maternal.tsv.gz")
cat("[load] tyrmi cols:", paste(names(ty), collapse = ","), "\n")
col_v <- function(pat) names(ty)[grepl(pat, names(ty), ignore.case = TRUE)]
c_chr <- col_v("^chromosome$")[1]
c_pos <- col_v("base_pair_location|^position$")[1]
c_ea  <- col_v("^effect_allele$")[1]; c_oa <- col_v("other_allele")[1]
c_ef  <- col_v("effect_allele_frequency")[1]
c_b   <- col_v("^beta$")[1]; c_se <- col_v("standard_error|^se$")[1]; c_p <- col_v("^p_value$")[1]
stopifnot(!is.na(c_chr), !is.na(c_pos), !is.na(c_ea), !is.na(c_oa),
          !is.na(c_ef), !is.na(c_b), !is.na(c_se), !is.na(c_p))
setnames(ty, c(c_chr, c_pos, c_ea, c_oa, c_ef, c_b, c_se, c_p),
         c("chrom","pos","ea","oa","ef","beta_gw","se_gw","p_gw"))
ty <- ty[, .(chr = as.character(chrom), pos = as.integer(pos), ea, oa, ef, beta_gw, se_gw, p_gw)]
ty <- ty[chr %in% c("3","6","17")]   # 只留目标染色体 (全表 ~8M 行, 减内存提速)
cat(sprintf("  tyrmi 目标染色体行数: %d\n", nrow(ty)))

## ============ 3. 主循环 ============
wald_rows <- list(); ivw_rows <- list(); loo_rows <- list(); heidi_rows <- list()
for (L in loci) {
  cat(sprintf("\n===== region %s chr%s (%s, %s) =====\n", L$region_id, L$chr, L$symbol, L$direction))
  sub <- mq[cpg %in% L$cpgs & chr == L$chr]
  cat(sprintf("  IV 池: %d 行 / %d CpG\n", nrow(sub), uniqueN(sub$cpg)))
  lf <- liftmap[chr == L$chr, .(chr, pos, pos38)]
  if (nrow(lf) == 0) { cat("  回映文件无该染色体, 跳过\n"); next }
  cat(sprintf("  回映 GRCh38: %d 位点\n", nrow(lf)))

  for (cg in L$cpgs) {
    sc <- sub[cpg == cg]
    sc38 <- merge(sc, lf, by = c("chr","pos"), all = FALSE)
    cat(sprintf("  --- %s: IV 池 %d, 回映 %d ---\n", cg, nrow(sc), nrow(sc38)))
    if (nrow(sc) < 2) next

    for (ph in names(GWAS_N)) {
      if (ph %in% c("PE_FinnGen","GH_FinnGen")) {
        fgx <- if (ph == "PE_FinnGen") fg_pe[chr == L$chr] else fg_gh[chr == L$chr]
        fgx <- fgx[, .(pos38, rsid, p_gw, beta_gw, se_gw, af, ref, alt)]
        j <- merge(sc38, fgx, by = "pos38", all = FALSE, allow.cartesian = TRUE)
        if (nrow(j) == 0) next
        sg <- harmonize(j$allele1, j$allele2, j$freq_a1, j$alt, j$ref, j$af)
        j[, sign_m := sg]
        j <- j[!is.na(sign_m)]
        j[, `:=`(beta_m = beta_m0 * sign_m, rsid = sub(",.*$", "", rsid))]
      } else {
        ## Tyrmi (GRCh37): 必须限定染色体, 否则同一 pos 跨染色体错误匹配
        j <- merge(sc, ty[chr == L$chr], by = "pos", all = FALSE, allow.cartesian = TRUE)
        if (nrow(j) == 0) next
        sg <- harmonize(j$allele1, j$allele2, j$freq_a1, j$ea, j$oa, j$ef)
        j[, sign_m := sg]
        j <- j[!is.na(sign_m)]
        j[, `:=`(beta_m = beta_m0 * sign_m, rsid = NA_character_, pos38 = NA_integer_)]
      }
      j <- j[!is.na(beta_m) & !is.na(se_gw) & se_gw > 0 & !is.na(se_m0) & se_m0 > 0]
      j <- j[j[, .I[which.min(p_gw)], by = pos]$V1]
      if (nrow(j) < 2) { cat(sprintf("    [%s] 定向后 SNP 不足\n", ph)); next }

      ## ---- LD clumping ----
      in1k <- ldpos[chr == L$chr, pos]
      iv <- clump_ivs(copy(j), ldt, L$chr, in1k)
      if (nrow(iv) < 1) next
      iv[, `:=`(b_zy = beta_gw / beta_m,
                se_zy = sqrt(se_gw^2 / beta_m^2 + beta_gw^2 * se_m0^2 / beta_m^4),
                F = (beta_m / se_m0)^2)]
      cat(sprintf("    [%s] 定向 %d / clump %d IV (meanF=%.1f)\n", ph, nrow(j), nrow(iv), mean(iv$F)))

      ## ---- MR 估计 ----
      est <- mr_estimates(iv$b_zy, iv$se_zy, iv$beta_m)
      est[, `:=`(region_id = L$region_id, symbol = L$symbol, direction = L$direction,
                 cpg = cg, pheno = ph, k_iv = nrow(iv),
                 mean_F = mean(iv$F), min_F = min(iv$F))]
      ivw_rows[[length(ivw_rows) + 1]] <- est
      wald_rows[[length(wald_rows) + 1]] <- iv[, .(region_id = L$region_id, symbol = L$symbol, cpg = cg, pheno = ph,
                                                   rsid, pos37 = pos, pos38,
                                                   beta_m, se_m = se_m0, p_m, F, beta_gw, se_gw, p_gw, b_zy, se_zy)]

      ## ---- leave-one-out ----
      for (i2 in seq_len(nrow(iv))) {
        w <- 1 / iv$se_zy[-i2]^2
        b <- sum(w * iv$b_zy[-i2]) / sum(w); se2 <- sqrt(1 / sum(w))
        loo_rows[[length(loo_rows) + 1]] <- data.table(
          region_id = L$region_id, cpg = cg, pheno = ph, left_out_rsid = iv$rsid[i2],
          n_left = nrow(iv) - 1, beta = b, se = se2, p = 2 * pnorm(-abs(b / se2)))
      }

      ## ---- HEIDI 近似: top IV 的 LD 伙伴 (r2 0.05-0.9) 的 Wald 比值同质性 ----
      top <- iv[which.min(p_m)]
      tp <- top$pos
      part <- ldt[chr == L$chr & (pos1 == tp | pos2 == tp)]
      partner_pos <- unique(fifelse(part$pos1 == tp, part$pos2, part$pos1)[part$r2 >= R2_HEIDI_LO & part$r2 <= R2_HEIDI_HI])
      hs <- j[pos %in% partner_pos]
      if (nrow(hs) >= 5) {
        hs[, `:=`(b_zy = beta_gw / beta_m,
                  se_zy = sqrt(se_gw^2 / beta_m^2 + beta_gw^2 * se_m0^2 / beta_m^4))]
        d <- hs$b_zy - top$b_zy
        v <- hs$se_zy^2 + top$se_zy^2
        chi <- sum(d^2 / v); df <- length(d)
        heidi_rows[[length(heidi_rows) + 1]] <- data.table(
          region_id = L$region_id, symbol = L$symbol, cpg = cg, pheno = ph,
          top_rsid = top$rsid, n_ld_snps = length(d),
          chi2 = chi, df = df, p_heidi = pchisq(chi, df, lower.tail = FALSE))
        cat(sprintf("    [%s] HEIDI-近似: %d LD SNP, chi2=%.1f df=%d p=%.3f\n",
                    ph, length(d), chi, df, pchisq(chi, df, lower.tail = FALSE)))
      }
    }
  }
}

## ============ 4. 汇总输出 ============
if (length(ivw_rows) == 0) stop("无任何 MR 分析完成 — 检查输入 (LD 表 / 回映 / GWAS 合并)")
ivw <- rbindlist(ivw_rows)
ivw[, p_bh := p.adjust(p, "BH"), by = pheno]
ivw[, OR := exp(beta)]
ivw[, OR_lo := exp(beta - 1.96 * se)]
ivw[, OR_hi := exp(beta + 1.96 * se)]
setorder(ivw, region_id, pheno, cpg)
fwrite(ivw, "results/GSE282512_mr_ivw_results.csv")
fwrite(rbindlist(wald_rows), "results/GSE282512_mr_wald_snps.csv")
if (length(loo_rows)) fwrite(rbindlist(loo_rows), "results/GSE282512_mr_leaveoneout.csv")
if (length(heidi_rows)) fwrite(rbindlist(heidi_rows), "results/GSE282512_mr_heidi.csv")

## ---- 森林图 (base R) ----
png("figures/GSE282512_mr_forest.png", width = 1600, height = 2000, res = 150)
par(mar = c(5, 14, 3, 2), cex = 0.9)
ph_lab <- c(PE_FinnGen = "FinnGen PE", GH_FinnGen = "FinnGen GH", PE_Tyrmi = "Tyrmi meta PE")
ord <- ivw[order(region_id, pheno, cpg)]
xr <- range(c(0.5, ord$OR_lo, ord$OR_hi), na.rm = TRUE)
plot(NA, xlim = xr, ylim = c(0.5, nrow(ord) + 0.5), yaxt = "n",
     xlab = "OR per 1-SD methylation", ylab = "", main = "cis-meQTL MR (IVW, LD-clumped)")
abline(v = 1, lty = 2, col = "grey50")
yy <- nrow(ord):1
for (i in seq_len(nrow(ord))) {
  col <- if (ord$p[i] < 0.05) "firebrick" else "grey30"
  segments(ord$OR_lo[i], yy[i], ord$OR_hi[i], yy[i], col = col, lwd = 2)
  points(ord$OR[i], yy[i], pch = 18, cex = 1.6, col = col)
}
axis(2, at = yy, labels = sprintf("%s | %s", ord$cpg, ph_lab[ord$pheno]), las = 1)
mtext(sprintf("beta=%+.3f p=%.2g k=%d F=%.0f", ord$beta, ord$p, ord$k_iv, ord$mean_F),
      side = 2, line = -0.4, cex = 0.5, las = 1, at = yy)
dev.off()

## ---- 汇总 txt ----
sink("results/GSE282512_mr_summary.txt")
cat("== MR 补强: 多 IV IVW + 敏感性 + 替代结局 ==\n\n")
cat(sprintf("暴露: GoDMC cis-meQTL (GRCh37, allele1=effect, N=%d); IV 池 p<%.0e; LD clump r2<%.1f (1000G phase3 EUR 本地计算, 58a)\n",
            N_MQTL, PVAL_IV, R2_CLUMP))
cat("结局: FinnGen R10 PE (7965/211852) 与 GH (9535/211957); Tyrmi 2023 JAMA Cardiol PE meta-GWAS (GCST90269903, 16743/280081)\n")
cat("方法: IVW (Q 显著则 DL 随机效应) + Wald ratio + MR-Egger 截距 + 加权中位数 (bootstrap CI) + leave-one-out + F 统计量\n")
cat("HEIDI 替代: (1) coloc PP4 (46 号脚本); (2) top IV LD 伙伴 (r2 0.05-0.9) Wald 比值同质性 chi2 (本脚本)\n\n")
cat("---- IVW 主结果 (OR per 1-SD methylation) ----\n")
print(ivw[, .(region = region_id, cpg, pheno, k_iv, OR = round(OR, 3),
              lo = round(OR_lo, 3), hi = round(OR_hi, 3), p = signif(p, 3),
              p_bh = signif(p_bh, 3), Qp = signif(Qp, 2), egger_p = signif(egger_p, 2),
              wmed_p = signif(wmed_p, 2), mean_F = round(mean_F, 1))])
cat("\n---- 说明 ----\n")
cat("- 与 46 号单 lead-SNP SMR 对比: SLC17A1 cg19490609/cg25753631 与 84604 cg15445000 在 FinnGen PE 中曾有名义显著 (p~0.002-0.014)\n")
cat("- 等位基因定向后 Wald 比值方向首次可判读 (此前 GoDMC 精简文件无等位基因列, 仅 T 统计量无方向)\n")
cat("- Tyrmi meta 含 FinnGen R6 (4285 例), 与 FinnGen R10 结局部分样本重叠; 判读为三角验证而非完全独立复制\n")
cat("- Egger 截距 p<0.05 提示定向多效性; Qp<0.05 提示 IV 间异质性; wmed_p 为加权中位数稳健估计\n")
cat("- HEIDI 近似 p<0.05 提示 top 变异外存在异质信号 (多因果变异或多效性)\n")
sink()

cat("\n[完成] 输出: GSE282512_mr_ivw_results.csv / _mr_wald_snps.csv / _mr_leaveoneout.csv / _mr_heidi.csv / _mr_summary.txt / figures/GSE282512_mr_forest.png\n")
