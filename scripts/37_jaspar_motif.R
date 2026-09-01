# 37_jaspar_motif.R —— 机制层 ④: JASPAR2024 CORE 人 TF motif 富集
# 候选166 DMR (含 hyper/hypo 分层) vs 长度+GC 匹配背景
# 双链 vcountPWM (min.score 80%), Fisher + BH
# 输出: results/GSE282512_jaspar_motif_enrichment.csv, GSE282512_jaspar_summary.txt,
#        figures/jaspar_top_tf.png

suppressPackageStartupMessages({
  library(data.table); library(Biostrings)
})
set.seed(1)

MIN_SCORE <- "80%"

## ============ 1. 读取序列 ============
read_seqs <- function(file_pref, ids) {
  fs <- sprintf("results/_tmp_seq/%s_%s.txt", file_pref, ids)
  ok <- file.exists(fs)
  seqs <- character(sum(ok))
  seqs <- vapply(fs[ok], function(f) gsub("\\s", "", readChar(f, file.info(f)$size)), character(1))
  list(ids = ids[ok], seqs = DNAStringSet(toupper(seqs)))
}
cand_ids <- fread("results/_tmp_cand_regions.tsv", header = TRUE)$region_id
bg_ids <- fread("results/_tmp_bg_regions.tsv", header = TRUE)$region_id
CA <- read_seqs("cand", as.character(cand_ids))
BG <- read_seqs("bg", as.character(bg_ids))

## 超长区域居中截断至 5kb (总序列 105MB 不可全扫; 截断后富集密度不被超长稀释)
MAXW <- 5000
cap_seq <- function(ss) {
  w <- width(ss)
  too <- w > MAXW
  if (any(too)) {
    st <- (w[too] - MAXW) %/% 2 + 1
    ss[too] <- narrow(ss[too], start = st, width = MAXW)
  }
  ss
}
CA$seqs <- cap_seq(CA$seqs); BG$seqs <- cap_seq(BG$seqs)
cat(sprintf("序列读取: 候选 %d (%.1f Mb), 背景 %d (%.1f Mb)\n",
            length(CA$seqs), sum(width(CA$seqs)) / 1e6,
            length(BG$seqs), sum(width(BG$seqs)) / 1e6))

gc_c <- letterFrequency(CA$seqs, "GC", as.prob = TRUE)[, 1]
gc_b <- letterFrequency(BG$seqs, "GC", as.prob = TRUE)[, 1]
cat(sprintf("GC: 候选 mean %.3f [%.3f-%.3f], 背景 mean %.3f\n",
            mean(gc_c), min(gc_c), max(gc_c), mean(gc_b)))

## GC 十分位匹配背景 (按候选 GC 分布逐层抽样)
qg <- quantile(gc_c, probs = seq(0, 1, 0.1), na.rm = TRUE)
st_c <- cut(gc_c, breaks = unique(qg), include.lowest = TRUE)
st_b <- cut(gc_b, breaks = unique(qg), include.lowest = TRUE)
keep <- rep(FALSE, length(gc_b))
for (lv in levels(st_c)) {
  idx <- which(st_b == lv)
  nc <- sum(st_c == lv, na.rm = TRUE)
  if (length(idx) > 0 && nc > 0) keep[sample(idx, min(length(idx), nc * 15))] <- TRUE
}
BG <- list(ids = BG$ids[keep], seqs = BG$seqs[keep])
gc_b <- gc_b[keep]
cat(sprintf("GC 匹配后背景: %d (mean GC %.3f)\n", length(BG$seqs), mean(gc_b)))

## ============ 2. 解析 JASPAR ============
jl <- readLines("data/annot/jaspar2024_human_core.jaspar")
hdr_idx <- grep("^>", jl)
motifs <- vector("list", length(hdr_idx))
for (i in seq_along(hdr_idx)) {
  h <- jl[hdr_idx[i]]
  mid <- sub("^>", "", sub("\\t.*", "", h))
  tf  <- sub("^[^\t]*\t", "", h)
  body <- jl[(hdr_idx[i] + 1):(hdr_idx[i] + 4)]
  rows <- c("A", "C", "G", "T")
  lst <- lapply(body, function(l) {
    v <- unlist(strsplit(gsub("[^0-9 ]", " ", l), "\\s+"))
    as.numeric(v[v != ""])
  })
  ncol_ <- unique(lengths(lst))
  if (length(ncol_) != 1 || ncol_ < 4 || any(lengths(lst) != ncol_)) next
  m <- matrix(unlist(lst), nrow = 4, byrow = TRUE,
              dimnames = list(rows, NULL))
  if (any(is.na(m))) next
  ## 958/1661 SELEX 来源 PFM 列和不齐, PWM() 拒绝 -> 手工建 PSSM (列内归一 + 0.01 伪计)
  pp <- t(t(m + 0.01) / (colSums(m) + 0.04))
  motifs[[i]] <- list(id = mid, tf = tf, pwm = log2(pp / 0.25))
}
motifs <- Filter(Negate(is.null), motifs)
cat(sprintf("JASPAR 有效 PWM: %d\n", length(motifs)))

## ============ 3. 双链扫描 (拼接单串 + N padding + 坐标回映射) ============
v2 <- fread("results/GSE282512_dmr_final_v2.csv",
            select = c("region_id", "direction", "symbol"))
hypo_ids <- as.character(v2[direction == "hypo"]$region_id)
hyper_ids <- as.character(v2[direction == "hyper"]$region_id)
is_hypo <- CA$ids %in% hypo_ids
is_hyper <- CA$ids %in% hyper_ids

all_seq <- c(CA$seqs, BG$seqs)
n_ca <- length(CA$seqs); n_all <- length(all_seq)
PAD <- 40L
comb <- DNAString(paste(as.character(all_seq), collapse = strrep("N", PAD)))
rc <- reverseComplement(comb)
L_tot <- length(comb)   # DNAString 单串用 length() 而非 width()
starts <- c(1L, head(cumsum(width(all_seq) + PAD), -1) + 1L)
starts <- as.integer(starts)
## 注意 starts[1]=1, starts[i] 为第 i 条序列起点
cat(sprintf("拼接总长: %.1f Mb × 2 链\n", L_tot / 1e6))

scan_both <- function(pwm) {
  vf <- matchPWM(pwm, comb, min.score = MIN_SCORE)
  vr <- matchPWM(pwm, rc,   min.score = MIN_SCORE)
  pf <- start(vf)
  pr <- L_tot - start(vr) + 1L   # 回映到正向坐标
  idx <- findInterval(c(pf, pr), starts)
  tabulate(idx, nbins = n_all)  # idx 从 1 开始恰对应序列序号
}

## 一次扫描, 两张表共用
cnts <- lapply(seq_along(motifs), function(i) {
  if (i %% 200 == 0) cat(sprintf("  scanning %d/%d...\n", i, length(motifs)))
  scan_both(motifs[[i]]$pwm)
})
res <- rbindlist(lapply(seq_along(motifs), function(i) {
  mo <- motifs[[i]]
  cnt <- cnts[[i]]
  ch <- cnt[seq_len(n_ca)]; bh <- cnt[(n_ca + 1):n_all]
  a <- sum(ch > 0); b <- sum(bh > 0)
  data.table(motif_id = mo$id, tf = mo$tf,
             cand_hit = a, cand_n = length(ch),
             hypo_hit = sum(ch[is_hypo] > 0), hypo_n = sum(is_hypo),
             hyper_hit = sum(ch[is_hyper] > 0), hyper_n = sum(is_hyper),
             bg_hit = b, bg_n = length(bh))
}))
cat(sprintf("全候选扫描完成: %d motifs\n", nrow(res)))
res[, or := (cand_hit / (cand_n - cand_hit)) / (bg_hit / (bg_n - bg_hit))]
res[, p_fisher := mapply(function(a, n1, b, n2) {
  fisher.test(matrix(c(a, n1 - a, b, n2 - b), nrow = 2))$p.value
}, cand_hit, cand_n, bg_hit, bg_n)]
res[, p_bh := p.adjust(p_fisher, "BH")]
setorder(res, p_fisher)
fwrite(res, "results/GSE282512_jaspar_motif_enrichment.csv")

## hypo 单独检验 (hypo DMR = 来源细胞中处于活跃/开放状态的元件)
res_h <- rbindlist(lapply(seq_along(motifs), function(i) {
  mo <- motifs[[i]]
  cnt <- cnts[[i]]
  ch <- cnt[seq_len(n_ca)][is_hypo]; bh <- cnt[(n_ca + 1):n_all]
  a <- sum(ch > 0); b <- sum(bh > 0)
  data.table(motif_id = mo$id, tf = mo$tf, hypo_hit = a, hypo_n = length(ch),
             bg_hit = b, bg_n = length(bh))
}))
cat(sprintf("hypo 子集表完成\n"))
res_h[, or := (hypo_hit / (hypo_n - hypo_hit)) / (bg_hit / (bg_n - bg_hit))]
res_h[, p_fisher := mapply(function(a, n1, b, n2) {
  fisher.test(matrix(c(a, n1 - a, b, n2 - b), nrow = 2))$p.value
}, hypo_hit, hypo_n, bg_hit, bg_n)]
res_h[, p_bh := p.adjust(p_fisher, "BH")]
setorder(res_h, p_fisher)
fwrite(res_h, "results/GSE282512_jaspar_hypo.csv")

## ============ 4. 汇总输出 ============
sink("results/GSE282512_jaspar_summary.txt")
cat("===== 机制层 ④: JASPAR2024 CORE 人 TF motif 富集 =====\n\n")
cat(sprintf("候选 DMR %d (hyper %d / hypo %d) vs GC+长度匹配背景 %d\n",
            length(CA$seqs), sum(is_hyper), sum(is_hypo), length(BG$seqs)))
cat(sprintf("扫描: %d 个 PWM, 双链, min.score=%s\n\n", length(motifs), MIN_SCORE))
cat("---- 全部候选 Top 20 (Fisher p) ----\n")
print(head(res[, .(motif_id, tf, cand_hit, cand_n, bg_hit, bg_n,
                   or = round(or, 2), p_fisher, p_bh)], 20))
cat("\n---- BH<0.05 的 TF ----\n")
sig <- res[p_bh < 0.05]
if (nrow(sig) == 0) cat("无\n") else print(sig[, .(motif_id, tf, cand_hit,
  cand_n, bg_hit, bg_n, or = round(or, 2), p_fisher, p_bh)])
cat("\n---- hypo 子集 Top 20 ----\n")
print(head(res_h[, .(motif_id, tf, hypo_hit, hypo_n, bg_hit, bg_n,
                     or = round(or, 2), p_fisher, p_bh)], 20))
cat("\n---- hypo BH<0.05 ----\n")
sigh <- res_h[p_bh < 0.05]
if (nrow(sigh) == 0) cat("无\n") else print(sigh[, .(motif_id, tf, hypo_hit,
  hypo_n, bg_hit, bg_n, or = round(or, 2), p_fisher, p_bh)])
cat("\n谱系参考: SPI1/CEBP/RUNX/IRF=髓系; GATA1/TAL1/KLF1=红系; SPIB/ETS1=淋巴;\n")
cat("GATA2/GATA3/TEAD4/TP63/TFAP2C/AP-2gamma=滋养层/胎盘; FLI1/ERG=内皮/巨核.\n")
sink()

## ============ 5. 图 ============
png("figures/jaspar_top_tf.png", width = 2400, height = 2000, res = 300)
par(mar = c(4, 8, 3, 2))
tp <- head(res[or > 1][order(p_fisher)], 20)
if (nrow(tp) < 5) tp <- head(res[order(p_fisher)], 20)
bp <- barplot(-log10(tp$p_fisher), horiz = TRUE, las = 1,
              names.arg = sprintf("%s (%.1fx)", tp$tf, tp$or),
              col = ifelse(tp$p_bh < 0.05, "#c0392b", "#7f8c8d"),
              xlab = "-log10(Fisher p)", main = "JASPAR motif enrichment (166 DMR vs GC-matched bg)")
abline(v = -log10(0.05), lty = 2, col = "grey40")
legend("bottomright", c("BH<0.05", "n.s.(BH)"), fill = c("#c0392b", "#7f8c8d"), bty = "n")
dev.off()
cat("DONE\n")
