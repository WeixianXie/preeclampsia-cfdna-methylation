# 22_select_subcohort.R ---------------------------------------------------
# GSE282512 cfDNA WGBS DMR 分析子队列筛选
# 设计原则:
#   1) QC 过滤: n_cpg_std >= 1e6 且 cov5x_pct >= 50
#   2) 防泄漏: 主对比每患者只取 1 个样本 (取采血孕周最晚者, 贴近病例发病时点)
#   3) 孕周匹配: 病例-对照按采样孕周 ±3 周贪心配对 (GA 是血浆 cfDNA 甲基化最强混杂)
#   4) 结构均衡: 记录 tube_type / batch 分布, 供 DMRcate 模型纳入协变量
# 输出: results/GSE282512_subcohort.csv (含 keep 标志与入选原因)
#       results/GSE282512_subcohort_summary.txt
# -------------------------------------------------------------------------

suppressWarnings(suppressMessages(library(data.table)))

# proj_dir: 改用相对路径(Windows R 无法 fread 中文绝对路径)
qc   <- fread("results/GSE282512_qc_summary.csv")
meta <- fread("results/GSE282512_samples_clean.csv")
rep  <- "results/GSE282512_subcohort_summary.txt"

MIN_CPG  <- 1e6    # 标准染色体 CpG 数下限
MIN_COV5 <- 50     # cov>=5x 位点百分比下限
GA_TOL   <- 3      # 孕周匹配容差
GA_MIN   <- 16     # 分析窗口下限 (周): PE 诊断大多 >=20w, 排除早孕对照
GA_MAX   <- 37

## 1. QC 过滤 ---------------------------------------------------------------
dt <- merge(meta, qc[, .(gsm, n_cpg_std, mean_depth, cov5x_pct, cov10x_pct,
                         mean_beta_std)],
            by = "gsm", all.x = TRUE, suffixes = c("", ".qc"))
dt[, qc_pass := !is.na(n_cpg_std) & n_cpg_std >= MIN_CPG & cov5x_pct >= MIN_COV5]
dt[, reason := fcase(
  is.na(n_cpg_std),       "qc_missing",
  n_cpg_std < MIN_CPG,    "qc_low_cpg",
  cov5x_pct < MIN_COV5,   "qc_low_cov5",
  default = "qc_pass")]

## 2. 病例选择: PE 患者, 每人取孕周最晚样本 --------------------------------
pe <- dt[category == "PE" & qc_pass == TRUE]
pe[, ga_weeks := as.integer(ga_weeks)]
pe <- pe[ga_weeks >= GA_MIN & ga_weeks <= GA_MAX]
pe_pick <- pe[order(patient, -ga_weeks, -ga_days)][, .SD[1], by = patient]

## 3. 对照选择: Control 患者, 每人取孕周最晚样本, 再按 GA ±3 周配对 -------
ctrl <- dt[category == "Control" & qc_pass == TRUE]
ctrl[, ga_weeks := as.integer(ga_weeks)]
ctrl <- ctrl[ga_weeks >= GA_MIN & ga_weeks <= GA_MAX]
ctrl_pick <- ctrl[order(patient, -ga_weeks, -ga_days)][, .SD[1], by = patient]

# 贪心配对: 每个病例找 GA 差最小且未被用过的对照
ctrl_avail <- copy(ctrl_pick)
res <- lapply(seq_len(nrow(pe_pick)), function(i) {
  if (!nrow(ctrl_avail)) return(NULL)
  ga_i <- pe_pick$ga_weeks[i]
  j <- which.min(abs(ctrl_avail$ga_weeks - ga_i))
  if (abs(ctrl_avail$ga_weeks[j] - ga_i) > GA_TOL) return(NULL)
  out <- copy(pe_pick[i]); out[, matched_ga_ctrl := ctrl_avail$ga_weeks[j]]
  ctl <- copy(ctrl_avail[j])
  ctrl_avail <<- ctrl_avail[-j]
  list(case = out, ctrl = ctl)
})
res <- res[!vapply(res, is.null, logical(1))]
cases <- rbindlist(lapply(res, `[[`, "case"))
ctrls <- rbindlist(lapply(res, `[[`, "ctrl"))

## 4. 汇总输出 ---------------------------------------------------------------
sink(rep, split = TRUE)
cat("===== GSE282512 DMR 子队列筛选汇总 =====\n\n")
cat(sprintf("QC 阈值: n_cpg_std>=%g, cov5x_pct>=%g%%\n", MIN_CPG, MIN_COV5))
cat(sprintf("QC: 通过 %d / 剔除 %d (低CpG %d, 低cov5x %d, 缺失 %d)\n",
  sum(dt$qc_pass), sum(!dt$qc_pass),
  sum(dt$reason == "qc_low_cpg"), sum(dt$reason == "qc_low_cov5"),
  sum(dt$reason == "qc_missing")))
cat(sprintf("分析孕周窗口: %d–%d 周, 匹配容差 ±%d 周\n\n", GA_MIN, GA_MAX, GA_TOL))
cat(sprintf("病例 (PE, 每患者1样本): %d 例\n", nrow(cases)))
cat(sprintf("  severity: %s\n", paste(sprintf("%s=%d",
  names(table(cases$severity)), table(cases$severity)), collapse = ", ")))
cat(sprintf("  onset:    %s\n", paste(sprintf("%s=%d",
  names(table(cases$onset)), table(cases$onset)), collapse = ", ")))
cat(sprintf("  GA: 中位 %d 周 (%d–%d)\n", median(cases$ga_weeks),
  min(cases$ga_weeks), max(cases$ga_weeks)))
cat(sprintf("配对成功对照: %d 例 (对照池 %d)\n", nrow(ctrls), nrow(ctrl_pick)))
cat(sprintf("  GA: 中位 %d 周 (%d–%d)\n", median(ctrls$ga_weeks),
  min(ctrls$ga_weeks), max(ctrls$ga_weeks)))
cat(sprintf("GA 配对差绝对值: 中位 %.1f 周, 最大 %d 周\n",
  median(abs(cases$ga_weeks - cases$matched_ga_ctrl)),
  max(abs(cases$ga_weeks - cases$matched_ga_ctrl))))
cat(sprintf("\ntube_type (病例/对照): PAXgene %d/%d, EDTA %d/%d\n",
  sum(cases$tube_type == "PAXgene DNA"), sum(ctrls$tube_type == "PAXgene DNA"),
  sum(cases$tube_type == "EDTA"), sum(ctrls$tube_type == "EDTA")))
cat("batch 分布 (病例): ", paste(sprintf("b%s=%d",
  names(table(cases$batch)), table(cases$batch)), collapse = ", "), "\n")
cat("batch 分布 (对照): ", paste(sprintf("b%s=%d",
  names(table(ctrls$batch)), table(ctrls$batch)), collapse = ", "), "\n")
cat(sprintf("\nEOPE 敏感性分析子集: 病例 %d 例\n", sum(cases$onset == "EOPE")))
cat("[防泄漏] 全部按 patient 去重, 无同一患者跨组/重复样本。\n")
sink()

## 5. 输出主表 ---------------------------------------------------------------
sub <- rbind(
  cbind(cases[, .(gsm, sample_id, patient, group = "PE", category, severity,
                  onset, ga_weeks, tube_type, batch, n_cpg_std, mean_depth,
                  cov5x_pct, cov10x_pct, mean_beta_std)],
        matched_ga = cases$matched_ga_ctrl),
  cbind(ctrls[, .(gsm, sample_id, patient, group = "Control", category,
                  severity, onset, ga_weeks, tube_type, batch, n_cpg_std,
                  mean_depth, cov5x_pct, cov10x_pct, mean_beta_std)],
        matched_ga = NA_integer_))
out_csv <- "results/GSE282512_subcohort.csv"
fwrite(sub, out_csv)
message(sprintf("输出: %s (%d 行: PE %d + Control %d)", out_csv, nrow(sub),
  sum(sub$group == "PE"), sum(sub$group == "Control")))
message(sprintf("汇总: %s", rep))
