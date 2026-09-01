# 20_build_wgbs_metadata.R ------------------------------------------------
# GSE282512 (plasma cfDNA WGBS, 纵向队列) 样本元数据整合
# 输入: data/geo_methylation/GSE282512_sample_annot.csv.gz      (369 行)
#       data/geo_methylation/GSE282512_series_matrix.txt.gz     (GSM 映射)
# 输出: results/GSE282512_samples_clean.csv
#       控制台/日志: 队列结构汇总 + 患者级防泄漏提示
# -------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  library(data.table)
}))

proj_dir <- "E:/妊娠期高血压甲基化研究方案/hdp-methylation-project"
out_csv  <- file.path(proj_dir, "results", "GSE282512_samples_clean.csv")

## 1. 样本注释表 -----------------------------------------------------------
annot <- fread(file.path(proj_dir, "data/geo_methylation",
                         "GSE282512_sample_annot.csv.gz"))
message(sprintf("注释表: %d 样本, %d 患者", nrow(annot), uniqueN(annot$Patient)))

## 2. series_matrix 提取 GSM <-> Sample_id 映射 ----------------------------
sm_path <- file.path(proj_dir, "data/geo_methylation",
                     "GSE282512_series_matrix.txt.gz")
title_line <- geo_line <- NULL
con <- gzfile(sm_path, "rt")
repeat {
  ln <- readLines(con, n = 1, warn = FALSE)
  if (length(ln) == 0) break
  if (grepl("^!Sample_title", ln))            title_line <- ln
  if (grepl("^!Sample_geo_accession", ln))    geo_line   <- ln
  if (grepl("^!Series_matrix_table_begin", ln)) break
}
close(con)

parse_fields <- function(line) {
  s <- sub("^[^\t]+\t", "", line)
  unlist(strsplit(s, "\t"))
}
titles <- gsub('"', "", parse_fields(title_line))
gsms   <- gsub('"', "", parse_fields(geo_line))
stopifnot(length(titles) == length(gsms), length(gsms) == nrow(annot))

# title 形如 "cell-free DNA, Healthy pregnancy DNA031134" / "Preeclampsia pregnancy DNA03xxxx"
sample_ids <- sub(".*\\s(DNA\\d+)\\s*$", "\\1", titles)
sm_disease <- ifelse(grepl("Healthy", titles), "Healthy", "Preeclampsia")
if (any(!grepl("^DNA", sample_ids))) {
  warning("有 title 未解析出 DNA 编号, 前3个: ",
          paste(head(titles[!grepl("^DNA", sample_ids)], 3), collapse = " | "))
}

## 3. 合并元数据表 ---------------------------------------------------------
meta <- data.table(
  gsm           = gsms,
  sample_id     = sample_ids,
  sm_disease    = sm_disease,
  cov_file      = paste0(gsms, ".cov.gz"),
  patient       = annot$Patient[match(sample_ids, annot$Sample_id)],
  batch         = annot$Batch[match(sample_ids, annot$Sample_id)],
  tube_type     = annot$Tube_type[match(sample_ids, annot$Sample_id)],
  category      = annot$Category[match(sample_ids, annot$Sample_id)],
  severity      = annot$Severity[match(sample_ids, annot$Sample_id)],
  onset         = annot$Onset[match(sample_ids, annot$Sample_id)],
  ga_weeks      = annot$Gest_age_weeks[match(sample_ids, annot$Sample_id)],
  ga_days       = annot$Gest_age_days[match(sample_ids, annot$Sample_id)],
  ga_sampling   = annot$Gest_age_weeks[match(sample_ids, annot$Sample_id)] +
                  annot$Gest_age_days[match(sample_ids, annot$Sample_id)] / 7,
  ga_birth      = annot$Gest_age_birth[match(sample_ids, annot$Sample_id)],
  replicate     = annot$Replicate[match(sample_ids, annot$Sample_id)],
  aspirin       = annot$Aspirin[match(sample_ids, annot$Sample_id)],
  gravida       = annot$Gravida[match(sample_ids, annot$Sample_id)],
  fertility     = annot$Fertility_treatment[match(sample_ids, annot$Sample_id)],
  bmi           = as.numeric(annot$BMI[match(sample_ids, annot$Sample_id)]),
  race          = annot$Race[match(sample_ids, annot$Sample_id)],
  birth_term    = annot$Birth_term_category[match(sample_ids, annot$Sample_id)],
  fetal_sex     = annot$Sex[match(sample_ids, annot$Sample_id)],
  iugr          = annot$IUGR[match(sample_ids, annot$Sample_id)],
  group_status  = annot$Group[match(sample_ids, annot$Sample_id)]
)

stopifnot(!any(is.na(meta$patient)))   # 注释表全覆盖

## 4. 一致性校验: 注释表 Category vs series_matrix 疾病状态 ---------------
# 注: Healthy==Control, Preeclampsia==PE (命名差异非真不一致), 先归一再比对
norm_cat  <- ifelse(meta$sm_disease == "Healthy", "Control", "PE")
mism <- meta[tolower(norm_cat) != tolower(category)]
message(sprintf("疾病状态一致性: series_matrix vs 注释表 不一致 %d 条", nrow(mism)))
if (nrow(mism)) print(mism[, .(gsm, sample_id, category, sm_disease)])

## 5. 队列结构汇总 ---------------------------------------------------------
message("\n===== 队列结构汇总 =====")
message(sprintf("样本总数: %d | 患者总数: %d", nrow(meta), uniqueN(meta$patient)))
message(sprintf("疾病分组: %s", paste(
  sprintf("%s=%d", names(table(meta$category)), table(meta$category)),
  collapse = ", ")))
pe <- meta[category == "PE"]
message(sprintf("PE 细分: severity %s | onset %s",
  paste(sprintf("%s=%d", names(table(pe$severity)), table(pe$severity)), collapse = ", "),
  paste(sprintf("%s=%d", names(table(pe$onset)),   table(pe$onset)),   collapse = ", ")))
n_pt <- meta[, .(n_samples = .N), by = patient]
message(sprintf("每患者样本数: %d 名患者, 中位 %d 个样本/人 (范围 %d-%d)",
  nrow(n_pt), median(n_pt$n_samples), min(n_pt$n_samples), max(n_pt$n_samples)))
message(sprintf("采样管: %s", paste(
  sprintf("%s=%d", names(table(meta$tube_type)), table(meta$tube_type)),
  collapse = ", ")))
message(sprintf("批次: %s | 技术重复: %d", paste(
  sprintf("batch%s=%d", names(table(meta$batch)), table(meta$batch)), collapse = ", "),
  sum(meta$replicate == TRUE)))

# 病例-对照的可配对纵向采样窗口
message("\n病例 vs 对照的孕周分布 (采样时):")
print(meta[, .(n = .N, ga_min = min(ga_weeks), ga_max = max(ga_weeks),
               ga_med = median(ga_weeks)), by = category])

# 防泄漏核心提示: 患者级统计
both <- meta[, .(has_case = any(category == "PE"),
                 has_ctrl = any(category == "Control")), by = patient]
message(sprintf("\n[防泄漏] 同一患者跨疾病分组的: %d (应为 0, 纵向设计同一患者状态固定)",
  sum(both$has_case & both$has_ctrl)))
message("[防泄漏] CV / train-test 划分必须按 patient 列分层, 严禁按样本随机划分。")

## 6. 输出 -----------------------------------------------------------------
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
fwrite(meta, out_csv)
message(sprintf("\n输出: %s (%d 行 x %d 列)", out_csv, nrow(meta), ncol(meta)))
