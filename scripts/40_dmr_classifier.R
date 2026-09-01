# 40_dmr_classifier.R — B: 166 候选 DMR 区域甲基化分类器性能评估 (任务 #38)
# 设计: LOOCV 弹性网 logistic (fold 内 NA 插补 + fold 内 lambda 选择);
#      标签置换零分布 (n=20); top-k 精简 panel (fold 内按 |t| 选择, 无泄漏);
#      EOPE/LOPE 分层; tube/batch 协变量 sanity 模型
suppressPackageStartupMessages({
  library(data.table); library(glmnet); library(pROC); library(ggplot2)
})
set.seed(42)

## ============ 1. 数据准备 ============
sub <- fread("results/GSE282512_subcohort.csv")
breg <- fread("results/GSE282512_region_beta.csv.gz")
breg[, region_id := as.character(region_id)]
bmat <- as.matrix(breg[, -1]); rownames(bmat) <- breg$region_id
v2 <- fread("results/GSE282512_dmr_final_v2.csv", select = c("region_id","direction","candidate"))
v2[, region_id := as.character(region_id)]
cand <- v2[isTRUE(candidate) | candidate == "TRUE"]$region_id
cand <- cand[cand %in% rownames(bmat)]
gs <- sub$gsm
X <- bmat[cand, gs, drop = FALSE]          # 特征 x 样本
y <- as.integer(sub$group == "PE")
tube <- factor(sub$tube_type); batch <- factor(sub$batch)
cat(sprintf("特征: %d DMR 区域 x %d 样本 (PE %d / CT %d); 特征 NA 率 %.1f%%\n",
            nrow(X), ncol(X), sum(y==1), sum(y==0), 100*mean(is.na(X))))

## ============ 2. LOOCV 弹性网函数 (fold 内插补 + lambda) ============
loocv_enet <- function(X, y, alpha = 0.5) {
  p <- nrow(X); n <- ncol(X); pr <- numeric(n)
  for (i in 1:n) {
    tr <- setdiff(1:n, i)
    Xi <- X[, tr, drop = FALSE]; xi <- X[, i]
    med <- apply(Xi, 1, median, na.rm = TRUE)     # fold 内行中位数
    med[is.na(med)] <- 0.5
    Xi[is.na(Xi)] <- med[which(is.na(Xi), arr.ind = TRUE)[, "row"]]
    xi[is.na(xi)] <- med
    cvfit <- suppressWarnings(cv.glmnet(t(Xi), factor(y[tr]), family = "binomial",
                                        alpha = alpha, nfolds = 5, type.measure = "auc"))
    pr[i] <- as.numeric(predict(cvfit, t(xi), s = "lambda.min", type = "response"))
  }
  pr
}

## ============ 3. 主模型: 全 166 DMR LOOCV ============
pr_full <- loocv_enet(X, y)
roc_full <- roc(y, pr_full, quiet = TRUE)
ci_full <- ci.auc(roc_full, method = "delong")
th <- coords(roc_full, "best", ret = "threshold", best.method = "youden")[[1]]
pred <- as.integer(pr_full >= th)
perf_full <- data.table(
  model = "elastic_net_166DMR", n = length(y), auc = as.numeric(auc(roc_full)),
  auc_lo = as.numeric(ci_full[1]), auc_hi = as.numeric(ci_full[3]),
  sens = mean(pred[y == 1] == 1), spec = mean(pred[y == 0] == 0),
  acc = mean(pred == y))
cat(sprintf("全模型 LOOCV AUC = %.3f [%.3f, %.3f]\n", perf_full$auc, perf_full$auc_lo, perf_full$auc_hi))

## ============ 4. 置换零分布 (20 次, 已有缓存则复用) ============
if (file.exists("results/GSE282512_classifier_perm_null.csv")) {
  perm_auc <- fread("results/GSE282512_classifier_perm_null.csv")$auc
  cat(sprintf("复用置换缓存: %d 次\n", length(perm_auc)))
} else {
perm_auc <- replicate(20, {
  y2 <- sample(y)
  r <- roc(y2, loocv_enet(X, y2), quiet = TRUE)
  as.numeric(auc(r))
})
fwrite(data.table(perm = seq_along(perm_auc), auc = perm_auc),
       "results/GSE282512_classifier_perm_null.csv")
}
perm_p <- mean(c(perm_auc, Inf) >= perf_full$auc)  # 含观测自身, 下界 1/21
cat(sprintf("置换零分布 AUC: mean=%.3f max=%.3f; 置换 p=%.3f\n", mean(perm_auc), max(perm_auc), perm_p))

## ============ 5. top-k 精简 panel (fold 内 |t| 选择, 无泄漏) ============
topk_cv <- function(X, y, k, alpha = 0) {
  n <- ncol(X); pr <- numeric(n)
  for (i in 1:n) {
    tr <- setdiff(1:n, i)
    Xi <- X[, tr, drop = FALSE]; xi <- X[, i]
    med <- apply(Xi, 1, median, na.rm = TRUE); med[is.na(med)] <- 0.5
    Xi[is.na(Xi)] <- med[which(is.na(Xi), arr.ind = TRUE)[, "row"]]
    xi[is.na(xi)] <- med
    tt <- apply(Xi, 1, function(v) abs(suppressWarnings(t.test(v ~ factor(y[tr])))$statistic))
    tt[is.na(tt)] <- 0
    sel <- names(sort(tt, decreasing = TRUE))[seq_len(min(k, sum(tt > 0)))]
    if (length(sel) < 2) { pr[i] <- mean(y[tr]); next }
    fit <- suppressWarnings(glmnet(t(Xi[sel,,drop=FALSE]), factor(y[tr]), family = "binomial",
                                   alpha = alpha, lambda = 0.1))
    pr[i] <- as.numeric(predict(fit, t(xi[sel,drop=FALSE]), type = "response"))
  }
  pr
}
panel_rows <- list(); panel_sel_freq <- list(); roc_list <- list(roc_full)
for (k in c(10, 20)) {
  pr_k <- topk_cv(X, y, k)
  r <- roc(y, pr_k, quiet = TRUE); roc_list[[length(roc_list)+1]] <- r
  tk <- coords(r, "best", ret = "threshold", best.method = "youden")[[1]]
  pk <- as.integer(pr_k >= tk)
  panel_rows[[as.character(k)]] <- data.table(
    model = sprintf("top%d_foldselect", k), n = length(y), auc = as.numeric(auc(r)),
    auc_lo = as.numeric(ci.auc(r, method="delong")[1]), auc_hi = as.numeric(ci.auc(r, method="delong")[3]),
    sens = mean(pk[y==1]==1), spec = mean(pk[y==0]==0), acc = mean(pk==y))
  cat(sprintf("top-%d panel AUC = %.3f\n", k, as.numeric(auc(r))))
}

## ============ 6. EOPE / LOPE 分层 (子集内重新 LOOCV) ============
strat_rows <- list()
for (s in c("EOPE","LOPE")) {
  idx <- which(sub$group == "Control" | (sub$group == "PE" & sub$onset == s))
  Xs <- X[, idx, drop = FALSE]; ys <- y[idx]
  if (sum(ys) < 5 || sum(1-ys) < 5) next
  prs <- loocv_enet(Xs, ys)
  r <- roc(ys, prs, quiet = TRUE); roc_list[[length(roc_list)+1]] <- r
  strat_rows[[s]] <- data.table(
    model = sprintf("elastic_net_%s_vs_CT", s), n = length(idx),
    auc = as.numeric(auc(r)), auc_lo = as.numeric(ci.auc(r, method="delong")[1]),
    auc_hi = as.numeric(ci.auc(r, method="delong")[3]),
    sens = NA_real_, spec = NA_real_, acc = NA_real_)
  cat(sprintf("%s (n=%d, PE %d) AUC = %.3f\n", s, length(idx), sum(ys), as.numeric(auc(r))))
}

## ============ 7. 协变量 sanity: tube/batch 只有 ============
Xd <- t(model.matrix(~ tube + batch)[, -1, drop = FALSE])  # 特征 x 样本
storage.mode(Xd) <- "double"
pr_cov <- loocv_enet(Xd, y)
r_cov <- roc(y, pr_cov, quiet = TRUE)
cov_row <- data.table(model = "tube_batch_only", n = length(y), auc = as.numeric(auc(r_cov)),
                      auc_lo = NA_real_, auc_hi = NA_real_, sens = NA_real_, spec = NA_real_, acc = NA_real_)
cat(sprintf("tube/batch 协变量模型 AUC = %.3f (sanity, 应接近 0.5)\n", as.numeric(auc(r_cov))))

## ============ 8. 汇总输出 ============
perf <- rbindlist(c(list(perf_full), panel_rows, strat_rows, list(cov_row)), use.names = TRUE, fill = TRUE)
perf[, perm_p := NA_real_]; perf[model == "elastic_net_166DMR", perm_p := perm_p]
fwrite(perf, "results/GSE282512_classifier_performance.csv")
fwrite(data.table(perm = seq_along(perm_auc), auc = perm_auc),
       "results/GSE282512_classifier_perm_null.csv")

# 选择频率 (全模型 LOOCV 中非零系数区域)
sel_freq <- rbindlist(lapply(1:ncol(X), function(i) {
  tr <- setdiff(1:ncol(X), i)
  Xi <- X[, tr, drop = FALSE]
  med <- apply(Xi, 1, median, na.rm = TRUE); med[is.na(med)] <- 0.5
  Xi[is.na(Xi)] <- med[which(is.na(Xi), arr.ind = TRUE)[, "row"]]
  cvfit <- suppressWarnings(cv.glmnet(t(Xi), factor(y[tr]), family = "binomial",
                                      alpha = 0.5, nfolds = 5, type.measure = "auc"))
  cf <- as.matrix(coef(cvfit, s = "lambda.min"))[-1, 1]
  data.table(region = names(cf)[cf != 0], coef = cf[cf != 0])
}))
panel <- sel_freq[, .(sel_freq = .N, mean_coef = mean(coef), freq = .N / ncol(X)), by = region]
dir_map <- v2[!duplicated(region_id)][, .(region_id, direction)]
panel <- merge(panel, dir_map, by.x = "region", by.y = "region_id", all.x = TRUE)
panel <- panel[order(-sel_freq)]
fwrite(panel, "results/GSE282512_classifier_panel.csv")

sink("results/GSE282512_classifier_summary.txt")
cat("===== B: DMR 分类器 LOOCV 性能评估 =====\n\n")
cat(sprintf("特征: %d 候选 DMR 区域 | 样本: PE %d / Control %d (EOPE 21, LOPE 10, onset 缺失 1)\n\n",
            nrow(X), sum(y), sum(1-y)))
print(perf)
cat(sprintf("\n置换检验: 20 次标签置换零分布 AUC mean=%.3f, max=%.3f; 观测 AUC p=%.4f\n",
            mean(perm_auc), max(perm_auc), perm_p))
cat("\nLOOCV 选择频率 top-15 区域:\n"); print(head(panel, 15))
cat("\n注: top-k 与分层模型均在 fold 内完成特征选择/插补/lambda 调优, 无信息泄漏;\n")
cat("    sens/spec 阈值取 OOF 预测 Youden 点 (轻度乐观); 置换 p 含观测自身, 下界 1/21。\n")
sink()

## ============ 9. 图 ============
lab <- c("elastic_net_166DMR" = "全 166 DMR", "top10_foldselect" = "top-10 panel",
         "top20_foldselect" = "top-20 panel", "elastic_net_EOPE_vs_CT" = "EOPE vs CT",
         "elastic_net_LOPE_vs_CT" = "LOPE vs CT")
pd <- rbindlist(lapply(seq_along(roc_list), function(i) {
  r <- roc_list[[i]]; m <- perf$model[i]
  data.table(model = if (m %in% names(lab)) lab[[m]] else m,
             fpr = 1 - r$specificities, tpr = r$sensitivities)
}))
pd$auc_lab <- sprintf("%s (AUC=%.3f)", pd$model, perf$auc[match(pd$model, lab)])
lv <- sprintf("%s (AUC=%.3f)", lab, perf$auc[match(names(lab), perf$model)])
pd$auc_lab <- factor(pd$auc_lab, levels = lv[!is.na(perf$auc[match(names(lab), perf$model)])])
p1 <- ggplot(pd, aes(fpr, tpr, color = auc_lab)) + geom_line(linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60") +
  scale_color_brewer(palette = "Set1", name = NULL) +
  labs(x = "1 - 特异度", y = "灵敏度", title = "PE cfDNA 分类器 LOOCV ROC") +
  theme_classic(base_size = 11)
ggsave("figures/roc_loocv.png", p1, width = 6.5, height = 5, dpi = 300)
cat("完成: results/GSE282512_classifier_*.csv/txt + figures/roc_loocv.png\n")
