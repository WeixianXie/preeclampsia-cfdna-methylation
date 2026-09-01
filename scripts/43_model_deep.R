## ============================================================
## 43. Phase 8 模型深度报告（任务 #41）
## 在 40 号脚本全模型 LOOCV 基础上补充:
##   1) bootstrap 校准曲线 + 校准截距/斜率 (rms::val.prob)
##   2) Hosmer-Lemeshow 检验 (ResourceSelection)
##   3) DCA 决策曲线 (dcurves)
##   4) 对照模型: XGBoost / SVM / 单变量评分 LR (LOOCV)
##   5) SHAP 特征重要性 (xgboost predcontrib, 全数据模型)
##   6) 6 稳定区域列线图 (rms::lrm + nomogram)
## 输出: results/GSE282512_classifier_depth_*.csv/txt
##       figures/calibration.png / dca.png / shap_importance.png / nomogram.png
## ============================================================
suppressPackageStartupMessages({
  library(data.table); library(glmnet); library(pROC); library(ggplot2)
  library(xgboost); library(e1071); library(ResourceSelection)
})
set.seed(42)

## ---------- 数据（与 40 号完全一致） ----------
sub <- fread("results/GSE282512_subcohort.csv")
breg <- fread("results/GSE282512_region_beta.csv.gz")
breg[, region_id := as.character(region_id)]
bmat <- as.matrix(breg[, -1]); rownames(bmat) <- breg$region_id
v2 <- fread("results/GSE282512_dmr_final_v2.csv", select = c("region_id","direction","candidate"))
v2[, region_id := as.character(region_id)]
cand <- v2[isTRUE(candidate) | candidate == "TRUE"]$region_id
cand <- cand[cand %in% rownames(bmat)]
gs <- sub$gsm
X <- bmat[cand, gs, drop = FALSE]
y <- as.integer(sub$group == "PE")
n <- ncol(X)

## LOOCV 弹性网（同 40 号函数, 确定性复现 OOF 预测）
loocv_enet <- function(X, y, alpha = 0.5) {
  p <- nrow(X); n <- ncol(X); pr <- numeric(n)
  for (i in 1:n) {
    tr <- setdiff(1:n, i)
    Xi <- X[, tr, drop = FALSE]; xi <- X[, i]
    med <- apply(Xi, 1, median, na.rm = TRUE); med[is.na(med)] <- 0.5
    Xi[is.na(Xi)] <- med[which(is.na(Xi), arr.ind = TRUE)[, "row"]]
    xi[is.na(xi)] <- med
    cvfit <- suppressWarnings(cv.glmnet(t(Xi), factor(y[tr]), family = "binomial",
                                        alpha = alpha, nfolds = 5, type.measure = "auc"))
    pr[i] <- as.numeric(predict(cvfit, t(xi), s = "lambda.min", type = "response"))
  }
  pr
}
pr_enet <- loocv_enet(X, y)
fwrite(data.table(gsm = gs, group = sub$group, y, oof_pred = pr_enet),
       "results/GSE282512_classifier_oof_predictions.csv")

## ---------- 1. 校准 (bootstrap, rms::val.prob) ----------
suppressPackageStartupMessages(library(rms))
opr <- options(digits = 4)
cal <- val.prob(pr_enet, y, m = 16, smooth = FALSE, pl = FALSE)
sink("results/GSE282512_classifier_depth_summary.txt")
cat("===== Phase 8 模型深度报告（弹性网全模型 LOOCV, AUC=0.932）=====\n\n")
cat("1) 校准（rms::val.prob, 分组 m=16）:\n"); print(cal)
cat("\n  Slope = 校准斜率(理想=1), Intercept = 校准截距(理想=0),\n")
cat("  P-val 为校准直线/二次曲线不可靠性检验(H0: 校准良好)。\n")
sink()
png("figures/calibration.png", 2000, 1600, res = 240)
val.prob(pr_enet, y, m = 16, smooth = TRUE, xlab = "预测概率", ylab = "观测 PE 频率")
title("校准曲线 (LOOCV 预测)")
dev.off()

## ---------- 2. Hosmer-Lemeshow ----------
hl <- hoslem.test(y, pr_enet, g = 10)
cat(sprintf("\n2) Hosmer-Lemeshow (10 分组): chi2=%.2f, df=%d, p=%.4f\n",
            hl$statistic, hl$parameter, hl$p.value),
    file = "results/GSE282512_classifier_depth_summary.txt", append = TRUE)

## ---------- 3. DCA (dcurves) ----------
suppressPackageStartupMessages(library(dcurves))
dca_df <- data.frame(pe = factor(ifelse(y == 1, "yes", "no"), levels = c("no", "yes")),
                     model = pr_enet)
dca_res <- tryCatch(dca(pe ~ model, data = dca_df,
                        thresholds = seq(0.05, 0.95, by = 0.05),
                        label = list("DMR classifier" = "model")),
                    error = function(e) NULL)
if (!is.null(dca_res)) {
  nb <- as.data.table(dca_res$dca)
  fwrite(nb, "results/GSE282512_classifier_dca.csv")
  png("figures/dca.png", 2200, 1700, res = 240)
  print(plot(dca_res))
  dev.off()
  nbg <- nb[variable == "model" & threshold %in% c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)]
  cat("\n3) DCA: 分类器净获益（阈值 0.1-0.6）:\n",
      file = "results/GSE282512_classifier_depth_summary.txt", append = TRUE)
  cat(capture.output(print(nbg[, .(threshold, net_benefit)])),
      file = "results/GSE282512_classifier_depth_summary.txt", append = TRUE, sep = "\n")
}

## ---------- 4. 对照模型 LOOCV ----------
xgb_params <- list(objective = "binary:logistic", eval_metric = "auc",
                   eta = 0.05, max_depth = 2, subsample = 0.8, colsample_bytree = 0.6)
best_nr <- function(cv) {
  if (!is.null(cv$best_iteration)) return(as.integer(cv$best_iteration))
  el <- cv$evaluation_log
  nm <- grep("test_.*_mean", names(el), value = TRUE)[1]
  as.integer(which.max(el[[nm]]))
}
## 4a. XGBoost（fold 内 xgb.cv 定 nrounds, 其余固定）
loocv_xgb <- function(X, y) {
  n <- ncol(X); pr <- numeric(n)
  for (i in 1:n) {
    tr <- setdiff(1:n, i)
    Xi <- X[, tr, drop = FALSE]; xi <- X[, i]
    med <- apply(Xi, 1, median, na.rm = TRUE); med[is.na(med)] <- 0.5
    Xi[is.na(Xi)] <- med[which(is.na(Xi), arr.ind = TRUE)[, "row"]]
    xi[is.na(xi)] <- med
    dtr <- xgb.DMatrix(t(Xi), label = y[tr])
    cvf <- xgb.cv(params = xgb_params, data = dtr, nrounds = 200, nfold = 5,
                  early_stopping_rounds = 20, verbose = 0)
    bst <- xgb.train(params = xgb_params, data = dtr, nrounds = best_nr(cvf))
    pr[i] <- predict(bst, xgb.DMatrix(t(xi)))
  }
  pr
}
## 4b. SVM (radial, 概率输出)
loocv_svm <- function(X, y) {
  n <- ncol(X); pr <- numeric(n)
  for (i in 1:n) {
    tr <- setdiff(1:n, i)
    Xi <- X[, tr, drop = FALSE]; xi <- X[, i]
    med <- apply(Xi, 1, median, na.rm = TRUE); med[is.na(med)] <- 0.5
    Xi[is.na(Xi)] <- med[which(is.na(Xi), arr.ind = TRUE)[, "row"]]
    xi[is.na(xi)] <- med
    fit <- svm(x = data.frame(t(Xi)), y = factor(y[tr]), probability = TRUE)
    pp <- attr(predict(fit, newdata = data.frame(t(xi)), probability = TRUE), "probabilities")
    ci <- which(colnames(pp) == "1")
    if (length(ci) == 1) pr[i] <- pp[1, ci] else pr[i] <- 1 - pp[1, 1]
  }
  pr
}
## 4c. 单变量方向加权评分 + 单变量 LR（最简基线）
dir_w <- ifelse(v2[match(cand, region_id)]$direction == "hyper", 1, -1)
Xs <- sweep(X, 1, dir_w, `*`)                     # 方向统一后求均值
score <- colMeans(Xs, na.rm = TRUE)
pr_lr <- numeric(n)
for (i in 1:n) {
  tr <- setdiff(1:n, i)
  fit <- glm(y[tr] ~ score[tr], family = binomial)
  pr_lr[i] <- predict(fit, data.frame(score = score[i]), type = "response")
}

cmp <- rbindlist(lapply(list(
  list("elastic_net_166DMR" = pr_enet),
  list("xgboost_LOOCV" = NULL),
  list("svm_rbf_LOOCV" = NULL),
  list("univariate_score_LR" = pr_lr)), function(l) NULL))
pr_xgb <- tryCatch(loocv_xgb(X, y), error = function(e) {cat("xgb err:", conditionMessage(e), "\n"); NULL})
pr_svm <- tryCatch(loocv_svm(X, y), error = function(e) {cat("svm err:", conditionMessage(e), "\n"); NULL})
models <- list(`elastic_net_166DMR` = pr_enet, `univariate_score_LR` = pr_lr)
if (!is.null(pr_xgb)) models$`xgboost_LOOCV` <- pr_xgb
if (!is.null(pr_svm)) models$`svm_rbf_LOOCV` <- pr_svm
cmp <- rbindlist(lapply(names(models), function(m) {
  r <- roc(y, models[[m]], quiet = TRUE)
  data.table(model = m, auc = as.numeric(auc(r)),
             auc_lo = as.numeric(ci.auc(r, method = "delong")[1]),
             auc_hi = as.numeric(ci.auc(r, method = "delong")[3]))
}))
fwrite(cmp, "results/GSE282512_classifier_model_comparison.csv")
cat("\n4) 对照模型 LOOCV AUC:\n", file = "results/GSE282512_classifier_depth_summary.txt", append = TRUE)
cat(capture.output(print(cmp)),
    file = "results/GSE282512_classifier_depth_summary.txt", append = TRUE, sep = "\n")

## ---------- 5. SHAP（全数据 XGBoost, predcontrib 原生 SHAP） ----------
med <- apply(X, 1, median, na.rm = TRUE); med[is.na(med)] <- 0.5
Xf <- X; Xf[is.na(Xf)] <- med[which(is.na(Xf), arr.ind = TRUE)[, "row"]]
dall <- xgb.DMatrix(t(Xf), label = y)
cvf <- xgb.cv(params = xgb_params, data = dall, nrounds = 200, nfold = 5,
              early_stopping_rounds = 20, verbose = 0)
bst <- xgb.train(params = xgb_params, data = dall, nrounds = best_nr(cvf))
sh <- predict(bst, dall, predcontrib = TRUE)          # n x (p+1)
sh_mat <- sh[, 1:nrow(X)]
rownames(sh_mat) <- colnames(sh) <- NULL
sh_imp <- rowMeans(abs(sh_mat))
sh_dt <- data.table(region = cand, direction = v2[match(cand, region_id)]$direction,
                    mean_abs_shap = sh_imp)[order(-mean_abs_shap)]
fwrite(sh_dt, "results/GSE282512_classifier_shap_importance.csv")
png("figures/shap_importance.png", 2000, 2000, res = 240)
top20 <- sh_dt[1:min(20, nrow(sh_dt))]
top20$region <- factor(top20$region, levels = rev(top20$region))
ggplot(top20, aes(mean_abs_shap, region, fill = direction)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(hyper = "#c0392b", hypo = "#2471a3"), name = "DMR 方向") +
  labs(x = "平均 |SHAP|（对数几率贡献）", y = "DMR 区域",
       title = "XGBoost SHAP 特征重要性（top-20）") +
  theme_classic(base_size = 11)
dev.off()
cat(sprintf("\n5) SHAP: 全数据 XGBoost(nrounds=%d); top-5 区域 %s\n",
            best_nr(cvf), paste(sh_dt$region[1:5], collapse = ", ")),
    file = "results/GSE282512_classifier_depth_summary.txt", append = TRUE)

## ---------- 6. 列线图（6 稳定区域, rms lrm + nomogram） ----------
panel <- fread("results/GSE282512_classifier_panel.csv")
stable6 <- as.character(panel[sel_freq == n]$region)
if (length(stable6) >= 3) {
  dd <- data.frame(y = y, t(Xf[stable6, , drop = FALSE]))
  names(dd) <- c("y", paste0("R", stable6))
  ddist <- datadist(dd); options(datadist = "ddist")
  fit6 <- tryCatch(lrm(y ~ ., data = dd, x = TRUE, y = TRUE), error = function(e) NULL)
  if (!is.null(fit6)) {
    nom <- tryCatch(nomogram(fit6, fun = function(x) 1 / (1 + exp(-x)),
                             funlabel = "PE 概率", fun.at = c(0.1, 0.3, 0.5, 0.7, 0.9)),
                    error = function(e) NULL)
    if (!is.null(nom)) {
      png("figures/nomogram.png", 2400, 1500, res = 240)
      plot(nom, cex.axis = 0.65, cex.var = 0.7, main = NULL)
      dev.off()
      cat(sprintf("\n6) 列线图: %d 个 100%% 选择频率区域, lrm C-index=%.3f\n",
                  length(stable6), fit6$stats["C"]),
          file = "results/GSE282512_classifier_depth_summary.txt", append = TRUE)
    }
  }
}
cat("\n完成\n", file = "results/GSE282512_classifier_depth_summary.txt", append = TRUE)
cat("完成: results/GSE282512_classifier_depth_*.csv/txt + figures/{calibration,dca,shap_importance,nomogram}.png\n")
