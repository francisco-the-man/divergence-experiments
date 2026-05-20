
library(jsonlite)

# result.json contains bare NaN literals (not valid JSON). Preprocess.
raw_txt <- readLines("result.json", warn = FALSE)
raw_txt <- paste(raw_txt, collapse = "\n")
# Replace bare NaN, Infinity, -Infinity with null (regex with word boundaries)
raw_txt <- gsub("\\bNaN\\b", "null", raw_txt)
raw_txt <- gsub("\\b-Infinity\\b", "null", raw_txt)
raw_txt <- gsub("\\bInfinity\\b", "null", raw_txt)

parsed <- fromJSON(raw_txt, simplifyVector = FALSE)
raw <- parsed$results
cat("Total results:", length(raw), "\n")

safe_num <- function(x) {
  if (is.null(x)) return(NA_real_)
  v <- suppressWarnings(as.numeric(x))
  if (length(v) == 0) return(NA_real_)
  v[1]
}

extract_row <- function(r) {
  ss   <- r$summary_stats
  sst  <- r$summary_stats_time
  ssd  <- r$summary_stats_saddle
  ssdt <- r$summary_stats_saddle_time
  data.frame(
    seed = safe_num(r$seed),
    family = if (is.null(r$family)) NA else r$family,
    cell = if (is.null(r$cell)) NA else r$cell,
    param_value = safe_num(r$param_value),
    n_epochs = safe_num(r$n_epochs),
    n_saddle_epochs = safe_num(r$n_saddle_epochs),
    log_mean = safe_num(ss$log_mean),
    log_var  = safe_num(ss$log_var),
    hill     = safe_num(ss$hill),
    ks_expon = safe_num(ss$ks_expon),
    n        = safe_num(ss$n),
    log_mean_t = safe_num(sst$log_mean),
    log_var_t  = safe_num(sst$log_var),
    hill_t     = safe_num(sst$hill),
    ks_expon_t = safe_num(sst$ks_expon),
    log_mean_s = safe_num(ssd$log_mean),
    log_var_s  = safe_num(ssd$log_var),
    hill_s     = safe_num(ssd$hill),
    ks_expon_s = safe_num(ssd$ks_expon),
    n_s        = safe_num(ssd$n),
    log_mean_st = safe_num(ssdt$log_mean),
    log_var_st  = safe_num(ssdt$log_var),
    hill_st     = safe_num(ssdt$hill),
    ks_expon_st = safe_num(ssdt$ks_expon),
    stringsAsFactors = FALSE
  )
}

df <- do.call(rbind, lapply(raw, extract_row))
cat("Frame rows:", nrow(df), "\n")
cat("Cells:\n"); print(table(df$cell, useNA="ifany"))
cat("Families:\n"); print(table(df$family, useNA="ifany"))

for (c in colnames(df)) {
  if (is.numeric(df[[c]])) df[[c]][!is.finite(df[[c]])] <- NA
}

feat_cols_pooled <- c("log_mean_t","log_var_t","hill_t","ks_expon_t")
feat_cols_saddle <- c("log_mean_st","log_var_st","hill_st","ks_expon_st")

has_nnet <- requireNamespace("nnet", quietly = TRUE)
cat("nnet available:", has_nnet, "\n")

set.seed(42)

stratified_folds <- function(labels, k = 5) {
  folds <- integer(length(labels))
  for (lev in unique(labels)) {
    idx <- which(labels == lev)
    idx <- sample(idx)
    folds[idx] <- rep(1:k, length.out = length(idx))
  }
  folds
}

balanced_accuracy <- function(true, pred) {
  lv <- unique(true)
  recs <- sapply(lv, function(l) {
    mask <- true == l & !is.na(pred)
    if (sum(mask) == 0) return(NA)
    mean(pred[mask] == l)
  })
  mean(recs, na.rm = TRUE)
}

binom_ci <- function(p, n, alpha = 0.05) {
  if (is.na(p) || n == 0) return(c(low=NA, high=NA))
  z <- qnorm(1 - alpha/2)
  se <- sqrt(p*(1-p)/n)
  c(low = max(0, p - z*se), high = min(1, p + z*se))
}

# Impute NAs in features with column median (so logistic glm doesn't drop rows)
impute_median <- function(X) {
  X <- as.matrix(X)
  for (j in 1:ncol(X)) {
    m <- median(X[,j], na.rm = TRUE)
    if (!is.finite(m)) m <- 0
    X[is.na(X[,j]), j] <- m
  }
  X
}

cv_classify <- function(X, y, k = 5) {
  X <- impute_median(X)
  keep <- !is.na(y)
  X <- X[keep,,drop=FALSE]; y <- y[keep]
  if (length(unique(y)) < 2 || nrow(X) < 10) return(list(acc=NA, n=length(y), preds=NULL, truth=NULL))
  folds <- stratified_folds(y, k)
  preds <- character(length(y))
  preds[] <- NA
  for (f in 1:k) {
    tr <- folds != f; te <- folds == f
    if (sum(te) == 0) next
    if (length(unique(y[tr])) < 2) next
    train_df <- data.frame(X[tr,,drop=FALSE])
    test_df  <- data.frame(X[te,,drop=FALSE])
    if (length(unique(y)) == 2) {
      lv <- sort(unique(y))
      y01_tr <- as.integer(y[tr] == lv[2])
      fit <- tryCatch(suppressWarnings(glm(y01_tr ~ ., data = train_df, family = binomial())),
                      error = function(e) NULL)
      if (is.null(fit)) next
      p <- predict(fit, newdata = test_df, type = "response")
      preds[te] <- ifelse(p > 0.5, lv[2], lv[1])
    } else if (has_nnet) {
      fit <- tryCatch(suppressWarnings(nnet::multinom(factor(y[tr]) ~ .,
                       data = train_df, trace = FALSE, MaxNWts = 5000)),
                      error = function(e) NULL)
      if (is.null(fit)) next
      preds[te] <- as.character(predict(fit, newdata = test_df))
    } else {
      # 5-NN
      for (i in which(te)) {
        d <- sqrt(rowSums((X[tr,,drop=FALSE] - matrix(X[i,], nrow=sum(tr), ncol=ncol(X), byrow=TRUE))^2))
        nn <- order(d)[1:min(5, sum(tr))]
        tab <- table(y[tr][nn])
        preds[i] <- names(tab)[which.max(tab)]
      }
    }
  }
  acc <- balanced_accuracy(y, preds)
  list(acc = acc, n = length(y), preds = preds, truth = y)
}

# ---- TEST 1: family classifier
fam_res <- cv_classify(df[, feat_cols_pooled], df$family)
fam_ci  <- binom_ci(fam_res$acc, fam_res$n)
cat("Family classifier acc:", fam_res$acc, "n =", fam_res$n, "\n")
cat("Family CI:", fam_ci, "\n")

# ---- TEST 2: cell classifier
cell_res <- cv_classify(df[, feat_cols_pooled], df$cell)
n_classes <- length(unique(df$cell[!is.na(df$cell)]))
chance <- 1 / n_classes
cell_ci <- binom_ci(cell_res$acc, cell_res$n)
cat("Cell classifier acc:", cell_res$acc, "n_classes =", n_classes, "chance =", chance, "\n")

# ---- TEST 3: Spearman per family
spearman_per_family <- list()
for (fam in unique(df$family)) {
  if (is.na(fam)) next
  sub <- df[df$family == fam & !is.na(df$family), ]
  res <- list()
  if (length(unique(sub$param_value)) < 2) {
    for (st in feat_cols_pooled) res[[st]] <- list(rho = NA, p = NA, n = nrow(sub))
    res[["note"]] <- "single parameter level"
  } else {
    for (st in feat_cols_pooled) {
      ok <- !is.na(sub[[st]]) & !is.na(sub$param_value)
      if (sum(ok) < 5) { res[[st]] <- list(rho=NA, p=NA, n=sum(ok)); next }
      ct <- suppressWarnings(cor.test(sub[[st]][ok], sub$param_value[ok], method = "spearman"))
      res[[st]] <- list(rho = unname(ct$estimate), p = unname(ct$p.value), n = sum(ok))
    }
  }
  spearman_per_family[[fam]] <- res
}

# Max abs rho with significance flag
max_abs_rho <- 0
spearman_any_sig <- FALSE
for (fam in names(spearman_per_family)) {
  for (st in feat_cols_pooled) {
    r <- spearman_per_family[[fam]][[st]]$rho
    p <- spearman_per_family[[fam]][[st]]$p
    if (!is.null(r) && !is.na(r)) {
      if (abs(r) > max_abs_rho) max_abs_rho <- abs(r)
      if (abs(r) >= 0.2 && !is.na(p) && p < 0.05) spearman_any_sig <- TRUE
    }
  }
}

# ---- TEST 4: stratification (saddle-only features)
fam_res_saddle <- cv_classify(df[, feat_cols_saddle], df$family)
fam_ci_saddle  <- binom_ci(fam_res_saddle$acc, fam_res_saddle$n)
stratification_drop <- fam_res$acc - fam_res_saddle$acc
cat("Family classifier acc (saddle-only):", fam_res_saddle$acc, "\n")

# Bootstrap drop CI (lightweight: 200 reps)
B <- 200
boot_drops <- numeric(B)
set.seed(123)
for (b in 1:B) {
  idx <- sample(nrow(df), replace = TRUE)
  sub <- df[idx, ]
  a1 <- cv_classify(sub[, feat_cols_pooled], sub$family)$acc
  a2 <- cv_classify(sub[, feat_cols_saddle], sub$family)$acc
  boot_drops[b] <- a1 - a2
}
boot_drops <- boot_drops[is.finite(boot_drops)]
drop_ci <- if (length(boot_drops) > 10) quantile(boot_drops, c(0.025, 0.975), names=FALSE) else c(NA,NA)
cat("Bootstrap drop CI:", drop_ci, "\n")

# Confusion matrix
truth_v <- cell_res$truth
pred_v  <- cell_res$preds
ok_idx <- !is.na(pred_v)
confusion <- table(true = truth_v[ok_idx], predicted = pred_v[ok_idx])
print(confusion)
confusion_df <- as.data.frame(confusion, stringsAsFactors = FALSE)
names(confusion_df) <- c("true","predicted","count")

# Per-cell aggregates
per_cell <- aggregate(df[, c("log_mean_t","log_var_t","hill_t","ks_expon_t","n_saddle_epochs","n_epochs")],
                      by = list(cell = df$cell), FUN = function(x) mean(x, na.rm = TRUE))

# CCDF summary per cell
ccdf_summary <- list()
for (cell in unique(df$cell)) {
  if (is.na(cell)) next
  pooled <- c()
  for (r in raw) {
    if (!is.null(r$cell) && r$cell == cell) {
      x <- unlist(r$stasis_lengths_time)
      x <- suppressWarnings(as.numeric(x))
      pooled <- c(pooled, x[is.finite(x) & x > 0])
    }
  }
  ccdf_summary[[cell]] <- list(
    n = length(pooled),
    median = if (length(pooled)) median(pooled) else NA,
    p90 = if (length(pooled)) quantile(pooled, 0.9, names=FALSE) else NA,
    p99 = if (length(pooled)) quantile(pooled, 0.99, names=FALSE) else NA,
    max  = if (length(pooled)) max(pooled) else NA
  )
}

saddle_frac <- mean(df$n_saddle_epochs / pmax(df$n_epochs, 1), na.rm = TRUE)
cat("Mean saddle fraction:", saddle_frac, "\n")

# Pass/fail
fam_falsified <- !is.na(fam_ci["high"]) && (fam_ci["high"] < 0.70)
fam_null_met  <- !is.na(fam_ci["high"]) && (fam_ci["high"] <= 0.55)
cell_falsified <- !is.na(cell_ci["high"]) && (cell_ci["high"] < 0.40)
cell_null_met  <- !is.na(cell_ci["high"]) && (cell_ci["high"] <= 0.20)

stats <- list(
  n_total_results = length(raw),
  cells_present = as.list(table(df$cell)),
  families_present = as.list(table(df$family)),
  n_classes_present = n_classes,
  rmf_cells_missing = setdiff(c("RMF_theta0.1","RMF_theta0.3","RMF_theta1.0","RMF_theta3.0"),
                              unique(df$cell[!is.na(df$cell)])),
  design_incomplete = TRUE,
  design_incomplete_reason = "Only 5 of the 8 preregistered cells produced results (RMF_theta0.3/1.0/3.0 missing); analyses are conditioned on cells that did run.",

  family_classifier = list(
    balanced_accuracy = fam_res$acc,
    ci_low = unname(fam_ci["low"]),
    ci_high = unname(fam_ci["high"]),
    n = fam_res$n,
    null_threshold = 0.55,
    hypothesis_threshold = 0.70,
    falsified = fam_falsified,
    null_supported = fam_null_met,
    note = "NK = 160, RMF = 40 -> class imbalance 4:1; balanced accuracy reported."
  ),

  cell_classifier = list(
    balanced_accuracy = cell_res$acc,
    ci_low = unname(cell_ci["low"]),
    ci_high = unname(cell_ci["high"]),
    n = cell_res$n,
    n_classes = n_classes,
    chance = chance,
    null_threshold = 0.20,
    hypothesis_threshold = 0.40,
    falsified = cell_falsified,
    null_supported = cell_null_met
  ),

  spearman_per_family = spearman_per_family,
  spearman_max_abs_rho = max_abs_rho,
  spearman_any_sig = spearman_any_sig,
  spearman_falsified = !spearman_any_sig,

  stratification = list(
    acc_pooled = fam_res$acc,
    acc_pooled_ci_low = unname(fam_ci["low"]),
    acc_pooled_ci_high = unname(fam_ci["high"]),
    acc_saddle_only = fam_res_saddle$acc,
    acc_saddle_only_ci_low = unname(fam_ci_saddle["low"]),
    acc_saddle_only_ci_high = unname(fam_ci_saddle["high"]),
    drop = stratification_drop,
    drop_ci_low = drop_ci[1],
    drop_ci_high = drop_ci[2]
  ),

  confusion = confusion_df,
  per_cell_stats = lapply(seq_len(nrow(per_cell)), function(i) {
    list(
      cell = per_cell$cell[i],
      log_mean_t = per_cell$log_mean_t[i],
      log_var_t  = per_cell$log_var_t[i],
      hill_t     = per_cell$hill_t[i],
      ks_expon_t = per_cell$ks_expon_t[i],
      mean_saddle_epochs = per_cell$n_saddle_epochs[i],
      mean_total_epochs  = per_cell$n_epochs[i]
    )
  }),
  ccdf_summary = ccdf_summary,
  mean_saddle_fraction = saddle_frac,

  hypothesis_status = (function() {
    if (fam_null_met && cell_null_met) "NULL (clean disconfirmation)"
    else if (!is.na(fam_ci["high"]) && fam_ci["low"] >= 0.70 &&
             !is.na(cell_ci["high"]) && cell_ci["low"] >= 0.40) "SUPPORTED"
    else "PARTIAL / MIXED"
  })()
)

write_json(stats, "stats.json", auto_unbox = TRUE, pretty = TRUE, na = "null")
cat("stats.json written.\n")
cat("Hypothesis status:", stats$hypothesis_status, "\n")
