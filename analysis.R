
suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(tidyr)
})

set.seed(42)

# Read raw file and sanitize JSON-illegal NaN / Infinity tokens
raw_txt <- paste(readLines("result.json", warn = FALSE), collapse = "\n")
raw_txt <- gsub("\\bNaN\\b", "null", raw_txt)
raw_txt <- gsub("\\b-Infinity\\b", "null", raw_txt)
raw_txt <- gsub("\\bInfinity\\b", "null", raw_txt)
raw <- fromJSON(raw_txt, simplifyVector = FALSE)$results
N <- length(raw)

feat_names <- c("log_mean","log_var","hill","ks_expon",
                paste0("ccdf_b", sprintf("%02d", 0:19)),
                "jump_log_mean","jump_log_var","stasis_jump_rankcorr","lag1_autocorr_logwt")

extract_feat <- function(r) {
  f <- r$features
  sapply(feat_names, function(nm) {
    v <- f[[nm]]
    if (is.null(v)) return(NA_real_)
    v <- suppressWarnings(as.numeric(v))
    if (length(v) == 0) return(NA_real_)
    if (is.nan(v) || is.infinite(v)) return(NA_real_)
    v
  })
}

X_full <- t(sapply(raw, extract_feat))
colnames(X_full) <- feat_names
cells     <- sapply(raw, function(r) r$cell)
families  <- sapply(raw, function(r) r$family)
param_val <- sapply(raw, function(r) as.numeric(r$param_value))
n_epochs  <- sapply(raw, function(r) as.integer(r$n_epochs))

cell_counts   <- as.list(table(cells))
family_counts <- as.list(table(families))
unique_cells <- sort(unique(cells))
K_classes <- length(unique_cells)
chance_acc <- 1 / K_classes

impute_median <- function(M) {
  for (j in seq_len(ncol(M))) {
    col <- M[, j]
    med <- median(col, na.rm = TRUE)
    if (is.na(med)) med <- 0
    col[is.na(col)] <- med
    M[, j] <- col
  }
  M
}
X_imp <- impute_median(X_full)

orig_cols  <- c("log_mean","log_var","hill","ks_expon")
ccdf_cols  <- paste0("ccdf_b", sprintf("%02d", 0:19))
jump_cols  <- c("jump_log_mean","jump_log_var","stasis_jump_rankcorr")
autoc_cols <- "lag1_autocorr_logwt"
set_4scalar <- orig_cols
set_pcc     <- c(orig_cols, ccdf_cols)
set_full    <- c(orig_cols, ccdf_cols, jump_cols, autoc_cols)

standardize_train_test <- function(Xtr, Xte) {
  mu <- colMeans(Xtr)
  sigma <- apply(Xtr, 2, sd)
  sigma[sigma < 1e-9] <- 1
  list(tr = sweep(sweep(Xtr, 2, mu, "-"), 2, sigma, "/"),
       te = sweep(sweep(Xte, 2, mu, "-"), 2, sigma, "/"))
}

stratified_folds <- function(y, k = 5, seed = 1) {
  set.seed(seed)
  folds <- integer(length(y))
  for (lev in unique(y)) {
    idx <- which(y == lev); idx <- sample(idx)
    folds[idx] <- (seq_along(idx) - 1) %% k + 1
  }
  folds
}

balanced_accuracy <- function(true, pred) {
  lv <- sort(unique(true))
  recalls <- sapply(lv, function(l) {
    m <- true == l
    if (sum(m) == 0) return(NA)
    mean(pred[m] == l)
  })
  mean(recalls, na.rm = TRUE)
}

knn_predict <- function(Xtr, ytr, Xte, k = 7) {
  yhat <- character(nrow(Xte))
  for (i in seq_len(nrow(Xte))) {
    d <- sqrt(rowSums((Xtr - matrix(Xte[i, ], nrow(Xtr), ncol(Xtr), byrow = TRUE))^2))
    ord <- order(d)
    nn  <- ytr[ord[1:k]]
    tab <- table(nn); mx <- max(tab)
    cand <- names(tab)[tab == mx]
    if (length(cand) == 1) { yhat[i] <- cand } else {
      for (j in ord) if (ytr[j] %in% cand) { yhat[i] <- ytr[j]; break }
    }
  }
  yhat
}

knn_cv_balacc <- function(X, y, k_nn = 7, k_folds = 5, seed = 1) {
  folds <- stratified_folds(y, k = k_folds, seed = seed)
  preds <- character(length(y))
  for (fold in seq_len(k_folds)) {
    tr <- folds != fold; te <- folds == fold
    sc <- standardize_train_test(X[tr, , drop = FALSE], X[te, , drop = FALSE])
    preds[te] <- knn_predict(sc$tr, y[tr], sc$te, k = k_nn)
  }
  list(bal_acc = balanced_accuracy(y, preds), preds = preds, truth = y)
}

boot_balacc_ci <- function(truth, preds, B = 800, seed = 7) {
  set.seed(seed); n <- length(truth)
  vals <- replicate(B, {
    ix <- sample.int(n, n, replace = TRUE)
    balanced_accuracy(truth[ix], preds[ix])
  })
  c(quantile(vals, 0.025, na.rm = TRUE), quantile(vals, 0.975, na.rm = TRUE))
}

K_NN <- 7

res_full <- knn_cv_balacc(X_imp[, set_full], cells, k_nn = K_NN, seed = 11)
ci_full  <- boot_balacc_ci(res_full$truth, res_full$preds, seed = 11)
res_4    <- knn_cv_balacc(X_imp[, set_4scalar], cells, k_nn = K_NN, seed = 12)
ci_4     <- boot_balacc_ci(res_4$truth, res_4$preds, seed = 12)
res_pcc  <- knn_cv_balacc(X_imp[, set_pcc], cells, k_nn = K_NN, seed = 13)
ci_pcc   <- boot_balacc_ci(res_pcc$truth, res_pcc$preds, seed = 13)

boot_delta_ci <- function(t1, p1, t2, p2, B = 800, seed = 21) {
  set.seed(seed); n <- length(t1)
  vals <- replicate(B, {
    ix <- sample.int(n, n, replace = TRUE)
    balanced_accuracy(t2[ix], p2[ix]) - balanced_accuracy(t1[ix], p1[ix])
  })
  c(mean(vals), quantile(vals, 0.025, na.rm = TRUE), quantile(vals, 0.975, na.rm = TRUE))
}

delta_pcc_4    <- boot_delta_ci(res_4$truth, res_4$preds, res_pcc$truth, res_pcc$preds, seed = 31)
delta_full_4   <- boot_delta_ci(res_4$truth, res_4$preds, res_full$truth, res_full$preds, seed = 32)
delta_full_pcc <- boot_delta_ci(res_pcc$truth, res_pcc$preds, res_full$truth, res_full$preds, seed = 33)

res_fam   <- knn_cv_balacc(X_imp[, set_full], families, k_nn = K_NN, seed = 41)
ci_fam    <- boot_balacc_ci(res_fam$truth, res_fam$preds, seed = 41)
res_fam_4 <- knn_cv_balacc(X_imp[, set_4scalar], families, k_nn = K_NN, seed = 42)
ci_fam_4  <- boot_balacc_ci(res_fam_4$truth, res_fam_4$preds, seed = 42)

# Permutation importance
set.seed(99); idx <- sample.int(N); te <- idx[1:floor(0.25*N)]; tr <- setdiff(seq_len(N), te)
sc <- standardize_train_test(X_imp[tr, set_full], X_imp[te, set_full])
base_pred <- knn_predict(sc$tr, cells[tr], sc$te, k = K_NN)
base_acc  <- balanced_accuracy(cells[te], base_pred)

perm_imp <- numeric(length(set_full)); names(perm_imp) <- set_full
n_repeats <- 15
set.seed(123)
Xtr_s <- sc$tr; Xte_s_base <- sc$te
for (j in seq_along(set_full)) {
  drops <- numeric(n_repeats)
  for (r in seq_len(n_repeats)) {
    Xte_s <- Xte_s_base
    Xte_s[, j] <- sample(Xte_s[, j])
    pp <- knn_predict(Xtr_s, cells[tr], Xte_s, k = K_NN)
    drops[r] <- base_acc - balanced_accuracy(cells[te], pp)
  }
  perm_imp[j] <- mean(drops)
}

feature_group_map <- c(
  setNames(rep("original", length(orig_cols)), orig_cols),
  setNames(rep("ccdf",     length(ccdf_cols)), ccdf_cols),
  setNames(rep("jump",     length(jump_cols)), jump_cols),
  setNames(rep("autocorr", length(autoc_cols)), autoc_cols)
)
perm_imp_df <- data.frame(
  feature = set_full,
  importance = unname(perm_imp),
  group = unname(feature_group_map[set_full]),
  stringsAsFactors = FALSE
)
perm_imp_df <- perm_imp_df[order(-perm_imp_df$importance), ]
rownames(perm_imp_df) <- NULL

# Spearman per (family, feature) -- only NK has param variation
spearman_results <- list()
for (fam in unique(families)) {
  rows <- which(families == fam)
  if (length(unique(param_val[rows])) < 2) next
  for (fn in set_full) {
    x <- X_full[rows, fn]; p <- param_val[rows]
    ok <- !is.na(x) & !is.na(p)
    if (sum(ok) < 10 || length(unique(x[ok])) < 2) {
      rho <- NA_real_; pv <- NA_real_
    } else {
      ct <- suppressWarnings(cor.test(x[ok], p[ok], method = "spearman"))
      rho <- as.numeric(ct$estimate); pv <- as.numeric(ct$p.value)
    }
    spearman_results[[length(spearman_results) + 1]] <- data.frame(
      family = fam, feature = fn, rho = rho, p = pv, n = sum(ok), stringsAsFactors = FALSE
    )
  }
}
spearman_no_variation_RMF <- !("RMF" %in% sapply(spearman_results, function(d) d$family[1]))
if (length(spearman_results) > 0) {
  sp_df <- do.call(rbind, spearman_results)
  sp_df$p_bh <- p.adjust(sp_df$p, method = "BH")
  sp_df$abs_rho <- abs(sp_df$rho)
  n_sig_strong <- sum(!is.na(sp_df$rho) & sp_df$abs_rho >= 0.2 & sp_df$p_bh < 0.05)
  max_abs_rho <- max(sp_df$abs_rho, na.rm = TRUE)
  top_feat_idx <- which.max(sp_df$abs_rho)
  top_feat <- sp_df[top_feat_idx, ]
  sp_top <- list(feature = top_feat$feature, family = top_feat$family,
                 rho = top_feat$rho, p_bh = top_feat$p_bh)
} else {
  sp_df <- data.frame(family=character(), feature=character(), rho=numeric(),
                      p=numeric(), n=integer(), p_bh=numeric(), abs_rho=numeric())
  n_sig_strong <- 0L; max_abs_rho <- NA_real_
  sp_top <- list(feature = NA, family = NA, rho = NA, p_bh = NA)
}

# Stratification drop
recompute_features_saddle <- function(r) {
  sl_strat <- r$stasis_lengths_stratified
  sl_saddle <- if (!is.null(sl_strat$saddle)) as.numeric(unlist(sl_strat$saddle)) else numeric(0)
  out <- list(log_mean = NA_real_, log_var = NA_real_, hill = NA_real_, ks_expon = NA_real_)
  if (length(sl_saddle) >= 4) {
    pos <- sl_saddle[sl_saddle > 0]
    if (length(pos) >= 4) {
      L <- log(pos)
      out$log_mean <- mean(L); out$log_var <- var(L)
    }
    if (length(pos) >= 8) {
      xs <- sort(pos); k <- max(3, floor(0.25 * length(xs)))
      if (k < length(xs)) {
        top <- xs[(length(xs)-k):length(xs)]; xk <- top[1]
        if (xk > 0) out$hill <- mean(log(top[-1] / xk))
      }
      rate <- 1 / mean(sl_saddle)
      if (rate > 0) {
        ks <- suppressWarnings(ks.test(sl_saddle, "pexp", rate))
        out$ks_expon <- as.numeric(ks$statistic)
      }
    }
  }
  out
}
X_saddle <- t(sapply(raw, function(r) {
  f <- recompute_features_saddle(r)
  sapply(orig_cols, function(nm) {
    v <- f[[nm]]
    if (is.null(v) || is.na(v) || is.nan(v) || is.infinite(v)) NA_real_ else v
  })
}))
colnames(X_saddle) <- orig_cols
X_saddle_imp <- impute_median(X_saddle)

res_fam_pool <- knn_cv_balacc(X_imp[, set_4scalar], families, k_nn = K_NN, seed = 51)
res_fam_sad  <- knn_cv_balacc(X_saddle_imp, families, k_nn = K_NN, seed = 52)
strat_drop <- res_fam_pool$bal_acc - res_fam_sad$bal_acc
set.seed(61); B <- 800
drops <- replicate(B, {
  ix <- sample.int(length(res_fam_pool$truth), length(res_fam_pool$truth), replace = TRUE)
  balanced_accuracy(res_fam_pool$truth[ix], res_fam_pool$preds[ix]) -
    balanced_accuracy(res_fam_sad$truth[ix], res_fam_sad$preds[ix])
})
strat_ci <- as.numeric(quantile(drops, c(0.025, 0.975), na.rm = TRUE))

conf_levels <- sort(unique(c(res_full$truth, res_full$preds)))
conf_mat <- table(true = factor(res_full$truth, levels = conf_levels),
                  pred = factor(res_full$preds, levels = conf_levels))
conf_df <- as.data.frame(conf_mat, stringsAsFactors = FALSE)
names(conf_df) <- c("true", "pred", "count")
class_recall <- sapply(conf_levels, function(l) {
  m <- res_full$truth == l
  if (sum(m) == 0) NA else mean(res_full$preds[m] == l)
})
names(class_recall) <- conf_levels

epoch_summary_by_cell <- lapply(split(n_epochs, cells), function(v)
  list(mean = mean(v), median = median(v), min = min(v), max = max(v)))

stats <- list(
  n_tasks = N,
  n_classes = K_classes,
  unique_cells = unique_cells,
  cell_counts = cell_counts,
  family_counts = family_counts,
  chance_accuracy = chance_acc,
  classifier_type = "k-NN (k=7) on standardized features, stratified 5-fold CV (randomForest unavailable in env)",
  hypothesis_threshold_bin = 0.40,
  hypothesis_threshold_fam = 0.70,
  parent_bin_acc = 0.2844,
  parent_fam_acc = 0.7062,

  bin_classifier_full_acc     = unname(res_full$bal_acc),
  bin_classifier_full_ci_low  = unname(ci_full[1]),
  bin_classifier_full_ci_high = unname(ci_full[2]),
  bin_classifier_full_pass    = unname(ci_full[2] >= 0.40),

  bin_classifier_4_acc     = unname(res_4$bal_acc),
  bin_classifier_4_ci_low  = unname(ci_4[1]),
  bin_classifier_4_ci_high = unname(ci_4[2]),

  bin_classifier_pcc_acc     = unname(res_pcc$bal_acc),
  bin_classifier_pcc_ci_low  = unname(ci_pcc[1]),
  bin_classifier_pcc_ci_high = unname(ci_pcc[2]),

  delta_pcc_minus_4_mean    = unname(delta_pcc_4[1]),
  delta_pcc_minus_4_ci_low  = unname(delta_pcc_4[2]),
  delta_pcc_minus_4_ci_high = unname(delta_pcc_4[3]),

  delta_full_minus_4_mean    = unname(delta_full_4[1]),
  delta_full_minus_4_ci_low  = unname(delta_full_4[2]),
  delta_full_minus_4_ci_high = unname(delta_full_4[3]),

  delta_full_minus_pcc_mean    = unname(delta_full_pcc[1]),
  delta_full_minus_pcc_ci_low  = unname(delta_full_pcc[2]),
  delta_full_minus_pcc_ci_high = unname(delta_full_pcc[3]),

  family_classifier_full_acc     = unname(res_fam$bal_acc),
  family_classifier_full_ci_low  = unname(ci_fam[1]),
  family_classifier_full_ci_high = unname(ci_fam[2]),
  family_classifier_full_pass    = unname(ci_fam[2] >= 0.70),

  family_classifier_4_acc     = unname(res_fam_4$bal_acc),
  family_classifier_4_ci_low  = unname(ci_fam_4[1]),
  family_classifier_4_ci_high = unname(ci_fam_4[2]),

  permutation_importance = perm_imp_df,
  top5_features          = perm_imp_df[1:min(5, nrow(perm_imp_df)), ],

  spearman_no_variation_RMF = spearman_no_variation_RMF,
  spearman_max_abs_rho      = max_abs_rho,
  spearman_top_feature      = sp_top$feature,
  spearman_top_family       = sp_top$family,
  spearman_top_rho          = sp_top$rho,
  spearman_top_p_bh         = sp_top$p_bh,
  spearman_n_sig_strong     = n_sig_strong,
  spearman_table            = sp_df,

  strat_acc_pooled  = unname(res_fam_pool$bal_acc),
  strat_acc_saddle  = unname(res_fam_sad$bal_acc),
  strat_drop        = unname(strat_drop),
  strat_drop_ci_low = unname(strat_ci[1]),
  strat_drop_ci_high= unname(strat_ci[2]),

  confusion        = conf_df,
  per_class_recall = as.list(class_recall),

  n_epochs_mean     = mean(n_epochs),
  n_epochs_median   = median(n_epochs),
  n_epochs_min      = min(n_epochs),
  n_epochs_max      = max(n_epochs),
  n_epochs_pct_lt_4 = mean(n_epochs < 4),
  n_epochs_by_cell  = epoch_summary_by_cell,
  nan_rate_original_features = mean(is.na(X_full[, "log_mean"])),

  hypothesis_status   = ifelse(ci_full[2] >= 0.40, "SUSTAINED", "FALSIFIED"),
  notes_missing_cells = paste(
    "Only 5 of the planned 8 cells appear in result.json (NK_K2/K4/K8/K16, RMF_theta0.1).",
    "RMF_theta0.3/1.0/3.0 are absent. The bin classifier is therefore effectively 5-way (chance = 0.20), not 8-way (chance = 0.125).",
    "Only one RMF parameter value is present, so the per-family Spearman-vs-theta test cannot be estimated on the RMF side and is reported NK-only."
  )
)

write_json(stats, "stats.json", auto_unbox = TRUE, pretty = TRUE, na = "null", digits = 6)
cat("done\n")
