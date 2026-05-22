
library(jsonlite)
library(rpart)

# ---- Robust JSON load -------------------------------------------------------
raw <- readLines("result.json", warn = FALSE)
raw <- paste(raw, collapse = "\n")
raw <- gsub("\\bNaN\\b", "null", raw)
raw <- gsub("\\bInfinity\\b", "null", raw)
raw <- gsub("\\b-Infinity\\b", "null", raw)
dat <- fromJSON(raw, simplifyVector = FALSE)$results
rm(raw); gc()

n_total <- length(dat)
cat("n_total =", n_total, "\n")

safe_num <- function(x) {
  if (is.null(x)) return(NA_real_)
  v <- suppressWarnings(as.numeric(x))
  if (length(v) == 0) NA_real_ else v[1]
}

param_bin <- vapply(dat, function(r) r$param_bin_label, character(1))
family_label <- vapply(dat, function(r) r$family_label, character(1))
param_value <- vapply(dat, function(r) safe_num(r$param_value), numeric(1))
n_epochs <- vapply(dat, function(r) safe_num(r$n_epochs_primary), numeric(1))

fnames <- names(dat[[1]]$features_time_primary)
fg <- dat[[1]]$feature_groups
feat_group <- vapply(fnames, function(fn) {
  v <- fg[[fn]]; if (is.null(v)) "unknown" else as.character(v)
}, character(1))

build_mat <- function(feat_key) {
  M <- matrix(NA_real_, nrow = n_total, ncol = length(fnames),
              dimnames = list(NULL, fnames))
  for (i in seq_len(n_total)) {
    fv <- dat[[i]][[feat_key]]
    if (is.null(fv)) next
    for (fn in fnames) {
      v <- fv[[fn]]
      if (is.null(v)) next
      vv <- suppressWarnings(as.numeric(v))
      if (length(vv)) M[i, fn] <- vv[1]
    }
  }
  M
}

X_time <- build_mat("features_time_primary")
X_subs <- build_mat("features_subs_primary")

eps_keys <- c("0.001", "0.005", "0.02")
X_eps <- list()
for (ek in eps_keys) {
  M <- matrix(NA_real_, nrow = n_total, ncol = length(fnames),
              dimnames = list(NULL, fnames))
  for (i in seq_len(n_total)) {
    es <- dat[[i]]$eps_sensitivity
    if (is.null(es)) next
    blk <- es[[ek]]
    if (is.null(blk)) next
    fv <- blk$features
    if (is.null(fv)) next
    for (fn in fnames) {
      v <- fv[[fn]]
      if (is.null(v)) next
      vv <- suppressWarnings(as.numeric(v))
      if (length(vv)) M[i, fn] <- vv[1]
    }
  }
  X_eps[[ek]] <- M
}

scalar_feats <- fnames[feat_group == "scalar"]
ccdf_stasis_feats <- fnames[feat_group == "ccdf_stasis"]
ccdf_jump_feats <- fnames[feat_group == "ccdf_jump"]

impute_med <- function(M) {
  for (j in seq_len(ncol(M))) {
    col <- M[, j]
    bad <- is.na(col) | !is.finite(col)
    if (any(bad)) {
      v <- col[!bad]; m <- if (length(v)) median(v) else 0
      if (!is.finite(m)) m <- 0
      M[bad, j] <- m
    }
  }
  M
}

# ---- Bagged trees (rpart) as RandomForest surrogate -------------------------
# Bootstrap n_trees rpart fits on random feature subsets (mtry = sqrt(p)).
fit_bagged <- function(Xtr, ytr, n_trees = 100, seed = 1) {
  set.seed(seed)
  ytr_f <- factor(ytr)
  classes <- levels(ytr_f)
  n <- nrow(Xtr); p <- ncol(Xtr)
  mtry <- max(2, floor(sqrt(p)))
  trees <- vector("list", n_trees)
  feat_cols_list <- vector("list", n_trees)
  for (t in seq_len(n_trees)) {
    boot_idx <- sample.int(n, n, replace = TRUE)
    feat_cols <- sample.int(p, mtry)
    feat_cols_list[[t]] <- feat_cols
    df <- as.data.frame(Xtr[boot_idx, feat_cols, drop = FALSE])
    df$.y <- ytr_f[boot_idx]
    trees[[t]] <- rpart(.y ~ ., data = df,
                        control = rpart.control(minsplit = 5, cp = 0.001,
                                                maxdepth = 12, xval = 0))
  }
  list(trees = trees, feat_cols_list = feat_cols_list, classes = classes,
       feat_names = colnames(Xtr))
}

predict_bagged <- function(fit, Xte) {
  n <- nrow(Xte)
  votes <- matrix(0, nrow = n, ncol = length(fit$classes),
                  dimnames = list(NULL, fit$classes))
  for (t in seq_along(fit$trees)) {
    fc <- fit$feat_cols_list[[t]]
    df <- as.data.frame(Xte[, fc, drop = FALSE])
    pr <- predict(fit$trees[[t]], df, type = "class")
    for (i in seq_len(n)) {
      cl <- as.character(pr[i])
      votes[i, cl] <- votes[i, cl] + 1
    }
  }
  fit$classes[apply(votes, 1, which.max)]
}

make_folds <- function(y, k = 5, seed = 7) {
  set.seed(seed)
  folds <- integer(length(y))
  for (lev in unique(y)) {
    idx <- which(y == lev)
    idx <- sample(idx)
    folds[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  folds
}

bal_acc <- function(true, pred) {
  lv <- unique(true)
  rec <- vapply(lv, function(l) {
    sel <- true == l
    if (sum(sel) == 0) NA_real_ else mean(pred[sel] == l)
  }, numeric(1))
  mean(rec, na.rm = TRUE)
}

cv_bal_acc <- function(X, y, k = 5, seed = 7, n_trees = 100) {
  folds <- make_folds(y, k, seed)
  per_fold <- numeric(k)
  preds_all <- character(length(y))
  for (f in seq_len(k)) {
    te <- folds == f; tr <- !te
    Xt <- impute_med(X[tr, , drop = FALSE])
    meds <- apply(Xt, 2, function(c) {
      v <- c[is.finite(c)]; if (!length(v)) 0 else median(v)
    })
    Xte <- X[te, , drop = FALSE]
    for (j in seq_len(ncol(Xte))) {
      bad <- !is.finite(Xte[, j])
      if (any(bad)) Xte[bad, j] <- meds[j]
    }
    fit <- fit_bagged(Xt, y[tr], n_trees = n_trees, seed = seed + f)
    p <- predict_bagged(fit, Xte)
    preds_all[te] <- p
    per_fold[f] <- bal_acc(y[te], p)
  }
  list(per_fold = per_fold, mean = mean(per_fold), preds = preds_all, folds = folds)
}

binom_ci <- function(acc, n, level = 0.95) {
  if (n <= 0 || is.na(acc)) return(c(NA, NA))
  z <- qnorm(1 - (1 - level) / 2)
  se <- sqrt(acc * (1 - acc) / n)
  c(max(0, acc - z * se), min(1, acc + z * se))
}

y_bin <- param_bin
n_classes <- length(unique(y_bin))
cat("n bin classes =", n_classes, "\n")

# Use moderate n_trees to keep runtime sane
NT <- 80

set.seed(42)
cat("Running full time CV...\n")
res_full_time <- cv_bal_acc(X_time, y_bin, 5, 7, n_trees = NT)
acc_full_time <- res_full_time$mean
ci_full_time <- binom_ci(acc_full_time, n_total)
cat("acc_full_time =", acc_full_time, " CI =", ci_full_time, "\n")

cat("Running 4-scalar time CV...\n")
res_4_time <- cv_bal_acc(X_time[, scalar_feats, drop=FALSE], y_bin, 5, 7, n_trees = NT)
acc_4_time <- res_4_time$mean

ccdf_set <- c(scalar_feats, ccdf_stasis_feats)
cat("Running +CCDF time CV...\n")
res_ccdf_time <- cv_bal_acc(X_time[, ccdf_set, drop=FALSE], y_bin, 5, 7, n_trees = NT)
acc_ccdf_time <- res_ccdf_time$mean

cat("Running full subs CV...\n")
res_full_subs <- cv_bal_acc(X_subs, y_bin, 5, 7, n_trees = NT)
acc_full_subs <- res_full_subs$mean
ci_full_subs <- binom_ci(acc_full_subs, n_total)

cat("Running 4-scalar subs CV...\n")
res_4_subs <- cv_bal_acc(X_subs[, scalar_feats, drop=FALSE], y_bin, 5, 7, n_trees = NT)
acc_4_subs <- res_4_subs$mean

cat("Running +CCDF subs CV...\n")
res_ccdf_subs <- cv_bal_acc(X_subs[, ccdf_set, drop=FALSE], y_bin, 5, 7, n_trees = NT)
acc_ccdf_subs <- res_ccdf_subs$mean

# Bootstraps
pf_time <- res_full_time$per_fold
pf_subs <- res_full_subs$per_fold
delta_per_fold <- pf_time - pf_subs
set.seed(99)
boot_delta <- replicate(2000, mean(sample(delta_per_fold, replace = TRUE)))
ci_delta_time_subs <- quantile(boot_delta, c(0.025, 0.975))
mean_delta_time_subs <- mean(delta_per_fold)

boot_pair <- function(a, b, B = 2000) {
  d <- a - b
  bb <- replicate(B, mean(sample(d, replace = TRUE)))
  list(mean = mean(d), ci = quantile(bb, c(0.025, 0.975)))
}
abl_full_minus_4 <- boot_pair(res_full_time$per_fold, res_4_time$per_fold)
abl_ccdf_minus_4 <- boot_pair(res_ccdf_time$per_fold, res_4_time$per_fold)
abl_full_minus_ccdf <- boot_pair(res_full_time$per_fold, res_ccdf_time$per_fold)

# Permutation importance
cat("Permutation importance...\n")
set.seed(13)
trn_idx <- unlist(lapply(unique(y_bin), function(l) {
  ii <- which(y_bin == l); sample(ii, size = floor(0.8 * length(ii)))
}))
tst_idx <- setdiff(seq_len(n_total), trn_idx)
Xtr_imp <- impute_med(X_time[trn_idx, , drop = FALSE])
meds_imp <- apply(Xtr_imp, 2, function(c) {
  v <- c[is.finite(c)]; if (!length(v)) 0 else median(v)
})
Xte_imp <- X_time[tst_idx, , drop = FALSE]
for (j in seq_len(ncol(Xte_imp))) {
  bad <- !is.finite(Xte_imp[, j])
  if (any(bad)) Xte_imp[bad, j] <- meds_imp[j]
}
ytr_imp <- y_bin[trn_idx]; yte_imp <- y_bin[tst_idx]

fit_imp <- fit_bagged(Xtr_imp, ytr_imp, n_trees = 150, seed = 2025)
base_pred <- predict_bagged(fit_imp, Xte_imp)
base_acc <- bal_acc(yte_imp, base_pred)
cat("perm base acc =", base_acc, "\n")

n_rep <- 15  # reduce for runtime, still gives stable ranking
perm_importance <- numeric(length(fnames))
names(perm_importance) <- fnames
set.seed(2030)
for (jf in seq_along(fnames)) {
  drops <- numeric(n_rep)
  for (r in seq_len(n_rep)) {
    Xp <- Xte_imp
    Xp[, jf] <- sample(Xp[, jf])
    pp <- predict_bagged(fit_imp, Xp)
    drops[r] <- base_acc - bal_acc(yte_imp, pp)
  }
  perm_importance[jf] <- mean(drops)
}

imp_order <- order(perm_importance, decreasing = TRUE)
top10_idx <- imp_order[1:10]
top10_names <- fnames[top10_idx]
top10_groups <- feat_group[top10_idx]
top10_imp <- perm_importance[top10_idx]
n_jump_ccdf_in_top10 <- sum(top10_groups == "ccdf_jump")
cat("top10:\n"); print(data.frame(f=top10_names, g=top10_groups, i=round(top10_imp,4)))

# Eps robustness
cat("Eps sensitivity...\n")
eps_accs <- numeric(length(eps_keys)); names(eps_accs) <- eps_keys
eps_ci_low <- numeric(length(eps_keys)); names(eps_ci_low) <- eps_keys
eps_ci_high <- numeric(length(eps_keys)); names(eps_ci_high) <- eps_keys
eps_per_fold_list <- list()
for (ek in eps_keys) {
  cat("  eps =", ek, "...\n")
  r <- cv_bal_acc(X_eps[[ek]], y_bin, 5, 7, n_trees = NT)
  eps_accs[ek] <- r$mean
  ci <- binom_ci(r$mean, n_total)
  eps_ci_low[ek] <- ci[1]; eps_ci_high[ek] <- ci[2]
  eps_per_fold_list[[ek]] <- r$per_fold
}
eps_spread <- max(eps_accs) - min(eps_accs)
set.seed(123)
k_folds <- length(eps_per_fold_list[[1]])
boot_spread <- replicate(2000, {
  idx <- sample(seq_len(k_folds), replace = TRUE)
  mat <- sapply(eps_per_fold_list, function(v) mean(v[idx]))
  max(mat) - min(mat)
})
spread_ci <- quantile(boot_spread, c(0.025, 0.975))
cat("eps:", eps_accs, " spread =", eps_spread, "\n")

# Family
cat("Family classifier...\n")
res_fam <- cv_bal_acc(X_time, family_label, 5, 7, n_trees = NT)
acc_fam <- res_fam$mean
ci_fam <- binom_ci(acc_fam, n_total)
cat("family acc =", acc_fam, "\n")

# Spearman per family
families_uniq <- unique(family_label)
sp_rows <- list()
for (fm in families_uniq) {
  sel <- family_label == fm
  pv <- param_value[sel]
  if (length(unique(pv)) < 2) next
  for (fn in fnames) {
    v <- X_time[sel, fn]
    good <- is.finite(v) & is.finite(pv)
    if (sum(good) < 5) {
      sp_rows[[length(sp_rows)+1]] <- data.frame(
        family=fm, feature=fn, rho=NA, p=NA, n=sum(good),
        stringsAsFactors=FALSE)
      next
    }
    ct <- suppressWarnings(cor.test(v[good], pv[good], method = "spearman"))
    sp_rows[[length(sp_rows)+1]] <- data.frame(
      family=fm, feature=fn, rho=unname(ct$estimate), p=ct$p.value,
      n=sum(good), stringsAsFactors=FALSE)
  }
}
sp_df <- do.call(rbind, sp_rows)
valid <- !is.na(sp_df$p)
sp_df$p_bh <- NA
if (any(valid)) sp_df$p_bh[valid] <- p.adjust(sp_df$p[valid], method = "BH")
top_sp_list <- list()
for (fm in families_uniq) {
  s <- sp_df[sp_df$family == fm & !is.na(sp_df$rho), ]
  if (nrow(s) == 0) next
  s <- s[order(-abs(s$rho)), ]
  top_sp_list[[fm]] <- list(feature = s$feature[1], rho = s$rho[1],
                            p_bh = s$p_bh[1])
}

# Context
bins_uniq <- unique(param_bin)
ctx_levels <- c("saddle", "plateau", "peak", "mixed")
ctx_tbl <- matrix(0, nrow = length(bins_uniq), ncol = length(ctx_levels),
                  dimnames = list(bins_uniq, ctx_levels))
for (i in seq_len(n_total)) {
  ccs <- dat[[i]]$context_counts_primary
  if (is.null(ccs)) next
  for (cl in ctx_levels) {
    v <- ccs[[cl]]
    if (!is.null(v)) ctx_tbl[param_bin[i], cl] <- ctx_tbl[param_bin[i], cl] + as.numeric(v)
  }
}
ctx_frac <- ctx_tbl / pmax(rowSums(ctx_tbl), 1)
chi <- suppressWarnings(chisq.test(ctx_tbl + 1e-9))
chi_p <- chi$p.value; chi_stat <- unname(chi$statistic)
context_dominant_in_every_bin <- any(apply(ctx_frac, 2, function(col) all(col > 0.9)))

# NaN rates
scalars_for_nan <- c("log_mean", "log_var", "hill", "ks_expon")
nan_rows <- list()
for (obs in c("time", "subs")) {
  for (fn in scalars_for_nan) {
    v2_nans <- 0; v3_nans <- 0; cnt <- 0
    for (i in seq_len(n_total)) {
      nv <- if (obs == "time") dat[[i]]$nan_indicator_time else dat[[i]]$nan_indicator_subs
      if (is.null(nv) || is.null(nv[[fn]])) next
      cnt <- cnt + 1
      v2_nans <- v2_nans + isTRUE(nv[[fn]]$v2_rule_nan)
      v3_nans <- v3_nans + isTRUE(nv[[fn]]$v3_rule_nan)
    }
    nan_rows[[length(nan_rows)+1]] <- list(
      observable=obs, feature=fn,
      v2_rate = if (cnt) v2_nans/cnt else NA,
      v3_rate = if (cnt) v3_nans/cnt else NA)
  }
}
mc_b <- 0; mc_c <- 0
for (fn in scalars_for_nan) {
  for (i in seq_len(n_total)) {
    nv <- dat[[i]]$nan_indicator_time
    if (is.null(nv) || is.null(nv[[fn]])) next
    a <- isTRUE(nv[[fn]]$v2_rule_nan); b <- isTRUE(nv[[fn]]$v3_rule_nan)
    if (a && !b) mc_b <- mc_b + 1
    if (!a && b) mc_c <- mc_c + 1
  }
}
mcnemar_stat <- if ((mc_b + mc_c) > 0) ((abs(mc_b - mc_c) - 1)^2) / (mc_b + mc_c) else NA
mcnemar_p <- if (!is.na(mcnemar_stat)) pchisq(mcnemar_stat, df = 1, lower.tail = FALSE) else NA

# Confusion
classes <- sort(unique(y_bin))
conf_mat <- matrix(0, nrow = length(classes), ncol = length(classes),
                   dimnames = list(classes, classes))
for (i in seq_len(n_total)) {
  t <- y_bin[i]; p <- res_full_time$preds[i]
  if (is.na(p) || p == "") next
  if (!(p %in% classes)) next
  conf_mat[t, p] <- conf_mat[t, p] + 1
}
conf_row_norm <- conf_mat / pmax(rowSums(conf_mat), 1)

primary_pass <- ci_full_time[2] >= 0.40
secondary_pass <- ci_delta_time_subs[1] > 0
tertiary_pass <- n_jump_ccdf_in_top10 >= 1
robustness_pass <- eps_spread <= 0.05
family_pass <- ci_fam[2] >= 0.70
joint_null_pass <- (!primary_pass) && (!secondary_pass) && (!tertiary_pass)

cat("\n--- VERDICTS ---\n")
cat("primary:", primary_pass, " secondary:", secondary_pass,
    " tertiary:", tertiary_pass, " robustness:", robustness_pass,
    " family:", family_pass, "\n")

ctx_frac_list <- lapply(seq_len(nrow(ctx_frac)), function(i) {
  list(bin = rownames(ctx_frac)[i],
       saddle = unname(ctx_frac[i, "saddle"]),
       plateau = unname(ctx_frac[i, "plateau"]),
       peak = unname(ctx_frac[i, "peak"]),
       mixed = unname(ctx_frac[i, "mixed"]),
       n_epochs = sum(ctx_tbl[i, ]))
})

conf_list <- lapply(seq_len(nrow(conf_row_norm)), function(i) {
  l <- as.list(unname(conf_row_norm[i, ]))
  names(l) <- colnames(conf_row_norm)
  l$true_class <- rownames(conf_row_norm)[i]
  l
})

top10_list <- lapply(seq_along(top10_names), function(i) {
  list(rank = i, feature = top10_names[i],
       group = unname(top10_groups[i]),
       importance = unname(top10_imp[i]))
})

all_imp_list <- lapply(seq_along(fnames), function(i) {
  list(feature = fnames[i], group = unname(feat_group[i]),
       importance = unname(perm_importance[i]))
})

stats <- list(
  n_total = n_total,
  n_classes = n_classes,
  classes = classes,
  families = unique(family_label),
  classifier = "bagged-rpart (RF surrogate, 80 trees / fold, mtry=sqrt(p))",
  acc_full_time = acc_full_time,
  ci_full_time_low = unname(ci_full_time[1]),
  ci_full_time_high = unname(ci_full_time[2]),
  primary_threshold = 0.40,
  chance_level = 1 / n_classes,
  primary_pass = primary_pass,
  acc_full_subs = acc_full_subs,
  ci_full_subs_low = unname(ci_full_subs[1]),
  ci_full_subs_high = unname(ci_full_subs[2]),
  delta_time_minus_subs = mean_delta_time_subs,
  delta_time_minus_subs_ci_low = unname(ci_delta_time_subs[1]),
  delta_time_minus_subs_ci_high = unname(ci_delta_time_subs[2]),
  secondary_pass = secondary_pass,
  acc_4scalar_time = acc_4_time,
  acc_ccdf_time = acc_ccdf_time,
  acc_4scalar_subs = acc_4_subs,
  acc_ccdf_subs = acc_ccdf_subs,
  abl_full_minus_4 = abl_full_minus_4$mean,
  abl_full_minus_4_ci_low = unname(abl_full_minus_4$ci[1]),
  abl_full_minus_4_ci_high = unname(abl_full_minus_4$ci[2]),
  abl_ccdf_minus_4 = abl_ccdf_minus_4$mean,
  abl_ccdf_minus_4_ci_low = unname(abl_ccdf_minus_4$ci[1]),
  abl_ccdf_minus_4_ci_high = unname(abl_ccdf_minus_4$ci[2]),
  abl_full_minus_ccdf = abl_full_minus_ccdf$mean,
  abl_full_minus_ccdf_ci_low = unname(abl_full_minus_ccdf$ci[1]),
  abl_full_minus_ccdf_ci_high = unname(abl_full_minus_ccdf$ci[2]),
  n_jump_ccdf_in_top10 = n_jump_ccdf_in_top10,
  tertiary_pass = tertiary_pass,
  top10_features = top10_list,
  perm_importance_all = all_imp_list,
  perm_base_acc_holdout = base_acc,
  eps_acc_0_001 = unname(eps_accs["0.001"]),
  eps_acc_0_005 = unname(eps_accs["0.005"]),
  eps_acc_0_02 = unname(eps_accs["0.02"]),
  eps_ci_0_001_low = unname(eps_ci_low["0.001"]),
  eps_ci_0_001_high = unname(eps_ci_high["0.001"]),
  eps_ci_0_005_low = unname(eps_ci_low["0.005"]),
  eps_ci_0_005_high = unname(eps_ci_high["0.005"]),
  eps_ci_0_02_low = unname(eps_ci_low["0.02"]),
  eps_ci_0_02_high = unname(eps_ci_high["0.02"]),
  eps_spread = eps_spread,
  eps_spread_ci_low = unname(spread_ci[1]),
  eps_spread_ci_high = unname(spread_ci[2]),
  robustness_threshold = 0.05,
  robustness_pass = robustness_pass,
  acc_family = acc_fam,
  ci_family_low = unname(ci_fam[1]),
  ci_family_high = unname(ci_fam[2]),
  family_threshold = 0.70,
  family_pass = family_pass,
  top_spearman_per_family = top_sp_list,
  context_fractions = ctx_frac_list,
  chi_square_stat = chi_stat,
  chi_square_p = chi_p,
  context_dominant_in_every_bin = context_dominant_in_every_bin,
  nan_rates = nan_rows,
  mcnemar_stat = mcnemar_stat,
  mcnemar_p = mcnemar_p,
  mcnemar_b_v2only = mc_b,
  mcnemar_c_v3only = mc_c,
  confusion_row_normalized = conf_list,
  joint_null_strong = joint_null_pass,
  parent_acc_full = 0.3406,
  parent_acc_full_ci_low = 0.2979,
  parent_acc_full_ci_high = 0.3828,
  parent_family_acc = 0.7812
)

write_json(stats, "stats.json", auto_unbox = TRUE, pretty = TRUE, na = "null")
cat("stats.json written\n")
