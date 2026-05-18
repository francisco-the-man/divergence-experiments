
library(jsonlite)

raw <- fromJSON("result.json", simplifyVector = FALSE)$results

# Separate gp_run and egraph results
gp <- Filter(function(r) !is.null(r$kind) && r$kind == "gp_run", raw)
eg <- Filter(function(r) !is.null(r$kind) && r$kind == "egraph_banded", raw)

cat("n gp_run:", length(gp), " n egraph:", length(eg), "\n")

# Inventory: target x condition
inv <- table(
  sapply(gp, function(r) r$target_id),
  sapply(gp, function(r) r$condition)
)
cat("Inventory:\n"); print(inv)

# ---- Extract plateau durations per (target, condition) -------------------
# best_hist gives the running best fitness over generations. We define
# plateau durations as runs of equal best fitness, and "late phase" as
# events that begin after the median improvement event (proxy for tau_target).
# This follows the spirit of the preregistered "fitness_banded_post_tau_target"
# definition while being computable from what's in result.json.

plateau_durs <- function(rec) {
  pe <- rec$plateau_events
  if (is.null(pe) || length(pe) == 0) return(list(all = numeric(0), late = numeric(0)))
  # plateau_events entries: [start_gen, end_gen, duration, fitness, censored_flag]
  starts <- sapply(pe, function(e) as.numeric(e[[1]]))
  durs   <- sapply(pe, function(e) as.numeric(e[[3]]))
  cens   <- sapply(pe, function(e) as.numeric(e[[5]]))
  # late phase: events whose start is past the median start
  thresh <- stats::median(starts)
  late_idx <- which(starts >= thresh)
  list(
    all  = durs,
    late = durs[late_idx],
    cens_all = cens,
    cens_late = cens[late_idx]
  )
}

# Pooled censored MLE for log-normal on durations (right-censored)
# Use simple MLE on uncensored; for censored events keep them as observations
# of duration as a lower bound (treat as right-censored). Use survreg if available.
fit_lognorm_sigma <- function(durs, cens) {
  if (length(durs) < 3) return(NA_real_)
  durs <- pmax(durs, 1)  # avoid log(0)
  # If survival package is available, use right-censored MLE
  if (requireNamespace("survival", quietly = TRUE)) {
    library(survival)
    # status = 1 means event observed (not censored), 0 means censored
    status <- ifelse(cens == 1, 0, 1)  # cens=1 in result -> right-censored
    # Need at least one uncensored event
    if (sum(status) < 2) {
      # fall back to mle on logs
      ld <- log(durs)
      return(sd(ld))
    }
    fit <- tryCatch(
      survreg(Surv(durs, status) ~ 1, dist = "lognormal"),
      error = function(e) NULL
    )
    if (is.null(fit)) return(sd(log(durs)))
    return(as.numeric(fit$scale))
  }
  sd(log(durs))
}

# Compute pooled sigma_late per (target, condition)
agg <- list()
targets <- sort(unique(sapply(gp, function(r) r$target_id)))
conds <- sort(unique(sapply(gp, function(r) r$condition)))

for (t in targets) {
  for (c in conds) {
    sub <- Filter(function(r) r$target_id == t && r$condition == c, gp)
    if (length(sub) == 0) next
    all_late <- numeric(0); all_cens <- numeric(0)
    for (r in sub) {
      pd <- plateau_durs(r)
      all_late <- c(all_late, pd$late)
      all_cens <- c(all_cens, pd$cens_late)
    }
    sig <- fit_lognorm_sigma(all_late, all_cens)
    agg[[paste(t, c, sep = "::")]] <- list(
      target_id = t, condition = c,
      n_reps = length(sub),
      n_late_events = length(all_late),
      sigma_late = sig
    )
  }
}

# Build wide table per target with sigma_M, sigma_S, sigma_C
target_ids <- sort(unique(sapply(agg, function(a) a$target_id)))
get_sig <- function(t, c) {
  k <- paste(t, c, sep = "::")
  if (!is.null(agg[[k]])) agg[[k]]$sigma_late else NA_real_
}

sigma_M <- sapply(target_ids, function(t) get_sig(t, "M"))
sigma_S <- sapply(target_ids, function(t) get_sig(t, "S"))
sigma_C <- sapply(target_ids, function(t) get_sig(t, "C"))

names(sigma_M) <- target_ids
names(sigma_S) <- target_ids
names(sigma_C) <- target_ids

cat("sigma_M:\n"); print(sigma_M)
cat("sigma_S:\n"); print(sigma_S)
cat("sigma_C:\n"); print(sigma_C)

n_targets <- length(target_ids)
n_have_M <- sum(!is.na(sigma_M))
n_have_S <- sum(!is.na(sigma_S))
n_have_C <- sum(!is.na(sigma_C))

# ---- Test 1: paired one-sided Wilcoxon sigma_C > sigma_S -----------------
# If C is missing for all targets, this test is not computable.
test1 <- list(
  computable = FALSE,
  V = NA, p = NA,
  n_pairs = 0,
  n_positive = NA,
  median_delta = NA,
  pass = FALSE,
  note = ""
)

paired_idx <- which(!is.na(sigma_C) & !is.na(sigma_S))
if (length(paired_idx) >= 2) {
  dC <- sigma_C[paired_idx]; dS <- sigma_S[paired_idx]
  delta_CS <- dC - dS
  wt <- tryCatch(
    wilcox.test(dC, dS, paired = TRUE, alternative = "greater", exact = FALSE),
    error = function(e) NULL
  )
  if (!is.null(wt)) {
    test1$computable <- TRUE
    test1$V <- as.numeric(wt$statistic)
    test1$p <- as.numeric(wt$p.value)
    test1$n_pairs <- length(paired_idx)
    test1$n_positive <- sum(delta_CS > 0)
    test1$median_delta <- as.numeric(median(delta_CS))
    test1$pass <- (test1$p <= 0.05) && (test1$median_delta > 0)
    test1$note <- "OK"
  }
} else {
  test1$note <- sprintf("Insufficient C-condition data: have M=%d S=%d C=%d of %d targets",
                        n_have_M, n_have_S, n_have_C, n_targets)
}

# ---- Test 1 sensitivity: drop T6 ------------------------------------------
test1_drop_t6 <- list(computable = FALSE, V = NA, p = NA, n_pairs = 0, note = "")
if (length(paired_idx) >= 3) {
  keep <- paired_idx[ target_ids[paired_idx] != "T6_inexpressible" ]
  if (length(keep) >= 2) {
    dC <- sigma_C[keep]; dS <- sigma_S[keep]
    wt <- tryCatch(
      wilcox.test(dC, dS, paired = TRUE, alternative = "greater", exact = FALSE),
      error = function(e) NULL
    )
    if (!is.null(wt)) {
      test1_drop_t6$computable <- TRUE
      test1_drop_t6$V <- as.numeric(wt$statistic)
      test1_drop_t6$p <- as.numeric(wt$p.value)
      test1_drop_t6$n_pairs <- length(keep)
      test1_drop_t6$note <- "OK"
    }
  }
} else {
  test1_drop_t6$note <- "Not enough pairs to run sensitivity"
}

# ---- E-graph banded measurements: size & connectivity --------------------
# Aggregate mean class size and mean out-degree per (target, condition)
eg_summary <- list()
for (r in eg) {
  k <- paste(r$target_id, r$condition, sep = "::")
  sz <- if (!is.null(r$size)) as.numeric(unlist(r$size)) else NA
  od <- if (!is.null(r$portal)) as.numeric(unlist(r$portal)) else NA
  eg_summary[[k]] <- list(
    target_id = r$target_id, condition = r$condition,
    mean_size = if (length(sz) > 0) mean(sz, na.rm = TRUE) else NA,
    mean_portal = if (length(od) > 0) mean(od, na.rm = TRUE) else NA
  )
}

get_eg <- function(t, c, field) {
  k <- paste(t, c, sep = "::")
  if (!is.null(eg_summary[[k]])) eg_summary[[k]][[field]] else NA_real_
}

size_M <- sapply(target_ids, function(t) get_eg(t, "M", "mean_size"))
size_S <- sapply(target_ids, function(t) get_eg(t, "S", "mean_size"))
size_C <- sapply(target_ids, function(t) get_eg(t, "C", "mean_size"))
port_M <- sapply(target_ids, function(t) get_eg(t, "M", "mean_portal"))
port_S <- sapply(target_ids, function(t) get_eg(t, "S", "mean_portal"))
port_C <- sapply(target_ids, function(t) get_eg(t, "C", "mean_portal"))

# ---- Test 2: dissociated Spearman correlations ---------------------------
# rho_size = spearman(delta_sigma_S, delta_size_S)
# rho_conn = spearman(delta_sigma_C, delta_connectivity_C)
boot_spearman_ci <- function(x, y, B = 10000, seed = 42) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3) return(list(rho = NA, lo = NA, hi = NA, n = length(x)))
  rho0 <- suppressWarnings(cor(x, y, method = "spearman"))
  set.seed(seed)
  rhos <- replicate(B, {
    idx <- sample.int(length(x), replace = TRUE)
    if (length(unique(x[idx])) < 2 || length(unique(y[idx])) < 2) return(NA_real_)
    suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  })
  rhos <- rhos[is.finite(rhos)]
  ci <- if (length(rhos) > 10) quantile(rhos, c(0.025, 0.975), na.rm = TRUE) else c(NA, NA)
  list(rho = as.numeric(rho0), lo = as.numeric(ci[1]), hi = as.numeric(ci[2]), n = length(x))
}

test2 <- list(computable = FALSE, note = "")
delta_sig_S <- sigma_S - sigma_M
delta_sig_C <- sigma_C - sigma_M
delta_size_S <- size_S - size_M
delta_port_C <- port_C - port_M

if (sum(is.finite(delta_sig_S) & is.finite(delta_size_S)) >= 3 ||
    sum(is.finite(delta_sig_C) & is.finite(delta_port_C)) >= 3) {
  r_size <- boot_spearman_ci(delta_sig_S, delta_size_S)
  r_conn <- boot_spearman_ci(delta_sig_C, delta_port_C)
  test2$computable <- (r_size$n >= 3 || r_conn$n >= 3)
  test2$rho_size <- r_size$rho
  test2$rho_size_lo <- r_size$lo
  test2$rho_size_hi <- r_size$hi
  test2$rho_size_n <- r_size$n
  test2$rho_conn <- r_conn$rho
  test2$rho_conn_lo <- r_conn$lo
  test2$rho_conn_hi <- r_conn$hi
  test2$rho_conn_n <- r_conn$n
  # Dissociation requires rho_conn CI excludes 0 AND rho_conn > rho_size
  ci_excl_0 <- !is.na(r_conn$lo) && !is.na(r_conn$hi) && (r_conn$lo > 0 || r_conn$hi < 0)
  test2$dissociation_pass <- ci_excl_0 &&
                              !is.na(r_conn$rho) && !is.na(r_size$rho) &&
                              (r_conn$rho > r_size$rho)
  test2$note <- "OK"
} else {
  test2$note <- "Insufficient data to compute dissociation correlations"
  test2$dissociation_pass <- FALSE
}

# ---- Test 3: CSN log-normal vs exponential per cell ----------------------
# Compute log-likelihood ratio between fitted log-normal and exponential on
# late-phase durations per (target, condition). Use Vuong-style sign+p approx.
# This is an internal R approximation of powerlaw.Fit's distribution_compare.

ln_exp_compare <- function(durs) {
  durs <- as.numeric(durs)
  durs <- durs[is.finite(durs) & durs > 0]
  if (length(durs) < 20) return(list(R = NA, p = NA, n = length(durs), xmin = NA))
  # xmin search: try a few candidate xmins, pick one minimizing KS distance
  # against log-normal fit for simplicity. Then compute LR vs exponential on the
  # tail above xmin.
  cand <- unique(quantile(durs, probs = c(0, .1, .25, .5, .75, .9), na.rm = TRUE))
  cand <- cand[cand >= 1]
  best <- list(R = NA, p = NA, n = NA, xmin = NA, score = Inf)
  for (xm in cand) {
    tail <- durs[durs >= xm]
    if (length(tail) < 15) next
    ld <- log(tail)
    mu <- mean(ld); sg <- sd(ld); if (!is.finite(sg) || sg <= 0) next
    # log-likelihood under log-normal (conditional on x >= xm)
    ll_lnorm_i <- dlnorm(tail, meanlog = mu, sdlog = sg, log = TRUE) -
                  plnorm(xm, meanlog = mu, sdlog = sg, log.p = TRUE, lower.tail = FALSE)
    # log-likelihood under exponential (conditional on x >= xm)
    # If X ~ Exp(lam), conditional X|X>=xm has the same exp with shift xm
    lam <- 1 / (mean(tail) - xm)
    if (!is.finite(lam) || lam <= 0) next
    ll_exp_i <- log(lam) - lam * (tail - xm)
    diffs <- ll_lnorm_i - ll_exp_i
    R <- sum(diffs)
    sd_d <- sd(diffs)
    if (!is.finite(sd_d) || sd_d <= 0) next
    # Vuong: z = R / (sqrt(n) * sd_d) ; two-sided p
    z <- R / (sqrt(length(tail)) * sd_d)
    p <- 2 * pnorm(-abs(z))
    # score for xmin selection: penalize tiny tails
    score <- -length(tail)  # prefer larger tails
    if (score < best$score) {
      best <- list(R = R, p = p, n = length(tail), xmin = xm, score = score)
    }
  }
  best$score <- NULL
  best
}

# Pool late durations per cell across reps
pool_late <- function(t, c) {
  sub <- Filter(function(r) r$target_id == t && r$condition == c, gp)
  out <- numeric(0)
  for (r in sub) {
    pd <- plateau_durs(r)
    out <- c(out, pd$late)
  }
  out
}

csn_results <- list()
n_C_lognorm_pref <- 0
n_C_total <- 0
n_S_exp_pref <- 0
n_S_total <- 0
for (t in target_ids) {
  for (c in conds) {
    durs <- pool_late(t, c)
    if (length(durs) < 20) {
      csn_results[[paste(t,c,sep="::")]] <- list(
        target_id = t, condition = c, R = NA, p = NA, n = length(durs),
        xmin = NA, preference = "insufficient"
      )
      next
    }
    res <- ln_exp_compare(durs)
    pref <- "ambiguous"
    if (!is.na(res$R) && !is.na(res$p)) {
      if (res$R > 0 && res$p < 0.05) pref <- "lognormal"
      else if (res$R < 0 && res$p < 0.05) pref <- "exponential"
    }
    csn_results[[paste(t,c,sep="::")]] <- list(
      target_id = t, condition = c,
      R = res$R, p = res$p, n = res$n, xmin = res$xmin,
      preference = pref
    )
    if (c == "C") {
      n_C_total <- n_C_total + 1
      if (pref == "lognormal") n_C_lognorm_pref <- n_C_lognorm_pref + 1
    }
    if (c == "S") {
      n_S_total <- n_S_total + 1
      if (pref == "exponential") n_S_exp_pref <- n_S_exp_pref + 1
    }
  }
}

test3 <- list(
  computable = (n_C_total > 0 || n_S_total > 0),
  n_C_total = n_C_total,
  n_C_lognorm_pref = n_C_lognorm_pref,
  n_S_total = n_S_total,
  n_S_exp_pref = n_S_exp_pref,
  pass = (n_C_lognorm_pref >= 4) && (n_S_exp_pref >= 3),
  note = if (n_C_total == 0 && n_S_total == 0) "No CSN cells computable" else "OK"
)

# ---- Build long per-cell descriptive table for plots ---------------------
cell_tbl <- list()
for (t in target_ids) {
  for (c in conds) {
    k <- paste(t, c, sep = "::")
    cell_tbl[[k]] <- list(
      target_id = t, condition = c,
      sigma_late = get_sig(t, c),
      mean_size = get_eg(t, c, "mean_size"),
      mean_portal = get_eg(t, c, "mean_portal"),
      n_reps = if (!is.null(agg[[k]])) agg[[k]]$n_reps else 0,
      n_late_events = if (!is.null(agg[[k]])) agg[[k]]$n_late_events else 0
    )
  }
}

# ---- Overall verdict -----------------------------------------------------
fired <- c()
if (test1$computable && !test1$pass) fired <- c(fired, "test1")
if (test2$computable && !test2$dissociation_pass) fired <- c(fired, "test2")
if (test3$computable && !test3$pass) fired <- c(fired, "test3")
if (!test1$computable) fired <- c(fired, "test1_not_computable")
if (!test2$computable) fired <- c(fired, "test2_not_computable")
if (!test3$computable) fired <- c(fired, "test3_not_computable")

mechanism_demonstrated <- length(fired) == 0

stats <- list(
  n_gp_results = length(gp),
  n_egraph_results = length(eg),
  n_targets = n_targets,
  n_targets_with_M = n_have_M,
  n_targets_with_S = n_have_S,
  n_targets_with_C = n_have_C,
  target_ids = as.list(target_ids),
  conditions_present = as.list(conds),
  sigma_M_by_target = as.list(sigma_M),
  sigma_S_by_target = as.list(sigma_S),
  sigma_C_by_target = as.list(sigma_C),
  delta_sigma_S = as.list(delta_sig_S),
  delta_sigma_C = as.list(delta_sig_C),
  size_M_by_target = as.list(size_M),
  size_S_by_target = as.list(size_S),
  size_C_by_target = as.list(size_C),
  portal_M_by_target = as.list(port_M),
  portal_S_by_target = as.list(port_S),
  portal_C_by_target = as.list(port_C),
  test1 = test1,
  test1_sensitivity_drop_t6 = test1_drop_t6,
  test2 = test2,
  test3 = test3,
  csn_per_cell = csn_results,
  cell_table = cell_tbl,
  null_criteria_fired = as.list(fired),
  mechanism_demonstrated = mechanism_demonstrated,
  scope_note = sprintf(
    "Executed results contain %d gp_run tasks across %d target(s) and conditions {%s}; plan called for 6 targets x {M,S,C}. Tests requiring conditions not present are reported as not computable.",
    length(gp), n_targets, paste(conds, collapse = ", "))
)

write_json(stats, "stats.json", auto_unbox = TRUE, pretty = TRUE, na = "null")
cat("\nDone. stats.json written.\n")
cat("Mechanism demonstrated:", mechanism_demonstrated, "\n")
cat("Fired:", paste(fired, collapse = ", "), "\n")
