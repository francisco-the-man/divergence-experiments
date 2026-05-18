library(jsonlite)
library(ggplot2)

raw <- fromJSON("result.json", simplifyVector = FALSE)$results

# Build per-(target, condition) records
recs <- lapply(raw, function(r) {
  list(
    target_id = r$target_id,
    primitives = r$primitives,
    sigma_late = r$pooled_sigma_late,
    sigma_early = r$pooled_sigma_early,
    mu_late = r$pooled_mu_late,
    n_late = r$n_late_events,
    nn_size = r$neutral_proxy$mean_class_size_weighted,
    durations_late = unlist(r$pooled_late_durations_sample),
    csn_R = r$csn_late$R,
    csn_p = r$csn_late$p
  )
})

targets <- unique(sapply(recs, function(x) x$target_id))

# Build paired sigma vectors
sigma_min <- numeric(length(targets))
sigma_rich <- numeric(length(targets))
nn_min <- numeric(length(targets))
nn_rich <- numeric(length(targets))
csn_R_rich <- numeric(length(targets))
csn_p_rich <- numeric(length(targets))
names(sigma_min) <- targets
names(sigma_rich) <- targets
names(nn_min) <- targets
names(nn_rich) <- targets
names(csn_R_rich) <- targets
names(csn_p_rich) <- targets

for (r in recs) {
  if (r$primitives == "minimal") {
    sigma_min[r$target_id] <- r$sigma_late
    nn_min[r$target_id] <- r$nn_size
  } else {
    sigma_rich[r$target_id] <- r$sigma_late
    nn_rich[r$target_id] <- r$nn_size
    csn_R_rich[r$target_id] <- r$csn_R
    csn_p_rich[r$target_id] <- r$csn_p
  }
}

delta_sigma <- sigma_rich - sigma_min
nn_ratio <- nn_rich / nn_min

# Test 1: paired Wilcoxon, one-sided greater
w <- wilcox.test(sigma_rich, sigma_min, paired = TRUE, alternative = "greater")

# Test 2: Spearman rho with bootstrap CI
sp <- suppressWarnings(cor.test(delta_sigma, nn_ratio, method = "spearman"))

set.seed(42)
n_boot <- 5000
n_t <- length(targets)
boot_rhos <- numeric(n_boot)
for (i in seq_len(n_boot)) {
  idx <- sample(seq_len(n_t), n_t, replace = TRUE)
  if (length(unique(idx)) < 2) {
    boot_rhos[i] <- NA
    next
  }
  boot_rhos[i] <- suppressWarnings(cor(delta_sigma[idx], nn_ratio[idx], method = "spearman"))
}
boot_rhos <- boot_rhos[is.finite(boot_rhos)]
ci <- quantile(boot_rhos, c(0.025, 0.975), na.rm = TRUE)

# Test 3: CSN per-target rich condition (already computed in result.json as csn_late R,p)
csn_per_target <- lapply(targets, function(t) {
  list(
    target_id = t,
    R = csn_R_rich[[t]],
    p = csn_p_rich[[t]],
    significant = csn_p_rich[[t]] < 0.05,
    lognormal_preferred = (csn_R_rich[[t]] < 0) && (csn_p_rich[[t]] < 0.05)
  )
})

# Pass/fail per null_result_criteria
crit_i_pass <- (w$p.value < 0.05) && all(delta_sigma > 0) == FALSE  
# direction: hypothesis says rich>minimal. Mixed direction = fail.
direction_consistent <- all(delta_sigma > 0)
# null criterion (i) fails if p>0.05 OR mixed direction
crit_i_fail <- (w$p.value > 0.05) || !direction_consistent
crit_ii_fail <- (ci[1] <= 0) && (ci[2] >= 0)
# crit iii: distinguishable in any target = at least one significant. Fails if NONE significant.
any_csn_sig <- any(csn_p_rich < 0.05)
crit_iii_fail <- !any_csn_sig

mechanism_demonstrated <- !crit_i_fail && !crit_ii_fail && !crit_iii_fail

# Per-target table
per_target_table <- lapply(targets, function(t) {
  list(
    target_id = t,
    sigma_minimal = sigma_min[[t]],
    sigma_rich = sigma_rich[[t]],
    delta_sigma = delta_sigma[[t]],
    nn_minimal = nn_min[[t]],
    nn_rich = nn_rich[[t]],
    nn_ratio = nn_ratio[[t]]
  )
})

stats <- list(
  targets = targets,
  sigma_minimal = as.list(sigma_min),
  sigma_rich = as.list(sigma_rich),
  delta_sigma = as.list(delta_sigma),
  nn_minimal = as.list(nn_min),
  nn_rich = as.list(nn_rich),
  nn_ratio = as.list(nn_ratio),
  per_target_table = per_target_table,
  wilcoxon_V = unname(w$statistic),
  wilcoxon_p = w$p.value,
  wilcoxon_alternative = "greater",
  direction_consistent = direction_consistent,
  n_targets_rich_gt_minimal = sum(delta_sigma > 0),
  n_targets_total = length(targets),
  spearman_rho = unname(sp$estimate),
  spearman_p = sp$p.value,
  spearman_ci_low = unname(ci[1]),
  spearman_ci_high = unname(ci[2]),
  spearman_ci_crosses_zero = unname((ci[1] <= 0) && (ci[2] >= 0)),
  csn_per_target = csn_per_target,
  n_targets_csn_significant = sum(csn_p_rich < 0.05),
  any_csn_significant = any_csn_sig,
  crit_i_fail = crit_i_fail,
  crit_ii_fail = crit_ii_fail,
  crit_iii_fail = crit_iii_fail,
  mechanism_demonstrated = mechanism_demonstrated,
  wilcoxon_pass = !crit_i_fail,
  spearman_pass = !crit_ii_fail,
  csn_pass = !crit_iii_fail
)

write_json(stats, "stats.json", auto_unbox = TRUE, pretty = TRUE)
cat("done\n")
