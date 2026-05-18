library(jsonlite)
library(dplyr)

raw <- fromJSON("result.json", simplifyVector = FALSE)$results

# Build per-(target,condition) table
rows <- lapply(raw, function(r) {
  data.frame(
    target_id = r$target_id,
    primitives = r$primitives,
    sigma_late = r$pooled_sigma_late,
    mu_late = r$pooled_mu_late,
    n_late = r$n_late_events,
    sigma_early = r$pooled_sigma_early,
    mu_early = r$pooled_mu_early,
    n_early = r$n_early_events,
    csn_R = r$csn_late$R,
    csn_p = r$csn_late$p,
    csn_n_tail = r$csn_late$n_tail,
    nn_mean_class = r$neutral_proxy$mean_class_size_weighted,
    nn_max_class = r$neutral_proxy$max_class_size,
    stringsAsFactors = FALSE
  )
})
df <- do.call(rbind, rows)

# Reshape: per target, both conditions
targets <- unique(df$target_id)
per_target <- lapply(targets, function(t) {
  m <- df[df$target_id == t & df$primitives == "minimal", ]
  r <- df[df$target_id == t & df$primitives == "rich", ]
  list(
    target_id = t,
    sigma_minimal = m$sigma_late,
    sigma_rich = r$sigma_late,
    delta_sigma = r$sigma_late - m$sigma_late,
    nn_minimal = m$nn_mean_class,
    nn_rich = r$nn_mean_class,
    nn_ratio = r$nn_mean_class / m$nn_mean_class,
    csn_R_rich = r$csn_R,
    csn_p_rich = r$csn_p,
    csn_n_tail_rich = r$csn_n_tail
  )
})

sigma_minimal <- sapply(per_target, function(x) x$sigma_minimal)
sigma_rich <- sapply(per_target, function(x) x$sigma_rich)
delta_sigma <- sapply(per_target, function(x) x$delta_sigma)
nn_ratio <- sapply(per_target, function(x) x$nn_ratio)
target_ids <- sapply(per_target, function(x) x$target_id)

# Test 1: paired Wilcoxon (one-sided rich > minimal)
wilc <- wilcox.test(sigma_rich, sigma_minimal, paired = TRUE, alternative = "greater", exact = FALSE)

# Test 2: Spearman + bootstrap CI
cor_res <- cor.test(delta_sigma, nn_ratio, method = "spearman", exact = FALSE)
set.seed(42)
B <- 5000
boot_rho <- replicate(B, {
  idx <- sample(seq_along(delta_sigma), replace = TRUE)
  if (length(unique(delta_sigma[idx])) < 2 || length(unique(nn_ratio[idx])) < 2) return(NA_real_)
  suppressWarnings(cor(delta_sigma[idx], nn_ratio[idx], method = "spearman"))
})
boot_rho <- boot_rho[is.finite(boot_rho)]
ci_low <- as.numeric(quantile(boot_rho, 0.025))
ci_high <- as.numeric(quantile(boot_rho, 0.975))

# Test 3: CSN log-normal vs exponential — already done in coder, R<0 favors log-normal
# Report per-target rich condition R, p, n_tail
csn_table <- lapply(per_target, function(x) {
  list(
    target_id = x$target_id,
    R = x$csn_R_rich,
    p = x$csn_p_rich,
    n_tail = x$csn_n_tail_rich,
    favors_lognormal = (x$csn_R_rich < 0) && (x$csn_p_rich < 0.05),
    significant_at_05 = (x$csn_p_rich < 0.05)
  )
})

# Pass/fail per null criterion
# (i) Wilcoxon p>0.05 OR direction mixed => null
direction_positive_count <- sum(delta_sigma > 0)
direction_mixed <- !(all(delta_sigma > 0) || all(delta_sigma < 0))
wilc_pass <- (wilc$p.value <= 0.05) && !direction_mixed && (mean(sigma_rich) > mean(sigma_minimal) || median(delta_sigma) > 0)

# (ii) Spearman CI crosses 0 => null
ci_crosses_zero <- (ci_low <= 0) && (ci_high >= 0)
spearman_pass <- !ci_crosses_zero && (cor_res$estimate > 0)

# (iii) CSN: requires significant log-normal vs exponential in EVERY target (rich)
# "not significantly distinguishable in any target" = null
# So pass = at least one target has p<0.05 AND favors log-normal? Plan: "not significantly distinguishable in any target (log-likelihood ratio not significant at α=0.05)" — read as: if NO target shows significance, null. Pass = at least one target significant. Stronger reading: every target shows significance. Use strict: all targets.
csn_p_vec <- sapply(csn_table, function(x) x$p)
csn_R_vec <- sapply(csn_table, function(x) x$R)
csn_all_significant <- all(csn_p_vec < 0.05)
csn_any_significant <- any(csn_p_vec < 0.05)
csn_pass <- csn_all_significant  # strict reading

# Overall mechanism demonstrated requires all three to pass
mechanism_demonstrated <- wilc_pass && spearman_pass && csn_pass

stats <- list(
  # per-target table
  target_ids = as.list(target_ids),
  sigma_minimal_by_target = as.list(sigma_minimal),
  sigma_rich_by_target = as.list(sigma_rich),
  delta_sigma_by_target = as.list(delta_sigma),
  nn_minimal_by_target = as.list(sapply(per_target, function(x) x$nn_minimal)),
  nn_rich_by_target = as.list(sapply(per_target, function(x) x$nn_rich)),
  nn_ratio_by_target = as.list(nn_ratio),

  # Wilcoxon
  wilcoxon_V = as.numeric(wilc$statistic),
  wilcoxon_p = as.numeric(wilc$p.value),
  wilcoxon_alternative = "greater (rich > minimal)",
  mean_sigma_minimal = mean(sigma_minimal),
  mean_sigma_rich = mean(sigma_rich),
  median_delta_sigma = median(delta_sigma),
  n_targets_positive_shift = direction_positive_count,
  n_targets_total = length(delta_sigma),
  direction_mixed = direction_mixed,
  wilcoxon_pass = wilc_pass,

  # Spearman
  spearman_rho = as.numeric(cor_res$estimate),
  spearman_p = as.numeric(cor_res$p.value),
  spearman_ci_low = ci_low,
  spearman_ci_high = ci_high,
  spearman_ci_crosses_zero = ci_crosses_zero,
  spearman_n_boot = length(boot_rho),
  spearman_pass = spearman_pass,

  # CSN per target (rich)
  csn_per_target = csn_table,
  csn_all_significant = csn_all_significant,
  csn_any_significant = csn_any_significant,
  csn_n_significant = sum(csn_p_vec < 0.05),
  csn_pass = csn_pass,

  # Overall
  mechanism_demonstrated = mechanism_demonstrated,

  # Ranges for plots
  nn_ratio_range = c(min(nn_ratio), max(nn_ratio)),
  delta_sigma_range = c(min(delta_sigma), max(delta_sigma))
)

write_json(stats, "stats.json", auto_unbox = TRUE, pretty = TRUE, digits = 6)
cat("WROTE stats.json\n")
print(stats)
