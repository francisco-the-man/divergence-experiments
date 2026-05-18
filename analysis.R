library(jsonlite)
library(dplyr)

raw <- fromJSON("result.json", simplifyVector = FALSE)$results

# Build a tidy frame of per-task summary
rows <- lapply(raw, function(r) {
  data.frame(
    target_id = r$target_id,
    primitives = r$primitives,
    pooled_sigma_late = r$pooled_sigma_late,
    pooled_mu_late = r$pooled_mu_late,
    n_late_events = r$n_late_events,
    pooled_sigma_early = r$pooled_sigma_early,
    csn_R = r$csn_late$R,
    csn_p = r$csn_late$p,
    csn_n_tail = r$csn_late$n_tail,
    nn_mean_class = r$neutral_proxy$mean_class_size_weighted,
    nn_n_classes = r$neutral_proxy$n_classes,
    nn_max_class = r$neutral_proxy$max_class_size,
    stringsAsFactors = FALSE
  )
})
df <- do.call(rbind, rows)

# Per-target paired
targets <- unique(df$target_id)
paired <- data.frame(target_id = targets, stringsAsFactors = FALSE)
paired$sigma_min  <- sapply(targets, function(t) df$pooled_sigma_late[df$target_id==t & df$primitives=="minimal"])
paired$sigma_rich <- sapply(targets, function(t) df$pooled_sigma_late[df$target_id==t & df$primitives=="rich"])
paired$delta_sigma <- paired$sigma_rich - paired$sigma_min
paired$nn_min  <- sapply(targets, function(t) df$nn_mean_class[df$target_id==t & df$primitives=="minimal"])
paired$nn_rich <- sapply(targets, function(t) df$nn_mean_class[df$target_id==t & df$primitives=="rich"])
paired$nn_ratio <- paired$nn_rich / paired$nn_min
paired$csn_R_rich <- sapply(targets, function(t) df$csn_R[df$target_id==t & df$primitives=="rich"])
paired$csn_p_rich <- sapply(targets, function(t) df$csn_p[df$target_id==t & df$primitives=="rich"])
paired$n_tail_rich <- sapply(targets, function(t) df$csn_n_tail[df$target_id==t & df$primitives=="rich"])

# Test 1: paired Wilcoxon one-sided (rich > minimal) on pooled_sigma_late
w <- wilcox.test(paired$sigma_rich, paired$sigma_min, paired = TRUE, alternative = "greater", exact = TRUE)
wilcoxon_V <- as.numeric(w$statistic)
wilcoxon_p <- as.numeric(w$p.value)
wilcoxon_pass <- (wilcoxon_p < 0.05) && all(paired$delta_sigma > 0)  # plan: p<=.05 and consistent direction
# Direction sign: count positive deltas
n_positive <- sum(paired$delta_sigma > 0)
mean_delta <- mean(paired$delta_sigma)

# Test 2: Spearman rho between delta_sigma and nn_ratio + bootstrap CI
sp <- suppressWarnings(cor.test(paired$delta_sigma, paired$nn_ratio, method = "spearman"))
spearman_rho <- as.numeric(sp$estimate)
spearman_p <- as.numeric(sp$p.value)

set.seed(20240101)
B <- 5000
boot_rhos <- replicate(B, {
  idx <- sample(seq_len(nrow(paired)), replace = TRUE)
  if (length(unique(paired$delta_sigma[idx])) < 2 || length(unique(paired$nn_ratio[idx])) < 2) return(NA_real_)
  suppressWarnings(cor(paired$delta_sigma[idx], paired$nn_ratio[idx], method = "spearman"))
})
boot_rhos <- boot_rhos[is.finite(boot_rhos)]
spearman_ci_low  <- as.numeric(quantile(boot_rhos, 0.025))
spearman_ci_high <- as.numeric(quantile(boot_rhos, 0.975))
spearman_ci_crosses_zero <- (spearman_ci_low <= 0) && (spearman_ci_high >= 0)

# Test 3: CSN log-normal vs exponential per (target, rich). Provided in result.json as csn_late.R and csn_late.p
# R > 0 => log-normal preferred; significance at alpha=0.05
csn_per_target <- lapply(targets, function(t) {
  i <- which(df$target_id == t & df$primitives == "rich")
  list(
    target_id = t,
    R = df$csn_R[i],
    p = df$csn_p[i],
    n_tail = df$csn_n_tail[i],
    lognormal_preferred = df$csn_R[i] > 0,
    significant = df$csn_p[i] < 0.05
  )
})
names(csn_per_target) <- targets
n_csn_sig_lognormal <- sum(sapply(csn_per_target, function(x) x$lognormal_preferred && x$significant))
n_csn_sig_any <- sum(sapply(csn_per_target, function(x) x$significant))
csn_pass_all_significant <- all(sapply(csn_per_target, function(x) x$significant))

# Per-target shifts and ratios as named list
delta_sigma_per_target <- as.list(setNames(paired$delta_sigma, paired$target_id))
nn_ratio_per_target <- as.list(setNames(paired$nn_ratio, paired$target_id))
sigma_min_per_target <- as.list(setNames(paired$sigma_min, paired$target_id))
sigma_rich_per_target <- as.list(setNames(paired$sigma_rich, paired$target_id))
nn_min_per_target <- as.list(setNames(paired$nn_min, paired$target_id))
nn_rich_per_target <- as.list(setNames(paired$nn_rich, paired$target_id))

# Aggregate verdict
mechanism_supported <- wilcoxon_pass &&
                       (!spearman_ci_crosses_zero && spearman_rho > 0) &&
                       csn_pass_all_significant

stats <- list(
  n_targets = length(targets),
  targets = as.list(targets),

  # Test 1
  wilcoxon_V = wilcoxon_V,
  wilcoxon_p = wilcoxon_p,
  wilcoxon_alpha = 0.05,
  wilcoxon_alternative = "greater (rich > minimal)",
  wilcoxon_pass = wilcoxon_pass,
  n_positive_deltas = n_positive,
  mean_delta_sigma = mean_delta,

  # Test 2
  spearman_rho = spearman_rho,
  spearman_p = spearman_p,
  spearman_ci_low = spearman_ci_low,
  spearman_ci_high = spearman_ci_high,
  spearman_ci_crosses_zero = spearman_ci_crosses_zero,
  spearman_pass = !spearman_ci_crosses_zero && spearman_rho > 0,
  n_bootstrap = length(boot_rhos),

  # Test 3
  csn_per_target = csn_per_target,
  n_csn_significant = n_csn_sig_any,
  n_csn_significant_lognormal_preferred = n_csn_sig_lognormal,
  csn_pass_all_significant = csn_pass_all_significant,

  # Per-target quantities
  delta_sigma_per_target = delta_sigma_per_target,
  nn_ratio_per_target = nn_ratio_per_target,
  sigma_min_per_target = sigma_min_per_target,
  sigma_rich_per_target = sigma_rich_per_target,
  nn_min_per_target = nn_min_per_target,
  nn_rich_per_target = nn_rich_per_target,

  # Aggregate
  mechanism_supported = mechanism_supported
)

write_json(stats, "stats.json", auto_unbox = TRUE, pretty = TRUE, digits = 6)
cat("Wrote stats.json\n")
print(stats)
