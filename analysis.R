library(jsonlite)
library(dplyr)

res <- fromJSON("result.json", simplifyVector = FALSE)$results

# Build per-(target, condition) table
rows <- lapply(res, function(r) {
  data.frame(
    target_id = r$target_id,
    primitives = r$primitives,
    sigma_late = r$pooled_sigma_late,
    mu_late = r$pooled_mu_late,
    n_late = r$n_late_events,
    sigma_early = r$pooled_sigma_early,
    nn_mean_class = r$neutral_proxy$mean_class_size_weighted,
    csn_R = r$csn_late$R,
    csn_p = r$csn_late$p,
    csn_n_tail = r$csn_late$n_tail,
    stringsAsFactors = FALSE
  )
})
df <- do.call(rbind, rows)

# Pivot to wide form per target
targets <- unique(df$target_id)
wide <- data.frame(
  target_id = targets,
  sigma_late_minimal = sapply(targets, function(t) df$sigma_late[df$target_id == t & df$primitives == "minimal"]),
  sigma_late_rich    = sapply(targets, function(t) df$sigma_late[df$target_id == t & df$primitives == "rich"]),
  nn_minimal         = sapply(targets, function(t) df$nn_mean_class[df$target_id == t & df$primitives == "minimal"]),
  nn_rich            = sapply(targets, function(t) df$nn_mean_class[df$target_id == t & df$primitives == "rich"]),
  stringsAsFactors = FALSE
)
wide$delta_sigma <- wide$sigma_late_rich - wide$sigma_late_minimal
wide$nn_ratio <- wide$nn_rich / wide$nn_minimal

# Test 1: paired Wilcoxon, one-sided greater (rich > minimal)
wt <- wilcox.test(wide$sigma_late_rich, wide$sigma_late_minimal,
                  paired = TRUE, alternative = "greater", exact = FALSE)

# Test 2: Spearman correlation + bootstrap CI
ct <- cor.test(wide$delta_sigma, wide$nn_ratio, method = "spearman", exact = FALSE)

set.seed(42)
B <- 5000
n <- nrow(wide)
boot_rho <- numeric(B)
for (i in seq_len(B)) {
  idx <- sample.int(n, n, replace = TRUE)
  if (length(unique(idx)) < 3) { boot_rho[i] <- NA; next }
  x <- wide$delta_sigma[idx]; y <- wide$nn_ratio[idx]
  if (sd(x) == 0 || sd(y) == 0) { boot_rho[i] <- NA; next }
  boot_rho[i] <- suppressWarnings(cor(x, y, method = "spearman"))
}
boot_rho <- boot_rho[is.finite(boot_rho)]
ci <- quantile(boot_rho, c(0.025, 0.975))

# Test 3: csn_late R and p per target (rich)
csn_rich <- df[df$primitives == "rich", c("target_id", "csn_R", "csn_p", "csn_n_tail")]
csn_rich_list <- setNames(lapply(seq_len(nrow(csn_rich)), function(i) {
  list(R = csn_rich$csn_R[i], p = csn_rich$csn_p[i],
       n_tail = csn_rich$csn_n_tail[i],
       significant = csn_rich$csn_p[i] < 0.05)
}), csn_rich$target_id)

# Pass/fail evaluations per null criteria
wilcox_pass <- (wt$p.value < 0.05) && (median(wide$delta_sigma) > 0 || mean(wide$delta_sigma) > 0) &&
               all(sign(wide$delta_sigma) >= 0) == FALSE  # direction check below
direction_all_positive <- all(wide$delta_sigma > 0)
direction_mixed <- any(wide$delta_sigma > 0) && any(wide$delta_sigma < 0)
wilcox_pass <- (wt$p.value < 0.05) && !direction_mixed

spearman_ci_crosses_zero <- (ci[1] <= 0) && (ci[2] >= 0)
spearman_pass <- !spearman_ci_crosses_zero && (ct$estimate > 0)

csn_all_significant <- all(sapply(csn_rich_list, function(x) x$significant))
csn_pass <- csn_all_significant

overall_pass <- wilcox_pass && spearman_pass && csn_pass

# Per-target table for narrative
per_target <- setNames(lapply(seq_len(nrow(wide)), function(i) {
  list(
    sigma_late_minimal = wide$sigma_late_minimal[i],
    sigma_late_rich = wide$sigma_late_rich[i],
    delta_sigma = wide$delta_sigma[i],
    nn_minimal = wide$nn_minimal[i],
    nn_rich = wide$nn_rich[i],
    nn_ratio = wide$nn_ratio[i]
  )
}), wide$target_id)

stats <- list(
  n_targets = nrow(wide),
  wilcoxon_V = unname(wt$statistic),
  wilcoxon_p = wt$p.value,
  wilcoxon_alternative = "rich > minimal (one-sided)",
  wilcoxon_pass = wilcox_pass,
  direction_all_positive = direction_all_positive,
  direction_mixed = direction_mixed,
  n_positive_deltas = sum(wide$delta_sigma > 0),
  n_negative_deltas = sum(wide$delta_sigma < 0),
  median_delta_sigma = median(wide$delta_sigma),
  mean_delta_sigma = mean(wide$delta_sigma),

  spearman_rho = unname(ct$estimate),
  spearman_p = ct$p.value,
  spearman_ci_low = unname(ci[1]),
  spearman_ci_high = unname(ci[2]),
  spearman_ci_crosses_zero = spearman_ci_crosses_zero,
  spearman_pass = spearman_pass,
  n_boot = length(boot_rho),

  csn_per_target = csn_rich_list,
  csn_all_significant = csn_all_significant,
  csn_n_significant = sum(sapply(csn_rich_list, function(x) x$significant)),
  csn_pass = csn_pass,

  per_target = per_target,
  overall_pass = overall_pass,

  nn_ratio_min = min(wide$nn_ratio),
  nn_ratio_max = max(wide$nn_ratio),
  nn_ratio_median = median(wide$nn_ratio)
)

write_json(stats, "stats.json", auto_unbox = TRUE, pretty = TRUE, digits = 6)
cat("done\n")
