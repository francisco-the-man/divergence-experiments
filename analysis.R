library(jsonlite)
library(ggplot2)

res <- fromJSON("result.json", simplifyVector = FALSE)$results

# Build per-(target, condition) frame
rows <- lapply(res, function(r) {
  data.frame(
    target_id = r$target_id,
    primitives = r$primitives,
    pooled_sigma_late = r$pooled_sigma_late,
    pooled_mu_late = r$pooled_mu_late,
    n_late_events = r$n_late_events,
    pooled_sigma_early = r$pooled_sigma_early,
    pooled_mu_early = r$pooled_mu_early,
    n_early_events = r$n_early_events,
    csn_R = r$csn_late$R,
    csn_p = r$csn_late$p,
    csn_n_tail = r$csn_late$n_tail,
    nn_mean_class_size = r$neutral_proxy$mean_class_size_weighted,
    nn_n_classes = r$neutral_proxy$n_classes,
    stringsAsFactors = FALSE
  )
})
df <- do.call(rbind, rows)
print(df)

# Pivot to wide on primitives
targets <- sort(unique(df$target_id))
sigma_minimal <- sapply(targets, function(t) df$pooled_sigma_late[df$target_id == t & df$primitives == "minimal"])
sigma_rich    <- sapply(targets, function(t) df$pooled_sigma_late[df$target_id == t & df$primitives == "rich"])
nn_minimal    <- sapply(targets, function(t) df$nn_mean_class_size[df$target_id == t & df$primitives == "minimal"])
nn_rich       <- sapply(targets, function(t) df$nn_mean_class_size[df$target_id == t & df$primitives == "rich"])
delta_sigma   <- sigma_rich - sigma_minimal
nn_ratio      <- nn_rich / nn_minimal

cat("\nPer-target table:\n")
per_target_tbl <- data.frame(
  target_id = targets,
  sigma_minimal = sigma_minimal,
  sigma_rich = sigma_rich,
  delta_sigma = delta_sigma,
  nn_minimal = nn_minimal,
  nn_rich = nn_rich,
  nn_ratio = nn_ratio
)
print(per_target_tbl)

# Test 1: Paired Wilcoxon one-sided (rich > minimal)
w <- wilcox.test(sigma_rich, sigma_minimal, paired = TRUE, alternative = "greater", exact = TRUE)
cat("\nWilcoxon:\n"); print(w)

# Test 2: Spearman rho + bootstrap CI
sp <- suppressWarnings(cor.test(delta_sigma, nn_ratio, method = "spearman", exact = FALSE))
cat("\nSpearman:\n"); print(sp)

set.seed(42)
B <- 5000
n <- length(delta_sigma)
boot_rhos <- replicate(B, {
  idx <- sample.int(n, n, replace = TRUE)
  if (length(unique(idx)) < 2) return(NA_real_)
  suppressWarnings(cor(delta_sigma[idx], nn_ratio[idx], method = "spearman"))
})
boot_rhos <- boot_rhos[is.finite(boot_rhos)]
boot_ci <- quantile(boot_rhos, c(0.025, 0.975), na.rm = TRUE)
cat("\nBootstrap Spearman CI:\n"); print(boot_ci)

# Test 3: CSN log-normal vs exponential per (target, rich) — use the precomputed csn_late.R and .p
csn_rich <- lapply(targets, function(t) {
  r <- df[df$target_id == t & df$primitives == "rich", ]
  list(target_id = t, R = r$csn_R, p = r$csn_p, n_tail = r$csn_n_tail,
       significant_at_05 = (r$csn_p < 0.05))
})
names(csn_rich) <- targets

# Build stats
stats <- list(
  targets = targets,
  per_target = list(
    target_id = targets,
    sigma_minimal = unname(sigma_minimal),
    sigma_rich = unname(sigma_rich),
    delta_sigma = unname(delta_sigma),
    nn_minimal = unname(nn_minimal),
    nn_rich = unname(nn_rich),
    nn_ratio = unname(nn_ratio)
  ),
  wilcoxon = list(
    V = unname(w$statistic),
    p_value = w$p.value,
    alternative = "greater",
    n_pairs = length(delta_sigma),
    n_positive = sum(delta_sigma > 0),
    n_negative = sum(delta_sigma < 0),
    median_delta = median(delta_sigma),
    mean_delta = mean(delta_sigma),
    pass = (w$p.value < 0.05)
  ),
  spearman = list(
    rho = unname(sp$estimate),
    p_value = sp$p.value,
    boot_ci_low = unname(boot_ci[1]),
    boot_ci_high = unname(boot_ci[2]),
    boot_B = B,
    ci_crosses_zero = (boot_ci[1] <= 0) && (boot_ci[2] >= 0),
    pass = !((boot_ci[1] <= 0) && (boot_ci[2] >= 0))
  ),
  csn_rich = list(
    per_target = lapply(csn_rich, function(x) list(target_id = x$target_id, R = x$R, p = x$p, n_tail = x$n_tail, significant_at_05 = x$significant_at_05)),
    all_significant = all(sapply(csn_rich, function(x) x$significant_at_05)),
    n_significant = sum(sapply(csn_rich, function(x) x$significant_at_05)),
    n_total = length(csn_rich),
    pass = all(sapply(csn_rich, function(x) x$significant_at_05))
  ),
  overall_pass = NA  # filled below
)

stats$overall_pass <- stats$wilcoxon$pass && stats$spearman$pass && stats$csn_rich$pass

# Pre-compute CCDF data for plotting (pooled late durations per target × condition)
ccdf_data <- list()
for (r in res) {
  key <- paste(r$target_id, r$primitives, sep = "__")
  durs <- unlist(r$pooled_late_durations_sample)
  durs <- durs[durs > 0]
  if (length(durs) == 0) next
  s <- sort(durs)
  ccdf_data[[key]] <- list(
    target_id = r$target_id,
    primitives = r$primitives,
    x = s,
    ccdf = (length(s):1) / length(s)
  )
}
stats$ccdf_summary <- list(
  n_series = length(ccdf_data),
  max_duration = max(sapply(ccdf_data, function(d) max(d$x)))
)

write_json(stats, "stats.json", auto_unbox = TRUE, pretty = TRUE, digits = 8)
cat("\nstats.json written\n")
cat(readLines("stats.json"), sep = "\n")
