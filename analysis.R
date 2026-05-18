
library(jsonlite)

dat <- fromJSON("result.json", simplifyVector = FALSE)$results

# ---- Extract per-run late-phase plateau durations ----
# plateau_events: list of [start_gen, end_gen, length, fitness, censored_flag]
# Late phase: post-tau_target where tau_target = generation when best fitness
# drops below some threshold. We use: 50th percentile of best_hist as the
# "tau_target" cutoff (fitness-banded post-tau_target as plan says, simplified
# given the truncated results have no precomputed band).

extract_late_durations <- function(r) {
  bh <- unlist(r$best_hist)
  if (length(bh) < 10) return(list(durs = numeric(0), cens = integer(0)))
  # banded late phase: events whose plateau fitness <= median(best_hist)
  thresh <- median(bh, na.rm = TRUE)
  pe <- r$plateau_events
  durs <- c(); cens <- c()
  for (ev in pe) {
    ev <- unlist(ev)
    if (length(ev) < 5) next
    plateau_fit <- ev[4]
    plateau_len <- ev[3]
    plateau_cens <- ev[5]
    if (!is.na(plateau_fit) && plateau_fit <= thresh && plateau_len >= 1) {
      durs <- c(durs, plateau_len)
      cens <- c(cens, plateau_cens)
    }
  }
  list(durs = durs, cens = cens)
}

# ---- Censored MLE for log-normal sigma ----
# log(T) ~ Normal(mu, sigma^2); right-censored observations contribute
# survival prob 1 - Phi((log t - mu)/sigma).
censored_lognormal_sigma <- function(durs, cens) {
  if (length(durs) < 3) return(NA_real_)
  d <- durs[durs > 0]
  c <- cens[durs > 0]
  if (length(d) < 3) return(NA_real_)
  ld <- log(d)
  # quick uncensored estimate
  mu0 <- mean(ld); s0 <- sd(ld)
  if (!is.finite(s0) || s0 <= 0) return(NA_real_)
  if (sum(c == 0) < 2) return(s0)  # too few uncensored, fall back
  nll <- function(par) {
    mu <- par[1]; sig <- exp(par[2])
    z <- (ld - mu) / sig
    # uncensored: log density of log-normal in log-time = log(dnorm(z)/sig)
    ll_unc <- ifelse(c == 0, dnorm(z, log = TRUE) - log(sig), 0)
    # censored (right-censored): survival
    ll_cen <- ifelse(c == 1, pnorm(z, lower.tail = FALSE, log.p = TRUE), 0)
    -sum(ll_unc + ll_cen)
  }
  fit <- try(optim(c(mu0, log(s0)), nll, method = "Nelder-Mead",
                   control = list(maxit = 500)), silent = TRUE)
  if (inherits(fit, "try-error") || fit$convergence != 0) return(s0)
  exp(fit$par[2])
}

# ---- Pool late-phase durations per (target, condition) ----
groups <- list()
for (r in dat) {
  key <- paste(r$target_id, r$condition, sep = "|")
  if (is.null(groups[[key]])) groups[[key]] <- list(durs = c(), cens = c())
  ld <- extract_late_durations(r)
  groups[[key]]$durs <- c(groups[[key]]$durs, ld$durs)
  groups[[key]]$cens <- c(groups[[key]]$cens, ld$cens)
}

sigma_table <- list()
for (key in names(groups)) {
  g <- groups[[key]]
  parts <- strsplit(key, "\\|")[[1]]
  s <- censored_lognormal_sigma(g$durs, g$cens)
  sigma_table[[key]] <- list(
    target = parts[1], condition = parts[2],
    n_durs = length(g$durs),
    n_cens = sum(g$cens),
    sigma = s
  )
}

# ---- Identify which targets have all three conditions M, S, C ----
targets_all <- unique(sapply(sigma_table, function(x) x$target))
complete_targets <- c()
for (tg in targets_all) {
  conds <- sapply(sigma_table, function(x) if (x$target == tg) x$condition else NA)
  conds <- conds[!is.na(conds)]
  if (all(c("M", "S", "C") %in% conds)) {
    # require non-NA sigma in all three
    sM <- sigma_table[[paste(tg, "M", sep = "|")]]$sigma
    sS <- sigma_table[[paste(tg, "S", sep = "|")]]$sigma
    sC <- sigma_table[[paste(tg, "C", sep = "|")]]$sigma
    if (all(is.finite(c(sM, sS, sC)))) complete_targets <- c(complete_targets, tg)
  }
}
complete_targets <- sort(complete_targets)

# ---- Build paired sigma table ----
sigM <- sapply(complete_targets, function(tg) sigma_table[[paste(tg, "M", sep = "|")]]$sigma)
sigS <- sapply(complete_targets, function(tg) sigma_table[[paste(tg, "S", sep = "|")]]$sigma)
sigC <- sapply(complete_targets, function(tg) sigma_table[[paste(tg, "C", sep = "|")]]$sigma)
names(sigM) <- names(sigS) <- names(sigC) <- complete_targets

n_complete <- length(complete_targets)
n_planned_targets <- 6
n_observed_targets <- length(unique(sapply(dat, function(r) r$target_id)))

# ---- Test 1: paired Wilcoxon sigma_C > sigma_S ----
wilcoxon_V <- NA; wilcoxon_p <- NA
if (n_complete >= 2) {
  w <- try(wilcox.test(sigC, sigS, paired = TRUE, alternative = "greater"), silent = TRUE)
  if (!inherits(w, "try-error")) {
    wilcoxon_V <- unname(w$statistic); wilcoxon_p <- w$p.value
  }
}
delta_CS <- sigC - sigS
n_pos_CS <- sum(delta_CS > 0)

# ---- Test 1b: sensitivity drop T6 ----
# T6 isn't in data; sensitivity not applicable. Report NA.
wilcoxon_sensitivity_p <- NA
wilcoxon_sensitivity_note <- "T6 not present in data; sensitivity check N/A"

# ---- Test 2: Spearman correlations ----
# We lack e-graph measurement results entirely. Use proxy:
#   delta_size_S  = derived from (sigS - sigM) ranks against a structural proxy?
# Actually we have NO static observables. The plan required egraph-measurement
# tasks which produced no results in result.json. Report this fact.

# Best we can do: report bootstrap Spearman on (delta_sigma_S vs target rank
# placeholder) — but that's meaningless. Mark as unavailable.
spearman_available <- FALSE

# But we CAN compute the descriptive Spearmans of delta_sigma vs each other,
# and delta_sigma_C vs delta_sigma_S to see if conditions move together.
rho_C_vs_S <- NA; rho_C_vs_S_p <- NA
if (n_complete >= 3) {
  delta_S <- sigS - sigM
  delta_C <- sigC - sigM
  cr <- try(cor.test(delta_C, delta_S, method = "spearman", exact = FALSE), silent = TRUE)
  if (!inherits(cr, "try-error")) {
    rho_C_vs_S <- unname(cr$estimate); rho_C_vs_S_p <- cr$p.value
  }
}

# ---- Test 3: CSN log-normal vs exponential ----
# Likelihood-ratio R between log-normal and exponential fits to pooled late
# durations (xmin = 1, since plan said xmin search but we don't have powerlaw
# package guaranteed; use full distribution with continuous correction).
# We implement Vuong-style LR.
ln_vs_exp_LR <- function(durs) {
  d <- durs[durs > 0 & is.finite(durs)]
  if (length(d) < 20) return(list(R = NA, p = NA, n = length(d), xmin = NA))
  # xmin search: try several xmins, pick one maximizing tail goodness via
  # Kolmogorov-Smirnov to lognormal; simple: xmin = quantile 0.10.
  xmin <- max(1, quantile(d, 0.10))
  tail <- d[d >= xmin]
  if (length(tail) < 20) return(list(R = NA, p = NA, n = length(tail), xmin = xmin))
  lt <- log(tail)
  # log-normal MLE on tail (truncated)
  mu <- mean(lt); sig <- sd(lt)
  if (!is.finite(sig) || sig <= 0) return(list(R = NA, p = NA, n = length(tail), xmin = xmin))
  # truncated log-normal log-lik
  ll_ln_i <- dnorm(lt, mu, sig, log = TRUE) - log(tail) -
    pnorm(log(xmin), mu, sig, lower.tail = FALSE, log.p = TRUE)
  # exponential MLE on tail shifted by xmin: rate = 1/mean(t - xmin) using
  # truncation; for left-truncated exponential, rate = 1/mean(t - xmin).
  rate <- 1 / mean(tail - xmin)
  ll_ex_i <- log(rate) - rate * (tail - xmin)
  diff_i <- ll_ln_i - ll_ex_i
  R <- sum(diff_i)
  # Vuong normalized statistic
  sd_diff <- sd(diff_i)
  if (!is.finite(sd_diff) || sd_diff <= 0) return(list(R = R, p = NA, n = length(tail), xmin = xmin))
  z <- R / (sqrt(length(diff_i)) * sd_diff)
  p <- 2 * pnorm(-abs(z))
  list(R = R, p = p, n = length(tail), xmin = xmin, z = z)
}

csn_by_cell <- list()
for (key in names(groups)) {
  parts <- strsplit(key, "\\|")[[1]]
  tg <- parts[1]; cd <- parts[2]
  res <- ln_vs_exp_LR(groups[[key]]$durs)
  csn_by_cell[[key]] <- list(
    target = tg, condition = cd,
    R = res$R, p = res$p, n_tail = res$n, xmin = res$xmin
  )
}

# Count C-cells log-normal preferred (R>0, p<0.05) and S-cells exponential preferred
C_lognormal_sig <- 0; S_exponential_sig <- 0
n_C_cells <- 0; n_S_cells <- 0
for (key in names(csn_by_cell)) {
  cc <- csn_by_cell[[key]]
  if (cc$condition == "C" && cc$target %in% complete_targets) {
    n_C_cells <- n_C_cells + 1
    if (!is.na(cc$R) && !is.na(cc$p) && cc$R > 0 && cc$p < 0.05) C_lognormal_sig <- C_lognormal_sig + 1
  }
  if (cc$condition == "S" && cc$target %in% complete_targets) {
    n_S_cells <- n_S_cells + 1
    if (!is.na(cc$R) && !is.na(cc$p) && cc$R < 0 && cc$p < 0.05) S_exponential_sig <- S_exponential_sig + 1
  }
}

# ---- Verdict against null criteria ----
# (i) Wilcoxon: p > 0.05 OR median sign <= 0
crit_i_fires <- is.na(wilcoxon_p) || wilcoxon_p > 0.05 || median(delta_CS) <= 0
# (ii) Spearman dissociation NOT computable -> mark as fired (cannot demonstrate)
crit_ii_fires <- TRUE  # static observables unavailable
# (iii) CSN: need >=4/6 C log-normal AND >=3/6 S exponential
crit_iii_fires <- !(C_lognormal_sig >= 4 && S_exponential_sig >= 3)

mechanism_demonstrated <- !(crit_i_fires || crit_ii_fires || crit_iii_fires)

# ---- Build stats list ----
stats <- list(
  n_results = length(dat),
  n_planned_tasks = 540,
  n_planned_targets = n_planned_targets,
  n_observed_targets = n_observed_targets,
  observed_target_ids = unique(sapply(dat, function(r) r$target_id)),
  complete_targets = complete_targets,
  n_complete_targets = n_complete,

  sigma_M_by_target = as.list(round(sigM, 4)),
  sigma_S_by_target = as.list(round(sigS, 4)),
  sigma_C_by_target = as.list(round(sigC, 4)),

  delta_sigma_CS_by_target = as.list(round(sigC - sigS, 4)),
  delta_sigma_SM_by_target = as.list(round(sigS - sigM, 4)),
  delta_sigma_CM_by_target = as.list(round(sigC - sigM, 4)),

  wilcoxon_V = wilcoxon_V,
  wilcoxon_p = round(wilcoxon_p, 4),
  wilcoxon_n_pairs = n_complete,
  wilcoxon_n_positive_CS = n_pos_CS,
  wilcoxon_median_delta_CS = round(median(delta_CS), 4),

  wilcoxon_sensitivity_p = wilcoxon_sensitivity_p,
  wilcoxon_sensitivity_note = wilcoxon_sensitivity_note,

  spearman_static_available = spearman_available,
  spearman_static_note = "E-graph measurement tasks (kind=egraph_measure) absent from result.json; static observables (Δsize, Δconnectivity) cannot be computed. Test (ii) cannot be evaluated.",
  rho_deltaC_vs_deltaS = round(rho_C_vs_S, 4),
  rho_deltaC_vs_deltaS_p = round(rho_C_vs_S_p, 4),

  csn_by_cell = csn_by_cell,
  csn_n_C_cells = n_C_cells,
  csn_n_S_cells = n_S_cells,
  csn_C_lognormal_significant = C_lognormal_sig,
  csn_S_exponential_significant = S_exponential_sig,
  csn_C_lognormal_threshold = 4,
  csn_S_exponential_threshold = 3,

  crit_i_wilcoxon_fires = crit_i_fires,
  crit_ii_spearman_fires = crit_ii_fires,
  crit_iii_csn_fires = crit_iii_fires,
  mechanism_demonstrated = mechanism_demonstrated,

  data_completeness_note = sprintf(
    "Plan called for 540 tasks across 6 targets x 3 conditions x 30 reps; result.json contains %d gp_run results covering %d target(s): %s. Targets T4_deep_mul, T5_transcend, T6_inexpressible absent; T3_rational partial (20 reps).",
    length(dat), n_observed_targets, paste(unique(sapply(dat, function(r) r$target_id)), collapse = ", ")
  )
)

write_json(stats, "stats.json", auto_unbox = TRUE, pretty = TRUE, na = "null")
cat("Wrote stats.json\n")
cat("Complete targets:", paste(complete_targets, collapse = ", "), "\n")
cat("Wilcoxon V =", wilcoxon_V, "p =", wilcoxon_p, "\n")
cat("CSN C-lognormal sig:", C_lognormal_sig, "/", n_C_cells, "\n")
cat("CSN S-exponential sig:", S_exponential_sig, "/", n_S_cells, "\n")
