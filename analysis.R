
# analysis.R — tk_3b8ce4aa
suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
})

`%||%` <- function(a,b) if (is.null(a)) b else a

raw <- fromJSON("result.json", simplifyVector = FALSE)$results
cat("n_results:", length(raw), "\n")

kinds <- table(sapply(raw, function(r) r$kind %||% "NA"))
print(kinds)

gp     <- Filter(function(r) (r$kind %||% "") == "gp_run", raw)
egraph <- Filter(function(r) (r$kind %||% "") != "gp_run", raw)
cat("gp_run:", length(gp), "  egraph/other:", length(egraph), "\n")

target_ids <- sort(unique(sapply(gp, function(r) r$target_id)))
conds      <- sort(unique(sapply(gp, function(r) r$condition)))
cat("targets:", paste(target_ids, collapse=", "), "\n")
cat("conditions:", paste(conds, collapse=", "), "\n")

# --------------------------------------------------------------
# Plateau-duration extraction
# --------------------------------------------------------------
extract_plateaus <- function(r) {
  pe <- r$plateau_events
  if (is.null(pe) || length(pe) == 0) return(NULL)
  m <- do.call(rbind, lapply(pe, function(x) as.numeric(unlist(x))))
  if (ncol(m) < 5) return(NULL)
  colnames(m) <- c("start","end","dur","fit","obs")
  data.frame(
    target_id = r$target_id,
    condition = r$condition,
    seed      = r$seed %||% NA_integer_,
    start     = m[,"start"],
    end       = m[,"end"],
    duration  = m[,"dur"],
    fitness   = m[,"fit"],
    observed  = as.integer(m[,"obs"]),
    stringsAsFactors = FALSE
  )
}

pl <- do.call(rbind, lapply(gp, extract_plateaus))
cat("plateau events total:", nrow(pl), "\n")

pl <- pl %>%
  group_by(target_id, condition, seed) %>%
  mutate(fit_med_run = median(fitness),
         late_phase  = fitness <= fit_med_run) %>%
  ungroup() %>%
  as.data.frame()

late <- pl %>% filter(late_phase, duration > 0)
cat("late-phase plateau events:", nrow(late), "\n")

# --------------------------------------------------------------
# Censored log-normal MLE
# --------------------------------------------------------------
fit_lognormal_censored <- function(dur, obs) {
  ld <- log(dur)
  if (length(ld) < 5 || sd(ld) == 0) return(c(mu=NA_real_, sigma=NA_real_, n=length(ld)))
  nll <- function(par) {
    mu <- par[1]; s <- exp(par[2])
    if (!is.finite(s) || s <= 0) return(1e10)
    ll_obs <- sum(dnorm(ld[obs==1], mu, s, log=TRUE))
    surv   <- pnorm(ld[obs==0], mu, s, lower.tail=FALSE, log.p=TRUE)
    ll_cen <- sum(surv)
    -(ll_obs + ll_cen)
  }
  init <- c(mean(ld), log(sd(ld)))
  fit  <- tryCatch(optim(init, nll, method="Nelder-Mead",
                         control=list(reltol=1e-8, maxit=2000)),
                   error=function(e) NULL)
  if (is.null(fit) || fit$convergence != 0) return(c(mu=NA_real_, sigma=NA_real_, n=length(ld)))
  c(mu = fit$par[1], sigma = exp(fit$par[2]), n = length(ld))
}

sigma_tab <- late %>%
  group_by(target_id, condition) %>%
  group_modify(~ {
    f <- fit_lognormal_censored(.x$duration, .x$observed)
    data.frame(mu=f["mu"], sigma=f["sigma"], n=f["n"])
  }) %>%
  ungroup() %>%
  as.data.frame()

print(sigma_tab)

sigma_wide <- reshape(sigma_tab[,c("target_id","condition","sigma")],
                      idvar="target_id", timevar="condition",
                      direction="wide")
names(sigma_wide) <- sub("sigma\\.", "sigma_", names(sigma_wide))
print(sigma_wide)

have_C <- "sigma_C" %in% names(sigma_wide)
have_S <- "sigma_S" %in% names(sigma_wide)
have_M <- "sigma_M" %in% names(sigma_wide)

if (have_C && have_S) {
  paired <- sigma_wide[complete.cases(sigma_wide[,c("sigma_C","sigma_S")]),]
  wt <- tryCatch(
    wilcox.test(paired$sigma_C, paired$sigma_S, paired=TRUE,
                alternative="greater", exact=FALSE),
    error=function(e) list(statistic=c(V=NA), p.value=NA))
  wilcox_V <- as.numeric(wt$statistic)
  wilcox_p <- as.numeric(wt$p.value)
  n_pos    <- sum(paired$sigma_C > paired$sigma_S, na.rm=TRUE)
  n_pairs  <- nrow(paired)
  median_delta_CS <- median(paired$sigma_C - paired$sigma_S, na.rm=TRUE)
} else {
  wilcox_V <- NA; wilcox_p <- NA; n_pos <- NA; n_pairs <- 0; median_delta_CS <- NA
  paired <- data.frame()
}
cat("Wilcoxon C>S: V=", wilcox_V, "p=", wilcox_p, "\n")

t6_mask <- grepl("^T6", paired$target_id)
if (any(t6_mask)) {
  ps <- paired[!t6_mask,]
  wt2 <- tryCatch(
    wilcox.test(ps$sigma_C, ps$sigma_S, paired=TRUE,
                alternative="greater", exact=FALSE),
    error=function(e) list(statistic=c(V=NA), p.value=NA))
  wilcox_V_noT6 <- as.numeric(wt2$statistic)
  wilcox_p_noT6 <- as.numeric(wt2$p.value)
  n_pairs_noT6  <- nrow(ps)
  sensitivity_T6_present <- TRUE
} else {
  wilcox_V_noT6 <- NA; wilcox_p_noT6 <- NA; n_pairs_noT6 <- NA
  sensitivity_T6_present <- FALSE
}

# --------------------------------------------------------------
# Egraph banded
# --------------------------------------------------------------
to_num_vec <- function(x) {
  if (is.null(x)) return(numeric(0))
  v <- suppressWarnings(as.numeric(unlist(x)))
  v[is.finite(v)]
}

egraph_summary <- function(r) {
  size_v   <- to_num_vec(r$size)
  portal_v <- to_num_vec(r$portal)
  edges_v  <- to_num_vec(r$edges)
  conn_v <- if (length(portal_v) > 0) portal_v else edges_v
  data.frame(
    target_id   = r$target_id   %||% NA,
    condition   = r$condition   %||% NA,
    mean_size   = if (length(size_v))  mean(size_v)  else NA_real_,
    mean_portal = if (length(conn_v))  mean(conn_v)  else NA_real_,
    stringsAsFactors = FALSE
  )
}

if (length(egraph) > 0) {
  eg_tab <- do.call(rbind, lapply(egraph, egraph_summary))
} else {
  eg_tab <- data.frame(target_id=character(), condition=character(),
                       mean_size=numeric(), mean_portal=numeric())
}
print(eg_tab)

# Spearmans
delta_tab <- NULL
if (nrow(eg_tab) > 0 && all(c("sigma_M","sigma_S","sigma_C") %in% names(sigma_wide))) {
  eg_wide <- reshape(eg_tab, idvar="target_id", timevar="condition",
                     direction="wide")
  delta_tab <- merge(sigma_wide, eg_wide, by="target_id", all=FALSE)
  delta_tab$delta_sigma_S    <- delta_tab$sigma_S - delta_tab$sigma_M
  delta_tab$delta_sigma_C    <- delta_tab$sigma_C - delta_tab$sigma_M
  if ("mean_size.S" %in% names(delta_tab) && "mean_size.M" %in% names(delta_tab))
    delta_tab$delta_size_S   <- delta_tab$mean_size.S - delta_tab$mean_size.M
  if ("mean_portal.C" %in% names(delta_tab) && "mean_portal.M" %in% names(delta_tab))
    delta_tab$delta_conn_C   <- delta_tab$mean_portal.C - delta_tab$mean_portal.M
}

boot_spearman_ci <- function(x, y, B=10000, seed=42) {
  if (length(x) < 3 || any(!is.finite(x)) || any(!is.finite(y)))
    return(list(rho=NA, lo=NA, hi=NA))
  set.seed(seed)
  rho_hat <- suppressWarnings(cor(x, y, method="spearman"))
  n <- length(x)
  reps <- replicate(B, {
    i <- sample.int(n, n, replace=TRUE)
    if (length(unique(x[i]))<2 || length(unique(y[i]))<2) return(NA_real_)
    suppressWarnings(cor(x[i], y[i], method="spearman"))
  })
  reps <- reps[is.finite(reps)]
  ci <- quantile(reps, c(0.025, 0.975), na.rm=TRUE)
  list(rho = unname(rho_hat),
       lo  = unname(ci[1]),
       hi  = unname(ci[2]))
}

rho_size_res <- list(rho=NA, lo=NA, hi=NA, p=NA)
rho_conn_res <- list(rho=NA, lo=NA, hi=NA, p=NA)
rho_contrast <- list(diff=NA, lo=NA, hi=NA)

if (!is.null(delta_tab) &&
    all(c("delta_sigma_S","delta_size_S","delta_sigma_C","delta_conn_C") %in% names(delta_tab))) {

  good <- complete.cases(delta_tab[,c("delta_sigma_S","delta_size_S",
                                      "delta_sigma_C","delta_conn_C")])
  d <- delta_tab[good,]

  if (nrow(d) >= 3) {
    rs <- boot_spearman_ci(d$delta_sigma_S, d$delta_size_S, B=10000, seed=1)
    rc <- boot_spearman_ci(d$delta_sigma_C, d$delta_conn_C, B=10000, seed=2)
    ps <- suppressWarnings(cor.test(d$delta_sigma_S, d$delta_size_S, method="spearman", exact=FALSE))
    pc <- suppressWarnings(cor.test(d$delta_sigma_C, d$delta_conn_C, method="spearman", exact=FALSE))

    rho_size_res <- list(rho=rs$rho, lo=rs$lo, hi=rs$hi, p=as.numeric(ps$p.value))
    rho_conn_res <- list(rho=rc$rho, lo=rc$lo, hi=rc$hi, p=as.numeric(pc$p.value))

    set.seed(7)
    B <- 10000
    n <- nrow(d)
    diffs <- replicate(B, {
      i <- sample.int(n, n, replace=TRUE)
      a <- suppressWarnings(cor(d$delta_sigma_C[i], d$delta_conn_C[i], method="spearman"))
      b <- suppressWarnings(cor(d$delta_sigma_S[i], d$delta_size_S[i], method="spearman"))
      a - b
    })
    diffs <- diffs[is.finite(diffs)]
    rho_contrast <- list(diff = rc$rho - rs$rho,
                         lo   = unname(quantile(diffs, 0.025, na.rm=TRUE)),
                         hi   = unname(quantile(diffs, 0.975, na.rm=TRUE)))
  }
}

cat("rho_size:", unlist(rho_size_res), "\n")
cat("rho_conn:", unlist(rho_conn_res), "\n")
cat("rho_contrast:", unlist(rho_contrast), "\n")

# --------------------------------------------------------------
# CSN log-normal vs exponential with xmin search
# --------------------------------------------------------------
fit_tail_compare <- function(x, min_n = 50) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < min_n) return(list(xmin=NA, R=NA, p=NA, n=length(x)))
  cand <- sort(unique(x))
  cand <- cand[cand < quantile(x, 0.95)]
  if (length(cand) < 2) return(list(xmin=NA, R=NA, p=NA, n=length(x)))

  ks_at <- function(xm) {
    tail <- x[x >= xm]
    if (length(tail) < min_n) return(Inf)
    lt <- log(tail); mu <- mean(lt); s <- sd(lt)
    if (s <= 0) return(Inf)
    st <- sort(tail)
    ecdf_v <- seq_along(st)/length(st)
    th     <- plnorm(st, meanlog=mu, sdlog=s)
    # rescale theoretical to tail: F_tail(x) = (F(x)-F(xmin))/(1-F(xmin))
    F_xm <- plnorm(xm, meanlog=mu, sdlog=s)
    th_tail <- (th - F_xm) / (1 - F_xm)
    max(abs(ecdf_v - th_tail))
  }
  ks_vals <- sapply(cand, ks_at)
  if (all(!is.finite(ks_vals))) return(list(xmin=NA, R=NA, p=NA, n=length(x)))
  xmin <- cand[which.min(ks_vals)]
  tail <- x[x >= xmin]
  if (length(tail) < min_n) return(list(xmin=xmin, R=NA, p=NA, n=length(tail)))

  lt <- log(tail); mu <- mean(lt); s <- sd(lt)
  # log f_ln(x) conditioned on x>=xmin
  ln_pdf <- dlnorm(tail, meanlog=mu, sdlog=s, log=TRUE)
  ln_norm <- plnorm(xmin, meanlog=mu, sdlog=s, lower.tail=FALSE, log.p=TRUE)
  log_f_ln <- ln_pdf - ln_norm

  # exponential shifted: y = x - xmin
  y <- tail - xmin
  if (mean(y) <= 0) return(list(xmin=xmin, R=NA, p=NA, n=length(tail)))
  rate <- 1/mean(y)
  log_f_exp <- dexp(y, rate=rate, log=TRUE)

  diff_ <- log_f_ln - log_f_exp
  R <- sum(diff_)
  sd_d <- sd(diff_)
  if (!is.finite(sd_d) || sd_d == 0) return(list(xmin=xmin, R=R, p=NA, n=length(tail)))
  z <- R / (sqrt(length(diff_)) * sd_d)
  p <- 2 * pnorm(-abs(z))
  list(xmin=xmin, R=R, p=p, n=length(tail), z=z)
}

csn_rows <- list()
for (tg in unique(late$target_id)) {
  for (cd in unique(late$condition)) {
    sub <- late[late$target_id==tg & late$condition==cd & late$observed==1, ]
    fit <- fit_tail_compare(sub$duration, min_n = 50)
    csn_rows[[paste(tg,cd,sep="|")]] <- data.frame(
      target_id=tg, condition=cd,
      xmin=fit$xmin %||% NA, R=fit$R %||% NA, p=fit$p %||% NA, n_tail=fit$n %||% NA,
      stringsAsFactors=FALSE
    )
  }
}
csn_tab <- do.call(rbind, csn_rows)
rownames(csn_tab) <- NULL
print(csn_tab)

csn_C <- csn_tab[csn_tab$condition=="C", ]
csn_S <- csn_tab[csn_tab$condition=="S", ]

n_C_lognorm_sig <- sum(csn_C$R > 0 & csn_C$p < 0.05, na.rm=TRUE)
n_S_exp_sig     <- sum(csn_S$R < 0 & csn_S$p < 0.05, na.rm=TRUE)
n_targets_csn   <- nrow(csn_C)

csn_pass <- (n_C_lognorm_sig >= 4) && (n_S_exp_sig >= 3)

wilcoxon_pass <- !is.na(wilcox_p) && (wilcox_p <= 0.05) &&
                 !is.na(median_delta_CS) && (median_delta_CS > 0)

diss_pass <- !is.na(rho_conn_res$lo) && !is.na(rho_size_res$rho) &&
             (rho_conn_res$lo > 0) &&
             (rho_conn_res$rho > rho_size_res$rho)

overall_pass <- wilcoxon_pass && diss_pass && csn_pass

sigma_by_target_cond <- list()
for (i in seq_len(nrow(sigma_tab))) {
  key <- paste(sigma_tab$target_id[i], sigma_tab$condition[i], sep="_")
  sigma_by_target_cond[[key]] <- sigma_tab$sigma[i]
}

ccdf_rows <- list()
for (tg in unique(late$target_id)) {
  for (cd in unique(late$condition)) {
    sub <- late[late$target_id==tg & late$condition==cd & late$observed==1, ]
    d <- sort(sub$duration)
    if (length(d) < 2) next
    ccdf_rows[[paste(tg,cd,sep="|")]] <- data.frame(
      target_id=tg, condition=cd,
      duration = d,
      ccdf = 1 - (seq_along(d)-1)/length(d)
    )
  }
}
ccdf_df <- do.call(rbind, ccdf_rows)

stats <- list(
  n_results            = length(raw),
  n_gp_runs            = length(gp),
  n_egraph             = length(egraph),
  target_ids           = target_ids,
  conditions           = conds,
  n_targets            = length(target_ids),

  sigma_table          = sigma_tab,
  sigma_wide           = sigma_wide,
  sigma_by_target_cond = sigma_by_target_cond,

  wilcoxon_V           = wilcox_V,
  wilcoxon_p           = wilcox_p,
  n_pos_CS             = n_pos,
  n_pairs_CS           = n_pairs,
  median_delta_CS      = median_delta_CS,
  wilcoxon_pass        = wilcoxon_pass,

  sensitivity_T6_present = sensitivity_T6_present,
  wilcoxon_V_noT6      = wilcox_V_noT6,
  wilcoxon_p_noT6      = wilcox_p_noT6,
  n_pairs_noT6         = n_pairs_noT6,

  rho_size             = rho_size_res$rho,
  rho_size_lo          = rho_size_res$lo,
  rho_size_hi          = rho_size_res$hi,
  rho_size_p           = rho_size_res$p,
  rho_conn             = rho_conn_res$rho,
  rho_conn_lo          = rho_conn_res$lo,
  rho_conn_hi          = rho_conn_res$hi,
  rho_conn_p           = rho_conn_res$p,
  rho_contrast_diff    = rho_contrast$diff,
  rho_contrast_lo      = rho_contrast$lo,
  rho_contrast_hi      = rho_contrast$hi,
  dissociation_pass    = diss_pass,
  delta_table          = if (!is.null(delta_tab)) delta_tab else data.frame(),

  csn_table            = csn_tab,
  n_C_lognorm_sig      = n_C_lognorm_sig,
  n_S_exp_sig          = n_S_exp_sig,
  n_targets_csn        = n_targets_csn,
  csn_pass             = csn_pass,

  egraph_summary       = eg_tab,
  ccdf_df              = ccdf_df,

  overall_pass         = overall_pass
)

write_json(stats, "stats.json", auto_unbox = TRUE, pretty = TRUE, na = "null", digits = 6)
cat("\nWROTE stats.json\n")
cat("overall_pass:", overall_pass, "\n")
cat("wilcoxon_pass:", wilcoxon_pass, "\n")
cat("dissociation_pass:", diss_pass, "\n")
cat("csn_pass:", csn_pass, "\n")
