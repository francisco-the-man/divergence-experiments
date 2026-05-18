library(ggplot2); library(wesanderson); library(dplyr); library(jsonlite)
stats <- fromJSON("stats.json", simplifyVector = FALSE)
darj  <- wes_palette("Darjeeling1", 5)

cs <- stats$csn_per_target
targets <- names(cs)
df <- data.frame(
  target_id = targets,
  R = sapply(cs, function(x) x$R),
  p = sapply(cs, function(x) x$p),
  significant = sapply(cs, function(x) x$significant),
  stringsAsFactors = FALSE
)
df$label <- sprintf("R = %.2f\np = %.3f", df$R, df$p)

ann <- sprintf("%d of %d targets significant at alpha = 0.05  -  %s",
               stats$csn_n_significant, length(targets),
               ifelse(stats$csn_pass, "PASS", "FAIL (preregistered: all must be significant)"))

p <- ggplot(df, aes(x = target_id, y = R, fill = target_id)) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_text(aes(label = label, vjust = ifelse(R >= 0, -0.2, 1.1)), size = 3.2) +
  scale_fill_manual(values = darj, name = "Target equation", guide = "none") +
  labs(title = "CSN log-likelihood ratio: log-normal vs exponential\n(rich condition, late-phase plateaus)",
       x = "Target equation",
       y = "R (positive = log-normal preferred, negative = exponential preferred)",
       caption = ann) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.caption = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(angle = 15, hjust = 1)) +
  expand_limits(y = c(min(df$R) - 1.2, max(df$R) + 1.2))

ggsave("csn_lognormal_vs_exponential.png", p, width = 7.5, height = 5, dpi = 300, units = "in")
