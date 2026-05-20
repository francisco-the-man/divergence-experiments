library(ggplot2); library(wesanderson); library(jsonlite); library(dplyr)
data  <- fromJSON("result.json", simplifyVector = FALSE)
stats <- fromJSON("stats.json", simplifyVector = TRUE)
darj  <- wes_palette("Darjeeling1", 5)

rows <- list()
for (r in data$results) {
  v <- r$features$stasis_jump_rankcorr
  if (is.null(v)) v <- NA
  v <- suppressWarnings(as.numeric(v))
  if (length(v) == 0) v <- NA
  if (is.nan(v) || is.infinite(v)) v <- NA
  rows[[length(rows)+1]] <- data.frame(cell = r$cell, family = r$family, rho = v, n_epochs = r$n_epochs)
}
df <- do.call(rbind, rows)

df$cell <- factor(df$cell,
  levels = c("NK_K2","NK_K4","NK_K8","NK_K16","RMF_theta0.1","RMF_theta0.3","RMF_theta1.0","RMF_theta3.0"))

ann <- sprintf("Mean epochs per trajectory: %.2f (median %.0f)\n%.1f%% of trajectories have <4 epochs (rho undefined or noise-limited)\nstasis_jump_rankcorr permutation importance: %.3f",
  stats$n_epochs_mean, stats$n_epochs_median, 100*stats$n_epochs_pct_lt_4,
  stats$permutation_importance$importance[stats$permutation_importance$feature == "stasis_jump_rankcorr"])

p <- ggplot(df, aes(x = cell, y = rho, color = family)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.6) +
  geom_boxplot(outlier.size = 0.8, linewidth = 0.7, na.rm = TRUE) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.5, na.rm = TRUE) +
  scale_color_manual(values = darj[c(1,2)], name = "Family") +
  annotate("text", x = 0.6, y = -0.95, label = ann, size = 3.2, family = "mono", hjust = 0) +
  labs(title = "Per-Trajectory Stasis–Jump Rank Correlation by Cell",
       x = "Landscape Cell",
       y = expression("Per-trajectory Spearman "*rho*" (log-stasis vs log-jump)")) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("jump_vs_stasis.png", p, width = 9, height = 5.5, dpi = 300, units = "in")
