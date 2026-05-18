library(ggplot2); library(wesanderson); library(dplyr); library(jsonlite); library(ggrepel)

data  <- fromJSON("result.json", simplifyVector = TRUE)$results
stats <- fromJSON("stats.json", simplifyVector = FALSE)
darj  <- wes_palette("Darjeeling1", 5)

tab <- do.call(rbind, lapply(stats$per_target_table, function(r) {
  data.frame(target_id = r$target_id,
             delta_sigma = r$delta_sigma,
             nn_ratio = r$nn_ratio)
}))

ann <- sprintf("Spearman rho = %.2f, 95%% bootstrap CI [%.2f, %.2f]\n(CI crosses zero: %s)",
               stats$spearman_rho, stats$spearman_ci_low, stats$spearman_ci_high,
               ifelse(stats$spearman_ci_crosses_zero, "yes", "no"))

p <- ggplot(tab, aes(x = nn_ratio, y = delta_sigma, color = target_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 4) +
  geom_text_repel(aes(label = target_id), size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = darj, name = "Target equation") +
  labs(title = expression(paste("Per-target ", Delta, sigma, " vs. e-graph neutral-network ratio")),
       x = "Neutral-network size ratio (rich / minimal)",
       y = expression(paste(Delta, sigma[late], " = ", sigma[rich], " \u2212 ", sigma[minimal]))) +
  annotate("text", x = min(tab$nn_ratio), y = max(tab$delta_sigma) + 0.05,
           label = ann, size = 3.2, hjust = 0) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("sigma_shift_vs_nn_ratio.png", p, width = 7.5, height = 5.5, dpi = 300, units = "in")
