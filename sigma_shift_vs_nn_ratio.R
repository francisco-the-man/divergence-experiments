library(ggplot2); library(wesanderson); library(dplyr); library(jsonlite); library(ggrepel)

data  <- fromJSON("result.json", simplifyVector = FALSE)
stats <- fromJSON("stats.json", simplifyVector = FALSE)
darj  <- wes_palette("Darjeeling1", 5)

targets <- unlist(stats$targets)
df <- data.frame(
  target_id   = targets,
  delta_sigma = sapply(targets, function(t) stats$delta_sigma_per_target[[t]]),
  nn_ratio    = sapply(targets, function(t) stats$nn_ratio_per_target[[t]]),
  stringsAsFactors = FALSE
)

rho_label <- sprintf("Spearman rho = %.2f\n95%% CI [%.3f, %.3f]\n(n = %d targets)",
                     stats$spearman_rho, stats$spearman_ci_low, stats$spearman_ci_high,
                     stats$n_targets)

p <- ggplot(df, aes(x = nn_ratio, y = delta_sigma, color = target_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_point(size = 5) +
  geom_text_repel(aes(label = target_id), size = 3.5, show.legend = FALSE, box.padding = 0.8) +
  scale_color_manual(values = darj, name = "Target equation") +
  annotate("label", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.2,
           label = rho_label, size = 3.5, fill = "white", label.size = 0.3) +
  labs(title = "Plateau-tail Shift vs Neutral-Network Inflation",
       x = "Neutral-network size ratio (rich / minimal)",
       y = expression(Delta * sigma[late] * " = " * sigma[rich] - sigma[minimal])) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "right")

ggsave("sigma_shift_vs_nn_ratio.png", p, width = 8, height = 5.5, dpi = 300, units = "in")
