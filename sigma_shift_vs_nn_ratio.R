library(ggplot2); library(wesanderson); library(dplyr); library(jsonlite); library(ggrepel)
stats <- fromJSON("stats.json", simplifyVector = TRUE)
darj <- wes_palette("Darjeeling1", 5)

df <- data.frame(
  target_id  = unlist(stats$target_ids),
  nn_ratio   = unlist(stats$nn_ratio_by_target),
  delta_sigma = unlist(stats$delta_sigma_by_target)
)

ann <- sprintf("Spearman rho = %.2f\n95%% bootstrap CI [%.2f, %.2f]\nWilcoxon V = %d, p = %.3f",
               stats$spearman_rho, stats$spearman_ci_low, stats$spearman_ci_high,
               as.integer(stats$wilcoxon_V), stats$wilcoxon_p)

p <- ggplot(df, aes(x = nn_ratio, y = delta_sigma, color = target_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_point(size = 5) +
  geom_text_repel(aes(label = target_id), size = 3.8, show.legend = FALSE, box.padding = 0.6) +
  scale_color_manual(values = darj[1:5], name = "Target Equation") +
  annotate("text", x = min(df$nn_ratio), y = max(df$delta_sigma),
           label = ann, hjust = 0, vjust = 1, size = 3.6, family = "mono") +
  labs(title = "Plateau-Tail Shift vs. Neutral-Network Inflation",
       x = "E-graph neutral-network size ratio (rich / minimal)",
       y = expression(Delta*sigma[late]*" = "*sigma[rich]*" \u2212 "*sigma[minimal])) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("sigma_shift_vs_nn_ratio.png", p, width = 8, height = 5.5, dpi = 300, units = "in")
