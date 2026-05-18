library(ggplot2); library(wesanderson); library(dplyr); library(jsonlite); library(ggrepel)
stats <- fromJSON("stats.json", simplifyVector = TRUE)
darj <- wes_palette("Darjeeling1", 5)

pt <- stats$per_target
df <- data.frame(
  target_id = pt$target_id,
  nn_ratio = pt$nn_ratio,
  delta_sigma = pt$delta_sigma
)

ann <- sprintf("Spearman rho = %.2f, 95%% bootstrap CI [%.2f, %.2f] (n = %d targets, B = %d)",
               stats$spearman$rho, stats$spearman$boot_ci_low, stats$spearman$boot_ci_high,
               nrow(df), stats$spearman$boot_B)

p <- ggplot(df, aes(x = nn_ratio, y = delta_sigma, color = target_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 4) +
  geom_text_repel(aes(label = target_id), size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = darj, name = "Target equation") +
  labs(
    title = "Per-target Delta sigma vs. e-graph neutral-network ratio",
    x = "Neutral-network size ratio (rich / minimal, mean class size)",
    y = "Delta sigma_late = sigma_rich - sigma_minimal",
    subtitle = ann
  ) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9))

ggsave("sigma_shift_vs_nn_ratio.png", p, width = 7.5, height = 5.5, dpi = 300, units = "in")
