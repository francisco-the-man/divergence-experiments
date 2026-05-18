library(ggplot2); library(wesanderson); library(dplyr); library(jsonlite); library(ggrepel)
stats <- fromJSON("stats.json", simplifyVector = FALSE)
darj  <- wes_palette("Darjeeling1", 5)

pt <- stats$per_target
targets <- names(pt)
df <- data.frame(
  target_id = targets,
  delta_sigma = sapply(pt, function(x) x$delta_sigma),
  nn_ratio    = sapply(pt, function(x) x$nn_ratio),
  stringsAsFactors = FALSE
)

ann <- sprintf("Spearman rho = %.2f, 95%% bootstrap CI [%.2f, %.2f]  -  %s",
               stats$spearman_rho, stats$spearman_ci_low, stats$spearman_ci_high,
               ifelse(stats$spearman_pass, "PASS (CI excludes 0)", "FAIL"))

p <- ggplot(df, aes(x = nn_ratio, y = delta_sigma, color = target_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 4) +
  geom_text_repel(aes(label = target_id), size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = darj, name = "Target equation") +
  labs(title = "Per-target shift in plateau sigma\nvs e-graph neutral-network ratio",
       x = "Neutral-network size ratio (rich / minimal)",
       y = "Delta sigma_late (rich - minimal)",
       caption = ann) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.caption = element_text(hjust = 0.5, size = 10))

ggsave("sigma_shift_vs_nn_ratio.png", p, width = 7, height = 5, dpi = 300, units = "in")
