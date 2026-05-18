library(ggplot2); library(wesanderson); library(dplyr); library(jsonlite); library(tidyr)
stats <- fromJSON("stats.json", simplifyVector = TRUE)
darj <- wes_palette("Darjeeling1", 5)

pt <- stats$per_target
df <- data.frame(
  target_id = rep(pt$target_id, 2),
  primitives = rep(c("minimal", "rich"), each = length(pt$target_id)),
  sigma = c(pt$sigma_minimal, pt$sigma_rich)
)
df$primitives <- factor(df$primitives, levels = c("minimal", "rich"))

ann <- sprintf("Paired Wilcoxon (rich > minimal): V = %d, p = %.3f\n%d / %d targets shift in predicted direction; mean Delta sigma = %.3f",
               as.integer(stats$wilcoxon$V), stats$wilcoxon$p_value,
               stats$wilcoxon$n_positive, stats$wilcoxon$n_pairs,
               stats$wilcoxon$mean_delta)

p <- ggplot(df, aes(x = primitives, y = sigma, group = target_id, color = target_id)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = darj, name = "Target equation") +
  labs(
    title = "Late-phase plateau sigma: minimal vs. rich primitives (paired by target)",
    x = "Primitive set",
    y = "Pooled log-normal sigma on late-phase plateau durations"
  ) +
  annotate("text", x = 1.5, y = min(df$sigma) - 0.08, label = ann, size = 3.2, hjust = 0.5) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("sigma_paired_by_target.png", p, width = 8, height = 5.5, dpi = 300, units = "in")
