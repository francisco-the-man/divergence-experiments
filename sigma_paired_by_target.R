library(ggplot2); library(wesanderson); library(dplyr); library(jsonlite); library(tidyr)
stats <- fromJSON("stats.json", simplifyVector = TRUE)
darj <- wes_palette("Darjeeling1", 5)

targets <- unlist(stats$target_ids)
df <- data.frame(
  target_id  = rep(targets, 2),
  primitives = rep(c("minimal", "rich"), each = length(targets)),
  sigma_late = c(unlist(stats$sigma_minimal_by_target), unlist(stats$sigma_rich_by_target))
)
df$primitives <- factor(df$primitives, levels = c("minimal", "rich"))

ann <- sprintf("Wilcoxon V = %d, p = %.3f (one-sided rich > minimal)\nDirection: %d/%d targets positive",
               as.integer(stats$wilcoxon_V), stats$wilcoxon_p,
               as.integer(stats$n_targets_positive_shift), as.integer(stats$n_targets_total))

p <- ggplot(df, aes(x = primitives, y = sigma_late, group = target_id, color = target_id)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 4) +
  scale_color_manual(values = darj[1:5], name = "Target Equation") +
  annotate("text", x = 1.5, y = max(df$sigma_late) + 0.05,
           label = ann, size = 3.6, family = "mono") +
  labs(title = "Late-Phase Plateau-Duration Sigma: Paired by Target",
       x = "Primitive Set",
       y = expression("Pooled log-normal "*sigma*" (late-phase plateaus)")) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("sigma_paired_by_target.png", p, width = 7.5, height = 5, dpi = 300, units = "in")
