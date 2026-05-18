library(ggplot2); library(wesanderson); library(dplyr); library(jsonlite); library(tidyr)

data  <- fromJSON("result.json", simplifyVector = TRUE)$results
stats <- fromJSON("stats.json", simplifyVector = FALSE)
darj  <- wes_palette("Darjeeling1", 5)

df <- data.frame(
  target_id = data$target_id,
  primitives = data$primitives,
  sigma_late = data$pooled_sigma_late
)
df$primitives <- factor(df$primitives, levels = c("minimal", "rich"))

ann <- sprintf("Paired Wilcoxon (one-sided rich > minimal):\nV = %.0f, p = %.4f\n%d / %d targets in predicted direction",
               stats$wilcoxon_V, stats$wilcoxon_p,
               stats$n_targets_rich_gt_minimal, stats$n_targets_total)

p <- ggplot(df, aes(x = primitives, y = sigma_late, group = target_id, color = target_id)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = darj, name = "Target equation") +
  labs(title = "Late-phase plateau heavy-tailedness: minimal vs rich primitives",
       x = "Primitive set",
       y = expression(paste("Log-normal ", sigma, " of late-phase plateau durations"))) +
  annotate("text", x = 1.5, y = max(df$sigma_late) + 0.05,
           label = ann, size = 3.2, hjust = 0.5) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("sigma_paired_by_target.png", p, width = 7, height = 5, dpi = 300, units = "in")
