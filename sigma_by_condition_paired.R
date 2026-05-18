library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr)
stats <- fromJSON("stats.json", simplifyVector = TRUE)
darj  <- wes_palette("Darjeeling1", 5)

targets <- stats$complete_targets
sigM <- unlist(stats$sigma_M_by_target)[targets]
sigS <- unlist(stats$sigma_S_by_target)[targets]
sigC <- unlist(stats$sigma_C_by_target)[targets]

df <- data.frame(
  target_id  = rep(targets, 3),
  condition  = factor(rep(c("M", "S", "C"), each = length(targets)), levels = c("M", "S", "C")),
  sigma_late = c(sigM, sigS, sigC)
)

ann <- sprintf("Paired Wilcoxon sigma_C > sigma_S\nV = %d, p = %.4f (one-sided)\n%d/%d targets positive Delta sigma_CS\nMedian Delta sigma_CS = %.3f",
               as.integer(stats$wilcoxon_V), stats$wilcoxon_p,
               as.integer(stats$wilcoxon_n_positive_CS), as.integer(stats$wilcoxon_n_pairs),
               stats$wilcoxon_median_delta_CS)

p <- ggplot(df, aes(x = condition, y = sigma_late, group = target_id, color = target_id)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 4) +
  scale_color_manual(values = darj[1:length(targets)], name = "Target Equation") +
  annotate("text", x = 1.5, y = max(df$sigma_late) * 0.95,
           label = ann, size = 3.6, family = "mono", hjust = 0) +
  labs(title = "Late-Phase Plateau-Duration Sigma: M to S to C, Paired by Target",
       x = "Primitive Set (M minimal, S size-only-rich, C connectivity-rich)",
       y = expression("Pooled censored-MLE log-normal " * sigma * " (late-phase plateaus)")) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("sigma_by_condition_paired.png", p, width = 8.5, height = 5.5, dpi = 300, units = "in")
