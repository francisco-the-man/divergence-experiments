library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr)
stats <- fromJSON("stats.json", simplifyVector = TRUE)
darj  <- wes_palette("Darjeeling1", 5)

st <- stats$sigma_table
st$condition <- factor(st$condition, levels = c("M", "S", "C"))
targets <- sort(unique(st$target_id))
# Color palette: 6 targets, recycle Darjeeling1
pal <- c(darj, "#666666")[seq_along(targets)]

ann <- sprintf("Wilcoxon (C > S, paired one-sided): V = %d, p = %.4f\nMedian Delta sigma (C-S) = %.3f; 0/%d targets with C > S",
               as.integer(stats$wilcoxon_V), stats$wilcoxon_p,
               stats$median_delta_CS, as.integer(stats$n_pairs_CS))

p <- ggplot(st, aes(x = condition, y = sigma, group = target_id, color = target_id)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 4) +
  scale_color_manual(values = pal, name = "Target Equation") +
  annotate("text", x = 2, y = max(st$sigma) + 0.7,
           label = ann, size = 3.6, family = "mono", hjust = 0.5) +
  labs(title = "Late-Phase Plateau-Duration Sigma: Paired by Target across M, S, C",
       x = "Primitive Set",
       y = expression("Censored log-normal "*sigma*" (late-phase plateaus)")) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("sigma_by_condition_paired.png", p, width = 8, height = 5.5, dpi = 300, units = "in")
