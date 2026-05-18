library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr)
stats <- fromJSON("stats.json", simplifyVector = TRUE)
darj  <- wes_palette("Darjeeling1", 5)

sw <- stats$sigma_wide
# Build long form just for C vs S
long <- data.frame(
  target_id = rep(sw$target_id, 2),
  condition = rep(c("S", "C"), each = nrow(sw)),
  sigma     = c(sw$sigma_S, sw$sigma_C)
)
long$condition <- factor(long$condition, levels = c("S", "C"))
targets <- sort(unique(long$target_id))
pal <- c(darj, "#666666")[seq_along(targets)]

ann <- sprintf("Wilcoxon V = %d, p = %.4f (one-sided C > S)\nMedian Delta sigma_CS = %.3f, %d/%d positive\nSensitivity (drop T6): V = %d, p = %.4f",
               as.integer(stats$wilcoxon_V), stats$wilcoxon_p,
               stats$median_delta_CS,
               as.integer(stats$n_pos_CS), as.integer(stats$n_pairs_CS),
               as.integer(stats$wilcoxon_V_noT6), stats$wilcoxon_p_noT6)

p <- ggplot(long, aes(x = condition, y = sigma, group = target_id, color = target_id)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 4) +
  scale_color_manual(values = pal, name = "Target Equation") +
  annotate("text", x = 1.5, y = max(long$sigma) + 0.7,
           label = ann, size = 3.6, family = "mono", hjust = 0.5) +
  labs(title = "Late-Phase Sigma: Connectivity-Rich C vs Size-Only-Rich S",
       x = "Primitive Set",
       y = expression("Censored log-normal "*sigma*" (late-phase plateaus)")) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("sigma_paired_by_target.png", p, width = 7.5, height = 5, dpi = 300, units = "in")
