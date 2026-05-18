library(ggplot2); library(wesanderson); library(dplyr); library(jsonlite); library(tidyr)

data  <- fromJSON("result.json", simplifyVector = FALSE)
stats <- fromJSON("stats.json", simplifyVector = FALSE)
darj  <- wes_palette("Darjeeling1", 5)

targets <- unlist(stats$targets)
df <- data.frame(
  target_id  = rep(targets, each = 2),
  primitives = rep(c("minimal", "rich"), times = length(targets)),
  sigma_late = as.numeric(unlist(lapply(targets, function(t) c(stats$sigma_min_per_target[[t]], stats$sigma_rich_per_target[[t]])))),
  stringsAsFactors = FALSE
)
df$primitives <- factor(df$primitives, levels = c("minimal", "rich"))

lab <- sprintf("Paired Wilcoxon (one-sided, rich > minimal):\nV = %g, p = %.3f  (alpha = %.2f) -> %s\nPositive Delta_sigma in %d / %d targets; mean Delta = %.3f",
               stats$wilcoxon_V, stats$wilcoxon_p, stats$wilcoxon_alpha,
               ifelse(stats$wilcoxon_pass, "PASS", "FAIL"),
               stats$n_positive_deltas, stats$n_targets, stats$mean_delta_sigma)

p <- ggplot(df, aes(x = primitives, y = sigma_late, group = target_id, color = target_id)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 4) +
  scale_color_manual(values = darj, name = "Target equation") +
  annotate("label", x = 1.5, y = max(df$sigma_late) + 0.15,
           label = lab, size = 3.2, fill = "white", label.size = 0.3) +
  labs(title = "Pooled log-normal sigma (late phase): Minimal vs Rich primitives",
       x = "Primitive set",
       y = expression(sigma[late] * " (pooled log-normal fit)")) +
  ylim(min(df$sigma_late) - 0.1, max(df$sigma_late) + 0.5) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "right")

ggsave("sigma_paired_by_target.png", p, width = 8, height = 5.5, dpi = 300, units = "in")
