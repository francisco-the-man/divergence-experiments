library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr)
darj <- wes_palette("Darjeeling1", 5)

data <- fromJSON("result.json", simplifyVector = FALSE)
results <- data$results

df <- bind_rows(lapply(results, function(r) {
  data.frame(
    target_id = r$target_id,
    primitives = r$primitives,
    sigma_late = r$pooled_sigma_late,
    stringsAsFactors = FALSE
  )
}))

df$primitives <- factor(df$primitives, levels = c("minimal", "rich"))

# paired wilcoxon (one-sided, rich > minimal)
wide <- df %>% pivot_wider(names_from = primitives, values_from = sigma_late)
wt <- wilcox.test(wide$rich, wide$minimal, paired = TRUE, alternative = "greater", exact = FALSE)
cat(sprintf("PRINT wilcox_V: %.3f\n", wt$statistic))
cat(sprintf("PRINT wilcox_p: %.4f\n", wt$p.value))
cat(sprintf("PRINT n_targets_up: %d\n", sum(wide$rich > wide$minimal)))
cat(sprintf("PRINT n_targets_down: %d\n", sum(wide$rich < wide$minimal)))

annot <- sprintf("Paired Wilcoxon (one-sided, rich > minimal):\nV = %.0f, p = %.3f\nrich > minimal in %d / 5 targets",
                 wt$statistic, wt$p.value, sum(wide$rich > wide$minimal))

p <- ggplot(df, aes(x = primitives, y = sigma_late, group = target_id, color = target_id)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = darj[1:5], name = "Target") +
  annotate("text", x = 1.5, y = max(df$sigma_late) + 0.05, label = annot,
           size = 3.2, hjust = 0.5, vjust = 1) +
  labs(x = "Primitive set",
       y = expression(sigma~"(log-normal fit, late-phase plateau durations)"),
       title = "Plateau-duration heavy-tailedness: minimal vs rich primitive set") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("plot_sigma_paired_by_target.png", p, width = 7, height = 5, dpi = 300, units = "in")
