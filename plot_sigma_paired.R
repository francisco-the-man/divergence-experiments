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

# paired wilcoxon for annotation
wide <- df %>% pivot_wider(names_from = primitives, values_from = sigma_late)
wt <- wilcox.test(wide$rich, wide$minimal, paired = TRUE, alternative = "greater")
annot <- sprintf("Paired Wilcoxon (rich > minimal): V = %g, p = %.3f", wt$statistic, wt$p.value)
cat(sprintf("PRINT wilcoxon_V: %g\n", wt$statistic))
cat(sprintf("PRINT wilcoxon_p: %.3f\n", wt$p.value))
for (i in seq_len(nrow(wide))) {
  cat(sprintf("PRINT delta_sigma_%s: %.3f\n", wide$target_id[i], wide$rich[i] - wide$minimal[i]))
}

p <- ggplot(df, aes(x = primitives, y = sigma_late, group = target_id, color = target_id)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = darj[1:5], name = "Target") +
  labs(x = "Primitive set",
       y = expression(sigma~"(log-normal fit, late-phase plateau durations)"),
       title = "Within-target paired σ: minimal vs rich",
       subtitle = annot) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("plot_sigma_paired.png", p, width = 6.5, height = 4.5, dpi = 300, units = "in")
