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

# Paired Wilcoxon
wide <- df %>% pivot_wider(names_from = primitives, values_from = sigma_late)
wt <- wilcox.test(wide$rich, wide$minimal, paired = TRUE, alternative = "greater")
mean_delta <- mean(wide$rich - wide$minimal)
n_pos <- sum((wide$rich - wide$minimal) > 0)
n_neg <- sum((wide$rich - wide$minimal) < 0)

cat(sprintf("PRINT wilcoxon_V: %.3f\n", as.numeric(wt$statistic)))
cat(sprintf("PRINT wilcoxon_p_one_sided_greater: %.3f\n", wt$p.value))
cat(sprintf("PRINT mean_delta_sigma: %.3f\n", mean_delta))
cat(sprintf("PRINT n_positive_shifts: %d\n", n_pos))
cat(sprintf("PRINT n_negative_shifts: %d\n", n_neg))

annot <- sprintf("Paired Wilcoxon (rich > minimal): V = %d, p = %.3f\nMean Δσ = %+.3f; direction %d+ / %d−",
                 as.integer(wt$statistic), wt$p.value, mean_delta, n_pos, n_neg)

p <- ggplot(df, aes(x = primitives, y = sigma_late, group = target_id, color = target_id)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = darj[1:5], name = "Target") +
  labs(x = "Primitive set",
       y = expression(sigma~"(log-normal fit, late-phase plateau durations)"),
       title = "Paired plateau-tail σ by target: rich vs minimal",
       subtitle = annot) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9))

ggsave("plot_sigma_paired_by_target.png", p, width = 7, height = 5, dpi = 300, units = "in")
