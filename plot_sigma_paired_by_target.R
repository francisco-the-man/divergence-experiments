library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr)
darj <- wes_palette("Darjeeling1", 5)

data <- fromJSON("result.json", simplifyVector = FALSE)
res <- data$results

df <- bind_rows(lapply(res, function(r) {
  data.frame(
    target_id = r$target_id,
    primitives = r$primitives,
    sigma_late = r$pooled_sigma_late
  )
}))
df$primitives <- factor(df$primitives, levels = c("minimal", "rich"))

# Wilcoxon
wide <- df %>% pivot_wider(names_from = primitives, values_from = sigma_late)
wt <- wilcox.test(wide$rich, wide$minimal, paired = TRUE, alternative = "greater")
subtitle_txt <- sprintf("Paired Wilcoxon (rich > minimal): V = %.0f, p = %.3f; %d/5 targets shift positive",
                       wt$statistic, wt$p.value, sum(wide$rich > wide$minimal))

p <- ggplot(df, aes(x = primitives, y = sigma_late, group = target_id, color = target_id)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3.5) +
  scale_color_manual(values = darj, name = "Target equation") +
  labs(x = "Primitive set",
       y = expression(sigma~"(log-normal fit, late-phase plateau durations)"),
       title = "Plateau heavy-tailedness: rich vs minimal primitive set",
       subtitle = subtitle_txt) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9))

ggsave("plot_sigma_paired_by_target.png", p, width = 7, height = 5, dpi = 300, units = "in")
