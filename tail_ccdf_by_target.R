library(ggplot2); library(wesanderson); library(dplyr); library(jsonlite); library(tidyr)

raw   <- fromJSON("result.json", simplifyVector = FALSE)$results
stats <- fromJSON("stats.json", simplifyVector = FALSE)
darj  <- wes_palette("Darjeeling1", 5)

rows <- list()
for (r in raw) {
  dur <- unlist(r$pooled_late_durations_sample)
  rows[[length(rows) + 1]] <- data.frame(
    target_id = r$target_id,
    primitives = r$primitives,
    duration = dur
  )
}
df <- do.call(rbind, rows)
df$primitives <- factor(df$primitives, levels = c("minimal", "rich"))

ccdf <- df %>%
  group_by(target_id, primitives) %>%
  arrange(duration) %>%
  mutate(rank = row_number(),
         n = n(),
         ccdf = 1 - (rank - 1) / n) %>%
  ungroup() %>%
  filter(duration > 0)

# annotation strings per target
csn_labels <- sapply(stats$csn_per_target, function(x) {
  sprintf("%s: R=%.2f, p=%.3f", x$target_id, x$R, x$p)
})
ann_df <- data.frame(
  target_id = sapply(stats$csn_per_target, function(x) x$target_id),
  label = sapply(stats$csn_per_target, function(x) sprintf("CSN rich: R=%.2f\np=%.3f", x$R, x$p))
)

p <- ggplot(ccdf, aes(x = duration, y = ccdf, color = primitives, linetype = primitives)) +
  geom_step(linewidth = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = darj[c(1, 2)], name = "Primitive set") +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Primitive set") +
  facet_wrap(~ target_id, ncol = 3) +
  geom_text(data = ann_df, aes(x = 1.5, y = 0.02, label = label),
            inherit.aes = FALSE, size = 2.8, hjust = 0) +
  labs(title = "Late-phase plateau-duration CCDF (rich vs minimal, per target)",
       x = "Plateau duration (generations, log scale)",
       y = "P(X >= x) (log scale)") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.text = element_text(face = "bold"))

ggsave("tail_ccdf_by_target.png", p, width = 9, height = 6, dpi = 300, units = "in")
