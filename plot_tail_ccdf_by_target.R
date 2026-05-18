library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr)
darj <- wes_palette("Darjeeling1", 5)

data <- fromJSON("result.json", simplifyVector = FALSE)
res <- data$results

rows <- list()
for (r in res) {
  durs <- unlist(r$pooled_late_durations_sample)
  durs <- durs[is.finite(durs) & durs > 0]
  if (length(durs) == 0) next
  rows[[length(rows) + 1]] <- data.frame(
    target_id = r$target_id,
    primitives = r$primitives,
    duration = durs,
    stringsAsFactors = FALSE
  )
}
df <- bind_rows(rows)
df$primitives <- factor(df$primitives, levels = c("minimal", "rich"))

ccdf <- df %>%
  group_by(target_id, primitives) %>%
  arrange(duration) %>%
  mutate(
    rank = row_number(),
    n = n(),
    ccdf = 1 - (rank - 1) / n
  ) %>%
  ungroup()

p <- ggplot(ccdf, aes(x = duration, y = ccdf, color = primitives, linetype = primitives)) +
  geom_step(linewidth = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = darj[c(1, 2)], name = "Primitive set") +
  scale_linetype_manual(values = c("minimal" = "dashed", "rich" = "solid"), name = "Primitive set") +
  facet_wrap(~ target_id, nrow = 2) +
  labs(
    x = "Plateau duration (generations, log scale)",
    y = expression(P(X >= x)~"(log scale)"),
    title = "Pooled late-phase plateau-duration CCDFs by target"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", color = NA)
  )

ggsave("plot_tail_ccdf_by_target.png", p, width = 9, height = 6, dpi = 300, units = "in")
