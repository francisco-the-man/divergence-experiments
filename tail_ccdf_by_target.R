library(ggplot2); library(wesanderson); library(dplyr); library(jsonlite); library(tidyr)
stats <- fromJSON("stats.json", simplifyVector = TRUE)
data <- fromJSON("result.json", simplifyVector = FALSE)$results
darj <- wes_palette("Darjeeling1", 5)

rows <- list()
for (r in data) {
  durs <- unlist(r$pooled_late_durations_sample)
  durs <- durs[durs > 0]
  if (length(durs) == 0) next
  s <- sort(durs)
  ccdf <- (length(s):1) / length(s)
  rows[[length(rows) + 1]] <- data.frame(
    target_id = r$target_id,
    primitives = r$primitives,
    x = s,
    ccdf = ccdf
  )
}
df <- do.call(rbind, rows)
df$primitives <- factor(df$primitives, levels = c("minimal", "rich"))

# Per-target CSN annotation
csn <- stats$csn_rich$per_target
csn_lab <- data.frame(
  target_id = names(csn),
  label = sapply(csn, function(x) sprintf("rich CSN R=%.2f, p=%.3f", x$R, x$p))
)

p <- ggplot(df, aes(x = x, y = ccdf, color = primitives, linetype = primitives)) +
  geom_step(linewidth = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = darj[c(1, 2)], name = "Primitive set") +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Primitive set") +
  facet_wrap(~ target_id, ncol = 3) +
  geom_text(data = csn_lab, aes(x = 2, y = 0.015, label = label),
            inherit.aes = FALSE, size = 2.7, hjust = 0) +
  labs(
    title = "Late-phase plateau-duration CCDFs (log-log), by target and primitive set",
    x = "Plateau duration (generations, log scale)",
    y = "P(X >= x), log scale"
  ) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.text = element_text(face = "bold"))

ggsave("tail_ccdf_by_target.png", p, width = 10, height = 6.5, dpi = 300, units = "in")
