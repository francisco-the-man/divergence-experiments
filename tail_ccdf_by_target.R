library(ggplot2); library(wesanderson); library(dplyr); library(jsonlite); library(tidyr)
data  <- fromJSON("result.json", simplifyVector = FALSE)$results
stats <- fromJSON("stats.json", simplifyVector = TRUE)
darj  <- wes_palette("Darjeeling1", 5)

# Build CCDF data per (target, primitives)
ccdf_rows <- list()
for (r in data) {
  durs <- unlist(r$pooled_late_durations_sample)
  durs <- durs[durs > 0]
  if (length(durs) == 0) next
  durs_sorted <- sort(durs)
  n <- length(durs_sorted)
  ccdf <- 1 - (seq_len(n) - 1) / n
  ccdf_rows[[length(ccdf_rows) + 1]] <- data.frame(
    target_id  = r$target_id,
    primitives = r$primitives,
    x = durs_sorted,
    ccdf = ccdf,
    stringsAsFactors = FALSE
  )
}
df <- do.call(rbind, ccdf_rows)
df$primitives <- factor(df$primitives, levels = c("minimal", "rich"))

# Per-target CSN annotation from stats.json
csn <- stats$csn_per_target
csn_df <- data.frame(
  target_id = csn$target_id,
  R         = csn$R,
  p         = csn$p,
  stringsAsFactors = FALSE
)
csn_df$label <- sprintf("rich CSN: R=%.2f\np=%.4f", csn_df$R, csn_df$p)

p <- ggplot(df, aes(x = x, y = ccdf, color = primitives, linetype = primitives)) +
  geom_step(linewidth = 0.8) +
  geom_text(data = csn_df, aes(x = 2, y = 0.02, label = label),
            inherit.aes = FALSE, hjust = 0, size = 3, family = "mono") +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = darj[c(1, 2)], name = "Primitive Set") +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Primitive Set") +
  facet_wrap(~ target_id, ncol = 3) +
  labs(title = "Late-Phase Plateau-Duration CCDF by Target",
       x = "Plateau duration (generations, log scale)",
       y = expression("P(X \u2265 x), log scale")) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.background = element_rect(fill = "grey90", color = NA),
        strip.text = element_text(face = "bold"))

ggsave("tail_ccdf_by_target.png", p, width = 10, height = 6.5, dpi = 300, units = "in")
