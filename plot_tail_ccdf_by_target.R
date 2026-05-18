library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr)
darj <- wes_palette("Darjeeling1", 5)

data <- fromJSON("result.json", simplifyVector = FALSE)
res <- data$results

# Build CCDF data per (target, primitives) from pooled_late_durations_sample
ccdf_list <- list()
for (r in res) {
  durs <- unlist(r$pooled_late_durations_sample)
  durs <- durs[durs > 0]
  if (length(durs) == 0) next
  durs_sorted <- sort(durs)
  n <- length(durs_sorted)
  ccdf <- 1 - (seq_len(n) - 1) / n
  ccdf_list[[length(ccdf_list) + 1]] <- data.frame(
    target_id = r$target_id,
    primitives = r$primitives,
    x = durs_sorted,
    ccdf = ccdf,
    stringsAsFactors = FALSE
  )
}
df <- bind_rows(ccdf_list)
df$primitives <- factor(df$primitives, levels = c("minimal", "rich"))

# Annotation: CSN R and p per (target, rich)
ann <- bind_rows(lapply(res, function(r) {
  if (r$primitives != "rich") return(NULL)
  data.frame(
    target_id = r$target_id,
    label = sprintf("rich: R=%.2f, p=%.3f", r$csn_late$R, r$csn_late$p),
    stringsAsFactors = FALSE
  )
}))

p <- ggplot(df, aes(x = x, y = ccdf, color = primitives, linetype = primitives)) +
  geom_step(linewidth = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = darj[c(1, 2)], name = "Primitive set") +
  scale_linetype_manual(values = c("minimal" = "solid", "rich" = "dashed"), name = "Primitive set") +
  facet_wrap(~ target_id, ncol = 3) +
  geom_text(data = ann, aes(x = 2, y = 0.002, label = label), inherit.aes = FALSE, size = 2.8, hjust = 0) +
  labs(
    x = "Plateau duration (generations, log scale)",
    y = "P(X >= x) (log scale)",
    title = "Late-phase plateau-duration CCDFs by target and primitive set",
    subtitle = "CSN log-normal-vs-exponential R per target (rich condition); R<0 favors exponential"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.background = element_rect(fill = "grey95", color = NA),
    legend.position = "bottom"
  )

ggsave("plot_tail_ccdf_by_target.png", p, width = 9, height = 6, dpi = 300, units = "in")
