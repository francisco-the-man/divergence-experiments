library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr)
stats <- fromJSON("stats.json", simplifyVector = TRUE)
darj  <- wes_palette("Darjeeling1", 5)

cc <- stats$ccdf_df
cc <- cc[is.finite(cc$duration) & cc$duration > 0 & is.finite(cc$ccdf) & cc$ccdf > 0, ]
cc$condition <- factor(cc$condition, levels = c("M", "S", "C"))

csn <- stats$csn_table
# Build annotation per target_id summarising C-cell R / p
csn_C <- csn[csn$condition == "C", ]

ann_df <- data.frame(
  target_id = csn_C$target_id,
  label     = sprintf("C: R=%.1f, p=%.2g", csn_C$R, csn_C$p)
)

p <- ggplot(cc, aes(x = duration, y = ccdf, color = condition)) +
  geom_line(linewidth = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("M" = darj[1], "S" = darj[2], "C" = darj[3]),
                     name = "Primitive Set") +
  facet_wrap(~ target_id, ncol = 3, scales = "free") +
  geom_text(data = ann_df, aes(x = 2, y = 0.002, label = label),
            inherit.aes = FALSE, size = 3.0, family = "mono", hjust = 0) +
  labs(title = "Pooled Late-Phase Plateau-Duration CCDF: Log-Log by Target",
       x = "Plateau Duration (Generations)",
       y = expression(P(X >= x))) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.text = element_text(face = "bold"))

ggsave("tail_ccdf_by_condition.png", p, width = 11, height = 7, dpi = 300, units = "in")
