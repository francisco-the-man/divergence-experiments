library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr)
darj <- wes_palette("Darjeeling1", 5)

data <- fromJSON("result.json", simplifyVector = FALSE)
results <- data$results

df <- bind_rows(lapply(results, function(r) {
  data.frame(
    target_id = r$target_id,
    primitives = r$primitives,
    sigma_late = r$pooled_sigma_late,
    nn_size = r$neutral_proxy$mean_class_size_weighted,
    csn_R = r$csn_late$R,
    csn_p = r$csn_late$p,
    n_late = r$n_late_events,
    stringsAsFactors = FALSE
  )
}))
df$primitives <- factor(df$primitives, levels = c("minimal", "rich"))

# Print CSN per-target rich results
rich_df <- df[df$primitives == "rich", ]
for (i in seq_len(nrow(rich_df))) {
  cat(sprintf("PRINT csn_rich_%s_R: %.3f  p: %.4f  n_tail: %d\n",
              rich_df$target_id[i], rich_df$csn_R[i], rich_df$csn_p[i], rich_df$n_late[i]))
}

p.a <- ggplot(df, aes(x = target_id, y = sigma_late, fill = primitives)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("minimal" = darj[1], "rich" = darj[2]), name = "Primitive set") +
  labs(x = "Target", y = expression(sigma[late]~"(log-normal fit)"),
       title = "Pooled σ (late phase)") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1))

p.b <- ggplot(df, aes(x = target_id, y = nn_size, fill = primitives)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("minimal" = darj[1], "rich" = darj[2]), name = "Primitive set") +
  scale_y_log10() +
  labs(x = "Target", y = "Mean weighted e-graph class size (log scale)",
       title = "Independent NN-size proxy") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1))

plot <- p.a + p.b + plot_layout(widths = c(1, 1), guides = "collect") &
  theme(legend.position = "bottom")

ggsave("plot_sigma_and_nn.png", plot, width = 9, height = 4.5, dpi = 300, units = "in")
