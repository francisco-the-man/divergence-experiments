library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr)
darj <- wes_palette("Darjeeling1", 5)

data <- fromJSON("result.json", simplifyVector = FALSE)
results <- data$results

df <- bind_rows(lapply(results, function(r) {
  durs <- unlist(r$pooled_late_durations_sample)
  if (length(durs) == 0) return(NULL)
  durs <- durs[durs > 0]
  durs_sorted <- sort(durs)
  n <- length(durs_sorted)
  ccdf <- 1 - (seq_len(n) - 1) / n
  data.frame(
    target_id = r$target_id,
    primitives = r$primitives,
    x = durs_sorted,
    ccdf = ccdf,
    stringsAsFactors = FALSE
  )
}))
df$primitives <- factor(df$primitives, levels = c("minimal", "rich"))

p <- ggplot(df, aes(x = x, y = ccdf, color = primitives, linetype = primitives)) +
  geom_step(linewidth = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("minimal" = darj[1], "rich" = darj[2]), name = "Primitive set") +
  scale_linetype_manual(values = c("minimal" = "solid", "rich" = "dashed"), name = "Primitive set") +
  facet_wrap(~ target_id, ncol = 3) +
  labs(x = "Plateau duration (generations, log scale)",
       y = "P(X ≥ x)  (log scale)",
       title = "Late-phase plateau-duration CCDF, by target") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave("plot_ccdf_by_target.png", p, width = 8, height = 5.5, dpi = 300, units = "in")
