library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr)
darj <- wes_palette("Darjeeling1", 5)

data <- fromJSON("result.json", simplifyVector = FALSE)
res <- data$results

df <- bind_rows(lapply(res, function(r) {
  durs <- unlist(r$pooled_late_durations_sample)
  if (length(durs) == 0) return(NULL)
  data.frame(
    target_id = r$target_id,
    primitives = r$primitives,
    duration = as.numeric(durs),
    stringsAsFactors = FALSE
  )
}))

# Compute CCDF per (target, primitives)
ccdf <- df %>%
  group_by(target_id, primitives) %>%
  arrange(duration) %>%
  mutate(rank = row_number(),
         n = n(),
         ccdf = 1 - (rank - 1) / n) %>%
  ungroup()

ccdf$primitives <- factor(ccdf$primitives, levels = c("minimal", "rich"))

p <- ggplot(ccdf, aes(x = duration, y = ccdf, color = primitives, linetype = primitives)) +
  geom_step(linewidth = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = darj[c(1, 3)], name = "Primitive set") +
  scale_linetype_manual(values = c("minimal" = "dashed", "rich" = "solid"), name = "Primitive set") +
  facet_wrap(~ target_id, nrow = 2) +
  labs(x = "Plateau duration in generations (log scale)",
       y = "P(X ≥ x), survival function (log scale)",
       title = "Late-phase plateau-duration CCDFs by target") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.background = element_rect(fill = "grey95", color = NA))

ggsave("plot_tail_ccdf_by_target.png", p, width = 8.5, height = 5.5, dpi = 300, units = "in")
