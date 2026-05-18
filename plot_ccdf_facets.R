library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr)
darj <- wes_palette("Darjeeling1", 5)

data <- fromJSON("result.json", simplifyVector = FALSE)
res <- data$results

ccdf_df <- bind_rows(lapply(res, function(r) {
  d <- unlist(r$pooled_late_durations_sample)
  d <- d[d > 0]
  if (length(d) == 0) return(NULL)
  ds <- sort(d)
  n <- length(ds)
  data.frame(
    target_id = r$target_id,
    primitives = r$primitives,
    x = ds,
    ccdf = (n - seq_len(n) + 1) / n,
    stringsAsFactors = FALSE
  )
}))
ccdf_df$primitives <- factor(ccdf_df$primitives, levels = c("minimal", "rich"))

p <- ggplot(ccdf_df, aes(x = x, y = ccdf, color = primitives, linetype = primitives)) +
  geom_step(linewidth = 0.8) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = darj[c(1, 2)], name = "Primitive set") +
  scale_linetype_manual(values = c("minimal" = "dashed", "rich" = "solid"),
                        name = "Primitive set") +
  facet_wrap(~ target_id, ncol = 3) +
  labs(x = "Plateau duration (generations, log scale)",
       y = expression(P(X >= x)~"(log scale)"),
       title = "Pooled late-phase plateau-duration CCDFs by target") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave("plot_ccdf_facets.png", p, width = 9, height = 6, dpi = 300, units = "in")
