library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr)
darj <- wes_palette("Darjeeling1", 5)

data <- fromJSON("result.json", simplifyVector = FALSE)
res <- data$results

build_ccdf <- function(x) {
  x <- sort(x[x > 0])
  n <- length(x)
  data.frame(x = x, ccdf = (n - seq_len(n) + 1) / n)
}

df <- bind_rows(lapply(res, function(r) {
  durs <- unlist(r$pooled_late_durations_sample)
  c <- build_ccdf(durs)
  c$target_id <- r$target_id
  c$primitives <- r$primitives
  c
}))
df$primitives <- factor(df$primitives, levels = c("minimal", "rich"))

p <- ggplot(df, aes(x = x, y = ccdf, color = primitives, linetype = primitives)) +
  geom_step(linewidth = 0.8) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = darj[c(1, 2)], name = "Primitive set") +
  scale_linetype_manual(values = c("minimal" = "dashed", "rich" = "solid"), name = "Primitive set") +
  facet_wrap(~target_id, ncol = 3) +
  labs(x = "Plateau duration (generations, log scale)",
       y = expression(P(X >= x)~"(log scale)"),
       title = "Late-phase plateau-duration CCDFs, by target and primitive set") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.background = element_rect(fill = "grey90", color = NA))

ggsave("plot_tail_ccdf_by_target.png", p, width = 9, height = 6, dpi = 300, units = "in")
