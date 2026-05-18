library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr)
darj <- wes_palette("Darjeeling1", 5)

data <- fromJSON("result.json", simplifyVector = FALSE)
res <- data$results

df <- bind_rows(lapply(res, function(r) {
  data.frame(
    target_id = r$target_id,
    primitives = r$primitives,
    sigma_early = r$pooled_sigma_early
  )
}))
df$primitives <- factor(df$primitives, levels = c("minimal", "rich"))

p <- ggplot(df, aes(x = primitives, y = sigma_early, group = target_id, color = target_id)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3.5) +
  scale_color_manual(values = darj, name = "Target equation") +
  labs(x = "Primitive set",
       y = expression(sigma[early]~"(log-normal fit, early-phase plateau durations)"),
       title = "Exploratory: early-phase sigma — also fails to support the mechanism",
       subtitle = "Not in pre-registered analysis_plan; shown to verify the null isn't a windowing artifact") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9))

ggsave("plot_sigma_early.png", p, width = 7, height = 5, dpi = 300, units = "in")
