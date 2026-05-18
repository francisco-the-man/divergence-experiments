library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr)
darj <- wes_palette("Darjeeling1", 5)

data <- fromJSON("result.json", simplifyVector = FALSE)
res <- data$results

df <- bind_rows(lapply(res, function(r) {
  data.frame(
    target_id = r$target_id,
    primitives = r$primitives,
    sigma_late = r$pooled_sigma_late,
    stringsAsFactors = FALSE
  )
}))
df$primitives <- factor(df$primitives, levels = c("minimal", "rich"))

p <- ggplot(df, aes(x = primitives, y = sigma_late, group = target_id, color = target_id)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = darj[1:5], name = "Target equation") +
  labs(x = "Primitive set",
       y = expression(sigma~"(log-normal fit, late-phase plateaus)"),
       title = "Within-target paired shift in plateau heavy-tailedness",
       subtitle = "Mechanism predicts all lines slope up; 3 of 5 slope down") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("plot_sigma_paired.png", p, width = 7, height = 5, dpi = 300, units = "in")
