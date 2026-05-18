library(ggplot2); library(wesanderson); library(dplyr); library(jsonlite); library(tidyr)
data  <- fromJSON("result.json", simplifyVector = TRUE)$results
stats <- fromJSON("stats.json", simplifyVector = FALSE)
darj  <- wes_palette("Darjeeling1", 5)

df <- data.frame(
  target_id = data$target_id,
  primitives = data$primitives,
  sigma_late = data$pooled_sigma_late,
  stringsAsFactors = FALSE
)
df$primitives <- factor(df$primitives, levels = c("minimal", "rich"))

ann <- sprintf("Paired Wilcoxon (one-sided, rich > minimal):\nV = %g, p = %.3f  -  %s",
               stats$wilcoxon_V, stats$wilcoxon_p,
               ifelse(stats$wilcoxon_pass, "PASS", "FAIL (direction mixed)"))

p <- ggplot(df, aes(x = primitives, y = sigma_late, group = target_id, color = target_id)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = darj, name = "Target equation") +
  labs(title = "Late-phase plateau log-normal sigma:\nminimal vs rich primitive set",
       x = "Primitive set",
       y = "Pooled log-normal sigma (late phase)",
       caption = ann) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.caption = element_text(hjust = 0.5, size = 10))

ggsave("sigma_paired_by_target.png", p, width = 7, height = 5, dpi = 300, units = "in")
