library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr)
darj <- wes_palette("Darjeeling1", 5)

data <- fromJSON("result.json", simplifyVector = FALSE)
df <- bind_rows(lapply(data$results, function(r) {
  data.frame(n_samples = r$n_samples, seed = factor(r$seed), sample_mean = r$sample_mean)
}))

med_df <- df %>% group_by(n_samples) %>% summarise(med = median(sample_mean), .groups = "drop")

plot <- ggplot(df, aes(x = n_samples, y = sample_mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = c(-0.05, 0.05), linetype = "dotted", color = "firebrick") +
  geom_line(data = med_df, aes(y = med), color = darj[5], linewidth = 0.8) +
  geom_point(aes(color = seed), size = 3, alpha = 0.85) +
  scale_x_log10() +
  scale_color_manual(values = darj[1:3], name = "Seed") +
  labs(x = "Number of samples (log scale)", y = "Sample mean",
       title = "Mean convergence — red dotted = ±0.05 preregistered tolerance") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("plot_mean.png", plot, width = 6, height = 4, dpi = 300, units = "in")
