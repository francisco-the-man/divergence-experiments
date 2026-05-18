library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr)
darj <- wes_palette("Darjeeling1", 5)

data <- fromJSON("result.json", simplifyVector = FALSE)
df <- bind_rows(lapply(data$results, function(r) {
  data.frame(
    n_samples = r$n_samples,
    seed = factor(r$seed),
    sample_mean = r$sample_mean,
    sample_stddev = r$sample_stddev
  )
}))

med_mean <- df %>% group_by(n_samples) %>% summarise(m = median(sample_mean), .groups="drop")
med_sd   <- df %>% group_by(n_samples) %>% summarise(m = median(sample_stddev), .groups="drop")

p1 <- ggplot(df, aes(x = n_samples, y = sample_mean, color = seed)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = c(-0.05, 0.05), linetype = "dashed", color = "grey60") +
  geom_line(data = med_mean, aes(x = n_samples, y = m), inherit.aes = FALSE,
            linewidth = 0.8, color = "black") +
  geom_point(size = 3, alpha = 0.85) +
  scale_x_log10() +
  scale_color_manual(values = darj[1:3], name = "Seed") +
  labs(x = "Number of samples (log scale)",
       y = "Sample mean",
       title = "Sample mean converges to 0") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

p2 <- ggplot(df, aes(x = n_samples, y = sample_stddev, color = seed)) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = c(0.95, 1.05), linetype = "dashed", color = "grey60") +
  geom_line(data = med_sd, aes(x = n_samples, y = m), inherit.aes = FALSE,
            linewidth = 0.8, color = "black") +
  geom_point(size = 3, alpha = 0.85) +
  scale_x_log10() +
  scale_color_manual(values = darj[1:3], name = "Seed") +
  labs(x = "Number of samples (log scale)",
       y = "Sample standard deviation",
       title = "Sample stddev converges to 1") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

plot <- p1 + p2 + plot_layout(widths = c(1, 1))

ggsave("plot_convergence.png", plot, width = 10, height = 4, dpi = 300, units = "in")
