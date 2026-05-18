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

ns <- 10^seq(log10(80), log10(15000), length.out = 100)
env_mean <- data.frame(n = ns, upper = 1/sqrt(ns), lower = -1/sqrt(ns))
env_sd <- data.frame(n = ns, upper = 1 + 1/sqrt(2*ns), lower = 1 - 1/sqrt(2*ns))

p.a <- ggplot(df, aes(x = n_samples, y = sample_mean, color = seed)) +
  geom_ribbon(data = env_mean, aes(x = n, ymin = lower, ymax = upper), inherit.aes = FALSE,
              fill = "grey80", alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 3, alpha = 0.85) +
  scale_x_log10() +
  scale_color_manual(values = darj[1:3], name = "Seed") +
  labs(x = "Number of samples (log scale)", y = "Sample mean",
       title = "Sample mean converges to 0") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

p.b <- ggplot(df, aes(x = n_samples, y = sample_stddev, color = seed)) +
  geom_ribbon(data = env_sd, aes(x = n, ymin = lower, ymax = upper), inherit.aes = FALSE,
              fill = "grey80", alpha = 0.4) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_point(size = 3, alpha = 0.85) +
  scale_x_log10() +
  scale_color_manual(values = darj[1:3], name = "Seed") +
  labs(x = "Number of samples (log scale)", y = "Sample standard deviation",
       title = "Sample stddev converges to 1") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

plot <- p.a + p.b + plot_layout(widths = c(1, 1))
ggsave("plot_convergence.png", plot, width = 10, height = 4, dpi = 300, units = "in")
