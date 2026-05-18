library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr); library(ggrepel)
darj <- wes_palette("Darjeeling1", 5)

data <- fromJSON("result.json", simplifyVector = FALSE)
res <- data$results

df <- bind_rows(lapply(res, function(r) {
  data.frame(
    target_id = r$target_id,
    primitives = r$primitives,
    sigma_late = r$pooled_sigma_late,
    nn_size = r$neutral_proxy$mean_class_size_weighted
  )
}))

wide <- df %>%
  pivot_wider(names_from = primitives, values_from = c(sigma_late, nn_size)) %>%
  mutate(delta_sigma = sigma_late_rich - sigma_late_minimal,
         nn_ratio = nn_size_rich / nn_size_minimal)

# Spearman + bootstrap CI
set.seed(1)
B <- 2000
boot_rho <- replicate(B, {
  idx <- sample(seq_len(nrow(wide)), replace = TRUE)
  if (length(unique(idx)) < 2) return(NA)
  suppressWarnings(cor(wide$delta_sigma[idx], wide$nn_ratio[idx], method = "spearman"))
})
rho_obs <- cor(wide$delta_sigma, wide$nn_ratio, method = "spearman")
ci <- quantile(boot_rho, c(0.025, 0.975), na.rm = TRUE)
subtitle_txt <- sprintf("Spearman rho = %.2f, bootstrap 95%% CI [%.2f, %.2f] — crosses 0",
                       rho_obs, ci[1], ci[2])

p <- ggplot(wide, aes(x = nn_ratio, y = delta_sigma, color = target_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 5) +
  geom_text_repel(aes(label = target_id), size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = darj, name = "Target equation") +
  labs(x = "Neutral-network size ratio (rich / minimal, e-graph mean class size)",
       y = expression(Delta*sigma[late]~"="~sigma[rich]~-~sigma[minimal]),
       title = "Per-target sigma shift vs independent neutral-network ratio",
       subtitle = subtitle_txt) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9))

ggsave("plot_sigma_shift_vs_nn_ratio.png", p, width = 7, height = 5, dpi = 300, units = "in")
