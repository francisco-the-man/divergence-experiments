library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr); library(ggrepel)
darj <- wes_palette("Darjeeling1", 5)

data <- fromJSON("result.json", simplifyVector = FALSE)
res <- data$results

df <- bind_rows(lapply(res, function(r) {
  data.frame(
    target_id = r$target_id,
    primitives = r$primitives,
    sigma_late = r$pooled_sigma_late,
    nn_mean = r$neutral_proxy$mean_class_size_weighted,
    stringsAsFactors = FALSE
  )
}))

wide <- df %>%
  pivot_wider(id_cols = target_id, names_from = primitives,
              values_from = c(sigma_late, nn_mean)) %>%
  mutate(delta_sigma = sigma_late_rich - sigma_late_minimal,
         nn_ratio = nn_mean_rich / nn_mean_minimal)

rho <- suppressWarnings(cor(wide$delta_sigma, wide$nn_ratio, method = "spearman"))

# Bootstrap CI on n=5 (informational only)
set.seed(1)
boot_rho <- replicate(2000, {
  i <- sample(seq_len(nrow(wide)), replace = TRUE)
  suppressWarnings(cor(wide$delta_sigma[i], wide$nn_ratio[i], method = "spearman"))
})
ci <- quantile(boot_rho, c(0.025, 0.975), na.rm = TRUE)

lab <- sprintf("Spearman rho = %.2f\n95%% bootstrap CI [%.2f, %.2f]",
               rho, ci[1], ci[2])

p <- ggplot(wide, aes(x = nn_ratio, y = delta_sigma, color = target_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 4) +
  geom_text_repel(aes(label = target_id), size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = darj[1:5], name = "Target equation") +
  annotate("label", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.3,
           label = lab, size = 3.2, label.size = 0) +
  labs(x = "Neutral-network-size ratio (rich / minimal)",
       y = expression(Delta*sigma[late]~"(rich - minimal)"),
       title = "Per-target sigma shift vs. e-graph NN-size ratio",
       subtitle = "Mechanism predicts positive slope; observed slope is negative") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("plot_shift_vs_nn.png", p, width = 7, height = 5, dpi = 300, units = "in")
