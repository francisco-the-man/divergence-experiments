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
  mutate(
    delta_sigma = sigma_late_rich - sigma_late_minimal,
    nn_ratio = nn_mean_rich / nn_mean_minimal
  )

# bootstrap Spearman CI across the 5 targets
set.seed(1)
B <- 2000
n <- nrow(wide)
rhos <- replicate(B, {
  idx <- sample(seq_len(n), n, replace = TRUE)
  if (length(unique(idx)) < 2) return(NA_real_)
  suppressWarnings(cor(wide$delta_sigma[idx], wide$nn_ratio[idx], method = "spearman"))
})
rhos <- rhos[is.finite(rhos)]
rho_obs <- suppressWarnings(cor(wide$delta_sigma, wide$nn_ratio, method = "spearman"))
ci <- quantile(rhos, c(0.025, 0.975), na.rm = TRUE)
label_txt <- sprintf("Spearman rho = %.2f\n95%% bootstrap CI: [%.2f, %.2f]", rho_obs, ci[1], ci[2])

p <- ggplot(wide, aes(x = nn_ratio, y = delta_sigma, color = target_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 4) +
  geom_text_repel(aes(label = target_id), size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = darj[1:5], name = "Target") +
  annotate("text", x = Inf, y = Inf, label = label_txt,
           hjust = 1.05, vjust = 1.3, size = 3.5) +
  labs(
    x = "Neutral-network size ratio (rich / minimal, e-graph proxy)",
    y = expression(Delta*sigma[late]~"= "*sigma[rich]~"− "*sigma[minimal]),
    title = "Per-target sigma shift vs neutral-network size ratio"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

ggsave("plot_sigma_shift_vs_nn_ratio.png", p, width = 7.5, height = 5, dpi = 300, units = "in")
