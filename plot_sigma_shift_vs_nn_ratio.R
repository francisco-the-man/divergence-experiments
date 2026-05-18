library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr); library(ggrepel)
darj <- wes_palette("Darjeeling1", 5)

data <- fromJSON("result.json", simplifyVector = FALSE)
results <- data$results

df <- bind_rows(lapply(results, function(r) {
  data.frame(
    target_id = r$target_id,
    primitives = r$primitives,
    sigma_late = r$pooled_sigma_late,
    nn_mean = r$neutral_proxy$mean_class_size_weighted,
    stringsAsFactors = FALSE
  )
}))

agg <- df %>%
  pivot_wider(names_from = primitives, values_from = c(sigma_late, nn_mean)) %>%
  mutate(delta_sigma = sigma_late_rich - sigma_late_minimal,
         nn_ratio    = nn_mean_rich / nn_mean_minimal)

print(agg)

# Spearman + bootstrap CI
set.seed(1)
rho <- cor(agg$delta_sigma, agg$nn_ratio, method = "spearman")
B <- 5000
boot_rho <- numeric(B)
n <- nrow(agg)
for (i in seq_len(B)) {
  idx <- sample.int(n, n, replace = TRUE)
  if (length(unique(idx)) < 2) { boot_rho[i] <- NA; next }
  r_ <- suppressWarnings(cor(agg$delta_sigma[idx], agg$nn_ratio[idx], method = "spearman"))
  boot_rho[i] <- r_
}
ci <- quantile(boot_rho, c(0.025, 0.975), na.rm = TRUE)

cat(sprintf("PRINT spearman_rho_delta_sigma_vs_nn_ratio: %.3f\n", rho))
cat(sprintf("PRINT spearman_ci_low: %.3f\n", ci[1]))
cat(sprintf("PRINT spearman_ci_high: %.3f\n", ci[2]))

annot <- sprintf("Spearman ρ = %.3f, bootstrap 95%% CI [%.3f, %.3f]  (n = 5 targets)",
                  rho, ci[1], ci[2])

p <- ggplot(agg, aes(x = nn_ratio, y = delta_sigma, color = target_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = target_id), size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = darj[1:5], name = "Target") +
  labs(x = "Neutral-network size ratio (rich / minimal), from e-graph proxy",
       y = expression(Delta*sigma[late]~"= "*sigma[rich]~-~sigma[minimal]),
       title = "Per-target plateau-σ shift vs neutral-network-size ratio",
       subtitle = annot) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9))

ggsave("plot_sigma_shift_vs_nn_ratio.png", p, width = 7, height = 5, dpi = 300, units = "in")
