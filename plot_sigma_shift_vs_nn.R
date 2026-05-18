library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr); library(ggrepel)
darj <- wes_palette("Darjeeling1", 5)

data <- fromJSON("result.json", simplifyVector = FALSE)
results <- data$results

df <- bind_rows(lapply(results, function(r) {
  data.frame(
    target_id = r$target_id,
    primitives = r$primitives,
    sigma_late = r$pooled_sigma_late,
    nn_size = r$neutral_proxy$mean_class_size_weighted,
    stringsAsFactors = FALSE
  )
}))

wide <- df %>%
  pivot_wider(names_from = primitives, values_from = c(sigma_late, nn_size)) %>%
  mutate(delta_sigma = sigma_late_rich - sigma_late_minimal,
         nn_ratio = nn_size_rich / nn_size_minimal)

rho_obj <- suppressWarnings(cor.test(wide$delta_sigma, wide$nn_ratio, method = "spearman"))
rho <- as.numeric(rho_obj$estimate)

set.seed(42)
B <- 5000
boot_rhos <- replicate(B, {
  idx <- sample(seq_len(nrow(wide)), replace = TRUE)
  if (length(unique(wide$delta_sigma[idx])) < 2 || length(unique(wide$nn_ratio[idx])) < 2) return(NA_real_)
  suppressWarnings(cor(wide$delta_sigma[idx], wide$nn_ratio[idx], method = "spearman"))
})
ci <- quantile(boot_rhos, c(0.025, 0.975), na.rm = TRUE)

cat(sprintf("PRINT spearman_rho: %.3f\n", rho))
cat(sprintf("PRINT spearman_ci_low: %.3f\n", ci[1]))
cat(sprintf("PRINT spearman_ci_high: %.3f\n", ci[2]))

annot <- sprintf("Spearman ρ = %.3f   95%% bootstrap CI [%.3f, %.3f]", rho, ci[1], ci[2])

p <- ggplot(wide, aes(x = nn_ratio, y = delta_sigma, color = target_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 4) +
  geom_text_repel(aes(label = target_id), size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = darj[1:5], name = "Target") +
  labs(x = "Neutral-network size ratio (rich / minimal, e-graph)",
       y = expression(Delta*sigma[late]~"= "~sigma[rich]~"−"~sigma[minimal]),
       title = "Per-target Δσ vs e-graph neutral-network ratio",
       subtitle = annot) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("plot_sigma_shift_vs_nn.png", p, width = 6.5, height = 4.5, dpi = 300, units = "in")
