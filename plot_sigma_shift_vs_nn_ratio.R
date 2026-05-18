library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr); library(ggrepel)
darj <- wes_palette("Darjeeling1", 5)

data <- fromJSON("result.json", simplifyVector = FALSE)
res <- data$results

df <- bind_rows(lapply(res, function(r) {
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

# Bootstrap Spearman CI
set.seed(1)
B <- 2000
rhos <- replicate(B, {
  idx <- sample(seq_len(nrow(wide)), replace = TRUE)
  suppressWarnings(cor(wide$delta_sigma[idx], wide$nn_ratio[idx], method = "spearman"))
})
ci <- quantile(rhos, c(0.025, 0.975), na.rm = TRUE)
rho_point <- cor(wide$delta_sigma, wide$nn_ratio, method = "spearman")
ann <- sprintf("Spearman rho = %.2f\n95%% bootstrap CI [%.2f, %.2f]", rho_point, ci[1], ci[2])

p <- ggplot(wide, aes(x = nn_ratio, y = delta_sigma, color = target_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 4) +
  geom_text_repel(aes(label = target_id), size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = darj[1:5], name = "Target equation") +
  annotate("text", x = Inf, y = Inf, label = ann, hjust = 1.05, vjust = 1.5, size = 3.5) +
  labs(x = "E-graph neutral-network size ratio (rich / minimal)",
       y = expression(Delta*sigma["late"]~"= "*sigma["rich"]~"−"~sigma["minimal"]),
       title = "Per-target σ shift vs. independently measured NN-size inflation") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("plot_sigma_shift_vs_nn_ratio.png", p, width = 7, height = 5, dpi = 300, units = "in")
