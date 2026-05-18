library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr); library(ggrepel)
darj <- wes_palette("Darjeeling1", 5)

data <- fromJSON("result.json", simplifyVector = FALSE)
res <- data$results

df <- bind_rows(lapply(res, function(r) {
  data.frame(
    target_id = r$target_id,
    primitives = r$primitives,
    sigma_late = r$pooled_sigma_late,
    nn_proxy = r$neutral_proxy$mean_class_size_weighted,
    stringsAsFactors = FALSE
  )
}))

wide <- df %>%
  pivot_wider(id_cols = target_id, names_from = primitives, values_from = c(sigma_late, nn_proxy)) %>%
  mutate(delta_sigma = sigma_late_rich - sigma_late_minimal,
         nn_ratio = nn_proxy_rich / nn_proxy_minimal)

rho <- suppressWarnings(cor(wide$delta_sigma, wide$nn_ratio, method = "spearman"))

p <- ggplot(wide, aes(x = nn_ratio, y = delta_sigma, color = target_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = target_id), size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = darj[1:5], name = "Target equation") +
  labs(x = "E-graph neutral-network size ratio (rich / minimal)",
       y = expression(Delta*sigma[late]~"= "*sigma[rich]~"−"~sigma[minimal]),
       title = "Per-target shift in tail heaviness vs. neutral-network inflation",
       subtitle = paste0("Spearman ρ = ", sprintf("%.2f", rho), " (n=5; preregistered, low-power)")) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("plot_sigma_shift_vs_nn_ratio.png", p, width = 6.5, height = 4.5, dpi = 300, units = "in")
