library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr); library(ggrepel)
stats <- fromJSON("stats.json", simplifyVector = TRUE)
darj  <- wes_palette("Darjeeling1", 5)

dt <- stats$delta_table
targets <- sort(unique(dt$target_id))
pal <- c(darj, "#666666")[seq_along(targets)]

ann_size <- sprintf("Spearman rho = %.2f\n95%% CI [%.2f, %.2f]\np = %.3f",
                    stats$rho_size, stats$rho_size_lo, stats$rho_size_hi, stats$rho_size_p)
ann_conn <- sprintf("Spearman rho = %.2f\n95%% CI [%.2f, %.2f]\np = %.3f",
                    stats$rho_conn, stats$rho_conn_lo, stats$rho_conn_hi, stats$rho_conn_p)

p.a <- ggplot(dt, aes(x = delta_size_S, y = delta_sigma_S, color = target_id)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = target_id), size = 3, show.legend = FALSE) +
  scale_color_manual(values = pal, name = "Target Equation") +
  annotate("text", x = min(dt$delta_size_S) + diff(range(dt$delta_size_S))*0.25,
           y = max(dt$delta_sigma_S) - diff(range(dt$delta_sigma_S))*0.05,
           label = ann_size, size = 3.6, family = "mono", hjust = 0.5) +
  labs(title = "Size Path: Delta-Sigma_S vs Delta-Size_S",
       x = expression(Delta*" Mean Banded Class Size (S - M)"),
       y = expression(Delta*sigma[S]*" (late-phase)")) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "none")

p.b <- ggplot(dt, aes(x = delta_conn_C, y = delta_sigma_C, color = target_id)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = target_id), size = 3, show.legend = FALSE) +
  scale_color_manual(values = pal, name = "Target Equation") +
  annotate("text", x = min(dt$delta_conn_C) + diff(range(dt$delta_conn_C))*0.25,
           y = max(dt$delta_sigma_C) - diff(range(dt$delta_sigma_C))*0.05,
           label = ann_conn, size = 3.6, family = "mono", hjust = 0.5) +
  labs(title = "Connectivity Path: Delta-Sigma_C vs Delta-Connectivity_C",
       x = expression(Delta*" Mean Banded Out-Degree (C - M)"),
       y = expression(Delta*sigma[C]*" (late-phase)")) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

p <- p.a + p.b + plot_layout(widths = c(1, 1))

ggsave("delta_sigma_vs_delta_static.png", p, width = 11, height = 5, dpi = 300, units = "in")
