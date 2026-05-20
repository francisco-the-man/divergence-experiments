library(ggplot2); library(wesanderson); library(jsonlite); library(dplyr)
data  <- fromJSON("result.json", simplifyVector = FALSE)
stats <- fromJSON("stats.json", simplifyVector = TRUE)
darj  <- wes_palette("Darjeeling1", 5)

rows <- list()
for (r in data$results) {
  sl <- as.numeric(unlist(r$stasis_lengths))
  if (length(sl) == 0) next
  rows[[length(rows)+1]] <- data.frame(cell = r$cell, family = r$family, stasis = sl)
}
df <- do.call(rbind, rows)

ccdf_df <- df %>%
  group_by(cell, family) %>%
  arrange(stasis) %>%
  mutate(rank = rank(stasis, ties.method = "first"),
         ccdf = 1 - (rank - 1)/n()) %>%
  ungroup()

ccdf_df$cell <- factor(ccdf_df$cell,
  levels = c("NK_K2","NK_K4","NK_K8","NK_K16","RMF_theta0.1","RMF_theta0.3","RMF_theta1.0","RMF_theta3.0"))

ann <- sprintf("Bin classifier full = %.3f [%.3f, %.3f]\nCCDF-bin permutation importance ~ 0 across all 20 bins",
  stats$bin_classifier_full_acc, stats$bin_classifier_full_ci_low, stats$bin_classifier_full_ci_high)

pal8 <- c(darj[1], darj[2], darj[3], darj[4], darj[5], "#7E2954", "#1F8FA6", "#4F8C00")

p <- ggplot(ccdf_df, aes(x = stasis, y = ccdf, color = cell, linetype = family)) +
  geom_step(linewidth = 0.8, alpha = 0.85) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = pal8, name = "Cell") +
  scale_linetype_manual(values = c(NK = "solid", RMF = "dashed"), name = "Family") +
  annotate("text", x = 1.05, y = 0.012, label = ann, size = 3.3, family = "mono", hjust = 0) +
  labs(title = "Stasis-Length CCDFs by Landscape Cell (Pooled Trajectories)",
       x = "Stasis Length (log)", y = "Empirical CCDF (log)") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("ccdf_by_class.png", p, width = 9, height = 6, dpi = 300, units = "in")
