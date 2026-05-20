library(ggplot2); library(wesanderson); library(jsonlite); library(dplyr)
stats <- fromJSON("stats.json", simplifyVector = TRUE)
darj  <- wes_palette("Darjeeling1", 5)

pi_df <- stats$permutation_importance
pi_df$feature <- factor(pi_df$feature, levels = rev(pi_df$feature[order(pi_df$importance)]))
pi_df$group_lab <- factor(pi_df$group, levels = c("original","ccdf","jump","autocorr"),
                          labels = c("Original 4 Scalars","CCDF Bins (20)","Jump Magnitude","Autocorrelation"))

ann <- sprintf("Top feature: %s (importance = %.3f)\nSpearman vs RMF theta: rho = %.3f, p_BH = %.1e",
  stats$spearman_top_feature, stats$top5_features$importance[1], stats$spearman_top_rho, stats$spearman_top_p_bh)

p <- ggplot(pi_df, aes(x = importance, y = feature, fill = group_lab)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.6) +
  scale_fill_manual(values = darj[1:4], name = "Feature Group") +
  annotate("text", x = 0.13, y = 4, label = ann, size = 3.4, family = "mono", hjust = 0) +
  labs(title = "Permutation Feature Importance: 8-Way Classifier, Full Set",
       x = "Drop in Balanced Accuracy When Feature Shuffled",
       y = "Feature") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.y = element_text(size = 7),
        legend.position = "right")

ggsave("feature_importance.png", p, width = 9, height = 7.5, dpi = 300, units = "in")
