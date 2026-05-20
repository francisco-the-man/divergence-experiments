library(ggplot2); library(wesanderson); library(jsonlite); library(dplyr)
stats <- fromJSON("stats.json", simplifyVector = TRUE)
darj  <- wes_palette("Darjeeling1", 5)

df <- data.frame(
  feature_set = factor(c("4-scalar","+CCDF","Full"), levels = c("4-scalar","+CCDF","Full")),
  acc = c(stats$bin_classifier_4_acc, stats$bin_classifier_pcc_acc, stats$bin_classifier_full_acc),
  lo  = c(stats$bin_classifier_4_ci_low, stats$bin_classifier_pcc_ci_low, stats$bin_classifier_full_ci_low),
  hi  = c(stats$bin_classifier_4_ci_high, stats$bin_classifier_pcc_ci_high, stats$bin_classifier_full_ci_high)
)

ann <- sprintf("Full vs 4-scalar: Delta = %.3f, 95%% CI [%.3f, %.3f]\n+CCDF vs 4-scalar: Delta = %.3f, 95%% CI [%.3f, %.3f]\nFull bin acc = %.3f, 95%% CI [%.3f, %.3f] — upper < 0.40 ⇒ FALSIFIED",
  stats$delta_full_minus_4_mean, stats$delta_full_minus_4_ci_low, stats$delta_full_minus_4_ci_high,
  stats$delta_pcc_minus_4_mean, stats$delta_pcc_minus_4_ci_low, stats$delta_pcc_minus_4_ci_high,
  stats$bin_classifier_full_acc, stats$bin_classifier_full_ci_low, stats$bin_classifier_full_ci_high)

p <- ggplot(df, aes(x = feature_set, y = acc, color = feature_set)) +
  geom_hline(yintercept = stats$hypothesis_threshold_bin, linetype = "dashed", color = "grey25", linewidth = 0.8) +
  geom_hline(yintercept = stats$chance_accuracy, linetype = "dotted", color = "grey45", linewidth = 0.8) +
  geom_hline(yintercept = stats$parent_bin_acc, linetype = "dotdash", color = "grey55", linewidth = 0.8) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.8) +
  geom_point(size = 5) +
  annotate("text", x = 0.7, y = stats$hypothesis_threshold_bin + 0.012, label = "Hypothesis threshold 0.40", hjust = 0, size = 3.3, color = "grey25") +
  annotate("text", x = 0.7, y = stats$chance_accuracy + 0.012, label = "Chance 0.125", hjust = 0, size = 3.3, color = "grey45") +
  annotate("text", x = 0.7, y = stats$parent_bin_acc + 0.012, label = sprintf("Parent v1 = %.3f", stats$parent_bin_acc), hjust = 0, size = 3.3, color = "grey55") +
  annotate("text", x = 2, y = 0.50, label = ann, size = 3.4, family = "mono") +
  scale_color_manual(values = darj[c(1,3,2)], name = "Feature Set") +
  scale_y_continuous(limits = c(0, 0.6)) +
  labs(title = "Ablation: 8-Way Bin Classifier Balanced Accuracy by Feature Set",
       x = "Feature Set", y = "Balanced Accuracy (5-Fold CV, 95% Bootstrap CI)") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "none")

ggsave("ablation_accuracy.png", p, width = 8, height = 5.5, dpi = 300, units = "in")
