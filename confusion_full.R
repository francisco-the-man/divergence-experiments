library(ggplot2); library(wesanderson); library(jsonlite); library(dplyr); library(tidyr)
stats <- fromJSON("stats.json", simplifyVector = TRUE)
darj  <- wes_palette("Darjeeling1", 5)

cm <- stats$confusion
# Order rows/cols: NK by K, RMF by theta
order_cells <- c("NK_K2","NK_K4","NK_K8","NK_K16","RMF_theta0.1","RMF_theta0.3","RMF_theta1.0","RMF_theta3.0")
order_cells <- order_cells[order_cells %in% unique(c(cm$true, cm$pred))]
cm$true <- factor(cm$true, levels = order_cells)
cm$pred <- factor(cm$pred, levels = order_cells)

# Row-normalize for visual (per-true recall)
cm <- cm %>% group_by(true) %>% mutate(frac = count / sum(count)) %>% ungroup()

ann <- sprintf("8-way bal acc = %.3f, 95%% CI [%.3f, %.3f]\nThreshold 0.40 not met — hypothesis FALSIFIED",
  stats$bin_classifier_full_acc, stats$bin_classifier_full_ci_low, stats$bin_classifier_full_ci_high)

p <- ggplot(cm, aes(x = pred, y = true, fill = frac)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%d", count)), size = 3.2, family = "mono") +
  scale_fill_gradient(low = "white", high = darj[1], name = "Row Fraction", limits = c(0, 1)) +
  labs(title = "Confusion Matrix: 8-Way Bin Classifier, Full Feature Set",
       subtitle = ann,
       x = "Predicted Cell", y = "True Cell") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, family = "mono", size = 9),
        axis.text.x = element_text(angle = 35, hjust = 1))

ggsave("confusion_full.png", p, width = 8, height = 6.5, dpi = 300, units = "in")
