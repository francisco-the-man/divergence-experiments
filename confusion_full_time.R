library(ggplot2); library(wesanderson); library(jsonlite); library(dplyr); library(tidyr)
stats <- fromJSON('stats.json', simplifyVector = FALSE)
darj <- wes_palette('Darjeeling1', 5)

conf_rows <- stats$confusion_row_normalized
bin_levels <- c('NK_K2','NK_K4','NK_K8','NK_K16','RMF_theta0.1','RMF_theta0.3','RMF_theta1.0','RMF_theta3.0')

df <- do.call(rbind, lapply(conf_rows, function(r) {
  vals <- vapply(bin_levels, function(b) as.numeric(r[[b]]), numeric(1))
  data.frame(true_class = r$true_class, pred_class = bin_levels, value = vals,
             stringsAsFactors = FALSE)
}))
df$true_class <- factor(df$true_class, levels = bin_levels)
df$pred_class <- factor(df$pred_class, levels = bin_levels)

ann <- sprintf('Full waiting-time RF\nbal acc = %.4f, 95%% CI [%.4f, %.4f]',
               stats$acc_full_time, stats$ci_full_time_low, stats$ci_full_time_high)

p <- ggplot(df, aes(x = pred_class, y = true_class, fill = value)) +
  geom_tile(color = 'white') +
  geom_text(aes(label = ifelse(value > 0.01, sprintf('%.2f', value), '')),
            size = 3.2, color = 'black', family = 'mono') +
  scale_fill_gradient(low = 'white', high = darj[1], name = 'Row-Norm Fraction',
                      limits = c(0, 1)) +
  scale_y_discrete(limits = rev(bin_levels)) +
  annotate('text', x = 8.6, y = 8, label = ann, hjust = 1, size = 3.4, family = 'mono') +
  labs(title = 'Confusion Matrix: 8-Way Bin Classifier (Full Waiting-Time Features)',
       x = 'Predicted Class', y = 'True Class') +
  theme_classic() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        axis.text.x = element_text(angle = 35, hjust = 1))

ggsave('confusion_full_time.png', p, width = 9, height = 7, dpi = 300, units = 'in')
