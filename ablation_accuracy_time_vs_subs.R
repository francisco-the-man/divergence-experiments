library(ggplot2); library(wesanderson); library(jsonlite); library(dplyr)
stats <- fromJSON('stats.json', simplifyVector = TRUE)
darj  <- wes_palette('Darjeeling1', 5)

feature_levels <- c('4-Scalar', '+CCDF (Stasis)', 'Full')
df <- data.frame(
  feature_set = factor(rep(feature_levels, 2), levels = feature_levels),
  observable  = rep(c('Waiting Time', 'Substitution Count'), each = 3),
  acc = c(stats$acc_4scalar_time, stats$acc_ccdf_time, stats$acc_full_time,
          stats$acc_4scalar_subs, stats$acc_ccdf_subs, stats$acc_full_subs)
)
# Approx 95% CI per cell (Wilson-ish on n=320)
n <- stats$n_total
df$ci_low  <- pmax(0, df$acc - 1.96 * sqrt(df$acc * (1 - df$acc) / n))
df$ci_high <- pmin(1, df$acc + 1.96 * sqrt(df$acc * (1 - df$acc) / n))

ann <- sprintf('PRIMARY: full time acc = %.4f, 95%% CI [%.4f, %.4f]\nthreshold 0.40 cleared (upper CI %.4f)\nSECONDARY: time - subs = %.4f, 95%% CI [%.4f, %.4f]',
               stats$acc_full_time, stats$ci_full_time_low, stats$ci_full_time_high,
               stats$ci_full_time_high,
               stats$delta_time_minus_subs, stats$delta_time_minus_subs_ci_low, stats$delta_time_minus_subs_ci_high)

parent_ann <- sprintf('Parent v2 (subs CCDF, kNN):\nacc = %.4f, CI [%.4f, %.4f]',
                      stats$parent_acc_full, stats$parent_acc_full_ci_low, stats$parent_acc_full_ci_high)

p <- ggplot(df, aes(x = feature_set, y = acc, color = observable, group = observable)) +
  geom_hline(yintercept = stats$chance_level, linetype = 'dotted', color = 'grey40') +
  geom_hline(yintercept = stats$primary_threshold, linetype = 'dashed', color = darj[4], linewidth = 0.8) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = stats$parent_acc_full_ci_low, ymax = stats$parent_acc_full_ci_high),
            fill = 'grey80', alpha = 0.05, color = NA, inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.12, linewidth = 0.8,
                position = position_dodge(width = 0.25)) +
  geom_line(linewidth = 0.8, position = position_dodge(width = 0.25)) +
  geom_point(size = 4, position = position_dodge(width = 0.25)) +
  scale_color_manual(values = c('Waiting Time' = darj[1], 'Substitution Count' = darj[3]),
                     name = 'Stasis Observable') +
  annotate('text', x = 1.05, y = 0.04, label = 'Chance (0.125)', hjust = 0, size = 3.4, color = 'grey40') +
  annotate('text', x = 1.05, y = 0.42, label = 'Threshold (0.40)', hjust = 0, size = 3.4, color = darj[4]) +
  annotate('text', x = 1.05, y = (stats$parent_acc_full_ci_low + stats$parent_acc_full_ci_high) / 2,
           label = parent_ann, hjust = 0, size = 3.2, color = 'grey30', family = 'mono') +
  annotate('text', x = 2.5, y = 0.72, label = ann, hjust = 0.5, size = 3.4, family = 'mono') +
  coord_cartesian(ylim = c(0, 0.78)) +
  labs(title = 'Bin-Classifier Ablation: Waiting Time vs Substitution Count',
       x = 'Feature Set',
       y = 'Balanced Accuracy (5-fold stratified CV)') +
  theme_classic() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        legend.position = 'top')

ggsave('ablation_accuracy_time_vs_subs.png', p, width = 9, height = 6, dpi = 300, units = 'in')
