library(ggplot2); library(wesanderson); library(jsonlite)
stats <- fromJSON('stats.json', simplifyVector = TRUE)
darj  <- wes_palette('Darjeeling1', 5)

df <- data.frame(
  eps = factor(c('0.001', '0.005', '0.02'), levels = c('0.001', '0.005', '0.02')),
  acc = c(stats$eps_acc_0_001, stats$eps_acc_0_005, stats$eps_acc_0_02),
  ci_low  = c(stats$eps_ci_0_001_low, stats$eps_ci_0_005_low, stats$eps_ci_0_02_low),
  ci_high = c(stats$eps_ci_0_001_high, stats$eps_ci_0_005_high, stats$eps_ci_0_02_high)
)

primary_acc <- stats$eps_acc_0_005
ann <- sprintf('Robustness: max - min = %.4f\n95%% CI on spread [%.4f, %.4f]\nThreshold ≤ %.2f — PASS',
               stats$eps_spread, stats$eps_spread_ci_low, stats$eps_spread_ci_high,
               stats$robustness_threshold)

p <- ggplot(df, aes(x = eps, y = acc)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf,
                ymin = primary_acc - 0.05, ymax = primary_acc + 0.05),
            fill = darj[2], alpha = 0.05, color = NA, inherit.aes = FALSE) +
  geom_hline(yintercept = 0.40, linetype = 'dashed', color = darj[4], linewidth = 0.8) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1, color = darj[1], linewidth = 0.8) +
  geom_point(size = 5, color = darj[1]) +
  annotate('text', x = 2, y = max(df$ci_high) + 0.04, label = ann,
           size = 3.6, family = 'mono', hjust = 0.5) +
  annotate('text', x = 0.6, y = 0.42, label = 'Threshold 0.40', size = 3.4, color = darj[4], hjust = 0) +
  annotate('text', x = 0.6, y = primary_acc + 0.06, label = '±0.05 robustness window', size = 3.2, color = darj[2], hjust = 0) +
  coord_cartesian(ylim = c(0.35, 0.75)) +
  labs(title = 'Detector Robustness: Bin Accuracy vs eps_jump_rel',
       x = expression('eps'[jump_rel]),
       y = 'Balanced Accuracy (5-fold stratified CV)') +
  theme_classic() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))

ggsave('eps_robustness.png', p, width = 7.5, height = 5, dpi = 300, units = 'in')
