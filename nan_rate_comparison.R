library(ggplot2); library(wesanderson); library(jsonlite); library(dplyr); library(tidyr)
stats <- fromJSON('stats.json', simplifyVector = FALSE)
darj <- wes_palette('Darjeeling1', 5)

nr <- stats$nan_rates
df <- do.call(rbind, lapply(nr, function(r) {
  data.frame(observable = r$observable, feature = r$feature,
             v2_rate = as.numeric(r$v2_rate), v3_rate = as.numeric(r$v3_rate),
             stringsAsFactors = FALSE)
}))
df_t <- df[df$observable == 'time', ]
df_long <- pivot_longer(df_t, cols = c('v2_rate','v3_rate'),
                        names_to = 'rule', values_to = 'rate')
df_long$rule <- factor(ifelse(df_long$rule == 'v2_rate', 'v2 rule (n≥8)', 'v3 rule (relaxed)'),
                       levels = c('v2 rule (n≥8)', 'v3 rule (relaxed)'))
df_long$feature <- factor(df_long$feature, levels = c('log_mean','log_var','hill','ks_expon'))

ann <- sprintf('McNemar paired test (time, all 4 scalars):\nstat = %.2f, p = %.2e\nv2-only NaN -> non-NaN: %d cells\nv3-only NaN -> non-NaN: %d cells',
               as.numeric(stats$mcnemar_stat), as.numeric(stats$mcnemar_p),
               as.integer(stats$mcnemar_b_v2only), as.integer(stats$mcnemar_c_v3only))

p <- ggplot(df_long, aes(x = feature, y = rate, fill = rule)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf('%.3f', rate)),
            position = position_dodge(width = 0.7), vjust = -0.3,
            size = 3.2, family = 'mono') +
  scale_fill_manual(values = c('v2 rule (n≥8)' = darj[3], 'v3 rule (relaxed)' = darj[1]),
                    name = 'NaN Rule') +
  annotate('text', x = 3.4, y = 0.85, label = ann, size = 3.4, family = 'mono', hjust = 0.5) +
  coord_cartesian(ylim = c(0, 0.95)) +
  labs(title = 'NaN Rate per Scalar Feature: v2 Rule vs v3 Relaxed Rule (Time Observable)',
       x = 'Scalar Feature', y = 'NaN Rate Across 320 Trajectories') +
  theme_classic() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        legend.position = 'top')

ggsave('nan_rate_comparison.png', p, width = 9, height = 5.5, dpi = 300, units = 'in')
