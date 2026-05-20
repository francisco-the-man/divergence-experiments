library(ggplot2); library(wesanderson); library(jsonlite)
stats <- fromJSON('stats.json', simplifyVector = TRUE)
darj  <- wes_palette('Darjeeling1', 5)

strat <- stats$stratification
df <- data.frame(
  feature_set = factor(c('Pooled (All Stasis)','Saddle-Only Stratum'), levels=c('Pooled (All Stasis)','Saddle-Only Stratum')),
  acc = c(strat$acc_pooled, strat$acc_saddle_only),
  lo  = c(strat$acc_pooled_ci_low, strat$acc_saddle_only_ci_low),
  hi  = c(strat$acc_pooled_ci_high, strat$acc_saddle_only_ci_high)
)
ann <- sprintf('pooled acc = %.4f [%.4f, %.4f]\nsaddle-only acc = %.4f [%.4f, %.4f]\ndrop = %.4f\nbootstrap 95%% CI [%.4f, %.4f]\nCI crosses zero: confound NOT confirmed',
               strat$acc_pooled, strat$acc_pooled_ci_low, strat$acc_pooled_ci_high,
               strat$acc_saddle_only, strat$acc_saddle_only_ci_low, strat$acc_saddle_only_ci_high,
               strat$drop, strat$drop_ci_low, strat$drop_ci_high)

p <- ggplot(df, aes(x = feature_set, y = acc, group = 1)) +
  geom_line(linewidth = 0.8, color = darj[2]) +
  geom_point(size = 5, color = darj[2]) +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0.12, linewidth=0.8, color=darj[2]) +
  geom_hline(yintercept = 0.70, linetype='dashed', color=darj[1], linewidth=0.8) +
  geom_hline(yintercept = 0.55, linetype='dotted', color='grey40', linewidth=0.6) +
  geom_hline(yintercept = 0.5, color='black', linewidth=0.3) +
  annotate('text', x = 1.5, y = 0.92, label = ann, family='mono', size=3.4, hjust=0.5) +
  scale_y_continuous(limits=c(0.4, 1.0)) +
  labs(title='Stratification Drop: Pooled vs Saddle-Only Feature Sets',
       x='Feature Set', y='Family Classifier Balanced Accuracy') +
  theme_classic() +
  theme(plot.title = element_text(face='bold', hjust=0.5))

ggsave('stratification.png', p, width=8, height=5.5, dpi=300, units='in')
