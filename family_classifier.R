library(ggplot2); library(wesanderson); library(jsonlite)
stats <- fromJSON('stats.json', simplifyVector = TRUE)
darj  <- wes_palette('Darjeeling1', 5)

fam <- stats$family_classifier
df <- data.frame(
  metric = c('Null upper bound\n(disconfirms)','Hypothesis threshold\n(confirms)','Observed lower CI','Observed point','Observed upper CI'),
  value  = c(fam$null_threshold, fam$hypothesis_threshold, fam$ci_low, fam$balanced_accuracy, fam$ci_high),
  kind   = c('threshold','threshold','obs','obs','obs')
)
df$metric <- factor(df$metric, levels=df$metric)
ann <- sprintf('Family classifier (NK vs RMF)\nbalanced acc = %.4f\n95%% CI [%.4f, %.4f]\nn = %d\nHypothesis threshold: %.2f\nNull threshold: %.2f', fam$balanced_accuracy, fam$ci_low, fam$ci_high, fam$n, fam$hypothesis_threshold, fam$null_threshold)

p <- ggplot(df, aes(x = metric, y = value, fill = kind)) +
  geom_col(width = 0.55) +
  geom_hline(yintercept = fam$hypothesis_threshold, linetype='dashed', color=darj[1], linewidth=0.9) +
  geom_hline(yintercept = fam$null_threshold, linetype='dotted', color='grey40', linewidth=0.7) +
  geom_hline(yintercept = 0.5, color='black', linewidth=0.3) +
  scale_fill_manual(values = c(threshold=darj[4], obs=darj[3]), name='Quantity') +
  annotate('text', x = 1, y = 0.95, label = ann, family='mono', size=3.4, hjust=0) +
  scale_y_continuous(limits=c(0,1)) +
  labs(title='Family Classifier: NK vs RMF Balanced Accuracy',
       x='Quantity', y='Balanced Accuracy') +
  theme_classic() +
  theme(plot.title = element_text(face='bold', hjust=0.5),
        axis.text.x = element_text(angle=20, hjust=1, size=8))

ggsave('family_classifier.png', p, width=8, height=5.5, dpi=300, units='in')
