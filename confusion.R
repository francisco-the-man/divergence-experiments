library(ggplot2); library(wesanderson); library(jsonlite); library(dplyr)
stats <- fromJSON('stats.json', simplifyVector = TRUE)
darj  <- wes_palette('Darjeeling1', 5)

conf <- stats$confusion
conf$true <- factor(conf$true, levels=c('NK_K2','NK_K4','NK_K8','NK_K16','RMF_theta0.1','RMF_theta0.3','RMF_theta1.0','RMF_theta3.0'))
conf$predicted <- factor(conf$predicted, levels=c('NK_K2','NK_K4','NK_K8','NK_K16','RMF_theta0.1','RMF_theta0.3','RMF_theta1.0','RMF_theta3.0'))
# Row-normalize
conf <- conf %>% group_by(true) %>% mutate(row_total = sum(count), frac = ifelse(row_total>0, count/row_total, 0)) %>% ungroup()

cc <- stats$cell_classifier
ann <- sprintf('8-Way Cell Classifier\nbalanced acc = %.4f\n95%% CI [%.4f, %.4f]\nchance = %.3f\nhypothesis threshold = %.2f', cc$balanced_accuracy, cc$ci_low, cc$ci_high, cc$chance, cc$hypothesis_threshold)

p <- ggplot(conf, aes(x = predicted, y = true, fill = frac)) +
  geom_tile(color='white', linewidth=0.4) +
  geom_text(aes(label = ifelse(count>0, count, '')), size=3, family='mono') +
  scale_fill_gradient(low='white', high=darj[5], name='Row\nFraction', limits=c(0,1)) +
  annotate('text', x = 1, y = 8.7, label = ann, family='mono', size=3.2, hjust=0) +
  labs(title='Confusion Matrix: 8-Way Cell Classifier',
       x='Predicted Class', y='True Class') +
  theme_classic() +
  theme(plot.title = element_text(face='bold', hjust=0.5),
        axis.text.x = element_text(angle=30, hjust=1, size=8),
        axis.text.y = element_text(size=8))

ggsave('confusion.png', p, width=8.5, height=7, dpi=300, units='in')
