library(ggplot2); library(wesanderson); library(jsonlite); library(dplyr)
stats <- fromJSON('stats.json', simplifyVector = FALSE)
darj  <- wes_palette('Darjeeling1', 5)

imp_list <- stats$perm_importance_all
feat <- vapply(imp_list, function(r) r$feature, character(1))
grp  <- vapply(imp_list, function(r) r$group, character(1))
imp  <- vapply(imp_list, function(r) as.numeric(r$importance), numeric(1))
df <- data.frame(feature = feat, group = grp, importance = imp,
                 stringsAsFactors = FALSE)
df <- df[order(-df$importance), ]
df$feature <- factor(df$feature, levels = rev(df$feature))

# Top-10 cutoff
cutoff <- sort(df$importance, decreasing = TRUE)[10]

grp_pal <- c('scalar' = darj[1], 'ccdf_stasis' = darj[2], 'jump_scalar' = darj[3],
             'ccdf_jump' = darj[4], 'corr' = darj[5], 'autocorr' = 'grey40')

n_jump <- stats$n_jump_ccdf_in_top10
ann <- sprintf('Tertiary hypothesis: ≥1 jump-CCDF bin in top-10\nObserved: %d jump-CCDF bins in top-10 — PASS\nBase holdout acc = %.4f',
               n_jump, stats$perm_base_acc_holdout)

p <- ggplot(df, aes(x = importance, y = feature, fill = group)) +
  geom_col() +
  geom_vline(xintercept = cutoff, linetype = 'dashed', color = 'grey40', linewidth = 0.8) +
  scale_fill_manual(values = grp_pal, name = 'Feature Group') +
  annotate('text', x = max(df$importance) * 0.6, y = 3,
           label = ann, size = 3.4, family = 'mono', hjust = 0.5) +
  annotate('text', x = cutoff, y = nrow(df) - 1, label = 'Top-10 cutoff',
           hjust = -0.1, size = 3.2, color = 'grey40') +
  labs(title = 'Permutation Importance: Full Waiting-Time Feature Set',
       x = 'Mean Decrease in Balanced Accuracy',
       y = 'Feature') +
  theme_classic() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        axis.text.y = element_text(size = 6),
        legend.position = 'right')

ggsave('permutation_importance_full.png', p, width = 9, height = 11, dpi = 300, units = 'in')
