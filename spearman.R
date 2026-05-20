library(ggplot2); library(wesanderson); library(jsonlite); library(dplyr)
stats <- fromJSON('stats.json', simplifyVector = TRUE)
darj  <- wes_palette('Darjeeling1', 5)

sp <- stats$spearman_per_family
rows <- list()
for (fam_n in names(sp)) {
  for (st in c('log_mean_t','log_var_t','hill_t','ks_expon_t')) {
    e <- sp[[fam_n]][[st]]
    if (!is.null(e$rho) && !is.na(e$rho)) {
      rows[[length(rows)+1]] <- data.frame(family=fam_n, stat=st, rho=e$rho, p=e$p, n=e$n,
                                           sig=ifelse(!is.na(e$p) && e$p<0.05, 'p<0.05','n.s.'),
                                           stringsAsFactors=FALSE)
    }
  }
}
df <- do.call(rbind, rows)
df$stat <- factor(df$stat, levels=c('log_mean_t','log_var_t','hill_t','ks_expon_t'),
                  labels=c('log-mean','log-var','Hill index','KS-to-exp'))
df$family <- factor(df$family, levels=c('NK','RMF'))

ann <- sprintf('max |rho| = %.4f\nany |rho| >= 0.2 with p < 0.05: %s\nnull (Spearman) falsified: %s',
               stats$spearman_max_abs_rho,
               ifelse(stats$spearman_any_sig,'YES','NO'),
               ifelse(stats$spearman_falsified,'YES','NO'))

p <- ggplot(df, aes(x = stat, y = rho, fill = family, alpha = sig)) +
  geom_col(position = position_dodge(width=0.8), width=0.7) +
  geom_text(aes(label = sprintf('%.3f', rho)), position=position_dodge(width=0.8), vjust=ifelse(df$rho<0, 1.3, -0.5), size=3, family='mono') +
  geom_hline(yintercept = c(-0.2, 0.2), linetype='dashed', color='grey40', linewidth=0.6) +
  geom_hline(yintercept = 0, color='black', linewidth=0.3) +
  scale_fill_manual(values = c(NK=darj[3], RMF=darj[5]), name='Family') +
  scale_alpha_manual(values = c('p<0.05'=1.0, 'n.s.'=0.35), name='Significance') +
  annotate('text', x = 4, y = -0.55, label = ann, family='mono', size=3.4, hjust=1) +
  scale_y_continuous(limits = c(-0.8, 0.6)) +
  labs(title=expression(bold('Spearman '*rho*': Summary Stat vs Ruggedness Parameter')),
       x='Summary Statistic',
       y=expression('Spearman '*rho)) +
  theme_classic() +
  theme(plot.title = element_text(face='bold', hjust=0.5))

ggsave('spearman.png', p, width=9, height=5.5, dpi=300, units='in')
