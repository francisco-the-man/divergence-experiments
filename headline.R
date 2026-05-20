library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr)

raw_txt <- readLines('result.json', warn = FALSE)
raw_txt <- paste(raw_txt, collapse='\n')
raw_txt <- gsub('\\bNaN\\b','null',raw_txt)
raw_txt <- gsub('\\b-Infinity\\b','null',raw_txt)
raw_txt <- gsub('\\bInfinity\\b','null',raw_txt)
data  <- fromJSON(raw_txt, simplifyVector = FALSE)$results
stats <- fromJSON('stats.json', simplifyVector = TRUE)
darj  <- wes_palette('Darjeeling1', 5)
zis   <- wes_palette('Zissou1', 5)

# Build per-task frame for headline panels
rows <- lapply(data, function(r) {
  data.frame(
    cell = r$cell,
    family = r$family,
    param_value = as.numeric(r$param_value),
    log_mean_t = ifelse(is.null(r$summary_stats_time$log_mean), NA, as.numeric(r$summary_stats_time$log_mean)),
    log_var_t  = ifelse(is.null(r$summary_stats_time$log_var),  NA, as.numeric(r$summary_stats_time$log_var)),
    hill_t     = ifelse(is.null(r$summary_stats_time$hill),     NA, as.numeric(r$summary_stats_time$hill)),
    ks_expon_t = ifelse(is.null(r$summary_stats_time$ks_expon), NA, as.numeric(r$summary_stats_time$ks_expon)),
    stringsAsFactors = FALSE
  )
})
df <- do.call(rbind, rows)
df <- df[is.finite(df$log_mean_t), ]
df$cell <- factor(df$cell, levels = c('NK_K2','NK_K4','NK_K8','NK_K16','RMF_theta0.1','RMF_theta0.3','RMF_theta1.0','RMF_theta3.0'))

cell_pal <- c('NK_K2'=darj[1],'NK_K4'=darj[2],'NK_K8'=darj[3],'NK_K16'=darj[4],
              'RMF_theta0.1'=zis[1],'RMF_theta0.3'=zis[2],'RMF_theta1.0'=zis[4],'RMF_theta3.0'=zis[5])

# Panel A: log_mean_t by cell (boxplot)
p.a <- ggplot(df, aes(x = cell, y = log_mean_t, fill = cell)) +
  geom_boxplot(outlier.size=0.6, linewidth=0.4) +
  scale_fill_manual(values = cell_pal, guide='none') +
  labs(title = 'Panel A: Log-Mean Stasis vs Cell',
       x = 'Landscape Cell',
       y = expression('Mean of log stasis length (time units)')) +
  theme_classic() +
  theme(plot.title = element_text(face='bold', hjust=0.5),
        axis.text.x = element_text(angle=30, hjust=1, size=8))

# Panel B: family classifier accuracy bar with CI
fam <- stats$family_classifier
cell_c <- stats$cell_classifier
class_df <- data.frame(
  test = c('Family\n(NK vs RMF)','Cell bin\n(8-way)'),
  acc  = c(fam$balanced_accuracy, cell_c$balanced_accuracy),
  lo   = c(fam$ci_low, cell_c$ci_low),
  hi   = c(fam$ci_high, cell_c$ci_high),
  thr  = c(fam$hypothesis_threshold, cell_c$hypothesis_threshold),
  chance = c(0.5, cell_c$chance)
)
ann_a <- sprintf('balanced acc = %.4f\n95%% CI [%.4f, %.4f]', fam$balanced_accuracy, fam$ci_low, fam$ci_high)
ann_b <- sprintf('balanced acc = %.4f\n95%% CI [%.4f, %.4f]\nchance = %.3f', cell_c$balanced_accuracy, cell_c$ci_low, cell_c$ci_high, cell_c$chance)

p.b <- ggplot(class_df, aes(x = test, y = acc)) +
  geom_col(fill = darj[5], width = 0.55) +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0.2, linewidth=0.8) +
  geom_hline(aes(yintercept = thr), linetype='dashed', color=darj[1], linewidth=0.8) +
  geom_hline(aes(yintercept = chance), linetype='dotted', color='grey40', linewidth=0.6) +
  annotate('text', x=1, y=0.95, label=ann_a, family='mono', size=3.2, hjust=0.5) +
  annotate('text', x=2, y=0.95, label=ann_b, family='mono', size=3.2, hjust=0.5) +
  scale_y_continuous(limits=c(0,1.0)) +
  labs(title='Panel B: Classifier Balanced Accuracy',
       x='Preregistered Test', y='Balanced Accuracy (5-Fold CV)') +
  theme_classic() +
  theme(plot.title = element_text(face='bold', hjust=0.5))

# Panel C: Spearman summary - bar chart of |rho| per (family, stat)
sp <- stats$spearman_per_family
sp_df <- data.frame()
for (fam_n in names(sp)) {
  for (st in c('log_mean_t','log_var_t','hill_t','ks_expon_t')) {
    r <- sp[[fam_n]][[st]]$rho
    p <- sp[[fam_n]][[st]]$p
    if (!is.null(r) && !is.na(r)) {
      sp_df <- rbind(sp_df, data.frame(family=fam_n, stat=st, rho=r,
                                       p=p, sig=ifelse(!is.na(p)&&p<0.05,'p<0.05','n.s.'),
                                       stringsAsFactors=FALSE))
    }
  }
}
sp_df$stat <- factor(sp_df$stat, levels=c('log_mean_t','log_var_t','hill_t','ks_expon_t'),
                     labels=c('log-mean','log-var','Hill','KS-exp'))
ann_c <- sprintf('max |rho| = %.4f', stats$spearman_max_abs_rho)
p.c <- ggplot(sp_df, aes(x=stat, y=rho, fill=family, alpha=sig)) +
  geom_col(position=position_dodge(width=0.8), width=0.7) +
  geom_hline(yintercept = c(-0.2, 0.2), linetype='dashed', color='grey40', linewidth=0.6) +
  geom_hline(yintercept = 0, color='black', linewidth=0.3) +
  scale_fill_manual(values = c(NK=darj[3], RMF=darj[5]), name='Family') +
  scale_alpha_manual(values = c('p<0.05'=1.0, 'n.s.'=0.35), name='Significance') +
  annotate('text', x=4, y=-0.55, label=ann_c, family='mono', size=3.2, hjust=1) +
  labs(title='Panel C: Spearman Rho vs Ruggedness Parameter',
       x='Summary Statistic', y=expression('Spearman '*rho*' (stat vs K or '*theta*')')) +
  theme_classic() +
  theme(plot.title = element_text(face='bold', hjust=0.5),
        legend.position='right')

# Panel D: stratification drop
strat <- stats$stratification
strat_df <- data.frame(
  feature_set = c('Pooled','Saddle-Only'),
  acc = c(strat$acc_pooled, strat$acc_saddle_only),
  lo  = c(strat$acc_pooled_ci_low, strat$acc_saddle_only_ci_low),
  hi  = c(strat$acc_pooled_ci_high, strat$acc_saddle_only_ci_high)
)
ann_d <- sprintf('drop = %.4f\nbootstrap CI [%.4f, %.4f]\n(CI crosses 0)', strat$drop, strat$drop_ci_low, strat$drop_ci_high)
p.d <- ggplot(strat_df, aes(x=feature_set, y=acc)) +
  geom_col(fill=darj[2], width=0.5) +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0.18, linewidth=0.8) +
  geom_hline(yintercept=0.70, linetype='dashed', color=darj[1], linewidth=0.8) +
  geom_hline(yintercept=0.55, linetype='dotted', color='grey40', linewidth=0.6) +
  annotate('text', x=1.5, y=0.92, label=ann_d, family='mono', size=3.2, hjust=0.5) +
  scale_y_continuous(limits=c(0,1.0)) +
  labs(title='Panel D: Stratification — Pooled vs Saddle-Only',
       x='Feature Set', y='Family Classifier Balanced Accuracy') +
  theme_classic() +
  theme(plot.title = element_text(face='bold', hjust=0.5))

final <- (p.a + p.b) / (p.c + p.d) + plot_annotation(
  title = 'Stasis-Length Distributions as Landscape Signature: Four Preregistered Tests',
  theme = theme(plot.title = element_text(face='bold', hjust=0.5, size=13))
)
ggsave('headline.png', final, width=12, height=9, dpi=300, units='in')
