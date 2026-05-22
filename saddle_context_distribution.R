library(ggplot2); library(wesanderson); library(jsonlite); library(dplyr); library(tidyr)
stats <- fromJSON('stats.json', simplifyVector = FALSE)
darj <- wes_palette('Darjeeling1', 5)

cf <- stats$context_fractions
df <- do.call(rbind, lapply(cf, function(r) {
  data.frame(bin = r$bin,
             saddle = as.numeric(r$saddle),
             plateau = as.numeric(r$plateau),
             peak = as.numeric(r$peak),
             mixed = as.numeric(r$mixed),
             n_epochs = as.integer(r$n_epochs),
             stringsAsFactors = FALSE)
}))
bin_levels <- c('NK_K2','NK_K4','NK_K8','NK_K16','RMF_theta0.1','RMF_theta0.3','RMF_theta1.0','RMF_theta3.0')
df$bin <- factor(df$bin, levels = bin_levels)

df_long <- pivot_longer(df, cols = c('saddle','plateau','peak','mixed'),
                        names_to = 'context', values_to = 'fraction')
df_long$context <- factor(df_long$context, levels = c('saddle','plateau','peak','mixed'))

ctx_pal <- c('saddle' = darj[1], 'plateau' = darj[2], 'peak' = darj[3], 'mixed' = darj[5])

chi_p_val <- as.numeric(stats$chi_square_p)
chi_stat_val <- as.numeric(stats$chi_square_stat)
ann <- sprintf('Chi-square: stat = %.2f, p = %.4f\nAny category > 90%% in every bin? %s',
               chi_stat_val, chi_p_val,
               ifelse(isTRUE(stats$context_dominant_in_every_bin), 'YES', 'no'))

p <- ggplot(df_long, aes(x = bin, y = fraction, fill = context)) +
  geom_col() +
  scale_fill_manual(values = ctx_pal, name = 'Entry Context') +
  annotate('text', x = 4.5, y = 1.08, label = ann, size = 3.4, family = 'mono', hjust = 0.5) +
  coord_cartesian(ylim = c(0, 1.13), clip = 'off') +
  labs(title = 'Saddle-Aware Context Fractions per Landscape Bin (v3 Rule)',
       x = 'Landscape Bin', y = 'Fraction of Epochs') +
  theme_classic() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        axis.text.x = element_text(angle = 35, hjust = 1),
        plot.margin = margin(20, 10, 10, 10))

ggsave('saddle_context_distribution.png', p, width = 9, height = 5.5, dpi = 300, units = 'in')
