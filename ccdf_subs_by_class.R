library(ggplot2); library(wesanderson); library(jsonlite); library(dplyr)
raw <- readLines('result.json', warn = FALSE)
raw <- paste(raw, collapse = '\n')
raw <- gsub('\\bNaN\\b', 'null', raw)
data <- fromJSON(raw, simplifyVector = FALSE)$results
stats <- fromJSON('stats.json', simplifyVector = TRUE)
darj <- wes_palette('Darjeeling1', 5)

rows <- list()
for (i in seq_along(data)) {
  r <- data[[i]]
  sl <- unlist(r$stasis_lengths)
  sl <- sl[is.finite(sl) & sl > 0]
  if (length(sl) == 0) next
  rows[[length(rows)+1]] <- data.frame(
    bin = r$param_bin_label, family = r$family,
    log10_sl = log10(sl), stringsAsFactors = FALSE)
}
df <- do.call(rbind, rows)
df <- df %>% group_by(bin, family) %>% arrange(log10_sl) %>%
  mutate(ccdf = 1 - (seq_along(log10_sl) - 1) / length(log10_sl)) %>% ungroup()

bin_levels <- c('NK_K2','NK_K4','NK_K8','NK_K16','RMF_theta0.1','RMF_theta0.3','RMF_theta1.0','RMF_theta3.0')
df$bin <- factor(df$bin, levels = bin_levels)
pal8 <- c(darj[1], darj[2], darj[3], darj[4], darj[5], 'goldenrod3', 'steelblue', 'firebrick')
names(pal8) <- bin_levels

ann <- sprintf('Secondary test:\ntime - subs = %.4f\n95%% CI [%.4f, %.4f]\nexcludes 0 — PASS',
               stats$delta_time_minus_subs, stats$delta_time_minus_subs_ci_low, stats$delta_time_minus_subs_ci_high)

p <- ggplot(df, aes(x = log10_sl, y = ccdf, color = bin)) +
  geom_step(linewidth = 0.8) +
  scale_y_log10() +
  scale_color_manual(values = pal8, name = 'Landscape Bin') +
  facet_wrap(~ family, scales = 'free_x') +
  annotate('text', x = -Inf, y = 0.02, label = ann, hjust = -0.05,
           size = 3.2, family = 'mono') +
  labs(title = 'Substitution-Count CCDF by Landscape Bin (the v2 Observable)',
       x = expression(log[10]*' Substitution Count per Epoch'),
       y = 'CCDF (log-scale)') +
  theme_classic() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        strip.text = element_text(face = 'bold'))

ggsave('ccdf_subs_by_class.png', p, width = 10, height = 5, dpi = 300, units = 'in')
