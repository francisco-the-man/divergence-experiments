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
  jm <- unlist(r$jump_magnitudes)
  jm <- jm[is.finite(jm) & jm > 0]
  if (length(jm) == 0) next
  rows[[length(rows)+1]] <- data.frame(
    bin = r$param_bin_label, family = r$family,
    log_jm = log(jm), stringsAsFactors = FALSE)
}
df <- do.call(rbind, rows)
df <- df %>% group_by(bin, family) %>% arrange(log_jm) %>%
  mutate(ccdf = 1 - (seq_along(log_jm) - 1) / length(log_jm)) %>% ungroup()

bin_levels <- c('NK_K2','NK_K4','NK_K8','NK_K16','RMF_theta0.1','RMF_theta0.3','RMF_theta1.0','RMF_theta3.0')
df$bin <- factor(df$bin, levels = bin_levels)
pal8 <- c(darj[1], darj[2], darj[3], darj[4], darj[5], 'goldenrod3', 'steelblue', 'firebrick')
names(pal8) <- bin_levels

nk_top <- stats$top_spearman_per_family$NK
rmf_top <- stats$top_spearman_per_family$RMF
ann <- sprintf('Top Spearman (BH-corrected):\nNK: %s, rho=%.3f, p=%.2e\nRMF: %s, rho=%.3f, p=%.2e',
               nk_top$feature, nk_top$rho, nk_top$p_bh,
               rmf_top$feature, rmf_top$rho, rmf_top$p_bh)

p <- ggplot(df, aes(x = log_jm, y = ccdf, color = bin)) +
  geom_step(linewidth = 0.8) +
  scale_y_log10() +
  scale_color_manual(values = pal8, name = 'Landscape Bin') +
  facet_wrap(~ family, scales = 'free_x') +
  annotate('text', x = -Inf, y = 0.02, label = ann, hjust = -0.05,
           size = 3.0, family = 'mono') +
  labs(title = 'Jump-Magnitude CCDF by Landscape Bin',
       x = expression('log |'*Delta*'f| (per Fixation)'),
       y = 'CCDF (log-scale)') +
  theme_classic() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        strip.text = element_text(face = 'bold'))

ggsave('ccdf_jump_by_class.png', p, width = 10, height = 5, dpi = 300, units = 'in')
