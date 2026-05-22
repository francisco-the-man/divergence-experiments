library(ggplot2); library(wesanderson); library(jsonlite); library(dplyr)
# Robust JSON load (NaN -> null)
raw <- readLines('result.json', warn = FALSE)
raw <- paste(raw, collapse = '\n')
raw <- gsub('\\bNaN\\b', 'null', raw)
raw <- gsub('\\bInfinity\\b', 'null', raw)
raw <- gsub('\\b-Infinity\\b', 'null', raw)
data <- fromJSON(raw, simplifyVector = FALSE)$results
stats <- fromJSON('stats.json', simplifyVector = TRUE)
darj <- wes_palette('Darjeeling1', 5)

# Build long df
rows <- list()
for (i in seq_along(data)) {
  r <- data[[i]]
  wt <- unlist(r$waiting_times)
  wt <- wt[is.finite(wt) & wt > 0]
  if (length(wt) == 0) next
  rows[[length(rows)+1]] <- data.frame(
    bin = r$param_bin_label, family = r$family,
    log10_wt = log10(wt), stringsAsFactors = FALSE)
}
df <- do.call(rbind, rows)

# Per-bin CCDF
df <- df %>% group_by(bin, family) %>% arrange(log10_wt) %>%
  mutate(ccdf = 1 - (seq_along(log10_wt) - 1) / length(log10_wt)) %>% ungroup()

bin_levels <- c('NK_K2','NK_K4','NK_K8','NK_K16','RMF_theta0.1','RMF_theta0.3','RMF_theta1.0','RMF_theta3.0')
df$bin <- factor(df$bin, levels = bin_levels)

pal8 <- c(darj[1], darj[2], darj[3], darj[4], darj[5], 'goldenrod3', 'steelblue', 'firebrick')
names(pal8) <- bin_levels

p <- ggplot(df, aes(x = log10_wt, y = ccdf, color = bin)) +
  geom_step(linewidth = 0.8) +
  scale_y_log10() +
  scale_color_manual(values = pal8, name = 'Landscape Bin') +
  facet_wrap(~ family, scales = 'free_x') +
  labs(title = 'Waiting-Time CCDF by Landscape Bin',
       x = expression(log[10]*' Waiting Time per Epoch'),
       y = 'CCDF (log-scale)') +
  theme_classic() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        strip.text = element_text(face = 'bold'))

ggsave('ccdf_time_by_class.png', p, width = 10, height = 5, dpi = 300, units = 'in')
