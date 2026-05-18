library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr)
darj <- wes_palette("Darjeeling1", 5)

data <- fromJSON("result.json", simplifyVector = FALSE)
results <- data$results

df <- bind_rows(lapply(results, function(r) {
  data.frame(
    target_id = r$target_id,
    primitives = r$primitives,
    R = r$csn_late$R,
    p = r$csn_late$p,
    n_tail = r$csn_late$n_tail,
    stringsAsFactors = FALSE
  )
}))

df_rich <- df %>% filter(primitives == "rich") %>%
  mutate(sig = ifelse(p < 0.05, "p < 0.05", "n.s."),
         label = sprintf("R=%.2f\np=%.3f", R, p))

cat(sprintf("PRINT n_targets_csn_sig: %d\n", sum(df_rich$p < 0.05)))
for (i in seq_len(nrow(df_rich))) {
  cat(sprintf("PRINT csn_%s_R: %.3f, p: %.4f\n", df_rich$target_id[i], df_rich$R[i], df_rich$p[i]))
}

p <- ggplot(df_rich, aes(x = target_id, y = R, fill = target_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_col(width = 0.6) +
  geom_text(aes(label = label, vjust = ifelse(R >= 0, -0.2, 1.2)), size = 3) +
  scale_fill_manual(values = darj[1:5], name = "Target") +
  labs(x = "Target",
       y = "CSN log-likelihood ratio R\n(log-normal vs exponential; R>0 favors log-normal)",
       title = "Log-normal vs exponential fit on rich-condition late durations") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "none") +
  expand_limits(y = c(min(df_rich$R) - 1, max(df_rich$R) + 1))

ggsave("plot_csn_R.png", p, width = 7, height = 5, dpi = 300, units = "in")
