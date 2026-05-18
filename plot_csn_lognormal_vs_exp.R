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

dfr <- df %>% filter(primitives == "rich")
dfr$sig <- ifelse(dfr$p < 0.05,
                  ifelse(dfr$R > 0, "log-normal (p<.05)", "exponential (p<.05)"),
                  "indistinguishable")

for (i in seq_len(nrow(dfr))) {
  cat(sprintf("PRINT csn_%s_R: %.3f  p: %.4f  verdict: %s\n",
              dfr$target_id[i], dfr$R[i], dfr$p[i], dfr$sig[i]))
}
n_lognormal_sig <- sum(dfr$p < 0.05 & dfr$R > 0)
n_exponential_sig <- sum(dfr$p < 0.05 & dfr$R < 0)
cat(sprintf("PRINT csn_n_lognormal_sig: %d / 5\n", n_lognormal_sig))
cat(sprintf("PRINT csn_n_exponential_sig: %d / 5\n", n_exponential_sig))

dfr$target_id <- factor(dfr$target_id, levels = dfr$target_id[order(dfr$R)])

p <- ggplot(dfr, aes(x = target_id, y = R, fill = sig)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("p=%.3f", p)), vjust = ifelse(dfr$R > 0, -0.4, 1.3), size = 3.2) +
  scale_fill_manual(values = c("log-normal (p<.05)" = darj[2],
                                "exponential (p<.05)" = darj[1],
                                "indistinguishable" = "grey80"),
                    name = "Verdict (α=0.05)") +
  labs(x = "Target (rich condition, late phase)",
       y = "Vuong R  (positive ⇒ log-normal preferred,  negative ⇒ exponential preferred)",
       title = "CSN model comparison: log-normal vs exponential",
       subtitle = sprintf("%d / 5 targets significantly prefer log-normal;  %d / 5 significantly prefer exponential",
                          n_lognormal_sig, n_exponential_sig)) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        axis.text.x = element_text(angle = 20, hjust = 1))

ggsave("plot_csn_lognormal_vs_exp.png", p, width = 7, height = 5, dpi = 300, units = "in")
