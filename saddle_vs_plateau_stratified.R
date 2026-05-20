library(ggplot2); library(wesanderson); library(jsonlite); library(dplyr); library(tidyr)
data  <- fromJSON("result.json", simplifyVector = FALSE)
stats <- fromJSON("stats.json", simplifyVector = TRUE)
darj  <- wes_palette("Darjeeling1", 5)

rows <- list()
for (r in data$results) {
  ec <- as.character(unlist(r$entry_contexts))
  if (length(ec) == 0) next
  for (c in unique(ec)) {
    rows[[length(rows)+1]] <- data.frame(cell = r$cell, family = r$family, context = c, count = sum(ec == c))
  }
}
df <- do.call(rbind, rows)

agg <- df %>% group_by(cell, context) %>% summarize(n = sum(count), .groups = "drop") %>%
  group_by(cell) %>% mutate(frac = n / sum(n)) %>% ungroup()
agg$cell <- factor(agg$cell,
  levels = c("NK_K2","NK_K4","NK_K8","NK_K16","RMF_theta0.1","RMF_theta0.3","RMF_theta1.0","RMF_theta3.0"))
agg$context <- factor(agg$context, levels = c("saddle","plateau","peak"))

ann <- sprintf("Family acc pooled = %.3f, saddle-only = %.3f, drop = %.3f, 95%% CI [%.3f, %.3f]\nStratification confound check has no power: ~100%% saddle entries",
  stats$strat_acc_pooled, stats$strat_acc_saddle, stats$strat_drop, stats$strat_drop_ci_low, stats$strat_drop_ci_high)

p <- ggplot(agg, aes(x = cell, y = frac, fill = context)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = darj[c(1,3,2)], name = "Stasis Entry Context") +
  annotate("text", x = 0.6, y = -0.12, label = ann, size = 3.2, family = "mono", hjust = 0) +
  scale_y_continuous(limits = c(-0.18, 1.05), breaks = seq(0, 1, 0.25)) +
  labs(title = "Stasis Entry Context Composition by Cell",
       x = "Landscape Cell", y = "Fraction of Stasis Entries") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("saddle_vs_plateau_stratified.png", p, width = 9, height = 5.5, dpi = 300, units = "in")
