library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr)
stats <- fromJSON("stats.json", simplifyVector = TRUE)
darj  <- wes_palette("Darjeeling1", 5)

targets <- stats$complete_targets
dSM <- unlist(stats$delta_sigma_SM_by_target)[targets]
dCM <- unlist(stats$delta_sigma_CM_by_target)[targets]

df <- data.frame(
  target_id = rep(targets, 2),
  contrast  = factor(rep(c("Delta sigma (S - M): size-only effect",
                           "Delta sigma (C - M): connectivity-rich effect"),
                         each = length(targets)),
                     levels = c("Delta sigma (S - M): size-only effect",
                                "Delta sigma (C - M): connectivity-rich effect")),
  delta = c(dSM, dCM)
)

ann <- sprintf("Wilcoxon sigma_C > sigma_S\nV = %d, p = %.4f\nMedian Delta sigma_CS = %.3f\n\nCriterion (ii) Spearman: not evaluable\n(egraph_measure tasks absent)",
               as.integer(stats$wilcoxon_V), stats$wilcoxon_p,
               stats$wilcoxon_median_delta_CS)

p <- ggplot(df, aes(x = target_id, y = delta, fill = target_id)) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  facet_wrap(~ contrast, ncol = 2) +
  scale_fill_manual(values = darj[1:length(targets)], name = "Target Equation", guide = "none") +
  annotate("text", x = 1, y = max(df$delta) * 0.95, label = ann,
           size = 3.2, family = "mono", hjust = 0, vjust = 1) +
  labs(title = "Decomposition of Delta Sigma: Size-Only vs Connectivity-Rich",
       x = "Target Equation",
       y = expression("Delta " * sigma * " (condition - M)")) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"))

ggsave("delta_sigma_decomposition.png", p, width = 10, height = 5, dpi = 300, units = "in")
