library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr)
darj <- wes_palette("Darjeeling1", 5)

data <- fromJSON("result.json", simplifyVector = FALSE)
res <- data$results

df <- bind_rows(lapply(res, function(r) {
  data.frame(
    target_id = r$target_id,
    primitives = r$primitives,
    R = r$csn_late$R,
    p = r$csn_late$p,
    n_tail = r$csn_late$n_tail,
    stringsAsFactors = FALSE
  )
}))
df$primitives <- factor(df$primitives, levels = c("minimal", "rich"))
df$sig <- ifelse(df$p < 0.05, "p < 0.05", "n.s.")

p <- ggplot(df, aes(x = target_id, y = R, fill = primitives, alpha = sig)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "black") +
  geom_text(aes(label = sprintf("p=%.3f", p)),
            position = position_dodge(width = 0.7),
            vjust = ifelse(df$R >= 0, -0.4, 1.2), size = 3) +
  scale_fill_manual(values = darj[c(1, 2)], name = "Primitive set") +
  scale_alpha_manual(values = c("p < 0.05" = 1.0, "n.s." = 0.4), name = "Significance") +
  labs(x = "Target equation",
       y = "CSN log-likelihood ratio R (log-normal vs exponential)",
       title = "Log-normal vs. exponential fit on pooled late plateau durations",
       subtitle = "Positive R favours log-normal; null criterion fires only if every rich cell is n.s.") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("plot_csn_lr.png", p, width = 9, height = 5.5, dpi = 300, units = "in")
