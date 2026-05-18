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

df_rich <- df %>% filter(primitives == "rich") %>%
  mutate(verdict = case_when(
    p < 0.05 & R > 0 ~ "Favors log-normal",
    p < 0.05 & R < 0 ~ "Favors exponential",
    TRUE ~ "Indistinguishable"
  ),
  label = sprintf("R=%.2f\np=%.3f", R, p))

df_rich$verdict <- factor(df_rich$verdict,
  levels = c("Favors log-normal", "Indistinguishable", "Favors exponential"))

p <- ggplot(df_rich, aes(x = target_id, y = R, fill = verdict)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "black") +
  geom_text(aes(label = label, y = R + ifelse(R >= 0, 0.4, -0.4)), size = 3) +
  scale_fill_manual(values = c("Favors log-normal" = darj[2],
                                "Indistinguishable" = "grey75",
                                "Favors exponential" = darj[1]),
                    name = "CSN verdict (α=0.05)", drop = FALSE) +
  labs(x = "Target equation",
       y = "Log-likelihood ratio R (log-normal vs. exponential)",
       title = "CSN model comparison on rich-condition late-phase tails\n(positive R favors log-normal; mechanism requires all targets significantly positive)") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom")

ggsave("plot_csn_lr_by_target.png", p, width = 8, height = 5.5, dpi = 300, units = "in")
