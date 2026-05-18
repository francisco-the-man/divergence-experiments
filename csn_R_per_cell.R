library(ggplot2); library(wesanderson); library(patchwork); library(dplyr); library(jsonlite); library(tidyr)
stats <- fromJSON("stats.json", simplifyVector = FALSE)
darj  <- wes_palette("Darjeeling1", 5)

cells <- stats$csn_by_cell
rows <- lapply(cells, function(cc) {
  data.frame(
    target = cc$target,
    condition = cc$condition,
    R = if (is.null(cc$R)) NA_real_ else cc$R,
    p = if (is.null(cc$p)) NA_real_ else cc$p,
    n_tail = if (is.null(cc$n_tail)) NA_integer_ else cc$n_tail,
    stringsAsFactors = FALSE
  )
})
df <- do.call(rbind, rows)
df$condition <- factor(df$condition, levels = c("M", "S", "C"))
df$sig <- ifelse(!is.na(df$p) & df$p < 0.05, "*", "")
df$R_plot <- ifelse(is.na(df$R), 0, df$R)

ann <- sprintf("C-cells log-normal & p<0.05: %d/%d\nS-cells exponential & p<0.05: %d/%d\nPre-registered bar: >=4/6 C AND >=3/6 S",
               as.integer(stats$csn_C_lognormal_significant), as.integer(stats$csn_n_C_cells),
               as.integer(stats$csn_S_exponential_significant), as.integer(stats$csn_n_S_cells))

p <- ggplot(df, aes(x = target, y = R_plot, fill = target)) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_text(aes(label = sig, y = R_plot + sign(R_plot) * 2), size = 6, vjust = 0.5) +
  facet_wrap(~ condition, scales = "free_y") +
  scale_fill_manual(values = darj[1:6], name = "Target Equation", guide = "none") +
  annotate("text", x = 1, y = Inf, label = ann, hjust = 0, vjust = 1.2,
           size = 3.2, family = "mono") +
  labs(title = "CSN Tail Shape: Log-Normal vs Exponential Likelihood Ratio per Cell",
       x = "Target Equation",
       y = expression("R = log L"[lognormal] * " - log L"[exponential] * "  (positive = log-normal preferred)")) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "grey90"))

ggsave("csn_R_per_cell.png", p, width = 10, height = 5, dpi = 300, units = "in")
