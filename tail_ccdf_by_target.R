library(ggplot2); library(wesanderson); library(dplyr); library(jsonlite); library(tidyr)

data  <- fromJSON("result.json", simplifyVector = FALSE)
stats <- fromJSON("stats.json", simplifyVector = FALSE)
darj  <- wes_palette("Darjeeling1", 5)

results <- data$results
rows <- list()
for (r in results) {
  durs <- unlist(r$pooled_late_durations_sample)
  if (length(durs) == 0) next
  durs <- sort(durs[durs > 0])
  ccdf <- 1 - (seq_along(durs) - 1) / length(durs)
  rows[[length(rows) + 1]] <- data.frame(
    target_id  = r$target_id,
    primitives = r$primitives,
    duration   = durs,
    ccdf       = ccdf,
    stringsAsFactors = FALSE
  )
}
df <- do.call(rbind, rows)
df$primitives <- factor(df$primitives, levels = c("minimal", "rich"))

# CSN annotations per target (rich)
targets <- unlist(stats$targets)
ann <- data.frame(
  target_id = targets,
  label = sapply(targets, function(t) {
    s <- stats$csn_per_target[[t]]
    sprintf("R=%.2f, p=%.3f\n%s",
            s$R, s$p,
            ifelse(s$lognormal_preferred, "lognormal>exp", "exp>lognormal"))
  }),
  stringsAsFactors = FALSE
)

p <- ggplot(df, aes(x = duration, y = ccdf, color = primitives, linetype = primitives)) +
  geom_step(linewidth = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = darj[c(1, 2)], name = "Primitive set") +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Primitive set") +
  facet_wrap(~ target_id, ncol = 3) +
  geom_label(data = ann, aes(x = 2, y = 0.02, label = label),
             inherit.aes = FALSE, size = 2.8, hjust = 0, fill = "white", label.size = 0.2) +
  labs(title = "Pooled Late-Phase Plateau-Duration CCDFs (rich-condition CSN annotated)",
       x = "Plateau duration (generations, log scale)",
       y = "P(X >= x) (log scale)") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom",
        strip.text = element_text(face = "bold"))

ggsave("tail_ccdf_by_target.png", p, width = 10, height = 7, dpi = 300, units = "in")
