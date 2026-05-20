library(ggplot2); library(wesanderson); library(jsonlite); library(dplyr)
stats <- fromJSON("stats.json", simplifyVector = TRUE)
darj  <- wes_palette("Darjeeling1", 5)

sp <- stats$spearman_table
sp$sig <- !is.na(sp$p_bh) & sp$p_bh < 0.05 & abs(sp$rho) >= 0.2
sp$family <- factor(sp$family, levels = c("NK","RMF"))
# Order features by max |rho|
fo <- sp %>% group_by(feature) %>% summarize(m = max(abs(rho), na.rm = TRUE)) %>% arrange(desc(m))
sp$feature <- factor(sp$feature, levels = rev(fo$feature))

ann <- sprintf("Top: %s (%s) rho = %.3f, p_BH = %.1e\nFeatures with |rho|>=0.2 & p_BH<0.05: %d",
  stats$spearman_top_feature, stats$spearman_top_family, stats$spearman_top_rho, stats$spearman_top_p_bh, stats$spearman_n_sig_strong)

p <- ggplot(sp, aes(x = rho, y = feature, color = family, shape = sig)) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.6) +
  geom_vline(xintercept = c(-0.2, 0.2), color = "grey25", linewidth = 0.6, linetype = "dashed") +
  geom_point(size = 3) +
  scale_color_manual(values = darj[c(1,2)], name = "Family") +
  scale_shape_manual(values = c(`FALSE`=1, `TRUE`=19), name = "|rho|>=0.2 and p_BH<0.05") +
  annotate("text", x = -0.95, y = 5, label = ann, size = 3.4, family = "mono", hjust = 0) +
  labs(title = "Spearman Correlation: Feature vs Ruggedness Parameter",
       x = expression("Spearman "*rho*" vs ruggedness parameter (K for NK, "*theta*" for RMF)"),
       y = "Feature") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.y = element_text(size = 7))

ggsave("spearman_features.png", p, width = 9, height = 7.5, dpi = 300, units = "in")
