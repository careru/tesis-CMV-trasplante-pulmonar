# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.22. Anexo: Artículo: Figura P75 en blanco y negro
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Artículo: Figura P75 en blanco y negro ----
# --- Etiquetas y grupos en inglés ---
etiq <- det_largo %>%
  group_by(estabilidad) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = paste0(ifelse(estabilidad == "estable",
                               "Stable genotypes", "Fluctuating genotypes"),
                        "\n(n = ", n, ")"))

det_fig <- det_largo %>%
  mutate(group = factor(ifelse(estabilidad == "estable",
                               "Stable genotypes", "Fluctuating genotypes"),
                        levels = c("Stable genotypes", "Fluctuating genotypes")))

etiq_x <- etiq %>%
  mutate(group = factor(ifelse(estabilidad == "estable",
                               "Stable genotypes", "Fluctuating genotypes"),
                        levels = c("Stable genotypes", "Fluctuating genotypes")))

# --- Figura (blanco y negro) ---
figura_R22 <- ggplot(det_fig, aes(x = group, y = dist_p75)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "black", linewidth = 0.5) +
  annotate("text", x = 2.45, y = 0.9, label = "P75",
           size = 4, color = "black") +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.25,
              color = "grey35", show.legend = FALSE) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, width = 0.55,
               color = "black", linewidth = 0.6, alpha = 0.85) +
  scale_fill_manual(values = c("Stable genotypes"      = "grey55",
                               "Fluctuating genotypes"  = "grey85")) +
  # corchete de significación
  annotate("segment", x = 1, xend = 2, y = 9,   yend = 9,   linewidth = 0.5) +
  annotate("segment", x = 1, xend = 1, y = 8.6, yend = 9,   linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2, y = 8.6, yend = 9,   linewidth = 0.5) +
  annotate("text", x = 1.5, y = 9.8, label = "italic(p)~'< 0.001'",
           parse = TRUE, size = 4.2, color = "black") +
  scale_x_discrete(labels = setNames(etiq_x$label, etiq_x$group)) +
  scale_y_continuous(name = "Distance to genotype-specific P75 (Ct)",
                     breaks = seq(-6, 10, by = 2)) +
  coord_cartesian(ylim = c(-7, 11)) +
  labs(x = NULL) +
  theme_classic(base_size = 13, base_family = "sans") +
  theme(
    legend.position = "none",
    axis.line    = element_line(color = "black", linewidth = 0.5),
    axis.ticks   = element_line(color = "black", linewidth = 0.4),
    axis.text.x  = element_text(size = 12, color = "black"),
    axis.text.y  = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    plot.margin  = margin(8, 8, 8, 8)
  )

figura_R22

# Export 300 dpi para JID
ggsave("Figure_P75_stable_vs_fluctuating.tiff", figura_R22,
       width = 4.6, height = 4.8, dpi = 300, compression = "lzw")

