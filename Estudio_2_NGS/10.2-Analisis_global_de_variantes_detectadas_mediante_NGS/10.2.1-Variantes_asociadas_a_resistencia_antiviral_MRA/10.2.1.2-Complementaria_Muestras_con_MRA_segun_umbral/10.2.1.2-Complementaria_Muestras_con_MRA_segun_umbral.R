# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.1.2. Complementaria. Muestras con MRA según umbral
# Extracto del script completo (ver README.md del estudio)

## ---- Complementaria. Muestras con MRA según umbral ----
# ── Calcular mra_por_umbral ───────────────────────────────────────────────────
thresholds <- c(5, 10, 15, 20)

mra_por_umbral <- map_dfr(thresholds, function(thr) {
  all_data_txp %>%
    filter(class == "MRA", Prevalence >= thr) %>%
    group_by(Protein) %>%
    summarise(n_muestras = n_distinct(Npet), .groups = "drop") %>%
    mutate(umbral = paste0("≥", thr, "%"))
}) %>%
  filter(Protein %in% c("UL97", "UL54")) %>%
  mutate(
    umbral  = factor(umbral, levels = paste0("≥", thresholds, "%")),
    Protein = factor(Protein, levels = c("UL97", "UL54"))
  )

# ── Parámetros ────────────────────────────────────────────────────────────────
colores <- c("UL97" = "#4472C4", "UL54" = "#FFC000")

ref_sanger <- tibble(
  Protein     = c("UL97", "UL54"),
  ref_value   = c(8, 6),
  color_linea = c("#104E8B", "darkgoldenrod3")
)

mra_por_umbral <- mra_por_umbral %>%
  mutate(vjust_label = if_else(Protein == "UL54", 1.8, -1))

# ── Figura ────────────────────────────────────────────────────────────────────
fig_mra_umbrales <- ggplot(mra_por_umbral, aes(x = umbral, y = n_muestras,
                           color = Protein, group = Protein)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  geom_text(
    aes(label = n_muestras, vjust = vjust_label),
    size     = 13 / .pt,
    fontface = "bold",
    color    = "black",
    show.legend = FALSE
  ) +
  geom_hline(
    data = ref_sanger,
    aes(yintercept = ref_value),
    color     = alpha(ref_sanger$color_linea, 0.5),
    linetype  = "dashed",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  scale_color_manual(values = colores, name = NULL) +
  scale_y_continuous(
    limits = c(0, 13),
    breaks = seq(0, 12, by = 3),
    expand = c(0, 0)
  ) +
  labs(
    x = "Umbral de prevalencia",
    y = "N muestras ≥1 MRA"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position    = "top",
    legend.text        = element_text(face = "bold", size = 13),
    axis.title         = element_text(face = "bold", size = 13),
    axis.text          = element_text(color = "black", size = 13),
    axis.title.x       = element_text(face = "bold", size = 13, margin = margin(t = 10)),
    axis.title.y       = element_text(face = "bold", size = 13),
    panel.grid.major.y = element_line(color = "#EEEEEE", linewidth = 0.4)
  )

# [guardado antiguo] ggsave(
# [guardado antiguo]   filename = "~/Desktop/Figuras_tesis/figura_MRA_umbrales.png",
# [guardado antiguo]   plot     = last_plot(),
# [guardado antiguo]   width    = 10,
# [guardado antiguo]   height   = 10,
# [guardado antiguo]   units    = "cm",
# [guardado antiguo]   dpi      = 300
# [guardado antiguo] )

guardar_fig(fig_mra_umbrales, "Complementaria_MRA_umbrales", ancho = 10, alto = 10, dpi = 300)


