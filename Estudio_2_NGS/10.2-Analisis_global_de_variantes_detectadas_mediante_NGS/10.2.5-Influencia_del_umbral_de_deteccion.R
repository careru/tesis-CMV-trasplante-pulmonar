# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.5. Influencia del umbral de detección
# Extracto del script completo (ver README.md del estudio)

# ---- 10.2.5. Influencia del umbral de detección ----

## ---- Figura 48. Variantes únicas según el umbral de frecuencia de variante ----
# ── Paleta ──────────────────────────────────────────────────────────────────
pal_genes <- c("UL97" = "#4A90C4", "UL54" = "#E8A33D", "UL56" = "#D9534F")

# ── Datos ───────────────────────────────────────────────────────────────────
datos_umbral <- tribble(
  ~Gen,   ~Umbral,  ~Tipo,  ~Variantes,
  "UL97", "≥ 5%",   "MRA",   3,
  "UL97", "≥ 10%",  "MRA",   3,
  "UL97", "≥ 15%",  "MRA",   2,
  "UL97", "≥ 20%",  "MRA",   2,
  "UL54", "≥ 5%",   "MRA",   1,
  "UL54", "≥ 10%",  "MRA",   1,
  "UL54", "≥ 15%",  "MRA",   1,
  "UL54", "≥ 20%",  "MRA",   1,
  "UL56", "≥ 5%",   "MRA",   0,
  "UL56", "≥ 10%",  "MRA",   0,
  "UL56", "≥ 15%",  "MRA",   0,
  "UL56", "≥ 20%",  "MRA",   0,
  "UL97", "≥ 5%",   "PN",    9,      # incluye D605E (5,7%)
  "UL97", "≥ 10%",  "PN",    8,
  "UL97", "≥ 15%",  "PN",    8,
  "UL97", "≥ 20%",  "PN",    8,
  "UL54", "≥ 5%",   "PN",   12,
  "UL54", "≥ 10%",  "PN",   11,
  "UL54", "≥ 15%",  "PN",   11,
  "UL54", "≥ 20%",  "PN",   11,
  "UL56", "≥ 5%",   "PN",    3,
  "UL56", "≥ 10%",  "PN",    3,
  "UL56", "≥ 15%",  "PN",    3,
  "UL56", "≥ 20%",  "PN",    3,
  "UL97", "≥ 5%",   "VND",  25,
  "UL97", "≥ 10%",  "VND",   2,
  "UL97", "≥ 15%",  "VND",   2,
  "UL97", "≥ 20%",  "VND",   1,
  "UL54", "≥ 5%",   "VND",  67,
  "UL54", "≥ 10%",  "VND",  17,
  "UL54", "≥ 15%",  "VND",  13,
  "UL54", "≥ 20%",  "VND",   6,
  "UL56", "≥ 5%",   "VND",   9,
  "UL56", "≥ 10%",  "VND",   1,
  "UL56", "≥ 15%",  "VND",   1,
  "UL56", "≥ 20%",  "VND",   1
) %>%
  mutate(
    Umbral = factor(Umbral, levels = c("≥ 20%", "≥ 15%", "≥ 10%", "≥ 5%")),
    Tipo   = factor(Tipo,   levels = c("MRA", "PN", "VND")),
    Gen    = factor(Gen,    levels = c("UL97", "UL54", "UL56"))
  )

# ── Figura ──────────────────────────────────────────────────────────────────
p_umbral <- ggplot(datos_umbral, aes(x = Umbral, y = Variantes, fill = Gen)) +
  geom_col(
    position  = position_dodge(width = 0.8),
    width     = 0.7,
    color     = "white",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = if_else(Variantes > 0, as.character(Variantes), ""),
        color = Gen),
    position = position_dodge(width = 0.8),
    vjust    = -0.6,
    size     = 3.9,
    fontface = "bold",
    show.legend = FALSE
  ) +
  facet_wrap(~ Tipo, scales = "free_y", nrow = 1) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = pal_genes) +
  scale_color_manual(values = pal_genes) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.22))) +
  labs(
    x    = "Umbral de frecuencia de variante",
    y    = "Número de variantes únicas",
    fill = "Gen"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 11) +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = 12),
    axis.title.x       = element_text(face = "bold", size = 11,
                                      margin = margin(t = 8)),
    axis.title.y       = element_text(face = "bold", size = 11,
                                      margin = margin(r = 8)),
    axis.text          = element_text(color = "black", size = 11),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
    panel.spacing      = unit(0.7, "cm"),
    legend.position    = "bottom",
    legend.title       = element_text(face = "bold", size = 11),
    legend.text        = element_text(size = 11),
    legend.margin      = margin(t = 0),
    plot.margin        = margin(4, 6, 4, 4)
  )

p_umbral

# [guardado antiguo] ggsave(
# [guardado antiguo]   filename = "~/Desktop/Figuras_tesis/figura_47_variantes_umbral_NGS.png",
# [guardado antiguo]   plot     = p_umbral,
# [guardado antiguo]   width    = 16,
# [guardado antiguo]   height   = 8,
# [guardado antiguo]   units    = "cm",
# [guardado antiguo]   dpi      = 600
# [guardado antiguo] )

# [Nuevo] Los recuentos de la figura están escritos a mano. Se comprueba que coinciden
# con los datos (variantes únicas por gen, clase y umbral en all_data_txp).
recuento_datos <- expand_grid(Gen = c("UL97", "UL54", "UL56"), thr = c(5, 10, 15, 20),
                              Tipo = c("MRA", "PN", "VND")) %>%
  mutate(Variantes_datos = pmap_int(list(Gen, thr, Tipo), function(g, t, cl) {
    clase <- c(MRA = "MRA", PN = "POL", VND = "NEW")[[cl]]
    n_distinct(filter(all_data_txp, Protein == g, class == clase, Prevalence >= t)$Mutation)
  }),
  Umbral = paste0("≥ ", thr, "%"))
comprobacion_fig48 <- datos_umbral %>%
  mutate(across(c(Gen, Umbral, Tipo), as.character)) %>%
  left_join(recuento_datos, by = c("Gen", "Umbral", "Tipo"))
stopifnot(all(comprobacion_fig48$Variantes == comprobacion_fig48$Variantes_datos))
cat("Figura 48: los recuentos escritos a mano coinciden con los datos.\n")

guardar_fig(p_umbral, "Figura_48", ancho = 16, alto = 8, dpi = 600)


