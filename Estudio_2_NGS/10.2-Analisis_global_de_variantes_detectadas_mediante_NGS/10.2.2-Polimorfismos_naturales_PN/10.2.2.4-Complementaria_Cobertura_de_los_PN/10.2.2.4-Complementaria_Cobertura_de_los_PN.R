# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.2.4. Complementaria. Cobertura de los PN
# Extracto del script completo (ver README.md del estudio)

## ---- Complementaria. Cobertura de los PN ----
# ── Datos ─────────────────────────────────────────────────────────────────────
pol_cob <- all_data_txp %>%
  filter(class == "POL", Prevalence >= 5) %>%
  mutate(
    Position = as.integer(str_extract(Mutation, "\\d+")),
    Protein  = factor(Protein, levels = c("UL97", "UL54", "UL56"))
  ) %>%
  group_by(Protein, Mutation, Position) %>%
  summarise(
    mediana  = median(cobertura_total, na.rm = TRUE),
    p25      = quantile(cobertura_total, 0.25, na.rm = TRUE),
    p75      = quantile(cobertura_total, 0.75, na.rm = TRUE),
    n        = n_distinct(Npet),
    prev_max = max(Prevalence, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  arrange(Protein, Position) %>%
  mutate(
    Mutation  = factor(Mutation, levels = unique(Mutation)),
    unica     = n == 1,
    baja_prev = prev_max < 20
  )

# ── Colores por gen ───────────────────────────────────────────────────────────
col_gen <- c("UL97" = "#4472C4", "UL54" = "#C8880A", "UL56" = "#C0392B")

# ── Shapes: rombo (18) si prev_max < 20%, círculo (16/21) si no ──────────────
# Combinamos unica y baja_prev en un indicador de shape
pol_cob <- pol_cob %>%
  mutate(
    shape_tipo = case_when(
      baja_prev & unica  ~ "baja_unica",
      baja_prev & !unica ~ "baja_multi",
      !baja_prev & unica ~ "norm_unica",
      TRUE               ~ "norm_multi"
    )
  )

shape_map <- c(
  "norm_multi" = 16,   # círculo relleno
  "norm_unica" = 21,   # círculo hueco
  "baja_multi" = 18,   # rombo relleno
  "baja_unica" = 23    # rombo hueco
)

# ── Plot ──────────────────────────────────────────────────────────────────────
p_cob_pol <- ggplot(pol_cob,
                    aes(x = Mutation, y = mediana, color = Protein)) +
  
  # RIC (solo si n > 1)
  geom_errorbar(
    data = pol_cob %>% filter(!unica),
    aes(ymin = p25, ymax = p75),
    width = 0.35, linewidth = 0.6, alpha = 0.7
  ) +
  
  # Puntos con shape según tipo
  geom_point(
    aes(shape = shape_tipo,
        fill  = ifelse(unica | shape_tipo == "baja_unica", "white", Protein)),
    size = 3.2, stroke = 0.8
  ) +
  scale_shape_manual(
    values = shape_map,
    labels = c(
      "norm_multi" = "Frecuencia variante \u226520%",
      "norm_unica" = "Frecuencia variante \u226520%, muestra \u00fanica",
      "baja_multi" = "Frecuencia variante <20%",
      "baja_unica" = "Frecuencia variante <20%, muestra \u00fanica"
    ),
    name = NULL
  ) +
  scale_fill_manual(values = c(col_gen, "white" = "white"), guide = "none") +
  scale_color_manual(values = col_gen, guide = "none") +
  
  scale_y_continuous(
    labels  = function(x) format(x, big.mark = ".", decimal.mark = ","),
    expand  = expansion(mult = c(0.05, 0.12))
  ) +
  
  facet_wrap(~ Protein, scales = "free_x", nrow = 1) +
  
  labs(
    x       = NULL,
    y       = "Cobertura de posici\u00f3n (x)",
    caption = paste0(
    )
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1,
                                    size = 12, face = "bold"),
    axis.text.y      = element_text(size = 10),
    strip.text       = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "#EEF4FB", color = NA),
    panel.spacing    = unit(1.2, "cm"),
    legend.position  = "bottom",
    legend.text      = element_text(size = 10),
    plot.caption     = element_text(size = 10, color = "grey40",
                                    hjust = 0, margin = margin(t = 8)),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin      = margin(t = 10, r = 15, b = 10, l = 10)
  )

p_cob_pol

# [guardado antiguo] ggsave(
# [guardado antiguo]   "~/Desktop/Figuras_tesis/figura_cobertura_PN.png",
# [guardado antiguo]   p_cob_pol,
# [guardado antiguo]   width = 32, height = 14, units = "cm", dpi = 400,
# [guardado antiguo]   bg = "white"
# [guardado antiguo] )

guardar_fig(p_cob_pol, "Complementaria_cobertura_PN", ancho = 32, alto = 14, dpi = 400)


