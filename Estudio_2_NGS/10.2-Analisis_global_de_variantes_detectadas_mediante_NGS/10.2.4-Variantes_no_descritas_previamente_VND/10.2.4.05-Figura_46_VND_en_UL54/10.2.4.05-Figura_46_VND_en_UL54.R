# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.4.05. Figura 46. VND en UL54
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 46. VND en UL54 ----
# ── Dominios UL54 ─────────────────────────────────────────────────────────────
dominios_ul54 <- data.frame(
  start = c(290, 696, 1150),
  end   = c(550, 990, 1243),
  label = c("Exonucleasa", "Polimerasa", "NLS + unión ppUL44"),
  fill  = c("#82E0AA", "#7FB3D3", "#E8A0D0")
)

# ── Hotspots UL54 ─────────────────────────────────────────────────────────────
hotspots_ul54 <- data.frame(
  start = c(300, 408, 498, 985,   598, 802,   700,   750, 781, 945),
  end   = c(325, 420, 545, 990,   600, 825,   715,   760, 790, 955),
  fill  = c(rep("#E67E22", 4),
            rep("#8E44AD", 2),
            "#E74C3C",
            rep("#C0392B", 3))
)

# ── Datos UL54 con umbral_max ─────────────────────────────────────────────────
ul54_umbrales <- all_data_txp %>%
  filter(class == "NEW", Protein == "UL54") %>%
  group_by(Mutation, Position) %>%
  summarise(
    prev_max   = max(Prevalence, na.rm = TRUE),
    umbral_max = case_when(
      max(Prevalence) >= 20 ~ 20,
      max(Prevalence) >= 15 ~ 15,
      max(Prevalence) >= 10 ~ 10,
      TRUE                  ~ 5
    ),
    .groups = "drop"
  ) %>%
  mutate(Position = as.numeric(Position))

ul54 <- new_vars_dom %>%
  filter(Protein == "UL54") %>%
  left_join(ul54_umbrales %>% select(Mutation, umbral_max), by = "Mutation")

# ── Regiones de resistencia ───────────────────────────────────────────────────
en_region_resistencia <- function(pos) {
  (pos >= 300 & pos <= 325) |
    (pos >= 408 & pos <= 420) |
    (pos >= 498 & pos <= 545) |
    (pos >= 598 & pos <= 600) |
    (pos >= 700 & pos <= 715) |
    (pos >= 750 & pos <= 760) |
    (pos >= 781 & pos <= 790) |
    (pos >= 802 & pos <= 825) |
    (pos >= 945 & pos <= 955) |
    (pos >= 985 & pos <= 990)
}

# ── Niveles fijos ─────────────────────────────────────────────────────────────
nivel <- c("5" = 0.6, "10" = 1.1, "15" = 1.6, "20" = 2.1)

df <- ul54 %>%
  mutate(
    en_resist = en_region_resistencia(as.numeric(Position)),
    y_nivel   = nivel[as.character(umbral_max)],
    y_point   = ifelse(en_resist, -y_nivel, y_nivel),
    color     = ifelse(en_resist, "#C0392B", "#2980B9"),
    label_mut = case_when(
      en_resist                         ~ Mutation,
      prev_max >= 10 | n_pacientes >= 4 ~ Mutation,
      TRUE                              ~ NA_character_
    )
  )

# ── Etiquetas con agrupación por proximidad ───────────────────────────────────
df_lab_azul <- df %>%
  filter(!is.na(label_mut), !en_resist) %>%
  arrange(Position) %>%
  mutate(prox_group = cumsum(c(TRUE, diff(Position) > 40))) %>%
  group_by(prox_group) %>%
  mutate(
    within_rank = row_number(),
    y_lab       = 3.0 + (within_rank - 1) * 0.60
  ) %>%
  ungroup()

df_lab_rojo <- df %>%
  filter(!is.na(label_mut), en_resist) %>%
  arrange(Position) %>%
  mutate(
    escalon = case_when(
      Position < 817  ~ (row_number() - 1) %% 3,
      TRUE            ~ (row_number() - 1 - sum(Position < 817)) %% 3
    ),
    y_lab = -(3.0 + escalon * 0.60)
  )%>%
  ungroup()

ymax_plot <- max(df_lab_azul$y_lab, na.rm = TRUE) + 0.5
ymin_plot <- min(df_lab_rojo$y_lab, na.rm = TRUE) - 0.5

# ── Niveles de referencia ─────────────────────────────────────────────────────
niveles_ref <- data.frame(
  y     = c(0.6, 1.1, 1.6, 2.1, -0.6, -1.1, -1.6, -2.1),
  label = rep(c("\u22655%", "\u226510%", "\u226515%", "\u226520%"), 2),
  lado  = rep(c("arriba", "abajo"), each = 4)
)

# ── Plot UL54 ─────────────────────────────────────────────────────────────────
p_ul54 <- ggplot() +
  
  # Líneas de referencia — de -55 a aa 1243
  geom_segment(data = niveles_ref,
               aes(x = -55, xend = 1243, y = y, yend = y),
               color = "grey85", linewidth = 0.3, linetype = "dashed") +
  
  # Etiquetas de prevalencia
  geom_text(data = niveles_ref %>% filter(lado == "arriba"),
            aes(x = -70, y = y, label = label),
            hjust = 1, size = 4.0, color = "grey50", fontface = "bold") +
  geom_text(data = niveles_ref %>% filter(lado == "abajo"),
            aes(x = -70, y = y, label = label),
            hjust = 1, size = 4.0, color = "grey50", fontface = "bold") +
  
  # Dominios
  geom_rect(data = dominios_ul54,
            aes(xmin = start, xmax = end, ymin = 0.0, ymax = 0.40, fill = fill),
            alpha = 0.85, color = "white", linewidth = 0.25) +
  geom_text(data = dominios_ul54 %>% mutate(mid = (start + end) / 2),
            aes(x = mid, y = 0.23, label = label),
            size = 4.0, color = "grey15", fontface = "bold",
            vjust = 0.5, hjust = 0.5) +
  
  # Hotspots
  geom_rect(data = hotspots_ul54,
            aes(xmin = start, xmax = end, ymin = -0.40, ymax = 0.0, fill = fill),
            alpha = 0.9, color = "white", linewidth = 0.25) +
  
  scale_fill_identity() +
  
  # Backbone aa 1–1243 con tapas
  geom_segment(aes(x = 1,    xend = 1243, y = 0,     yend = 0),
               color = "grey10", linewidth = 1.0) +
  geom_segment(aes(x = 1,    xend = 1,    y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  geom_segment(aes(x = 1243, xend = 1243, y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  
  # Etiquetas extremos backbone
  annotate("text", x = 1,    y = -0.28, label = "1",    size = 4.0,
           color = "grey10", vjust = 1) +
  annotate("text", x = 1243, y = -0.28, label = "1243", size = 4.0,
           color = "grey10", vjust = 1) +
  
  # Stems
  geom_segment(data = df,
               aes(x = Position, xend = Position,
                   y = ifelse(en_resist, -0.42, 0.42),
                   yend = y_point, color = color),
               linewidth = 0.45, alpha = 0.8) +
  
  # Puntos
  geom_point(data = df,
             aes(x = Position, y = y_point, color = color, size = prev_max),
             alpha = 0.88) +
  scale_color_identity() +
  scale_size_continuous(range = c(1.2, 4.5), guide = "none") +
  
  # Guías azules
  geom_segment(data = df_lab_azul,
               aes(x = Position, xend = Position,
                   y = y_point + 0.08, yend = y_lab - 0.08),
               color = "#2980B9", linewidth = 0.3, alpha = 0.4,
               linetype = "dotted") +
  
  # Guías rojas
  geom_segment(data = df_lab_rojo,
               aes(x = Position, xend = Position,
                   y = y_point - 0.08, yend = y_lab + 0.08),
               color = "#C0392B", linewidth = 0.3, alpha = 0.4,
               linetype = "dotted") +
  
  # Etiquetas azules
  geom_text(data = df_lab_azul,
            aes(x = Position, y = y_lab, label = label_mut),
            color = "#2980B9", size = 4.0, fontface = "bold",
            hjust = 0.5, vjust = 0) +
  
  # Etiquetas rojas
  geom_text(data = df_lab_rojo,
            aes(x = Position, y = y_lab, label = label_mut),
            color = "#C0392B", size = 4.0, fontface = "bold",
            hjust = 0.5, vjust = 1) +
  
  scale_x_continuous(
    limits = c(-130, 1380),
    breaks = NULL,
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(limits = c(ymin_plot, ymax_plot)) +
  labs(title = NULL, x = "", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(t = 8, r = 12, b = 8, l = 50)
  )

# ── Leyenda hotspots ──────────────────────────────────────────────────────────
leyenda_hotspots_ul54 <- ggplot(data.frame(
  x     = c(1.0, 4.2, 7.4, 10.6),
  y     = 1,
  label = c("GCV/CDV", "GCV/CDV/FOS", "FOS > GCV", "FOS"),
  color = c("#E67E22", "#8E44AD", "#E74C3C", "#C0392B")
)) +
  geom_point(aes(x, y, color = color), size = 4, shape = 15) +
  geom_text(aes(x + 0.18, y, label = label),
            hjust = 0, size = 4.0, color = "grey20") +
  scale_color_identity() +
  scale_x_continuous(limits = c(0.5, 13)) +
  scale_y_continuous(limits = c(0.5, 1.5)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(t = 0, b = 6, l = 50)
  )

# ── Composición final ─────────────────────────────────────────────────────────
figura_ul54 <- p_ul54 / leyenda_hotspots_ul54 +
  plot_layout(heights = c(10, 1)) +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  )

figura_ul54

# [guardado antiguo] ggsave(
# [guardado antiguo]   "~/Desktop/Figuras_tesis/figura_lollipop_UL54.png",
# [guardado antiguo]   figura_ul54,
# [guardado antiguo]   width = 34, height = 18, units = "cm", dpi = 400
# [guardado antiguo] )

guardar_fig(figura_ul54, "Figura_46", ancho = 34, alto = 18, dpi = 400)

