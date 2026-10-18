# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.4. Variantes no descritas previamente (VND)
# Extracto del script completo (ver README.md del estudio)

# ---- 10.2.4. Variantes no descritas previamente (VND) ----

## ---- Datos de VND ----
# ── Resumen actualizado de variantes NEW ──────────────────────────────────────
new_vars <- all_data_txp %>%
  filter(class == "NEW", Prevalence >= 5) %>%
  mutate(Position = as.numeric(Position)) %>%
  group_by(Protein, Mutation, Position) %>%
  summarise(
    n_muestras  = n_distinct(Npet),
    n_pacientes = n_distinct(NHC),
    prev_media  = round(mean(Prevalence, na.rm = TRUE), 1),
    prev_max    = round(max(Prevalence, na.rm = TRUE), 1),
    .groups = "drop"
  )

cat("Total NEW por gen:\n")
new_vars %>% count(Protein) %>% print()

cat("\nNEW por gen y umbral:\n")
map_dfr(c(5, 10, 15, 20), function(thr) {
  all_data_txp %>%
    filter(class == "NEW", Prevalence >= thr) %>%
    group_by(Protein) %>%
    summarise(n_variantes = n_distinct(Mutation),
              n_pacientes = n_distinct(NHC), .groups = "drop") %>%
    mutate(umbral = paste0("≥", thr, "%"))
}) %>%
  arrange(Protein, umbral) %>%
  print(n = Inf)


new_vars_dom <- new_vars %>%
  mutate(
    dominio = case_when(
      Protein == "UL97" & Position >= 335 & Position <= 370 ~ "P-loop",
      Protein == "UL97" & Position >= 409 & Position <= 411 ~ "Gatekeeper",
      Protein == "UL97" & Position >= 440 & Position <= 480 ~ "Loop catalítico",
      Protein == "UL97" & Position >= 520 & Position <= 607 ~ "Región canónica GCV",
      Protein == "UL54" & Position >= 408 & Position <= 412 ~ "Región IV",
      Protein == "UL54" & Position >= 501 & Position <= 545 ~ "Delta-C",
      Protein == "UL54" & Position >= 700 & Position <= 715 ~ "Región II",
      Protein == "UL54" & Position >= 725 & Position <= 964 ~ "Unión ATP/dNTP",
      Protein == "UL54" & Position >= 802 & Position <= 821 ~ "Región III",
      Protein == "UL54" & Position >= 910 & Position <= 913 ~ "Región I",
      Protein == "UL54" & Position >= 962 & Position <= 987 ~ "Región V",
      Protein == "UL54" & Position >= 1153 & Position <= 1159 ~ "NLS-A",
      Protein == "UL54" & Position >= 1213 & Position <= 1242 ~ "UL44 binding",
      Protein == "UL56" & Position == 25 ~ "N-terminal",
      Protein == "UL56" & Position >= 202 & Position <= 220 ~ "LAGLIDADG",
      Protein == "UL56" & Position >= 231 & Position <= 369 ~ "Región LMV",
      Protein == "UL56" & Position >= 514 & Position <= 522 ~ "Interacción UL89",
      Protein == "UL56" & Position >= 671 & Position <= 680 ~ "Interacción UL89",
      Protein == "UL56" & Position >= 757 & Position <= 764 ~ "Interacción UL89",
      TRUE ~ "Fuera de dominio funcional"
    ),
    en_dominio = dominio != "Fuera de dominio funcional"
  ) %>%
  arrange(Protein, desc(en_dominio), desc(n_pacientes), desc(prev_max))

cat("── En dominio funcional ─────────────────────────\n")
new_vars_dom %>% filter(en_dominio) %>% print(n = Inf)

cat("\n── Fuera de dominio pero relevantes ─────────────\n")
new_vars_dom %>%
  filter(!en_dominio, n_pacientes >= 2 | prev_max >= 15) %>%
  print(n = Inf)

## ---- Figura 45. VND en UL97 ----
dominios_ul97 <- data.frame(
  start = c(230),
  end   = c(620),
  label = c("Dominio quinasa"),
  fill  = c("lightskyblue")
)

hotspots_ul97 <- data.frame(
  start = c(335, 397, 519, 590),
  end   = c(370, 480, 521, 607),
  label = c("", "", "", ""),
  fill  = c("#B03060", "#8E44AD", "#E67E22", "#E67E22")
)

ul97_umbrales <- all_data_txp %>%
  filter(class == "NEW", Protein == "UL97") %>%
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

ul97 <- new_vars_dom %>%
  filter(Protein == "UL97") %>%
  left_join(ul97_umbrales %>% select(Mutation, umbral_max), by = "Mutation")

nivel <- c("5" = 0.6, "10" = 1.1, "15" = 1.6, "20" = 2.1)
paso  <- 0.4

df <- ul97 %>%
  mutate(
    y_nivel    = nivel[as.character(umbral_max)],
    en_quinasa = Position >= 230 & Position <= 620,
    y_point    = ifelse(en_quinasa, -y_nivel, y_nivel),
    color      = ifelse(en_quinasa, "#C0392B", "#2980B9"),
    label_mut  = case_when(
      en_quinasa                        ~ Mutation,
      prev_max >= 10 | n_pacientes >= 4 ~ Mutation,
      TRUE                              ~ NA_character_
    )
  )

df_lab_azul <- df %>%
  filter(!is.na(label_mut), !en_quinasa) %>%
  arrange(Position) %>%
  mutate(prox_group = cumsum(c(TRUE, diff(Position) > 30))) %>%
  group_by(prox_group) %>%
  mutate(
    within_rank = row_number(),
    y_lab       = 3.0 + (within_rank - 1) * 0.60
  ) %>%
  ungroup()

df_lab_rojo <- df %>%
  filter(!is.na(label_mut), en_quinasa) %>%
  arrange(Position) %>%
  mutate(prox_group = cumsum(c(TRUE, diff(Position) > 30))) %>%
  group_by(prox_group) %>%
  mutate(
    within_rank = row_number(),
    y_lab       = -(3.0 + (within_rank - 1) * 0.60)
  ) %>%
  ungroup()

ymax_plot <- max(
  max(df_lab_azul$y_lab, na.rm = TRUE),
  max(abs(df_lab_rojo$y_lab), na.rm = TRUE)
) + 0.5

niveles_ref <- data.frame(
  y     = c(0.6, 1.1, 1.6, 2.1, -0.6, -1.1, -1.6, -2.1),
  label = rep(c("\u22655%", "\u226510%", "\u226515%", "\u226520%"), 2),
  lado  = rep(c("arriba", "abajo"), each = 4)
)

p_ul97 <- ggplot() +
  
  # Líneas de referencia — empiezan después de las etiquetas, acaban en aa 707
  geom_segment(data = niveles_ref,
               aes(x = -42, xend = 707, y = y, yend = y),
               color = "grey85", linewidth = 0.3, linetype = "dashed") +
  
  # Etiquetas de prevalencia
  geom_text(data = niveles_ref %>% filter(lado == "arriba"),
            aes(x = -55, y = y, label = label),
            hjust = 1, size = 4.0, color = "grey50", fontface = "bold") +
  geom_text(data = niveles_ref %>% filter(lado == "abajo"),
            aes(x = -55, y = y, label = label),
            hjust = 1, size = 4.0, color = "grey50", fontface = "bold") +
  
  # Dominios
  geom_rect(data = dominios_ul97,
            aes(xmin = start, xmax = end, ymin = 0.0, ymax = 0.40, fill = fill),
            alpha = 0.9, color = "white", linewidth = 0.25) +
  geom_text(data = dominios_ul97 %>% mutate(mid = (start + end) / 2),
            aes(x = mid, y = 0.20, label = label),
            size = 4.0, color = "black", fontface = "bold",
            vjust = 0.5, hjust = 0.5) +
  
  # Hotspots
  geom_rect(data = hotspots_ul97,
            aes(xmin = start, xmax = end, ymin = -0.40, ymax = 0.0, fill = fill),
            alpha = 0.9, color = "white", linewidth = 0.25) +
  geom_text(data = hotspots_ul97 %>% mutate(mid = (start + end) / 2),
            aes(x = mid, y = -0.20, label = label),
            size = 3, color = "white", fontface = "bold",
            vjust = 0.5, hjust = 0.5) +
  
  scale_fill_identity() +
  
  # Backbone aa 1–707 con tapas
  geom_segment(aes(x = 1,   xend = 707, y = 0,     yend = 0),
               color = "grey10", linewidth = 1.0) +
  geom_segment(aes(x = 1,   xend = 1,   y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  geom_segment(aes(x = 707, xend = 707, y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  
  # Etiquetas extremos backbone
  annotate("text", x = 1,   y = -0.28, label = "1",   size = 4.0,
           color = "grey10", vjust = 1) +
  annotate("text", x = 707, y = -0.28, label = "707", size = 4.0,
           color = "grey10", vjust = 1) +
  
  # Stems
  geom_segment(data = df,
               aes(x = Position, xend = Position,
                   y = ifelse(en_quinasa, -0.42, 0.42),
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
    limits = c(-100, 730),
    breaks = NULL,
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(limits = c(-ymax_plot, ymax_plot)) +
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
    plot.margin        = margin(t = 8, r = 12, b = 8, l = 35)
  )

leyenda_hotspots <- ggplot(data.frame(
  x     = c(1.5, 3.7, 6),
  y     = c(1, 1, 1),
  label = c("GCV/MBV-R (335-370)",
            "MBV-R (397-480)",
            "GCV-R cod\u00f3n 520 + (590-607)"),
  color = c("#B03060", "#8E44AD", "#E67E22")
)) +
  geom_point(aes(x, y, color = color), size = 4, shape = 15) +
  geom_text(aes(x + 0.15, y, label = label),
            hjust = 0, size = 4.0, color = "grey20") +
  scale_color_identity() +
  scale_x_continuous(limits = c(0.8, 9.5)) +
  scale_y_continuous(limits = c(0.5, 1.5)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(t = 0, b = 6, l = 35)
  )

figura_ul97 <- p_ul97 / leyenda_hotspots +
  plot_layout(heights = c(10, 1))

figura_ul97
# [guardado antiguo] ggsave(
# [guardado antiguo]   "~/Desktop/Figuras_tesis/figura_lollipop_UL97.png",
# [guardado antiguo]   figura_ul97,
# [guardado antiguo]   width = 30, height = 15, units = "cm", dpi = 400,
# [guardado antiguo]   bg = "white"
# [guardado antiguo] )

guardar_fig(figura_ul97, "Figura_45", ancho = 30, alto = 15, dpi = 400)

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

## ---- Figura 47. VND en UL56 ----
dominios_ul56 <- data.frame(
  start = c(4,   117, 280, 514, 671, 757),
  end   = c(20,  130, 310, 522, 680, 764),
  label = c("NLS", "endonucleasa", "Dedo de zinc", "Interacción pUL89", "Estab. 1", "Estab. 2"),
  fill  = c("cadetblue3", "cadetblue3", "cadetblue3", "#F0A500", "#F0A500", "#F0A500")
)

hotspots_ul56 <- data.frame(
  start = c(20,  230, 392),
  end   = c(30,  370, 398),
  label = c("",  "", ""),
  fill  = c("#B03060", "#B03060", "#B03060")
)

ul56_umbrales <- all_data_txp %>%
  filter(class == "NEW", Protein == "UL56") %>%
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

ul56 <- new_vars_dom %>%
  filter(Protein == "UL56") %>%
  left_join(ul56_umbrales %>% select(Mutation, umbral_max), by = "Mutation")

nivel <- c("5" = 0.6, "10" = 1.1, "15" = 1.6, "20" = 2.1)
paso  <- 0.4

df <- ul56 %>%
  mutate(
    y_nivel   = nivel[as.character(umbral_max)],
    en_lmv    = Position >= 231 & Position <= 412,
    y_point   = ifelse(en_lmv, -y_nivel, y_nivel),
    color     = ifelse(en_lmv, "#C0392B", "#2980B9"),
    label_mut = Mutation
  )

df_lab_azul_base <- df %>%
  filter(!is.na(label_mut), !en_lmv) %>%
  arrange(Position) %>%
  mutate(row_idx = row_number())

df_lab_azul <- df_lab_azul_base %>%
  mutate(
    y_lab_raw = 2.5 + (row_idx - 1) * paso,
    y_lab_raw = case_when(
      Position <= 30 ~ y_lab_raw + 0.8,
      TRUE           ~ y_lab_raw
    ),
    y_lab = case_when(
      Position >= 595 ~ y_lab_raw - 1.6,
      TRUE            ~ y_lab_raw
    )
  )

df_lab_rojo <- df %>%
  filter(!is.na(label_mut), en_lmv) %>%
  arrange(Position) %>%
  mutate(y_lab = -(0.9 + (row_number() - 1) * paso))

ymax_plot <- max(df_lab_azul$y_lab, na.rm = TRUE) + 0.5

niveles_ref <- data.frame(
  y     = c(0.6, 1.1, 1.6, 2.1, -0.6, -1.1, -1.6, -2.1),
  label = rep(c("\u22655%", "\u226510%", "\u226515%", "\u226520%"), 2),
  lado  = rep(c("arriba", "abajo"), each = 4)
)

p_ul56 <- ggplot() +
  geom_segment(data = niveles_ref,
               aes(x = -42, xend = 850, y = y, yend = y),
               color = "grey85", linewidth = 0.3, linetype = "dashed") +
  # Etiquetas prevalencia — más cerca del aa 1 (x = -55, límite -75)
  geom_text(data = niveles_ref %>% filter(lado == "arriba"),
            aes(x = -55, y = y, label = label),
            hjust = 1, size = 3.5, color = "grey50", fontface = "bold") +
  geom_text(data = niveles_ref %>% filter(lado == "abajo"),
            aes(x = -55, y = y, label = label),
            hjust = 1, size = 3.5, color = "grey50", fontface = "bold") +
  # Dominios
  geom_rect(data = dominios_ul56,
            aes(xmin = start, xmax = end, ymin = 0.0, ymax = 0.40, fill = fill),
            alpha = 0.9, color = "white", linewidth = 0.25) +
  geom_text(data = dominios_ul56 %>% mutate(mid = (start + end) / 2),
            aes(x = mid, y = 0.20, label = label),
            size = 3.5, color = "black", fontface = "bold",
            vjust = 0.5, hjust = 0.5, lineheight = 0.85) +
  # Hotspots
  geom_rect(data = hotspots_ul56,
            aes(xmin = start, xmax = end, ymin = -0.40, ymax = 0.0, fill = fill),
            alpha = 0.9, color = "white", linewidth = 0.25) +
  scale_fill_identity() +
  # Backbone aa 1–850
  geom_segment(aes(x = 1, xend = 850, y = 0, yend = 0),
               color = "grey10", linewidth = 1.0) +
  # Tapas extremos (palitos)
  geom_segment(aes(x = 1,   xend = 1,   y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  geom_segment(aes(x = 850, xend = 850, y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  # Stems
  geom_segment(data = df,
               aes(x = Position, xend = Position,
                   y = ifelse(en_lmv, -0.42, 0.42),
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
            color = "#2980B9", size = 3.5, fontface = "bold",
            hjust = 0.5, vjust = 0) +
  # Etiquetas rojas
  geom_text(data = df_lab_rojo,
            aes(x = Position, y = y_lab, label = label_mut),
            color = "#C0392B", size = 3.5, fontface = "bold",
            hjust = 0.5, vjust = 1) +
  annotate("text", x = 1,   y = -0.28, label = "1",   size = 3.5, color = "grey10", vjust = 1) +
  annotate("text", x = 850, y = -0.28, label = "850", size = 3.5, color = "grey10", vjust = 1) +
  scale_x_continuous(
    limits = c(-100, 870),
    breaks = NULL,
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(limits = c(-2.5, ymax_plot)) +
  # Sin título de panel — el nombre "UL56" va en el pie de figura
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
    plot.margin        = margin(t = 8, r = 12, b = 8, l = 35)
  )

leyenda_hotspots_ul56 <- ggplot(data.frame(
  x     = 1.5,
  y     = 1,
  label = "LMV-R",
  color = "#B03060"
)) +
  geom_point(aes(x, y, color = color), size = 4, shape = 15) +
  geom_text(aes(x + 0.15, y, label = label),
            hjust = 0, size = 3.5, color = "grey20") +
  scale_color_identity() +
  scale_x_continuous(limits = c(0.8, 4)) +
  scale_y_continuous(limits = c(0.5, 1.5)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(t = 0, b = 6, l = 35)
  )

figura_ul56 <- p_ul56 / leyenda_hotspots_ul56 +
  plot_layout(heights = c(10, 1))

figura_ul56

# [guardado antiguo] ggsave(
# [guardado antiguo]   "~/Desktop/Figuras_tesis/figura_lollipop_UL56.png",
# [guardado antiguo]   figura_ul56,
# [guardado antiguo]   width = 28, height = 15, units = "cm", dpi = 400,
# [guardado antiguo]   bg = "white"
# [guardado antiguo] )

guardar_fig(figura_ul56, "Figura_47", ancho = 28, alto = 15, dpi = 400)

## ---- Anexo II. VND con región funcional y evidencia técnica ----
# El bloque guarda el Word apaisado (Anexo_II.docx) en la carpeta de salidas.
## Requiere all_data_txp (con class y columnas DeepChek: cobertura_total, mut_fw, mut_rv, q_score).
## Una fila por VND; frecuencia y cobertura resumen todas las muestras en que se detectó.


fmt_n <- function(x) format(round(x), big.mark = ".", decimal.mark = ",", trim = TRUE)
fmt_1 <- function(x) sub("\\.", ",", sprintf("%.1f", x))

# ── Regiones funcionales y hotspots (editar límites según Material y métodos) ──
## Se asigna la primera fila que contenga la posición: los hotspots van antes que los dominios.
regiones_funcionales <- tribble(
  ~Protein, ~ini,  ~fin,  ~etiqueta,
  # UL97
  "UL97",   335,   370,   "GCV/MBV-R (aa 335\u2013370)",
  "UL97",   397,   480,   "MBV-R (aa 397\u2013480)",
  "UL97",   590,   607,   "GCV-R (aa 590\u2013607)",
  "UL97",   335,   665,   "GCV-R (aa 335\u2013665)",
  "UL97",   337,   651,   "pUL97-KD",
  # UL54
  "UL54",   300,   325,   "GCV/CDV-R (aa 300\u2013325)",
  "UL54",   498,   545,   "GCV/CDV-R (aa 498\u2013545)",
  "UL54",   700,   715,   "FOS/GCV-R (aa 700\u2013715)",
  "UL54",   750,   760,   "FOS-R (aa 750\u2013760)",
  "UL54",   802,   825,   "GCV/CDV/FOS-R (aa 802\u2013825)",
  "UL54",   985,   990,   "GCV/CDV-R (aa 985\u2013990)",
  "UL54",   300,   545,   "pUL54-Exo",
  "UL54",   696,   988,   "pUL54-Pol",
  "UL54",  1153,  1242,   "pUL54-Cter/NLS",
  # UL56
  "UL56",   230,   370,   "LMV-R (aa 230\u2013370)"
)
asignar_region <- function(prot, pos) {
  r <- regiones_funcionales %>% filter(Protein == prot, pos >= ini, pos <= fin)
  if (nrow(r) == 0) "" else r$etiqueta[1]
}

# ── Denominadores por gen ────────────────────────────────────────────────────
denom <- all_data_txp %>% group_by(Protein) %>%
  summarise(N_m = n_distinct(Npet), N_p = n_distinct(RTP), .groups = "drop")

# ── Tabla ────────────────────────────────────────────────────────────────────
anexo2 <- all_data_txp %>%
  filter(class == "NEW", Prevalence >= 5) %>%
  mutate(mut_tot = mut_fw + mut_rv,
         balance = pmin(mut_fw, mut_rv) / mut_tot) %>%
  group_by(Protein, Mutation) %>%
  summarise(
    n_m       = n_distinct(Npet),
    n_p       = n_distinct(RTP),
    pacientes = paste(unique(RTP[order(as.numeric(sub("RTP-", "", RTP)))]), collapse = ", "),
    prev_mean = mean(Prevalence), prev_sd = sd(Prevalence),
    cob_med   = median(cobertura_total),
    cob_q1    = quantile(cobertura_total, .25), cob_q3 = quantile(cobertura_total, .75),
    lect_min  = min(mut_tot),
    fw_rv_min = paste0(mut_fw[which.min(mut_tot)], "/", mut_rv[which.min(mut_tot)]),
    bal_min   = min(balance),
    q_min     = min(q_score),
    .groups = "drop") %>%
  left_join(denom, by = "Protein") %>%
  mutate(
    Posicion  = as.numeric(gsub("[^0-9]", "", Mutation)),
    Region    = mapply(asignar_region, Protein, Posicion),
    Gen       = factor(Protein, levels = c("UL97", "UL54", "UL56"),
                       labels = c(sprintf("pUL97 (%d muestras, %d pacientes)", denom$N_m[denom$Protein == "UL97"], denom$N_p[denom$Protein == "UL97"]),
                                  sprintf("pUL54 (%d muestras, %d pacientes)", denom$N_m[denom$Protein == "UL54"], denom$N_p[denom$Protein == "UL54"]),
                                  sprintf("pUL56 (%d muestras, %d pacientes)", denom$N_m[denom$Protein == "UL56"], denom$N_p[denom$Protein == "UL56"]))),
    Muestras   = sprintf("%d (%s)", n_m, fmt_1(100 * n_m / N_m)),
    Pacientes  = sprintf("%d (%s)", n_p, fmt_1(100 * n_p / N_p)),
    Frecuencia = ifelse(n_m == 1, fmt_1(prev_mean),
                        sprintf("%s \u00b1 %s", fmt_1(prev_mean), fmt_1(prev_sd))),
    Cobertura  = ifelse(n_m == 1, paste0(fmt_n(cob_med), " \u00d7"),
                        sprintf("%s \u00d7 (%s\u2013%s)", fmt_n(cob_med), fmt_n(cob_q1), fmt_n(cob_q3))),
    Lecturas   = sprintf("%d (%s)", lect_min, fw_rv_min),
    Hebra      = fmt_1(100 * bal_min),
    Q          = as.character(q_min)
  ) %>%
  arrange(Gen, Posicion) %>%
  select(Gen, Mutation, Muestras, Pacientes, pacientes, Frecuencia, Cobertura, Lecturas, Hebra, Q, Region)

cat("VND por gen:\n"); print(table(anexo2$Gen))

# ── Flextable con fila combinada por gen ─────────────────────────────────────
anexo2_grp <- as_grouped_data(anexo2, groups = "Gen")

ft_anexo2 <- as_flextable(anexo2_grp, hide_grouplabel = TRUE) %>%
  set_header_labels(
    Mutation = "Variante",
    Muestras = "Muestras, n (%)", Pacientes = "Pacientes, n (%)",
    pacientes = "Pacientes (Tabla X)",
    Frecuencia = "Frecuencia de variante (%), media \u00b1 DE",
    Cobertura = "Cobertura por posici\u00f3n, mediana (RIC)",
    Lecturas = "Lecturas mutadas, m\u00ednimo (fw/rv)",
    Hebra = "Hebra minoritaria, m\u00ednimo (%)",
    Q = "Q-score, m\u00ednimo",
    Region = "Regi\u00f3n funcional / hotspot de resistencia"
  ) %>%
  align(align = "center", part = "all") %>%
  align(j = c(1, 10), align = "left", part = "body") %>%
  bold(part = "header") %>%
  bold(i = ~ !is.na(Gen), bold = TRUE) %>%
  bg(i = ~ !is.na(Gen), bg = "#1F3B73") %>%
  color(i = ~ !is.na(Gen), color = "white") %>%
  bg(part = "header", bg = "#1F3B73") %>% color(part = "header", color = "white") %>%
  bold(j = "Mutation", bold = TRUE) %>%
  flextable::font(fontname = "Source Sans Pro", part = "all") %>%
  fontsize(size = 7.5, part = "all") %>%
  border_inner(border = fp_border(color = "grey70", width = 0.5)) %>%
  border_outer(border = fp_border(color = "grey40", width = 1)) %>%
  add_footer_lines(paste0(
    "VND: variantes no descritas previamente (no recogidas en CHARMD, versi\u00f3n actualizada a 20 de abril de 2026). ",
    "Muestras y pacientes sobre el total con resultado v\u00e1lido para cada gen. Frecuencia y cobertura resumen las muestras en ",
    "que se detect\u00f3 la variante (valor \u00fanico cuando se detect\u00f3 en una sola). Lecturas mutadas: n\u00famero m\u00ednimo de lecturas que ",
    "soportan la variante en cualquiera de las muestras, con su reparto entre hebras directa (fw) e inversa (rv). Hebra minoritaria: ",
    "proporci\u00f3n m\u00ednima de lecturas mutadas en la hebra menos representada. El filtro por Q-score del sistema experto estaba ",
    "desactivado. pUL97-KD: dominio quinasa de pUL97; Exo: dominio exonucleasa; Pol: dominio polimerasa; Cter/NLS: extremo ",
    "C-terminal y se\u00f1al de localizaci\u00f3n nuclear; -R: regi\u00f3n asociada a resistencia; GCV: ganciclovir; CDV: cidofovir; ",
    "FOS: foscarnet; MBV: maribavir; LMV: letermovir. * Cod\u00f3n de terminaci\u00f3n prematuro.")) %>%
  fontsize(size = 6.5, part = "footer") %>%
  set_table_properties(layout = "autofit", width = 1)

ft_anexo2

# ── Exportar a Word en horizontal ────────────────────────────────────────────
doc <- read_docx() %>%
  body_add_par("Anexo II. Variantes no descritas previamente identificadas mediante NGS", style = "heading 2") %>%
  body_add_flextable(ft_anexo2) %>%
  body_end_section_landscape()
print(doc, target = file.path(dir_salida, "Anexo_II.docx"))


