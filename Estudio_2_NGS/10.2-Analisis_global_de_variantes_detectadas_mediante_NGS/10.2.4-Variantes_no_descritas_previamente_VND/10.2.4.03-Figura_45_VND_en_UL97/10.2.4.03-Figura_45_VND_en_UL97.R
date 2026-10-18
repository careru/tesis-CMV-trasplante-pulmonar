# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.4.03. Figura 45. VND en UL97
# Extracto del script completo (ver README.md del estudio)

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

