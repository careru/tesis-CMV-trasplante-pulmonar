# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.4.07. Figura 47. VND en UL56
# Extracto del script completo (ver README.md del estudio)

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

