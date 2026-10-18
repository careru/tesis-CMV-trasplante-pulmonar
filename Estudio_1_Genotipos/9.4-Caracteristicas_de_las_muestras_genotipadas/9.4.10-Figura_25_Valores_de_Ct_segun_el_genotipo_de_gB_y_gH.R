# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.4.10. Figura 25. Valores de Ct según el genotipo de gB y gH
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 25. Valores de Ct según el genotipo de gB y gH ----
# Formato largo
df_ct_long <- muestras_sel %>%
  select(GB1CT, GB2CT, GB3CT, GB4CT, GH1CT, GH2CT) %>%
  pivot_longer(
    cols      = everything(),
    names_to  = "Genotipo",
    values_to = "CT"
  ) %>%
  mutate(
    Genotipo = factor(
      Genotipo,
      levels = c("GB1CT", "GB2CT", "GB3CT", "GB4CT", "GH1CT", "GH2CT")
    ),
    Etiqueta = recode(
      Genotipo,
      "GB1CT" = "gB1", "GB2CT" = "gB2", "GB3CT" = "gB3", "GB4CT" = "gB4",
      "GH1CT" = "gH1", "GH2CT" = "gH2"
    ),
    Glucoproteina = if_else(grepl("^gB", Etiqueta), "", "")
  ) %>%
  filter(!is.na(CT))

# Paletas
cols_fill <- c(
  "GB1CT" = "#5B9BD5",
  "GB2CT" = "#5B9BD5",
  "GB3CT" = "#5B9BD5",
  "GB4CT" = "#5B9BD5",
  "GH1CT" = "#D98FD6",
  "GH2CT" = "#D98FD6"
)

cols_jitter <- c(
  "GB1CT" = "#1A5FA8",
  "GB2CT" = "#1A5FA8",
  "GB3CT" = "#1A5FA8",
  "GB4CT" = "#1A5FA8",
  "GH1CT" = "#A04EA0",
  "GH2CT" = "#A04EA0"
)

y_max_plot <- max(df_ct_long$CT, na.rm = TRUE) + 1.5

# Plot
p_ct <- ggplot(df_ct_long, aes(x = Etiqueta, y = CT, fill = Genotipo)) +
  geom_boxplot(
    width          = 0.50,
    color          = "grey25",
    linewidth      = 0.45,
    outlier.shape  = 21,
    outlier.size   = 1.8,
    outlier.fill   = "white",
    outlier.color  = "grey50",
    outlier.stroke = 0.4,
    outlier.alpha  = 0.7
  ) +
  geom_jitter(
    aes(color = Genotipo),
    width = 0.12, size = 1.0, alpha = 0.35, show.legend = FALSE
  ) +
  facet_grid(
    ~ Glucoproteina,
    scales = "free_x",
    space  = "free_x"
  ) +
  scale_fill_manual(values  = cols_fill)  +
  scale_color_manual(values = cols_jitter) +
  scale_y_continuous(
    name   = "Valor Ct",
    limits = c(24, y_max_plot),
    breaks = seq(25, 40, 5),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(x = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "none",
    panel.grid.major.y = element_line(color = "white", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_blank(),
    # Línea de eje solo en el panel izquierdo para Y, y en ambos para X
    axis.line.x        = element_line(color = "grey30", linewidth = 0.5),
    axis.line.y        = element_line(color = "grey30", linewidth = 0.5),
    axis.ticks         = element_line(color = "grey30", linewidth = 0.4),
    axis.ticks.length  = unit(3, "pt"),
    axis.text.x        = element_text(size = 11, color = "grey20", margin = margin(t = 4)),
    axis.text.y        = element_text(size = 10, color = "grey30"),
    axis.title.y       = element_text(size = 11, color = "grey20", margin = margin(r = 8)),
    # Strip
    strip.background   = element_blank(),
    strip.text         = element_text(size = 12, face = "bold", color = "grey20",
                                      margin = margin(b = 6)),
    panel.spacing.x    = unit(20, "pt"),
    plot.margin        = margin(10, 15, 8, 8)
  ) +
  coord_cartesian(clip = "off")

p_ct

guardar_fig(p_ct, "Figura_25", ancho = 16, alto = 10)

