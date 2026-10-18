# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.6.2. Datos para la Figura 30 (variabilidad intraepisodio; la figura de la tesis es una infografía)
# Extracto del script completo (ver README.md del estudio)

## ---- Datos para la Figura 30 (variabilidad intraepisodio; la figura de la tesis es una infografía) ----
df_fig <- variabilidad_ep %>%
  group_by(clase, variable) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(clase) %>%
  mutate(
    total = sum(n),
    pct   = 100 * n / total,
    label = paste0(n, " (", sub("\\.", ",", round(pct, 1)), "%)")
  ) %>%
  ungroup() %>%
  mutate(
    clase    = factor(clase, levels = c("Único", "Mixto"),
                      labels = c("Genotipo único gB/gH", "Genotipo mixto gB/gH")),
    variable = factor(variable, levels = c("Variable", "Estable"))
  )

figura_R21 <- ggplot(df_fig, aes(x = clase, y = pct, fill = variable)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4, fontface = "bold", color = "white"
  ) +
  scale_fill_manual(
    values = c("Estable"  = "#2E7D32",
               "Variable" = "#E65100"),
    labels = c("Estable"  = "Patrón dominante",
               "Variable" = "Genotipo fluctuante"),
    name   = ""
  ) +
  scale_y_continuous(
    name   = "% de episodios",
    limits = c(0, 105),
    breaks = seq(0, 100, 25),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "top",
    legend.text        = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
    axis.line.y        = element_line(color = "grey60", linewidth = 0.4),
    axis.text.x        = element_text(size = 12, color = "grey20"),
    axis.text.y        = element_text(size = 11, color = "grey40"),
    axis.title.y       = element_text(size = 11, color = "grey20"),
    plot.margin        = margin(8, 8, 8, 8)
  )

figura_R21

