# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.1.1. Figura 14. Proceso de selección de muestras y pacientes (datos del diagrama)
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 14. Proceso de selección de muestras y pacientes (datos del diagrama) ----
# La Figura 14 es un diagrama de flujo montado fuera de R; este bloque genera los recuentos.
datos_cv <- data.frame(
  categoria = c(
    "Alguna CV ≥ 1.000 UI/mL",
    "Alguna CV < 1.000 UI/mL",
    "Todas CV neg / <LDC",
    "Sin registros de CV"
  ),
  pct = c(47.0, 14.7, 33.2, 5.1)
) %>%
  mutate(categoria = factor(categoria, levels = categoria))

ggplot(datos_cv, aes(x = reorder(categoria, pct), y = pct, fill = categoria)) +
  geom_col(width = 0.6) +
  geom_text(
    aes(label = paste0(gsub("\\.", ",", format(pct, nsmall = 1)), "%")),
    hjust = -0.2, fontface = "bold", size = 4.5, color = "grey20"
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 55),
    expand = expansion(mult = c(0, 0)),
    labels = function(x) gsub("\\.", ",", x)
  ) +
  scale_fill_manual(
    values = c(
      "Alguna CV ≥ 1.000 UI/mL" = "#C0392B",
      "Alguna CV < 1.000 UI/mL" = "#E8A598",
      "Todas CV neg / <LDC"     = "#5B9BD5",
      "Sin registros de CV"     = "grey60"
    )
  ) +
  theme_classic(base_size = 15) +
  theme(
    legend.position  = "none",
    axis.text.x      = element_text(size = 13),
    axis.text.y      = element_text(size = 13),
    axis.title.x     = element_text(size = 14),
    axis.line        = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
    plot.margin      = margin(10, 20, 10, 10)
  ) +
  labs(x = "", y = "% pacientes")

