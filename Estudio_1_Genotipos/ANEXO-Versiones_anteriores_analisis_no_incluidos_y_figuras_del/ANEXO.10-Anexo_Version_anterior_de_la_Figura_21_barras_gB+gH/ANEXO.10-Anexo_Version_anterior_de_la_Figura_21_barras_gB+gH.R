# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.10. Anexo: Versión anterior de la Figura 21 (barras gB+gH)
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Versión anterior de la Figura 21 (barras gB+gH) ----
## 1) Función para limpiar texto
clean_combo_total <- function(x) {
  x <- gsub("GB", "gB", x)
  x <- gsub("GH", "gH", x)
  x <- gsub("\\+", " + ", x)
  x <- gsub("\\s+/\\s+", " / ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

## 2) Etiquetas de grupo
desglose_gBunico_gHunico$grupo <- "Genotipos únicos gB / gH"
desglose_gBmixto_gHunico$grupo <- "Genotipos mixtos gB / únicos gH"
desglose_gBunico_gHmixto$grupo <- "Genotipos únicos gB / mixtos gH"
desglose_gBmixto_gHmixto$grupo <- "Genotipos mixtos gB / gH"

## 3) Unir todo
desglose_all_gBgH <- bind_rows(
  desglose_gBunico_gHunico,
  desglose_gBmixto_gHunico,
  desglose_gBunico_gHmixto,
  desglose_gBmixto_gHmixto
) %>%
  mutate(
    grupo = factor(
      grupo,
      levels = c(
        "Genotipos únicos gB / gH",
        "Genotipos únicos gB / mixtos gH",
        "Genotipos mixtos gB / únicos gH",
        "Genotipos mixtos gB / gH"
      )
    ),
    combo_plot = clean_combo_total(combo_total)
  )

## 4) Cálculo de etiquetas, con % mixtos sobre total = 131
TOTAL_MIXTOS <- 131

desglose_all_gBgH <- desglose_all_gBgH %>%
  mutate(
    tipo_general = if_else(grupo == "Genotipos únicos gB / gH", "único", "mixto")
  ) %>%
  group_by(grupo) %>%
  mutate(
    pct_grupo = if_else(
      tipo_general == "mixto",
      round(100 * n / TOTAL_MIXTOS, 1),   # sobre total mixtos
      round(100 * n / sum(n), 1)          # únicos siguen igual
    ),
    pct_label = if_else(
      pct_grupo %% 1 == 0,
      as.character(as.integer(pct_grupo)),
      formatC(pct_grupo, digits = 1, format = "f")
    ),
    etiqueta = sprintf("%d (%s%%)", n, pct_label),
    total_grupo = sum(n)
  ) %>%
  ungroup()

## 5) Variable de facet
desglose_all_gBgH <- desglose_all_gBgH %>%
  mutate(
    grupo2 = if_else(grupo == "Genotipos únicos gB / gH",
                     "Genotipos únicos gB / gH",
                     "Genotipos mixtos gB / gH")
  ) %>%
  group_by(grupo2) %>%
  mutate(total_grupo2 = sum(n)) %>%
  ungroup()

## 6) Orden de combinaciones
orden_uni <- desglose_all_gBgH %>%
  filter(grupo2 == "Genotipos únicos gB / gH") %>%
  group_by(combo_plot) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  arrange(combo_plot)

orden_mix <- desglose_all_gBgH %>%
  filter(grupo2 == "Genotipos mixtos gB / gH") %>%
  group_by(combo_plot) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  arrange(desc(n), combo_plot)

levels_combo <- c(orden_uni$combo_plot, orden_mix$combo_plot)
desglose_all_gBgH$combo_plot <- factor(desglose_all_gBgH$combo_plot, levels = levels_combo)

## 7) Etiquetas de los facets
desglose_all_gBgH <- desglose_all_gBgH %>%
  mutate(
    strip_lab = if_else(
      grupo2 == "Genotipos únicos gB / gH",
      sprintf("Genotipos únicos gB / gH\n(n=%d)",
              unique(total_grupo2[grupo2 == "Genotipos únicos gB / gH"])),
      sprintf("Genotipos mixtos gB / gH\n(n=%d)",
              unique(total_grupo2[grupo2 == "Genotipos mixtos gB / gH"]))
    )
  ) %>%
  mutate(
    strip_lab = factor(
      strip_lab,
      levels = c(
        sprintf("Genotipos únicos gB / gH\n(n=%d)",
                unique(total_grupo2[grupo2 == "Genotipos únicos gB / gH"])),
        sprintf("Genotipos mixtos gB / gH\n(n=%d)",
                unique(total_grupo2[grupo2 == "Genotipos mixtos gB / gH"]))
      )
    )
  )

## 8) Gráfico final
distribuciones_gBgH <- ggplot(desglose_all_gBgH,
                              aes(x = combo_plot, y = n, fill = grupo)) +
  geom_col(width = 0.75) +
  geom_text(
    aes(label = etiqueta),
    vjust = -0.2,
    size  = 3,
    angle = 50,
    hjust = 0,
    lineheight = 0.9
  ) +
  facet_grid(~ strip_lab, scales = "free_x", space = "free_x") +
  scale_fill_manual(
    values = c(
      "Genotipos únicos gB / gH"         = "#3CB371",
      "Genotipos únicos gB / mixtos gH"  = "#E57373",
      "Genotipos mixtos gB / únicos gH"  = "#F4A2A2",
      "Genotipos mixtos gB / gH"         = "#B22222"
    ),
    breaks = c(
      "Genotipos únicos gB / mixtos gH",
      "Genotipos mixtos gB / únicos gH",
      "Genotipos mixtos gB / gH"
    ),
    labels = c(
      "Genotipos únicos gB / mixtos gH",
      "Genotipos mixtos gB / únicos gH",
      "Genotipos mixtos gB / gH"
    )
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(x = "", y = "Número de muestras", fill = "") +
  theme_classic(base_size = 13) +
  theme(
    plot.margin      = margin(t = 25, r = 15, b = 30, l = 10),
    legend.position  = "top",
    legend.key.size  = unit(8, "pt"),
    legend.direction = "horizontal",
    legend.box       = "horizontal",
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.line        = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks       = element_line(linewidth = 0.3, colour = "grey50"),
    axis.title.y     = element_text(size = 11),
    strip.text       = element_text(face = "bold"),
    panel.spacing.y  = unit(2, "lines"),
    panel.spacing.x  = unit(0.9, "lines")
  )

distribuciones_gBgH

