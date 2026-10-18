# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.09. Anexo: Versión anterior de la Figura 20 (barras gH)
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Versión anterior de la Figura 20 (barras gH) ----
# 1) Asignar nombres de grupo
desglose_1_gh$grupo <- "Genotipos gH únicos"
desglose_2_gh$grupo <- "Genotipos gH mixtos"

# 2) Unir, calcular porcentajes dentro de cada grupo y etiquetas
desglose_all_gh <- bind_rows(desglose_1_gh, desglose_2_gh) %>%
  mutate(
    grupo = factor(grupo,
                   levels = c("Genotipos gH únicos", "Genotipos gH mixtos")),
    # Renombrar combos para mejor lectura
    combo_clean = case_when(
      combo == "GH1" ~ "gH1",
      combo == "GH2" ~ "gH2",
      combo == "GH1+GH2" ~ "gH1 + gH2",
      TRUE ~ combo
    )
  ) %>%
  group_by(grupo) %>%
  mutate(
    pct_grupo   = round(100 * n / sum(n), 1),
    pct_label   = if_else(
      pct_grupo %% 1 == 0,
      as.character(as.integer(pct_grupo)),
      formatC(pct_grupo, format = "f", digits = 1)
    ),
    etiqueta    = sprintf("%d\n(%s%%)", n, pct_label),
    total_grupo = sum(n)
  ) %>%
  ungroup() %>%
  mutate(
    strip_lab = case_when(
      grupo == "Genotipos gH únicos" ~
        sprintf("Genotipos gH únicos\n(n=%d)",
                unique(total_grupo[grupo == "Genotipos gH únicos"])),
      grupo == "Genotipos gH mixtos" ~
        sprintf("Genotipos gH mixtos\n(n=%d)",
                unique(total_grupo[grupo == "Genotipos gH mixtos"]))
    ),
    strip_lab = factor(
      strip_lab,
      levels = c(
        sprintf("Genotipos gH únicos\n(n=%d)",
                unique(total_grupo[grupo == "Genotipos gH únicos"])),
        sprintf("Genotipos gH mixtos\n(n=%d)",
                unique(total_grupo[grupo == "Genotipos gH mixtos"]))
      )
    )
  )


# 3) Gráfico
distribuciones_gH <- ggplot(desglose_all_gh,
                            aes(x = combo_clean, y = n, fill = grupo)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = etiqueta),
            vjust = -0.4, size = 3.3, lineheight = 0.9) +
  facet_grid(~ strip_lab, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("Genotipos gH únicos" = "#DDA0DD",
                               "Genotipos gH mixtos" = "#CD69C9")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(x = "", y = "Número de muestras") +
  theme_classic(base_size = 13) +
  theme(
    plot.margin      = margin(t = 18, r = 15, b = 10, l = 10),
    legend.position  = "none",
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.line        = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks       = element_line(linewidth = 0.3, colour = "grey50"),
    axis.title.y     = element_text(size = 11),
    strip.text       = element_text(face = "bold"),
    panel.spacing.x  = unit(2.5, "lines")
  )

distribuciones_gH

