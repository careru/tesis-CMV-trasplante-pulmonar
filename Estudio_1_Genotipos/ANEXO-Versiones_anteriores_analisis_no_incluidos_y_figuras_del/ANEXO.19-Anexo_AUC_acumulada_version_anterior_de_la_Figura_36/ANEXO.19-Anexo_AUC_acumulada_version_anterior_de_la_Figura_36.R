# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.19. Anexo: AUC acumulada, versión anterior de la Figura 36
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: AUC acumulada, versión anterior de la Figura 36 ----
# Distribución del número total de episodios por paciente según mixto/único

# Asegurar variable gbghm coherente con la usada en el resto del análisis
pacientes %>%
  mutate(genotipo = factor(ifelse(gbghm == 1, "Mixto", "Único"),
                           levels = c("Único", "Mixto"))) %>%
  group_by(genotipo) %>%
  summarise(
    n_pac      = n(),
    `1`        = sum(nepi_tot == 1, na.rm = TRUE),
    `2`        = sum(nepi_tot == 2, na.rm = TRUE),
    `3`        = sum(nepi_tot == 3, na.rm = TRUE),
    `4-5`      = sum(nepi_tot %in% 4:5, na.rm = TRUE),
    `>=6`      = sum(nepi_tot >= 6, na.rm = TRUE),
    mediana    = median(nepi_tot, na.rm = TRUE),
    q1         = quantile(nepi_tot, 0.25, na.rm = TRUE),
    q3         = quantile(nepi_tot, 0.75, na.rm = TRUE),
    rango_min  = min(nepi_tot, na.rm = TRUE),
    rango_max  = max(nepi_tot, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  print(width = Inf)

