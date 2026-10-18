# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.05. Anexo: Porcentaje de positividad por año (no incluida)
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Porcentaje de positividad por año (no incluida) ----
DTA_CV_TXP <- read_dta(f_cargas)

# Años totales que quieres mostrar (sin 2008 ni 2023)
anios_todos <- DTA_CV_TXP %>%
  filter(!is.na(Freg)) %>%
  mutate(anio = as.numeric(format(as.Date(Freg), "%Y"))) %>%
  distinct(anio) %>%
  filter(!anio %in% c(2008, 2023)) %>%
  arrange(anio) %>%
  pull(anio)


DTA_CV_TXP_cv <- DTA_CV_TXP %>%
  mutate(
    Resultado_interpretado = na_if(Resultado_interpretado, ""),  # ← CLAVE
    anio = as.numeric(format(as.Date(Freg), "%Y"))
  ) %>%
  filter(
    !is.na(Freg),
    !is.na(Resultado_interpretado),   # ahora sí funciona
    !anio %in% c(2008, 2023)
  ) %>%
  mutate(
    Resultado_interpretado = factor(Resultado_interpretado,
                                    levels = c("NEGATIVO", "< LDC", "POSITIVO")
    ),
    anio = factor(anio, levels = anios_todos)
  )

# Conteo por año y categoría (rellenando con 0 donde falte)
df_counts <- DTA_CV_TXP_cv %>%
  group_by(anio, Resultado_interpretado) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(
    anio,
    Resultado_interpretado,
    fill = list(n = 0)
  )

# Cálculo del % de positividad por año
df_pct <- df_counts %>%
  group_by(anio) %>%
  summarise(
    total     = sum(n),
    positivos = sum(n[Resultado_interpretado == "POSITIVO"]),
    pct_pos   = ifelse(total > 0, 100 * positivos / total, NA_real_)
  ) %>%
  ungroup() %>%
  filter(!is.na(pct_pos))

max_muestras <- df_counts %>%
  group_by(anio) %>%
  summarise(total = sum(n)) %>%
  pull(total) %>%
  max()

pct_min <- 0
pct_max <- 25

df_pct <- df_pct %>%
  mutate(
    pct_scaled = (pct_pos - pct_min) / (pct_max - pct_min) * max_muestras
  )


graph_cv_anual_linea <- ggplot(df_counts,
                               aes(x = anio, y = n,
                                   fill = Resultado_interpretado)) +
  geom_bar(stat = "identity", width = 0.85) +
  
  geom_line(
    data = df_pct,
    aes(x = anio, y = pct_scaled, group = 1),
    inherit.aes = FALSE,
    color = "#CD0000",
    linewidth = 0.8,
    alpha = 0.4
  ) +
  geom_text(
    data = df_pct,
    aes(
      x = anio,
      y = pct_scaled,
      label = paste0(round(pct_pos, 1), "%")
    ),
    inherit.aes = FALSE,
    vjust = -1.1,
    color = "firebrick4",
    fontface = "bold",
    size = 3
  ) +
  
  scale_fill_manual(
    values = c(
      "POSITIVO" = "#CD0000",
      "< LDC"    = "#A2B5CD",
      "NEGATIVO" = "#CAE1FF"
    ),
    labels = c(
      "NEGATIVO" = "CV negativas",
      "< LDC"    = "CV < LDC",
      "POSITIVO" = "CV positivas cuantificables"
    ),
    name = ""
  ) +
  
  scale_y_continuous(
    name = "Número de muestras",
    expand = expansion(mult = c(0, 0.10)),
    sec.axis = sec_axis(
      ~ (. * (pct_max - pct_min) / max_muestras) + pct_min,
      name = "% positividad",
      breaks = seq(0, 25, 5)
    )
  ) +
  labs(x = "") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x     = element_text(angle = 45, hjust = 1),
    axis.line       = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks      = element_line(linewidth = 0.3, colour = "grey50"),
    legend.key.size = unit(9, "pt")
  )

graph_cv_anual_linea

# Conteo por año y categoría
df_counts <- DTA_CV_TXP_cv %>%
  group_by(anio, Resultado_interpretado) %>%
  summarise(n = n(), .groups = "drop")

# Cálculo del % de positividad por año
df_pct <- df_counts %>%
  group_by(anio) %>%
  summarise(
    total     = sum(n),
    positivos = sum(n[Resultado_interpretado == "POSITIVO"]),
    pct_pos   = 100 * positivos / total
  ) %>%
  ungroup() %>%
  filter(!is.na(pct_pos))   # excluir NAs de pct_pos


max_muestras <- df_counts %>%
  group_by(anio) %>%
  summarise(total = sum(n)) %>%
  pull(total) %>%
  max()

pct_min <- 0
pct_max <- 25   # nuevo rango 0–25%

df_pct <- df_pct %>%
  mutate(
    pct_scaled = (pct_pos - pct_min) / (pct_max - pct_min) * max_muestras
  )


graph_cv_anual_linea <- ggplot(df_counts,
                               aes(x = anio, y = n,
                                   fill = Resultado_interpretado)) +
  geom_bar(stat = "identity", width = 0.85) +
  
  # Línea roja oscura y semi-transparente
  geom_line(
    data = df_pct,
    aes(x = anio, y = pct_scaled, group = 1),
    inherit.aes = FALSE,
    color = "#CD0000",
    linewidth = 0.8,
    alpha = 0.4
  ) +
  # Puntos pequeños y semi-transparente
  geom_point(
    data = df_pct,
    aes(x = anio, y = pct_scaled),
    inherit.aes = FALSE,
    color = "firebrick4",
    size = 0.0,      # puntos más pequeños
    alpha = 0.0
  ) +
  # Texto del % más separado y en rojo más oscuro
  geom_text(
    data = df_pct,
    aes(
      x = anio,
      y = pct_scaled,
      label = paste0(round(pct_pos, 1), "%")
    ),
    inherit.aes = FALSE,
    vjust = -1.1,          # más espacio entre punto y número
    color = "firebrick4",  # rojo más oscuro
    fontface = "bold",
    size = 3
  ) +
  
  scale_fill_manual(
    values = c(
      "POSITIVO" = "#CD0000",
      "< LDC"    = "#A2B5CD",
      "NEGATIVO" = "#CAE1FF"
    ),
    labels = c(
      "NEGATIVO" = "CV negativas",
      "< LDC"    = "CV < LDC",
      "POSITIVO" = "CV positivas cuantificables"
    ),
    name = ""
  ) +
  
  scale_y_continuous(
    name = "Número de muestras",
    expand = expansion(mult = c(0, 0.10)),
    sec.axis = sec_axis(
      ~ (. * (pct_max - pct_min) / max_muestras) + pct_min,
      name = "% positividad",
      breaks = seq(0, 25, 5)
    )
  ) +
  labs(x = "") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x     = element_text(angle = 45, hjust = 1),
    axis.line       = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks      = element_line(linewidth = 0.3, colour = "grey50"),
    legend.key.size = unit(9, "pt")
  )

graph_cv_anual_linea

