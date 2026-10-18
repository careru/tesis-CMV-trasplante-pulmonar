# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.4.02. Figura 18. Muestras con DNAemia cuantificable por año y categoría de carga viral
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 18. Muestras con DNAemia cuantificable por año y categoría de carga viral ----
DTA_CV_TXP <- read_dta(f_cargas)

# Años válidos
anios_todos <- DTA_CV_TXP %>%
  filter(!is.na(Freg)) %>%
  mutate(anio = as.numeric(format(as.Date(Freg), "%Y"))) %>%
  distinct(anio) %>%
  filter(!anio %in% c(2008, 2023)) %>%
  arrange(anio) %>%
  pull(anio)

# Limpieza + categorías CV
DTA_CV_TXP_cv2 <- DTA_CV_TXP %>%
  mutate(
    Resultado_interpretado = na_if(Resultado_interpretado, ""),
    anio = as.numeric(format(as.Date(Freg), "%Y"))
  ) %>%
  filter(
    !is.na(Freg),
    !is.na(Resultado_interpretado),
    !anio %in% c(2008, 2023)
  ) %>%
  mutate(
    anio = factor(anio, levels = anios_todos),
    cat_cv = case_when(
      Resultado_interpretado == "< LDC"                 ~ "CV < LDC",
      Resultado_interpretado == "POSITIVO" & CV < 1000  ~ "CV positivas < 1.000 UI/mL",
      Resultado_interpretado == "POSITIVO" & CV >= 1000 ~ "CV positivas ≥ 1.000 UI/mL",
      TRUE                                              ~ NA_character_
    )
  ) %>%
  filter(!is.na(cat_cv))

# Conteo
df_cv_cat <- DTA_CV_TXP_cv2 %>%
  group_by(anio, cat_cv) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(cat_cv = factor(cat_cv,
                         levels = c("CV < LDC",
                                    "CV positivas < 1.000 UI/mL",
                                    "CV positivas ≥ 1.000 UI/mL")
  ))


graph_cv_cat <- ggplot(df_cv_cat,
                       aes(x = anio, y = n, fill = cat_cv)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  
  scale_fill_manual(
    values = c(
      "CV < LDC"                     = "#FFC0CB",
      "CV positivas < 1.000 UI/mL"   = "#CD5555",
      "CV positivas ≥ 1.000 UI/mL"   = "#8B0000"
    ),
    name = "",
    labels = c(
      "CV < LDC",
      "CV positivas < 1.000 UI/mL",
      "CV positivas ≥ 1.000 UI/mL"
    )
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(
    x = "",
    y = "Número de muestras"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x     = element_text(angle = 45, hjust = 1),
    axis.line       = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks      = element_line(linewidth = 0.3, colour = "grey50"),
    legend.key.size = unit(9, "pt")
  )

graph_cv_cat

# Tabla con N y %
tabla_cv <- df_cv_cat %>%
  group_by(anio) %>%
  mutate(
    total_anio = sum(n),
    pct = round(100 * n / total_anio, 1)
  ) %>%
  ungroup() %>%
  arrange(anio, cat_cv)

# Mostrar en consola
tabla_cv


tabla_cv_wide <- tabla_cv %>%
  mutate(N_pct = paste0(n, " (", pct, "%)")) %>%
  select(anio, cat_cv, N_pct) %>%
  tidyr::pivot_wider(names_from = cat_cv, values_from = N_pct)

# Mostrar en consola
tabla_cv_wide

# Preparar datos para el gráfico (en %)
df_cv_cat <- DTA_CV_TXP_cv2 %>%
  group_by(anio, cat_cv) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(anio) %>%
  mutate(
    total_anio = sum(n),
    pct = round(100 * n / total_anio, 0)
  ) %>%
  ungroup()

pos_dodge <- position_dodge(width = 0.85)

graph_cv_cat <- ggplot(df_cv_cat,
                       aes(x = anio, y = pct, fill = cat_cv)) +
  geom_col(position = pos_dodge, width = 0.8) +
  
  # Etiquetas SOLO para CV ≥ 1.000 UI/mL (las rojas oscuras, a la derecha)
  geom_text(
    aes(
      label = ifelse(cat_cv == "CV ≥ 1.000 UI/mL",
                     paste0(round(pct), "%"),
                     "")
    ),
    position = pos_dodge,
    vjust = -0.5,
    size = 4
  ) +
  
  scale_fill_manual(
    values = c(
      "CV < 1.000 UI/mL" = "#FFAEB9",
      "CV ≥ 1.000 UI/mL" = "#8B0000"
    ),
    name = "",
    labels = c(
      "CV < 1.000 UI/mL",
      "CV ≥ 1.000 UI/mL"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.10)),
    limits = c(0, 100),
    labels = scales::label_number(accuracy = 1)
  ) +
  labs(
    x = "",
    y = "% muestras CV detectable"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    legend.text     = element_text(size = 12),
    axis.text       = element_text(size = 13),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 13),
    axis.title      = element_text(size = 14),
    axis.line       = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks      = element_line(linewidth = 0.3, colour = "grey50"),
    legend.key.size = unit(9, "pt")
  )

graph_cv_cat

guardar_fig(graph_cv_cat, "Figura_18", ancho = 18, alto = 12)

