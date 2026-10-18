# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.06. Anexo: Cargas positivas por año apiladas (no incluida)
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Cargas positivas por año apiladas (no incluida) ----
DTA_CV_TXP <- read_dta(f_cargas)

# Años válidos (sin 2008 ni 2023)
anios_todos <- DTA_CV_TXP %>%
  filter(!is.na(Freg)) %>%
  mutate(anio = as.numeric(format(as.Date(Freg), "%Y"))) %>%
  distinct(anio) %>%
  filter(!anio %in% c(2008, 2023)) %>%
  arrange(anio) %>%
  pull(anio)

# Datos limpios + solo CV POSITIVAS y categorización por umbral 1.000 UI/mL
DTA_CV_TXP_cv2 <- DTA_CV_TXP %>%
  mutate(
    Resultado_interpretado = na_if(Resultado_interpretado, ""),
    anio = as.numeric(format(as.Date(Freg), "%Y"))
  ) %>%
  filter(
    !is.na(Freg),
    !is.na(Resultado_interpretado),
    Resultado_interpretado == "POSITIVO",   # solo CV positivas
    !anio %in% c(2008, 2023)
  ) %>%
  mutate(
    anio = factor(anio, levels = anios_todos),
    cat_cv = case_when(
      CV < 1000  ~ "CV positivas < 1.000 UI/mL",
      CV >= 1000 ~ "CV positivas ≥ 1.000 UI/mL"
    )
  )

# Conteo por año y categoría
df_cv_cat <- DTA_CV_TXP_cv2 %>%
  group_by(anio, cat_cv) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    cat_cv = factor(
      cat_cv,
      levels = c(
        "CV positivas < 1.000 UI/mL",
        "CV positivas ≥ 1.000 UI/mL"
      )
    )
  )


graph_cv_cat_stacked <- ggplot(df_cv_cat,
                               aes(x = anio, y = n, fill = cat_cv)) +
  geom_col(width = 0.85) +   # stacked por defecto
  
  scale_fill_manual(
    values = c(
      "CV positivas < 1.000 UI/mL" = "#F4A6A6",
      "CV positivas ≥ 1.000 UI/mL" = "#CD0000"
    ),
    name = "",
    labels = c(
      "CV positivas < 1.000 UI/mL",
      "CV positivas ≥ 1.000 UI/mL"
    )
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(
    x = "",
    y = "Número de muestras con CV positiva"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x     = element_text(angle = 45, hjust = 1),
    axis.line       = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks      = element_line(linewidth = 0.3, colour = "grey50"),
    legend.key.size = unit(9, "pt")
  )

graph_cv_cat_stacked

