# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.5.04. Datos para el texto (9.5.2): frecuencias de mixtos por seroestatus
# Extracto del script completo (ver README.md del estudio)

## ---- Datos para el texto (9.5.2): frecuencias de mixtos por seroestatus ----
# 1) Preparar datos (mismos 3 serDR que estabas usando)
df <- dnaemias %>%
  filter(serDR %in% c("D-/R+", "D+/R+", "D+/R-")) %>%
  mutate(
    gb_mixto = (gbmix == 1),
    
    # gH mixto: si existe 'ghmix' lo usa; si no, lo deja en NA
    gh_mixto = if ("ghmix" %in% names(.)) (ghmix == 1) else NA,
    
    # gB+gH mixto: si existe 'gbghmix' lo usa; si no, si existe ghmix calcula ambos mixtos
    gbgh_mixto = if ("gbghmix" %in% names(.)) {
      (gbghmix == 1)
    } else if ("ghmix" %in% names(.)) {
      (gbmix == 1 & ghmix == 1)
    } else {
      NA
    }
  )

# 2) Función resumen por grupo
resumen_por_grupo <- function(data) {
  data %>%
    summarise(
      episodios = n(),
      gb_mix_n  = sum(gb_mixto, na.rm = TRUE),
      gb_mix_p  = round(100 * gb_mix_n / episodios, 1),
      
      gh_mix_n  = if (all(is.na(gh_mixto))) NA_integer_ else sum(gh_mixto, na.rm = TRUE),
      gh_mix_p  = if (all(is.na(gh_mixto))) NA_real_    else round(100 * gh_mix_n / episodios, 1),
      
      gbgh_mix_n = if (all(is.na(gbgh_mixto))) NA_integer_ else sum(gbgh_mixto, na.rm = TRUE),
      gbgh_mix_p = if (all(is.na(gbgh_mixto))) NA_real_    else round(100 * gbgh_mix_n / episodios, 1)
    )
}

# 3) Tabla por serDR + fila Global
tab <- df %>%
  group_by(serDR) %>%
  group_modify(~ resumen_por_grupo(.x)) %>%
  ungroup()

tab_global <- resumen_por_grupo(df) %>%
  mutate(serDR = "Global")

tab_final <- bind_rows(tab, tab_global) %>%
  relocate(serDR)

# 4) Formato "n (%)" y flextable final
tab_final_fmt <- tab_final %>%
  transmute(
    serDR,
    `Episodios (n)` = episodios,
    `gB mixto`      = ifelse(is.na(gb_mix_n),  NA, sprintf("%d (%.1f%%)", gb_mix_n, gb_mix_p)),
    `gH mixto`      = ifelse(is.na(gh_mix_n),  NA, sprintf("%d (%.1f%%)", gh_mix_n, gh_mix_p)),
    `gB+gH mixto`   = ifelse(is.na(gbgh_mix_n), NA, sprintf("%d (%.1f%%)", gbgh_mix_n, gbgh_mix_p))
  )

ft <- flextable(tab_final_fmt) %>%
  autofit()

ft


# 1) Preparar datos
df <- dnaemias %>%
  filter(serDR %in% c("D-/R+", "D+/R+", "D+/R-")) %>%
  mutate(
    gb_mixto = (gbmix == 1),
    gh_mixto = if ("ghmix" %in% names(.)) (ghmix == 1) else NA,
    gbgh_mixto = if ("gbghmix" %in% names(.)) {
      (gbghmix == 1)
    } else if ("ghmix" %in% names(.)) {
      (gbmix == 1 & ghmix == 1)
    } else {
      NA
    }
  )

# 2) Calcular frecuencias por grupo
tab_fig <- df %>%
  group_by(serDR) %>%
  summarise(
    episodios = n(),
    `gB mixto` = round(100 * sum(gb_mixto, na.rm = TRUE) / episodios, 1),
    `gH mixto` = if (all(is.na(gh_mixto))) NA_real_ 
    else round(100 * sum(gh_mixto, na.rm = TRUE) / episodios, 1),
    `gB/gH mixto` = if (all(is.na(gbgh_mixto))) NA_real_ 
    else round(100 * sum(gbgh_mixto, na.rm = TRUE) / episodios, 1)
  ) %>%
  ungroup()

# 3) Convertir a formato largo para ggplot
tab_long <- tab_fig %>%
  pivot_longer(
    cols = c(`gB mixto`, `gH mixto`, `gB/gH mixto`),
    names_to = "tipo_genotipo",
    values_to = "porcentaje"
  ) %>%
  # Ordenar grupos serológicos
  mutate(
    serDR = factor(serDR, levels = c("D+/R-", "D+/R+", "D-/R+")),
    tipo_genotipo = factor(tipo_genotipo, levels = c("gB mixto", "gH mixto", "gB/gH mixto"))
  )

# 4) Crear figura
figura <- ggplot(tab_long, aes(x = serDR, y = porcentaje, fill = tipo_genotipo)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    aes(label = sprintf("%.1f%%", porcentaje)),
    position = position_dodge(width = 0.8),
    vjust = -0.5,
    size = 3.5
  ) +
  scale_fill_manual(
    values = c("gB mixto" = "#E69F00", "gH mixto" = "#56B4E9", "gB/gH mixto" = "#009E73"),
    name = "Tipo de genotipo"
  ) +
  labs(
    x = "Perfil serológico D/R",
    y = "Frecuencia (%)",
    title = "Frecuencia de genotipos mixtos por episodio según perfil serológico"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.position = "bottom",
    panel.grid.major.x = element_blank()
  ) +
  ylim(0, 45)

# Mostrar figura
print(figura)

# Opcional: guardar figura
# ggsave("figura_genotipos_mixtos.png", figura, width = 8, height = 6, dpi = 300)

