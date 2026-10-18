# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.7.05. Figura 34. Distribución del número de episodios por paciente según genotipos mixtos
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 34. Distribución del número de episodios por paciente según genotipos mixtos ----
# ── 1. Preparar datos ────────────────────────────────────────────────────────
df_dist <- pacientes %>%
  filter(!is.na(nepi_tot), !is.na(gbghm)) %>%
  mutate(
    genotipo = factor(ifelse(gbghm == 1, "Mixto", "Único"),
                      levels = c("Único", "Mixto")),
    cat_ep = factor(case_when(
      nepi_tot == 1 ~ "1",
      nepi_tot == 2 ~ "2",
      nepi_tot == 3 ~ "3",
      nepi_tot == 4 ~ "4",
      nepi_tot >= 5 ~ "\u22655"
    ), levels = c("1", "2", "3", "4", "\u22655"))
  )

# Conteos por grupo y categoría
df_pir <- df_dist %>%
  count(genotipo, cat_ep) %>%
  group_by(genotipo) %>%
  mutate(
    pct       = 100 * n / sum(n),
    pct_signo = ifelse(genotipo == "Único", -pct, pct),
    etiqueta  = sprintf("%d (%.1f%%)", n, pct) %>%
      sub("\\.", ",", .)
  ) %>%
  ungroup()

# Mostrar conteos por consola para confirmar
cat("--- Distribución por número de episodios ---\n")
print(df_pir)

# ── 2. Paleta coherente con el resto del documento ──────────────────────────
cols_pir <- c("Único" = "#4CAF7D", "Mixto" = "#E07070")

n_unico <- sum(df_dist$genotipo == "Único")
n_mixto <- sum(df_dist$genotipo == "Mixto")

# ── 3. Figura pirámide ──────────────────────────────────────────────────────
fig_pir <- ggplot(df_pir, aes(x = pct_signo, y = cat_ep, fill = genotipo)) +
  geom_col(width = 0.7, color = "grey25", linewidth = 0.3) +
  geom_text(
    aes(label = etiqueta,
        hjust = ifelse(genotipo == "Único", 1.05, -0.05)),
    size = 3.4, color = "grey15"
  ) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.4) +
  scale_x_continuous(
    name   = "Porcentaje de pacientes (%)",
    breaks = seq(-50, 50, 10),
    labels = function(x) paste0(abs(x), "%"),
    limits = c(-55, 55),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = cols_pir,
                    labels = c(sprintf("Genotipos \u00fanicos (n = %d)", n_unico),
                               sprintf("Genotipos mixtos (n = %d)", n_mixto))) +
  labs(
    y    = "N episodios por paciente",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "bottom",
    legend.text        = element_text(size = 11),
    legend.key.size    = unit(12, "pt"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
    axis.text.y        = element_text(size = 12, color = "grey15", face = "bold"),
    axis.text.x        = element_text(size = 10, color = "grey40"),
    axis.title.x       = element_text(size = 11, color = "grey20",
                                      margin = margin(t = 6)),
    axis.title.y       = element_text(size = 11, color = "grey20",
                                      margin = margin(r = 8)),
    plot.margin        = margin(10, 14, 8, 8)
  )

fig_pir

# [guardado antiguo] ggsave("~/Desktop/Figuras_tesis/figura_R35bis_distribucion_episodios.png",
# [guardado antiguo]        fig_pir, width = 6, height = 3, dpi = 320, bg = "white")


df_test <- pacientes %>%
  filter(!is.na(nepi_tot), !is.na(gbghm)) %>%
  mutate(genotipo = factor(ifelse(gbghm == 1, "Mixto", "Único"),
                           levels = c("Único", "Mixto")))

# Mann-Whitney
test_mw <- wilcox.test(nepi_tot ~ genotipo, data = df_test, exact = FALSE)
cat("Mann-Whitney p =", format.pval(test_mw$p.value, digits = 3), "\n")

# Chi-cuadrado de la distribución categórica (≥4 vs <4 episodios)
df_test <- df_test %>% mutate(cat_4 = ifelse(nepi_tot >= 4, "≥4", "<4"))
test_chi <- chisq.test(table(df_test$genotipo, df_test$cat_4))
cat("Chi-cuadrado p =", format.pval(test_chi$p.value, digits = 3), "\n")
print(table(df_test$genotipo, df_test$cat_4))


# ── 1. AUCCV acumulada por paciente ──────────────────────────────────────────
# Limpieza previa de columnas si existen de ejecuciones anteriores
pacientes <- pacientes %>%
  select(-any_of(c("AUCCV_AC", "n_episodios_AUC",
                   "AUCCV_AC.x", "AUCCV_AC.y",
                   "n_episodios_AUC.x", "n_episodios_AUC.y")))

auc_paciente <- dnaemias %>%
  group_by(NHC) %>%
  summarise(
    AUCCV_AC        = sum(AUCCV, na.rm = TRUE),
    n_episodios_AUC = sum(!is.na(AUCCV)),
    .groups         = "drop"
  )

pacientes <- pacientes %>%
  left_join(auc_paciente, by = "NHC")

# ── 2. Datos para comparaciones ──────────────────────────────────────────────
df_comp <- pacientes %>%
  filter(!is.na(n_episodios_AUC), n_episodios_AUC > 0) %>%
  mutate(
    gbm_lab   = factor(ifelse(gbm   == 1, "Mixto", "Único"),
                       levels = c("Único", "Mixto")),
    ghm_lab   = factor(ifelse(ghm   == 1, "Mixto", "Único"),
                       levels = c("Único", "Mixto")),
    gbghm_lab = factor(ifelse(gbghm == 1, "Mixto", "Único"),
                       levels = c("Único", "Mixto"))
  )

# ── 3. Paletas por glucoproteína ─────────────────────────────────────────────
# gB: azul claro / azul oscuro
cols_fill_gb <- c("Único" = "#AED6F1", "Mixto" = "#1A5276")
cols_jit_gb  <- c("Único" = "#2980B9", "Mixto" = "#0D2B45")

# gH: rosa / rosa oscuro
cols_fill_gh <- c("Único" = "pink",    "Mixto" = "hotpink3")
cols_jit_gh  <- c("Único" = "pink3",   "Mixto" = "palevioletred4")

# gB+gH: verde / rojo
cols_fill_gbgh <- c("Único" = "#4CAF7D", "Mixto" = "#E07070")
cols_jit_gbgh  <- c("Único" = "#2E7D52", "Mixto" = "#B04040")

# ── 4. Función para construir cada panel ─────────────────────────────────────
# Límite Y fijo en 2000; asteriscos en vez de p
make_panel_paciente <- function(data, var_grupo, titulo,
                                cols_fill, cols_jit,
                                ylim_max = 1000) {
  
  d <- data %>%
    select(grupo = all_of(var_grupo), AUCCV_AC) %>%
    filter(!is.na(AUCCV_AC), !is.na(grupo))
  
  # Test sobre datos completos (sin filtrar por cap visual)
  stat <- d %>%
    rstatix::wilcox_test(AUCCV_AC ~ grupo) %>%
    rstatix::add_significance("p") %>%
    mutate(y.position = ylim_max * 0.95)
  
  # Datos con cap visual aplicado
  d_vis <- d %>% filter(AUCCV_AC <= ylim_max)
  
  ggplot(d_vis, aes(x = grupo, y = AUCCV_AC, fill = grupo)) +
    geom_boxplot(
      width         = 0.55,
      color         = "grey25",
      linewidth     = 0.45,
      outlier.shape = 21,
      outlier.size  = 1.8,
      outlier.fill  = "white",
      outlier.color = "grey50",
      outlier.alpha = 0.7
    ) +
    geom_jitter(
      aes(color = grupo),
      width = 0.12, size = 1.1, alpha = 0.45,
      show.legend = FALSE
    ) +
    scale_fill_manual(values  = cols_fill) +
    scale_color_manual(values = cols_jit) +
    stat_pvalue_manual(
      stat,
      label        = "p.signif",
      tip.length   = 0.01,
      bracket.size = 0.45,
      size         = 5,
      inherit.aes  = FALSE
    ) +
    labs(x = NULL,
         y = "AUCCV acumulada (UI/mL\u00b7d\u00eda)",
         title = titulo) +
    scale_y_continuous(
      limits = c(NA, ylim_max),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position    = "none",
      panel.border       = element_blank(),
      axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
      axis.line.y        = element_line(color = "grey60", linewidth = 0.4),
      axis.ticks         = element_line(color = "grey60", linewidth = 0.3),
      axis.text.x        = element_text(size = 11, color = "grey25"),
      axis.text.y        = element_text(size = 10, color = "grey30"),
      axis.title.y       = element_text(size = 10, color = "grey20"),
      plot.title         = element_text(size = 13, color = "grey15",
                                        face = "bold", hjust = 0.5),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.margin        = margin(8, 10, 6, 6)
    )
}

# ── 5. Construir los 3 paneles ───────────────────────────────────────────────
p_gb <- make_panel_paciente(
  df_comp, "gbm_lab", "gB",
  cols_fill = cols_fill_gb, cols_jit = cols_jit_gb
)

p_gh <- make_panel_paciente(
  df_comp, "ghm_lab", "gH",
  cols_fill = cols_fill_gh, cols_jit = cols_jit_gh
)

p_gbgh <- make_panel_paciente(
  df_comp, "gbghm_lab", "gB/gH",
  cols_fill = cols_fill_gbgh, cols_jit = cols_jit_gbgh
)

# ── 6. Combinar y mostrar ────────────────────────────────────────────────────
figura_AUCCV_AC <- p_gb | p_gh | p_gbgh

figura_AUCCV_AC

# ── 7. Guardar ───────────────────────────────────────────────────────────────
# ggsave("~/Desktop/Figuras_tesis/figura_R22_AUCCV_AC_mixto.png",
#        figura_AUCCV_AC,
#        width = 10, height = 4, dpi = 320, bg = "white")

guardar_fig(fig_pir, "Figura_34", ancho = 15.2, alto = 7.6)

