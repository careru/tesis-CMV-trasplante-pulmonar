# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.23. Anexo: Versión anterior de la Figura 33 (tres glucoproteínas)
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Versión anterior de la Figura 33 (tres glucoproteínas) ----
# ── Análisis B: clasificación por patrón del paciente ─────────────────────
df_primer_B <- dnaemias %>%
  filter(primer_mil == 1) %>%
  mutate(
    Ftx              = as.Date(Ftx),
    Fini             = as.Date(Fini),
    dias_tx_a_inicio = as.numeric(Fini - Ftx),
    meses_tx_inicio  = dias_tx_a_inicio / 30.44
  ) %>%
  # Sustituir las variables de episodio por las de paciente
  select(-any_of(c("gbm", "ghm", "gbghm"))) %>%
  left_join(pacientes %>% select(NHC, gbm, ghm, gbghm), by = "NHC") %>%
  mutate(
    riesgo = case_when(
      serDR == "D+/R-"                ~ "Alto riesgo (D+/R-)",
      serDR %in% c("D+/R+", "D-/R+")  ~ "Riesgo intermedio (R+)",
      TRUE                             ~ NA_character_
    ),
    riesgo = factor(riesgo,
                    levels = c("Alto riesgo (D+/R-)",
                               "Riesgo intermedio (R+)"))
  ) %>%
  filter(!is.na(riesgo))

cat(sprintf("── N primeros episodios analizados: %d ──\n", nrow(df_primer_B)))
print(table(df_primer_B$riesgo))

# ── Función ────────────────────────────────────────────────────────────────
resumen_y_test_B <- function(data, var_pac, etiqueta_glic) {
  d <- data %>%
    filter(!is.na(.data[[var_pac]])) %>%
    mutate(genotipo = factor(ifelse(.data[[var_pac]] == 1, "Mixto", "Único"),
                             levels = c("Único", "Mixto")))
  
  resumen <- d %>%
    group_by(genotipo) %>%
    summarise(
      n       = n(),
      mediana = round(median(meses_tx_inicio, na.rm = TRUE), 1),
      q1      = round(quantile(meses_tx_inicio, 0.25, na.rm = TRUE), 1),
      q3      = round(quantile(meses_tx_inicio, 0.75, na.rm = TRUE), 1),
      .groups = "drop"
    )
  
  cat(sprintf("\n--- %s ---\n", etiqueta_glic))
  print(resumen)
  
  if (length(unique(d$genotipo)) == 2 && all(table(d$genotipo) > 0)) {
    test <- wilcox.test(meses_tx_inicio ~ genotipo, data = d, exact = FALSE)
    cat(sprintf("p = %.4f\n", test$p.value))
  } else {
    cat("Solo un grupo disponible, no se realiza test.\n")
  }
}

# ── 1. D+/R- ───────────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════")
cat("\n  GRUPO 1: ALTO RIESGO (D+/R-)")
cat("\n  Clasificación por PATRÓN DEL PACIENTE")
cat("\n══════════════════════════════════════════════\n")

df_alto_B <- df_primer_B %>% filter(riesgo == "Alto riesgo (D+/R-)")
cat(sprintf("N = %d primeros episodios\n", nrow(df_alto_B)))

resumen_y_test_B(df_alto_B, "gbghm", "gB+gH paciente (D+/R-)")
resumen_y_test_B(df_alto_B, "gbm",   "gB paciente (D+/R-)")
resumen_y_test_B(df_alto_B, "ghm",   "gH paciente (D+/R-)")

# ── 2. R+ ──────────────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════")
cat("\n  GRUPO 2: RIESGO INTERMEDIO (R+)")
cat("\n  Clasificación por PATRÓN DEL PACIENTE")
cat("\n══════════════════════════════════════════════\n")

df_inter_B <- df_primer_B %>% filter(riesgo == "Riesgo intermedio (R+)")
cat(sprintf("N = %d primeros episodios\n", nrow(df_inter_B)))

resumen_y_test_B(df_inter_B, "gbghm", "gB+gH paciente (R+)")
resumen_y_test_B(df_inter_B, "gbm",   "gB paciente (R+)")
resumen_y_test_B(df_inter_B, "ghm",   "gH paciente (R+)")


# ── 1. Preparar datos ──────────────────────────────────────────────────────
df_fig <- dnaemias %>%
  filter(primer_mil== 1) %>%
  mutate(
    Ftx              = as.Date(Ftx),
    Fini             = as.Date(Fini),
    meses_tx_inicio  = as.numeric(Fini - Ftx) / 30.44
  ) %>%
  select(-any_of(c("gbm", "ghm", "gbghm"))) %>%
  left_join(pacientes %>% select(NHC, gbm, ghm, gbghm), by = "NHC") %>%
  mutate(
    riesgo = case_when(
      serDR == "D+/R-"                ~ "Alto riesgo (D+/R-)",
      serDR %in% c("D+/R+", "D-/R+")  ~ "Riesgo intermedio (R+)",
      TRUE                             ~ NA_character_
    ),
    riesgo = factor(riesgo,
                    levels = c("Alto riesgo (D+/R-)",
                               "Riesgo intermedio (R+)"))
  ) %>%
  filter(!is.na(riesgo), !is.na(meses_tx_inicio))

# ── 2. Pasar a formato largo ───────────────────────────────────────────────
df_long <- df_fig %>%
  select(NHC, riesgo, meses_tx_inicio, gbm, ghm, gbghm) %>%
  pivot_longer(
    cols      = c(gbm, ghm, gbghm),
    names_to  = "glucoproteina",
    values_to = "mixto"
  ) %>%
  filter(!is.na(mixto)) %>%
  mutate(
    glucoproteina = factor(glucoproteina,
                           levels = c("gbm", "ghm", "gbghm"),
                           labels = c("gB", "gH", "gB+gH")),
    genotipo      = factor(ifelse(mixto == 1, "Mixto", "Único"),
                           levels = c("Único", "Mixto")),
    # Combinación glucoproteína × genotipo para asignar color
    cat_color = paste(glucoproteina, genotipo, sep = "_"),
    cat_color = factor(cat_color,
                       levels = c("gB_Único",   "gB_Mixto",
                                  "gH_Único",   "gH_Mixto",
                                  "gB+gH_Único","gB+gH_Mixto"))
  )

# ── 3. N por panel ─────────────────────────────────────────────────────────
n_alto  <- df_fig %>% filter(riesgo == "Alto riesgo (D+/R-)") %>% nrow()
n_inter <- df_fig %>% filter(riesgo == "Riesgo intermedio (R+)") %>% nrow()

etiquetas_riesgo <- c(
  "Alto riesgo (D+/R-)"     = sprintf("Alto riesgo (D+/R-)\n(n = %d)", n_alto),
  "Riesgo intermedio (R+)"  = sprintf("Riesgo intermedio (R+)\n(n = %d)", n_inter)
)

# ── 4. Líneas horizontales con duración de profilaxis ─────────────────────
df_prof <- tibble(
  riesgo = factor(c("Alto riesgo (D+/R-)", "Riesgo intermedio (R+)"),
                  levels = levels(df_fig$riesgo)),
  prof_meses = c(12, 6)
)

# ── 5. Tests: solo significativos ──────────────────────────────────────────
DODGE <- 0.7

stat_test <- df_long %>%
  group_by(riesgo, glucoproteina) %>%
  rstatix::wilcox_test(meses_tx_inicio ~ genotipo) %>%
  rstatix::add_significance("p") %>%
  ungroup() %>%
  filter(p.signif != "ns") %>%   # solo significativos
  mutate(
    x_pos      = as.numeric(glucoproteina),
    xmin       = x_pos - DODGE / 4,
    xmax       = x_pos + DODGE / 4,
    y.position = 22
  )

# ── 6. Paletas por glucoproteína ───────────────────────────────────────────
# gB: azules     | gH: rosas    | gB+gH: verde/rojo
cols_fill <- c(
  "gB_Único"     = "#AED6F1",  "gB_Mixto"     = "#4682B4",
  "gH_Único"     = "pink",     "gH_Mixto"     = "hotpink3",
  "gB+gH_Único"  = "#4CAF7D",  "gB+gH_Mixto"  = "#E07070"
)

cols_jit <- c(
  "gB_Único"     = "#2980B9",  "gB_Mixto"     = "#0D2B45",
  "gH_Único"     = "pink3",    "gH_Mixto"     = "palevioletred4",
  "gB+gH_Único"  = "#2E7D52",  "gB+gH_Mixto"  = "#B04040"
)

# ── 7. Figura ──────────────────────────────────────────────────────────────
fig_tiempo <- ggplot(df_long, aes(x = glucoproteina, y = meses_tx_inicio,
                                  fill = cat_color,
                                  group = interaction(glucoproteina, genotipo))) +
  geom_hline(data = df_prof,
             aes(yintercept = prof_meses),
             color = "grey55", linetype = "dashed", linewidth = 0.5) +
  geom_boxplot(
    width          = 0.55,
    position       = position_dodge(width = DODGE),
    color          = "grey25",
    linewidth      = 0.4,
    outlier.shape  = 21,
    outlier.size   = 1.4,
    outlier.fill   = "white",
    outlier.color  = "grey50",
    outlier.stroke = 0.3,
    outlier.alpha  = 0.7
  ) +
  geom_jitter(
    aes(color = cat_color),
    position    = position_jitterdodge(jitter.width = 0.10,
                                       dodge.width  = DODGE),
    size        = 0.9,
    alpha       = 0.35,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = cols_fill, guide = "none") +
  scale_color_manual(values = cols_jit, guide = "none") +
  facet_wrap(~ riesgo, labeller = labeller(riesgo = etiquetas_riesgo)) +
  {
    if (nrow(stat_test) > 0)
      stat_pvalue_manual(
        stat_test,
        label        = "p.signif",
        xmin         = "xmin",
        xmax         = "xmax",
        y.position   = "y.position",
        tip.length   = 0.01,
        bracket.size = 0.4,
        size         = 5,
        color        = "grey20",
        inherit.aes  = FALSE
      )
  } +
  labs(
    x = NULL,
    y = "Meses desde el trasplante al primer episodio"
  ) +
  scale_y_continuous(
    breaks = seq(0, 24, 6),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  coord_cartesian(ylim = c(0, 24)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "none",
    panel.border       = element_blank(),
    panel.spacing.x    = unit(15, "pt"),
    axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
    axis.line.y        = element_line(color = "grey60", linewidth = 0.4),
    axis.ticks         = element_line(color = "grey60", linewidth = 0.3),
    axis.text.x        = element_text(size = 11, color = "grey25"),
    axis.text.y        = element_text(size = 10, color = "grey30"),
    axis.title.y       = element_text(size = 11, color = "grey20"),
    strip.background   = element_rect(fill = "grey94",
                                      color = "grey70", linewidth = 0.4),
    strip.text         = element_text(size = 11, color = "grey15",
                                      face = "bold",
                                      margin = margin(t = 4, b = 4)),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.margin        = margin(10, 12, 8, 10)
  )

fig_tiempo

ggsave("~/Desktop/Figuras_tesis/figura_R20_tiempo_primer_episodio.png",
       fig_tiempo,
       width = 7, height = 4, dpi = 320, bg = "white")

