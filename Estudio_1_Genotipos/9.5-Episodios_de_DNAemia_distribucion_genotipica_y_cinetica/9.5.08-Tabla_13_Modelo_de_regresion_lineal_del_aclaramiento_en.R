# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.5.08. Tabla 13. Modelo de regresión lineal del aclaramiento en episodios con genotipo único de gB
# Extracto del script completo (ver README.md del estudio)

## ---- Tabla 13. Modelo de regresión lineal del aclaramiento en episodios con genotipo único de gB ----
# Preparación, análisis univariante (Figura R19, no incluida) y selección backward;
# la Tabla 13 corresponde al modelo clínico con variables elegidas a priori (modelo_clinico_gb3).
#  APARTADO 4.5.2 — IMPACTO DE LOS GENOTIPOS gB Y gH EN LA CINÉTICA VIRAL
#
#  Tesis CMV en RTP — Carmen, HU12O
#
#  Contenido del script:
#    1. Preparación de datos (df_gb3_uni)
#    2. Análisis descriptivo por genotipo  → Tabla R6, Figura R18
#    3. Análisis univariante (confusores)  → Figura R19
#    4. Análisis multivariante (backward)  → Tabla R7
#    5. Análisis de sensibilidad (modelo mixto)
#    6. Análisis exploratorio gB-gH (retaguardia)
#    7. Análisis de sensibilidad sin gB4 (retaguardia)
#
#  Outputs:
#    Tabla R6   → tabla_R6_cinetica_genotipos.docx (o pegada en Word)
#    Figura R18 → figura_R18_boxplots_cinetica_gB.png
#    Figura R19 → figura_R19_univariante_gB3.png
#    Tabla R7   → tabla_R7_multivariante_gB3.docx (o pegada en Word)
#
# ════════════════════════════════════════════════════════════════════════════


# ════════════════════════════════════════════════════════════════════════════
#  1. PREPARACIÓN DE DATOS
# ════════════════════════════════════════════════════════════════════════════

df_gb3_uni <- dnaemias %>%
  left_join(pacientes %>%
              select(NHC, edad_tx, tipotx, ebase,
                     SexF, edad_don),
            by = "NHC") %>%
  filter(gbsum == 1) %>%
  mutate(
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini),
    Tdesdetx_meses    = as.numeric(Fini - as.Date(Ftx)) / 30.44,
    sintomatico       = case_when(
      clinica == "VA"  ~ 0L,
      !is.na(clinica)  ~ 1L,
      TRUE             ~ NA_integer_
    ),
    risk_f       = factor(
      ifelse(serDR == "D+/R-", "Alto riesgo", "Riesgo intermedio"),
      levels = c("Riesgo intermedio", "Alto riesgo")
    ),
    gb3_f        = factor(gb3,         levels = c(0,1), labels = c("No gB3","gB3")),
    ingr_f       = factor(ingr,        levels = c(0,1), labels = c("No","Sí")),
    ing_cmv_f    = factor(ing_cmv,     levels = c(0,1), labels = c("No","Sí")),
    sintom_f     = factor(sintomatico, levels = c(0,1), labels = c("No","Sí")),
    profilaxis_f = factor(antiviral,   levels = c(0,1), labels = c("No","Sí")),
    tto_f        = factor(tto,         levels = c(0,1), labels = c("No","Sí")),
    tipotx_f     = factor(tipotx,      levels = c("UNI","BI")),
    ebase_f      = factor(ebase),
    sexo_f       = factor(SexF, levels = c(0,1), labels = c("Hombre","Mujer"))
  ) %>%
  filter(!is.na(aclaramiento_dias), !is.na(gb3_f))

cat(sprintf("\n── N episodios mono-gB: %d │ N pacientes: %d ──\n",
            nrow(df_gb3_uni), n_distinct(df_gb3_uni$NHC)))


# ════════════════════════════════════════════════════════════════════════════
#  2. ANÁLISIS DESCRIPTIVO POR GENOTIPO
#     → Tabla R6 (medianas e IQR por genotipo)
#     → Figura R18 (boxplots, código no incluido aquí)
# ════════════════════════════════════════════════════════════════════════════

# 2.1 Construir df ampliado con genotipo gB de 4 niveles ────────────────────
df_gb_4niv <- dnaemias %>%
  filter(gbsum == 1) %>%
  mutate(
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini),
    gb_tipo = case_when(
      gb1 == 1 ~ "gB1",
      gb2 == 1 ~ "gB2",
      gb3 == 1 ~ "gB3",
      gb4 == 1 ~ "gB4"
    )
  ) %>%
  filter(!is.na(gb_tipo))

# 2.2 Resumen Tabla R6: cinética por genotipo gB ────────────────────────────
cat("\n── Tabla R6: Parámetros cinéticos por genotipo gB ──\n")
df_gb_4niv %>%
  group_by(gb_tipo) %>%
  summarise(
    n_cv         = sum(!is.na(cvini)),
    cv_inicial   = sprintf("%.0f (%.0f–%.0f)",
                           median(cvini, na.rm=TRUE),
                           quantile(cvini, 0.25, na.rm=TRUE),
                           quantile(cvini, 0.75, na.rm=TRUE)),
    cv_pico      = sprintf("%.0f (%.0f–%.0f)",
                           median(cvpico, na.rm=TRUE),
                           quantile(cvpico, 0.25, na.rm=TRUE),
                           quantile(cvpico, 0.75, na.rm=TRUE)),
    n_ac         = sum(!is.na(aclaramiento_dias)),
    aclaramiento = sprintf("%.0f (%.0f–%.0f)",
                           median(aclaramiento_dias, na.rm=TRUE),
                           quantile(aclaramiento_dias, 0.25, na.rm=TRUE),
                           quantile(aclaramiento_dias, 0.75, na.rm=TRUE)),
    n_auc        = sum(!is.na(AUCCV)),
    auc_viral    = sprintf("%.0f (%.0f–%.0f)",
                           median(AUCCV, na.rm=TRUE),
                           quantile(AUCCV, 0.25, na.rm=TRUE),
                           quantile(AUCCV, 0.75, na.rm=TRUE)),
    .groups      = "drop"
  ) %>%
  print()

# 2.3 Tests globales Kruskal-Wallis (4 genotipos gB) ────────────────────────
cat("\n── Tests Kruskal-Wallis (4 genotipos gB) ──\n")
cat(sprintf("CV inicial:    p = %.3f\n",
            kruskal.test(cvini ~ gb_tipo, data = df_gb_4niv)$p.value))
cat(sprintf("CV pico:       p = %.3f\n",
            kruskal.test(cvpico ~ gb_tipo, data = df_gb_4niv)$p.value))
cat(sprintf("Aclaramiento:  p = %.3f\n",
            kruskal.test(aclaramiento_dias ~ gb_tipo, data = df_gb_4niv)$p.value))
cat(sprintf("AUC viral:     p = %.3f\n",
            kruskal.test(AUCCV ~ gb_tipo, data = df_gb_4niv)$p.value))

# 2.4 Post-hoc Dunn con corrección Bonferroni ───────────────────────────────
if (requireNamespace("FSA", quietly = TRUE)) {
  cat("\n── Post-hoc Dunn (Bonferroni) - aclaramiento ──\n")
  print(FSA::dunnTest(aclaramiento_dias ~ gb_tipo,
                      data = df_gb_4niv, method = "bonferroni"))
  
  cat("\n── Post-hoc Dunn (Bonferroni) - AUC ──\n")
  print(FSA::dunnTest(AUCCV ~ gb_tipo,
                      data = df_gb_4niv, method = "bonferroni"))
} else {
  message("Instalar FSA para post-hoc Dunn: install.packages('FSA')")
}


# ════════════════════════════════════════════════════════════════════════════
#  3. ANÁLISIS UNIVARIANTE DE CONFUSORES
#     → Figura R19 (forest plot tipo tabla con cajas por dominio)
# ════════════════════════════════════════════════════════════════════════════

# 3.1 Correlaciones de Spearman (variables continuas) ───────────────────────
cor_vars <- list(
  "Edad receptor (años)"     = "edad_tx",
  "Edad donante (años)"      = "edad_don",
  "CV inicial (UI/mL)"       = "cvini",
  "CV pico (UI/mL)"          = "cvpico",
  "Tiempo desde TXP (meses)" = "Tdesdetx_meses"
)

res_cor <- map_dfr(names(cor_vars), function(lab) {
  v <- cor_vars[[lab]]
  x <- as.numeric(df_gb3_uni[[v]])
  if (v %in% c("cvini","cvpico")) x <- log10(x + 1)
  idx <- !is.na(x) & !is.na(df_gb3_uni$aclaramiento_dias)
  test <- cor.test(x[idx], df_gb3_uni$aclaramiento_dias[idx],
                   method = "spearman", exact = FALSE)
  tibble(Variable = lab, p_val = test$p.value,
         label_est = sprintf("\u03C1 = %+.3f", test$estimate))
})

# 3.2 Wilcoxon (categóricas 2 niveles) ──────────────────────────────────────
cat_vars <- list(
  "Genotipo gB3"                      = "gb3_f",
  "Sexo (Mujer vs Hombre)"            = "sexo_f",
  "D+/R- (vs R+)"                     = "risk_f",
  "Hospitalización cualquier causa"   = "ingr_f",
  "Hospitalización por CMV"           = "ing_cmv_f",
  "Enfermedad sintomática"            = "sintom_f",
  "Profilaxis al inicio del episodio" = "profilaxis_f",
  "Tratamiento antiviral"             = "tto_f",
  "Tipo de trasplante (BI vs UNI)"    = "tipotx_f"
)

res_cat <- map_dfr(names(cat_vars), function(lab) {
  v  <- cat_vars[[lab]]
  g  <- df_gb3_uni[[v]]
  ac <- df_gb3_uni$aclaramiento_dias
  d0 <- ac[g == levels(g)[1] & !is.na(g)]
  d1 <- ac[g == levels(g)[2] & !is.na(g)]
  test <- wilcox.test(d1, d0, exact = FALSE)
  dif  <- median(d1, na.rm = TRUE) - median(d0, na.rm = TRUE)
  tibble(Variable = lab, p_val = test$p.value,
         label_est = sprintf("\u0394 = %+.0f d\u00edas", dif))
})

# 3.3 Kruskal-Wallis (categórica multinivel: enfermedad de base) ────────────
kw_eb <- kruskal.test(aclaramiento_dias ~ ebase_f, data = df_gb3_uni)
res_eb <- tibble(Variable = "Enfermedad de base",
                 p_val = kw_eb$p.value,
                 label_est = "(no calculable)")

# 3.4 Combinar y categorizar ────────────────────────────────────────────────
orden_vars <- c(
  # Demográficas
  "Edad receptor (años)",
  "Edad donante (años)",
  "Sexo (Mujer vs Hombre)",
  "Tipo de trasplante (BI vs UNI)",
  "Enfermedad de base",
  "D+/R- (vs R+)",
  # Cinéticas
  "Tiempo desde TXP (meses)",
  "CV inicial (UI/mL)",
  "CV pico (UI/mL)",
  "Genotipo gB3",
  # Clínicas
  "Enfermedad sintomática",
  "Profilaxis al inicio del episodio",
  "Tratamiento antiviral",
  "Hospitalización cualquier causa",
  "Hospitalización por CMV"
)

grupos_map <- c(
  "Edad receptor (años)"              = "Demográficas",
  "Edad donante (años)"               = "Demográficas",
  "Sexo (Mujer vs Hombre)"            = "Demográficas",
  "Tipo de trasplante (BI vs UNI)"    = "Demográficas",
  "Enfermedad de base"                = "Demográficas",
  "D+/R- (vs R+)"                     = "Demográficas",
  "Tiempo desde TXP (meses)"          = "Cinéticas",
  "CV inicial (UI/mL)"                = "Cinéticas",
  "CV pico (UI/mL)"                   = "Cinéticas",
  "Genotipo gB3"                      = "Cinéticas",
  "Enfermedad sintomática"            = "Clínicas",
  "Profilaxis al inicio del episodio" = "Clínicas",
  "Tratamiento antiviral"             = "Clínicas",
  "Hospitalización cualquier causa"   = "Clínicas",
  "Hospitalización por CMV"           = "Clínicas"
)

res_all <- bind_rows(res_cor, res_cat, res_eb) %>%
  mutate(
    sig   = p_val < 0.05,
    cand  = p_val < 0.20 & !sig,
    p_fmt = case_when(
      p_val < 0.001 ~ "<0,001",
      TRUE          ~ sub("\\.", ",", sprintf("%.3f", p_val))
    ),
    grupo    = factor(grupos_map[as.character(Variable)],
                      levels = c("Demográficas", "Cinéticas", "Clínicas")),
    Variable = factor(Variable, levels = rev(orden_vars)),
    color_cat = factor(case_when(
      sig  ~ "Significativa (p < 0,05)",
      cand ~ "Tendencia (p < 0,20)",
      TRUE ~ "No asociada"
    ), levels = c("Significativa (p < 0,05)",
                  "Tendencia (p < 0,20)",
                  "No asociada"))
  )

cat("\n── Univariante gB3: resultados ordenados por p ──\n")
print(res_all %>% arrange(p_val) %>%
        select(Variable, label_est, p_fmt), n = Inf)


# 3.5 Figura R19 ─────────────────────────────────────────────────────────────
X_EST_COL  <- -0.05
X_PLOT_MIN <-  0.00
X_PLOT_MAX <-  1.10
X_P_COL    <-  1.12
X_LEFT     <- -0.40
X_RIGHT    <-  1.32

cabeceras <- tibble(
  grupo = factor("Demográficas",
                 levels = c("Demográficas", "Cinéticas", "Clínicas")),
  Variable = factor("Edad receptor (años)",
                    levels = levels(res_all$Variable))
)

figura_R19_univariante_gb3 <- ggplot(res_all, aes(x = p_val, y = Variable)) +
  annotate("rect",
           xmin = X_PLOT_MIN, xmax = X_PLOT_MAX,
           ymin = -Inf, ymax = Inf,
           fill = "grey97", color = NA) +
  geom_vline(xintercept = 0.05, color = "firebrick3",
             linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  geom_vline(xintercept = 0.20, color = "grey55",
             linetype = "dashed", linewidth = 0.4, alpha = 0.6) +
  geom_text(aes(x = X_EST_COL, label = label_est),
            hjust = 1, size = 3.6, color = "grey25",
            show.legend = FALSE) +
  geom_point(aes(color = color_cat), size = 3.0) +
  geom_text(aes(x = X_P_COL, label = p_fmt, color = color_cat),
            hjust = 0, size = 3.6, fontface = "bold",
            show.legend = FALSE) +
  geom_text(data = cabeceras,
            aes(x = X_EST_COL, y = Inf, label = "Estimador"),
            hjust = 1, vjust = -0.5, size = 3.8,
            fontface = "bold", color = "grey15",
            inherit.aes = FALSE) +
  geom_text(data = cabeceras,
            aes(x = (X_PLOT_MIN + X_PLOT_MAX) / 2, y = Inf,
                label = "valor p"),
            hjust = 0.5, vjust = -0.5, size = 3.8,
            fontface = "bold", color = "grey15",
            inherit.aes = FALSE) +
  facet_grid(grupo ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  scale_x_continuous(
    name   = NULL,
    limits = c(X_LEFT, X_RIGHT),
    breaks = c(0.05, 0.20, 0.50, 1.00),
    labels = c("0,05", "0,20", "0,50", "1,00"),
    expand = c(0, 0)
  ) +
  scale_color_manual(
    name   = NULL,
    values = c(
      "Significativa (p < 0,05)" = "firebrick3",
      "Tendencia (p < 0,20)"     = "#E08E45",
      "No asociada"              = "grey55"
    )
  ) +
  coord_cartesian(clip = "off") +
  labs(y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position      = "bottom",
    legend.text          = element_text(size = 10),
    legend.key.size      = unit(10, "pt"),
    legend.margin        = margin(t = 4, b = 0),
    legend.box.spacing   = unit(2, "pt"),
    panel.grid.major.y   = element_line(color = "grey95", linewidth = 0.3),
    panel.grid.major.x   = element_blank(),
    panel.grid.minor     = element_blank(),
    axis.line.x          = element_line(color = "grey60", linewidth = 0.4),
    axis.ticks.x         = element_line(color = "grey60", linewidth = 0.3),
    axis.ticks.length.x  = unit(3, "pt"),
    axis.text.y          = element_text(size = 11, color = "grey15"),
    axis.text.x          = element_text(size = 10, color = "grey40"),
    strip.placement      = "outside",
    strip.text.y.left    = element_text(angle = 90, size = 11.5,
                                        face = "bold", color = "white",
                                        margin = margin(r = 5, l = 5)),
    strip.background.y   = element_rect(fill = "grey35", color = NA),
    panel.spacing.y      = unit(8, "pt"),
    panel.border         = element_rect(color = "grey80",
                                        fill = NA, linewidth = 0.4),
    plot.margin          = margin(28, 14, 8, 8)
  )

figura_R19_univariante_gb3

# [guardado antiguo] ggsave("~/Desktop/Figuras_tesis/figura_R19_univariante_gB3.png",
# [guardado antiguo]        figura_R19_univariante_gb3,
# [guardado antiguo]        width = 9, height = 6, dpi = 320, bg = "white")


# ════════════════════════════════════════════════════════════════════════════
#  4. ANÁLISIS MULTIVARIANTE (regresión lineal con backward)
#     → Tabla R7
# ════════════════════════════════════════════════════════════════════════════

# 4.1 Construir df con variables candidatas (p < 0,20) ──────────────────────
df_multi_gb3 <- df_gb3_uni %>%
  select(NHC, aclaramiento_dias,
         gb3_f, risk_f, cvini, cvpico, Tdesdetx_meses,
         sintom_f, profilaxis_f, ing_cmv_f, ingr_f, tto_f) %>%
  filter(complete.cases(.))

cat(sprintf("\n── N para multivariante (datos completos): %d episodios / %d pacientes ──\n",
            nrow(df_multi_gb3), n_distinct(df_multi_gb3$NHC)))

# 4.2 Modelo completo ────────────────────────────────────────────────────────
modelo_completo_gb3 <- lm(
  aclaramiento_dias ~ gb3_f + risk_f + log10(cvini + 1) + log10(cvpico + 1) +
    Tdesdetx_meses + sintom_f + profilaxis_f + ing_cmv_f + ingr_f,
  data = df_multi_gb3
)

cat("\n── Modelo completo (todas las candidatas) ──\n")
print(summary(modelo_completo_gb3))

# 4.3 Selección backward (Akaike) ───────────────────────────────────────────
modelo_backward_gb3 <- step(modelo_completo_gb3, direction = "backward",
                            trace = FALSE)

cat("\n══════════════════════════════════════════════\n")
cat("  MODELO FINAL (post-backward)\n")
cat("══════════════════════════════════════════════\n")
print(summary(modelo_backward_gb3))

cat("\n── IC 95% del modelo final ──\n")
print(confint(modelo_backward_gb3))

cat("\n── VIF (multicolinealidad) ──\n")
print(car::vif(modelo_backward_gb3))


# ════════════════════════════════════════════════════════════════════════════
#  5. ANÁLISIS DE SENSIBILIDAD: modelo lineal mixto (paciente aleatorio)
# ════════════════════════════════════════════════════════════════════════════

modelo_mixto_gb3 <- lmer(
  aclaramiento_dias ~ gb3_f + log10(cvini + 1) + log10(cvpico + 1) +
    profilaxis_f + ing_cmv_f + (1 | NHC),
  data = df_multi_gb3
)

cat("\n══════════════════════════════════════════════\n")
cat("  MODELO MIXTO: paciente como efecto aleatorio\n")
cat("══════════════════════════════════════════════\n")
print(summary(modelo_mixto_gb3))

cat("\n── IC 95% (Wald) ──\n")
print(confint(modelo_mixto_gb3, method = "Wald"))

vc  <- as.data.frame(VarCorr(modelo_mixto_gb3))
icc <- vc$vcov[1] / sum(vc$vcov)
cat(sprintf("\n── ICC (varianza atribuible al paciente): %.3f ──\n", icc))

# 5.1 Comparación lineal vs mixto ───────────────────────────────────────────
coef_lin <- coef(modelo_backward_gb3)["gb3_fgB3"]
ci_lin   <- confint(modelo_backward_gb3)["gb3_fgB3", ]
p_lin    <- summary(modelo_backward_gb3)$coefficients["gb3_fgB3", "Pr(>|t|)"]

coef_mix <- fixef(modelo_mixto_gb3)["gb3_fgB3"]
ci_mix   <- confint(modelo_mixto_gb3, method = "Wald")["gb3_fgB3", ]
p_mix    <- summary(modelo_mixto_gb3)$coefficients["gb3_fgB3", "Pr(>|t|)"]

cat("\n── Comparación LINEAL vs MIXTO (efecto gB3) ──\n")
cat(sprintf("  Lineal:  %+.1f días (IC 95%%: %+.1f a %+.1f), p = %.4f\n",
            coef_lin, ci_lin[1], ci_lin[2], p_lin))
cat(sprintf("  Mixto:   %+.1f días (IC 95%%: %+.1f a %+.1f), p = %.4f\n",
            coef_mix, ci_mix[1], ci_mix[2], p_mix))


# ════════════════════════════════════════════════════════════════════════════
#  6. ANÁLISIS EXPLORATORIO: combinaciones gB-gH (RETAGUARDIA)
# ════════════════════════════════════════════════════════════════════════════

df_gbgh_combo <- dnaemias %>%
  filter(gbsum == 1, ghsum == 1) %>%
  mutate(
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini),
    gb_tipo = case_when(
      gb1 == 1 ~ "gB1", gb2 == 1 ~ "gB2",
      gb3 == 1 ~ "gB3", gb4 == 1 ~ "gB4"
    ),
    gh_tipo = case_when(
      gh1 == 1 ~ "gH1", gh2 == 1 ~ "gH2"
    ),
    combo = paste0(gb_tipo, "-", gh_tipo)
  ) %>%
  filter(!is.na(aclaramiento_dias), !is.na(combo))

cat("\n── Combinaciones gB-gH (8 grupos) ──\n")
print(
  df_gbgh_combo %>%
    group_by(combo) %>%
    summarise(
      n       = n(),
      mediana = median(aclaramiento_dias),
      q1      = quantile(aclaramiento_dias, 0.25),
      q3      = quantile(aclaramiento_dias, 0.75),
      .groups = "drop"
    ) %>% arrange(desc(mediana))
)

cat("\n── Kruskal-Wallis global (8 combinaciones) ──\n")
print(kruskal.test(aclaramiento_dias ~ combo, data = df_gbgh_combo))

# 6.1 Foco en gB3: ¿efecto del gH? ──────────────────────────────────────────
df_gb3_combo <- df_gbgh_combo %>% filter(gb_tipo == "gB3")

cat("\n── Dentro de gB3: gH1 vs gH2 ──\n")
print(df_gb3_combo %>%
        group_by(combo) %>%
        summarise(n = n(),
                  mediana = median(aclaramiento_dias),
                  .groups = "drop"))

cat("\n── Wilcoxon gB3-gH1 vs gB3-gH2 ──\n")
print(wilcox.test(aclaramiento_dias ~ combo,
                  data = df_gb3_combo, exact = FALSE))


# ════════════════════════════════════════════════════════════════════════════
#  7. ANÁLISIS DE SENSIBILIDAD: excluyendo gB4 (RETAGUARDIA)
# ════════════════════════════════════════════════════════════════════════════

df_gb_no_gb4 <- df_gbgh_combo %>% filter(gb_tipo != "gB4")

cat("\n── Kruskal-Wallis sin gB4 (3 grupos: gB1, gB2, gB3) ──\n")
print(kruskal.test(aclaramiento_dias ~ gb_tipo, data = df_gb_no_gb4))


# ════════════════════════════════════════════════════════════════════════════
#  FIN DEL SCRIPT
# ════════════════════════════════════════════════════════════════════════════

# ── Modelo clínico gB3 (variables elegidas a priori) ──────────────────────────
df_clin_gb3 <- df_gb3_uni %>%
  select(NHC, aclaramiento_dias,
         gb3_f, risk_f, cvini, cvpico, profilaxis_f) %>%
  filter(complete.cases(.))

cat(sprintf("\n── N modelo clínico: %d episodios / %d pacientes ──\n",
            nrow(df_clin_gb3), n_distinct(df_clin_gb3$NHC)))

modelo_clinico_gb3 <- lm(
  aclaramiento_dias ~ gb3_f + risk_f + log10(cvini + 1) + log10(cvpico + 1) +
    profilaxis_f,
  data = df_clin_gb3
)

cat("\n══════════════════════════════════════════════\n")
cat("  MODELO CLÍNICO gB3 (variables a priori)\n")
cat("══════════════════════════════════════════════\n")
print(summary(modelo_clinico_gb3))

cat("\n── IC 95% ──\n")
print(confint(modelo_clinico_gb3))

cat("\n── VIF (multicolinealidad) ──\n")
print(car::vif(modelo_clinico_gb3))

# Efecto ajustado de gB3 en una línea
coef_gb3 <- coef(modelo_clinico_gb3)["gb3_fgB3"]
ci_gb3   <- confint(modelo_clinico_gb3)["gb3_fgB3", ]
p_gb3    <- summary(modelo_clinico_gb3)$coefficients["gb3_fgB3", "Pr(>|t|)"]

cat(sprintf("\n── Efecto ajustado gB3: %+.1f días (IC 95%%: %+.1f a %+.1f), p = %.4f ──\n",
            coef_gb3, ci_gb3[1], ci_gb3[2], p_gb3))


modelo_mixto_clinico <- lmer(
  aclaramiento_dias ~ gb3_f + risk_f + log10(cvini + 1) + log10(cvpico + 1) +
    profilaxis_f + (1 | NHC),
  data = df_clin_gb3
)

print(summary(modelo_mixto_clinico))
print(confint(modelo_mixto_clinico, method = "Wald")["gb3_fgB3", ])

vc  <- as.data.frame(VarCorr(modelo_mixto_clinico))
icc <- vc$vcov[1] / sum(vc$vcov)
cat(sprintf("\nICC: %.3f\n", icc))

guardar_modelo(modelo_clinico_gb3, "Tabla_13_modelo")

