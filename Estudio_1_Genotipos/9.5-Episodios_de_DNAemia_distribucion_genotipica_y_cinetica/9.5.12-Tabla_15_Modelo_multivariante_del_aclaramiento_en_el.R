# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.5.12. Tabla 15. Modelo multivariante del aclaramiento en el análisis de genotipos mixtos gB/gH
# Extracto del script completo (ver README.md del estudio)

## ---- Tabla 15. Modelo multivariante del aclaramiento en el análisis de genotipos mixtos gB/gH ----
# Análisis univariante de confusores y modelos previos; la Tabla 15 corresponde a modelo_clinico_mix.
# ── 1. Preparar datos ─────────────────────────────────────────────────────────
df_uni <- dnaemias %>%
  left_join(pacientes %>% select(NHC, edad_tx, tipotx, ebase), by = "NHC") %>%
  mutate(
    risk_f = factor(
      ifelse(serDR == "D+/R-", "Alto riesgo", "Riesgo intermedio"),
      levels = c("Riesgo intermedio", "Alto riesgo")
    ),
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini),
    Tdesdetx_meses    = as.numeric(Fini - as.Date(Ftx)) / 30.44,
    sintomatico       = case_when(
      clinica == "VA"  ~ 0L,
      !is.na(clinica)  ~ 1L,
      TRUE             ~ NA_integer_
    ),
    gbghmix_f    = factor(gbghmix,  levels = c(0,1), labels = c("Único","Mixto")),
    ingr_f       = factor(ingr,     levels = c(0,1), labels = c("No","Sí")),
    ing_cmv_f    = factor(ing_cmv,  levels = c(0,1), labels = c("No","Sí")),
    sintom_f     = factor(sintomatico, levels = c(0,1), labels = c("No","Sí")),
    profilaxis_f = factor(antiviral, levels = c(0,1), labels = c("No","Sí")),
    tto_f        = factor(tto,      levels = c(0,1), labels = c("No","Sí")),
    cmvr_f       = factor(CMV_R,    levels = c(0,1), labels = c("No","Sí")),
    tipotx_f     = factor(tipotx,   levels = c("UNI","BI")),
    ebase_f      = factor(ebase)
  ) %>%
  filter(!is.na(aclaramiento_dias))

# ── 2. Correlaciones de Spearman (variables continuas) ────────────────────────
cor_vars <- list(
  "Edad receptor (años)"        = "edad_tx",
  "CV inicial (UI/mL)"          = "cvini",
  "CV pico (UI/mL)"             = "cvpico",
  "Tiempo desde TXP (meses)"     = "Tdesdetx_meses"
)

res_cor <- map_dfr(names(cor_vars), function(lab) {
  v <- cor_vars[[lab]]
  x <- as.numeric(df_uni[[v]])
  if (v %in% c("cvini","cvpico")) x <- log10(x + 1)
  idx  <- !is.na(x) & !is.na(df_uni$aclaramiento_dias)
  x_ok <- x[idx]
  y_ok <- df_uni$aclaramiento_dias[idx]
  if (length(x_ok) < 5) return(tibble(Variable = lab, tipo = "continua",
                                      estimador = NA_real_, p_val = NA_real_,
                                      label_est = "n insuf."))
  test <- tryCatch(
    cor.test(x_ok, y_ok, method = "spearman", exact = FALSE),
    error = function(e) NULL
  )
  if (is.null(test)) return(tibble(Variable = lab, tipo = "continua",
                                   estimador = NA_real_, p_val = NA_real_,
                                   label_est = "error"))
  tibble(
    Variable  = lab,
    tipo      = "continua",
    estimador = round(test$estimate, 3),
    p_val     = test$p.value,
    label_est = sprintf("rho = %+.3f", test$estimate)
  )
})

# ── 3. Comparaciones de grupos (variables categóricas 2 grupos) ───────────────
cat_vars <- list(
  "Genotipos mixtos gB/gH"            = "gbghmix_f",
  "Perfil de riesgo serológico D/R"   = "risk_f",      
  "Hospitalización cualquier causa"   = "ingr_f",
  "Hospitalización por CMV"           = "ing_cmv_f",
  "Enfermedad sintomática"            = "sintom_f",
  "Profilaxis al inicio del episodio" = "profilaxis_f",
  "Tratamiento antiviral"             = "tto_f",
  "Tipo de trasplante (BI vs UNI)"    = "tipotx_f"
)

res_cat <- map_dfr(names(cat_vars), function(lab) {
  v  <- cat_vars[[lab]]
  g  <- df_uni[[v]]
  ac <- df_uni$aclaramiento_dias
  d0 <- ac[g == levels(g)[1] & !is.na(g)]
  d1 <- ac[g == levels(g)[2] & !is.na(g)]
  test <- wilcox.test(d1, d0, exact = FALSE)
  dif  <- median(d1, na.rm = TRUE) - median(d0, na.rm = TRUE)
  tibble(
    Variable  = lab,
    tipo      = "categorica",
    estimador = round(dif, 1),
    p_val     = test$p.value,
    label_est = sprintf("\u0394 = %+.1f d\u00edas", dif)
  )
})

# ── 4. Kruskal-Wallis (variables categóricas >2 grupos) ───────────────────────
kw_eb  <- kruskal.test(aclaramiento_dias ~ ebase_f, data = df_uni)

res_eb <- tibble(
  Variable  = "Enfermedad de base",
  tipo      = "categorica",
  estimador = NA_real_,
  p_val     = kw_eb$p.value,
  label_est = "Kruskal-Wallis"
)

res_eb <- tibble(
  Variable  = "Enfermedad de base",
  tipo      = "categorica",
  estimador = NA_real_,
  p_val     = kw_eb$p.value,
  label_est = "Kruskal-Wallis"
)

# ── 5. Combinar y preparar para figura ────────────────────────────────────────
orden_vars <- c(
  "Edad receptor (años)",
  "Tipo de trasplante (BI vs UNI)",
  "Enfermedad de base",
  "Perfil de riesgo serológico D/R",    # <-- corregido
  "Tiempo desde TXP (meses)",
  "CV inicial (UI/mL)",
  "CV pico (UI/mL)",
  "Genotipos mixtos gB/gH",
  "Enfermedad sintomática",
  "Profilaxis al inicio del episodio",
  "Tratamiento antiviral",
  "Hospitalización cualquier causa",
  "Hospitalización por CMV"
)

grupos_map <- c(
  "Edad receptor (años)"              = "Demográficas",
  "Tipo de trasplante (BI vs UNI)"    = "Demográficas",
  "Enfermedad de base"                = "Demográficas",
  "Perfil de riesgo serológico D/R"   = "Demográficas",   
  "Tiempo desde TXP (meses)"          = "Cinéticas",
  "CV inicial (UI/mL)"                = "Cinéticas",
  "CV pico (UI/mL)"                   = "Cinéticas",
  "Genotipos mixtos gB/gH"            = "Cinéticas",
  "Enfermedad sintomática"            = "Clínicas",
  "Profilaxis al inicio del episodio" = "Clínicas",
  "Tratamiento antiviral"             = "Clínicas",
  "Hospitalización cualquier causa"   = "Clínicas",
  "Hospitalización por CMV"           = "Clínicas"
)

res_all <- bind_rows(res_cor, res_cat, res_eb) %>%
  mutate(
    sig     = p_val < 0.05,
    p_fmt   = case_when(
      p_val < 0.001 ~ "<0,001",
      TRUE          ~ sub("\\.", ",", sprintf("%.3f", p_val))
    ),
    p_label = paste0("p = ", p_fmt),
    grupo   = factor(grupos_map[as.character(Variable)],
                     levels = c("Demográficas", "Cinéticas", "Clínicas")),
    Variable = factor(Variable, levels = rev(orden_vars))
  )

# ── 6. Figura ─────────────────────────────────────────────────────────────────
figura_univariante <- ggplot(res_all, aes(x = p_val, y = Variable)) +
  geom_vline(
    xintercept = 0.05,
    color = "firebrick3",
    linetype = "dashed",
    linewidth = 0.5,
    alpha = 0.5
  ) +
  geom_vline(xintercept = 0.20, color = "grey65",
             linetype = "dashed", linewidth = 0.4) +
  geom_point(aes(color = sig, size = sig)) +
  geom_text(aes(label = p_label, color = sig),
            hjust = -0.15, size = 4) +
  geom_text(aes(x = 1.1, label = label_est),
            hjust = 1, size = 4, color = "grey40") +
  facet_grid(grupo ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(
    name   = "valor p",
    limits = c(0, 1.15),
    breaks = c(0.05, 0.20, 0.50, 1.00),
    labels = c("0,05", "0,20", "0,50", "1,00")
  ) +
  scale_color_manual(values = c("TRUE" = "firebrick3", "FALSE" = "grey50")) +
  scale_size_manual(values  = c("TRUE" = 3.5,          "FALSE" = 2.5)) +
  labs(y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "none",
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.3),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
    axis.text.y        = element_text(size = 12, color = "grey20"),
    axis.text.x        = element_text(size = 12, color = "grey40"),
    axis.title.x       = element_text(size = 11, color = "grey20"),
    strip.text.y       = element_text(angle = -90, size = 12,
                                      face = "bold", color = "grey20",
                                      margin = margin(l = 6, r = 6)),
    strip.background   = element_rect(fill = "grey93", color = NA),
    panel.spacing.y    = unit(0.6, "lines"),
    plot.margin        = margin(8, 60, 8, 8)
  )

figura_univariante

# Variables candidatas (p < 0,20 en univariante)
df_multi <- df_uni %>%
  filter(
    !is.na(aclaramiento_dias),
    !is.na(gbghmix_f),
    !is.na(risk_f),
    !is.na(cvini),
    !is.na(cvpico),
    !is.na(Tdesdetx_meses),
    !is.na(sintom_f),
    !is.na(profilaxis_f),
    !is.na(ing_cmv_f)
  )

modelo_completo <- lm(
  aclaramiento_dias ~ gbghmix_f + risk_f + log10(cvini + 1) + log10(cvpico + 1) +
    Tdesdetx_meses + sintom_f + profilaxis_f + ing_cmv_f,
  data = df_multi
)

summary(modelo_completo)

# ── Selección backward ────────────────────────────────────────────────────────
modelo_backward <- step(modelo_completo, direction = "backward", trace = TRUE)
summary(modelo_backward)

# ── Intervalos de confianza y VIF ─────────────────────────────────────────────

confint(modelo_backward, level = 0.95)
vif(modelo_backward)


confint(modelo_backward, level = 0.95)

# ── 1. Datos completos para las variables del modelo clínico ──────────────────
df_clin_mix <- df_uni %>%
  select(NHC, aclaramiento_dias,
         gbghmix_f, risk_f, cvini, cvpico, profilaxis_f) %>%
  filter(complete.cases(.))

cat(sprintf("\n── N modelo clínico mixtos: %d episodios / %d pacientes ──\n",
            nrow(df_clin_mix), n_distinct(df_clin_mix$NHC)))

# ── 2. Modelo lineal múltiple ─────────────────────────────────────────────────
modelo_clinico_mix <- lm(
  aclaramiento_dias ~ gbghmix_f + risk_f + log10(cvini + 1) + log10(cvpico + 1) +
    profilaxis_f,
  data = df_clin_mix
)

cat("\n══════════════════════════════════════════════\n")
cat("  MODELO CLÍNICO MIXTOS gB+gH (variables a priori)\n")
cat("══════════════════════════════════════════════\n")
print(summary(modelo_clinico_mix))

cat("\n── IC 95% ──\n")
print(confint(modelo_clinico_mix))

cat("\n── VIF (multicolinealidad) ──\n")
print(car::vif(modelo_clinico_mix))

# ── 3. Efecto ajustado de los mixtos en una línea ─────────────────────────────
# El nombre del coeficiente depende del nivel: comprobamos cuál es
coef_name <- grep("gbghmix_f", names(coef(modelo_clinico_mix)), value = TRUE)
cat("\n── Coeficiente de interés:", coef_name, "──\n")

coef_mix <- coef(modelo_clinico_mix)[coef_name]
ci_mix   <- confint(modelo_clinico_mix)[coef_name, ]
p_mix    <- summary(modelo_clinico_mix)$coefficients[coef_name, "Pr(>|t|)"]

cat(sprintf("\n── Efecto ajustado mixtos gB+gH: %+.1f días (IC 95%%: %+.1f a %+.1f), p = %.4f ──\n",
            coef_mix, ci_mix[1], ci_mix[2], p_mix))

# ── 4. Análisis de sensibilidad: modelo mixto (paciente aleatorio) ────────────
modelo_mixto_mix <- lmer(
  aclaramiento_dias ~ gbghmix_f + risk_f + log10(cvini + 1) + log10(cvpico + 1) +
    profilaxis_f + (1 | NHC),
  data = df_clin_mix
)

cat("\n══════════════════════════════════════════════\n")
cat("  MODELO MIXTO: paciente como efecto aleatorio\n")
cat("══════════════════════════════════════════════\n")
print(summary(modelo_mixto_mix))

cat("\n── IC 95% (Wald) del efecto mixtos ──\n")
print(confint(modelo_mixto_mix, method = "Wald")[coef_name, ])

vc  <- as.data.frame(VarCorr(modelo_mixto_mix))
icc <- vc$vcov[1] / sum(vc$vcov)
cat(sprintf("\n── ICC (varianza atribuible al paciente): %.3f ──\n", icc))

guardar_modelo(modelo_clinico_mix, "Tabla_15_modelo")

