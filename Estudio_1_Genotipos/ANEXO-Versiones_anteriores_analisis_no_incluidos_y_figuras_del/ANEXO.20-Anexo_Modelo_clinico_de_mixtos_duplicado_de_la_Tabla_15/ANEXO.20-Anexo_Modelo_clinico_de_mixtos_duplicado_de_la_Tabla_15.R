# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.20. Anexo: Modelo clínico de mixtos (duplicado de la Tabla 15)
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Modelo clínico de mixtos (duplicado de la Tabla 15) ----
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

