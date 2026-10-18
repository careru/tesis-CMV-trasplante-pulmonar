# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.7.02. Tabla 19. Asociación entre genotipos mixtos y enfermedad sintomática por CMV
# Extracto del script completo (ver README.md del estudio)

## ---- Tabla 19. Asociación entre genotipos mixtos y enfermedad sintomática por CMV ----
# ── 1. Preparar datos a nivel paciente ────────────────────────────────────────
pacientes_sintom <- pacientes %>%
  mutate(
    sintomatico = as.integer(vira == 0),
    gbghm_f     = factor(gbghm, levels = c(0,1), labels = c("Único","Mixto")),
    gbm_f       = factor(gbm,   levels = c(0,1), labels = c("Único","Mixto")),
    ghm_f       = factor(ghm,   levels = c(0,1), labels = c("Único","Mixto")),
    serDR_f     = factor(serDR)
  )

cat("── Distribución sintomático vs grupo ──\n")
pacientes_sintom %>%
  count(sintomatico, gbghm_f) %>%
  group_by(gbghm_f) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()

# ── 2. Funciones ──────────────────────────────────────────────────────────────
calc_or <- function(var_mix, data, label) {
  m_crudo <- glm(sintomatico ~ get(var_mix),
                 data = data, family = binomial)
  or_c <- exp(coef(m_crudo)[2])
  ci_c <- exp(confint.default(m_crudo)[2,])
  p_c  <- summary(m_crudo)$coefficients[2,4]
  
  m_ajus <- glm(sintomatico ~ get(var_mix) + serDR_f,
                data = data, family = binomial)
  or_a <- exp(coef(m_ajus)[2])
  ci_a <- exp(confint.default(m_ajus)[2,])
  p_a  <- summary(m_ajus)$coefficients[2,4]
  
  cat(sprintf("\n── %s ──\n", label))
  cat(sprintf("OR cruda:    %.2f (%.2f–%.2f)  p = %.3f\n",
              or_c, ci_c[1], ci_c[2], p_c))
  cat(sprintf("OR ajustada: %.2f (%.2f–%.2f)  p = %.3f\n",
              or_a, ci_a[1], ci_a[2], p_a))
  
  tibble(
    genotipo    = label,
    or_cruda    = sprintf("%.2f (%.2f\u2013%.2f)", or_c, ci_c[1], ci_c[2]),
    p_cruda     = ifelse(p_c < 0.001, "<0,001",
                         sub("\\.", ",", sprintf("%.3f", p_c))),
    or_ajustada = sprintf("%.2f (%.2f\u2013%.2f)", or_a, ci_a[1], ci_a[2]),
    p_ajustada  = ifelse(p_a < 0.001, "<0,001",
                         sub("\\.", ",", sprintf("%.3f", p_a)))
  )
}

calc_or_mixto <- function(var_mix, data, label) {
  m <- tryCatch(
    glmer(sintomatico ~ get(var_mix) + (1|NHC),
          data = data, family = binomial,
          control = glmerControl(optimizer = "bobyqa")),
    error = function(e) NULL
  )
  if (is.null(m)) {
    cat(sprintf("\n── %s: error en modelo mixto ──\n", label))
    return(tibble(genotipo = label, or_mixto = "error", p_mixto = "error"))
  }
  coefs <- summary(m)$coefficients
  or    <- exp(coefs[2, 1])
  se    <- coefs[2, 2]
  ci_lo <- exp(coefs[2, 1] - 1.96 * se)
  ci_hi <- exp(coefs[2, 1] + 1.96 * se)
  p     <- coefs[2, 4]
  
  cat(sprintf("\n── %s (modelo mixto) ──\n", label))
  cat(sprintf("OR: %.2f (%.2f–%.2f)  p = %.3f\n", or, ci_lo, ci_hi, p))
  
  tibble(
    genotipo = label,
    or_mixto = sprintf("%.2f (%.2f\u2013%.2f)", or, ci_lo, ci_hi),
    p_mixto  = ifelse(p < 0.001, "<0,001",
                      sub("\\.", ",", sprintf("%.3f", p)))
  )
}

freq_gt <- function(var, label, data) {
  tab    <- table(data[[var]], data$sintomatico)
  test   <- fisher.test(tab, simulate.p.value = TRUE, B = 10000)
  n0     <- sum(data$sintomatico == 0)
  n1     <- sum(data$sintomatico == 1)
  n_var0 <- sum(data[[var]] == 1 & data$sintomatico == 0, na.rm = TRUE)
  n_var1 <- sum(data[[var]] == 1 & data$sintomatico == 1, na.rm = TRUE)
  tibble(
    genotipo = label,
    asintom  = sprintf("%d (%.1f%%)", n_var0, 100*n_var0/n0),
    sintom   = sprintf("%d (%.1f%%)", n_var1, 100*n_var1/n1),
    p        = ifelse(test$p.value < 0.001, "<0,001",
                      sub("\\.", ",", sprintf("%.3f", test$p.value)))
  )
}

# ── 3. OR cruda y ajustada a nivel paciente ───────────────────────────────────
res_or <- bind_rows(
  calc_or("gbm_f",   pacientes_sintom, "gB mixtos"),
  calc_or("ghm_f",   pacientes_sintom, "gH mixtos"),
  calc_or("gbghm_f", pacientes_sintom, "gB+gH mixtos")
)

cat("\n── Tabla OR resumen ──\n")
print(res_or)

# ── 4. Preparar datos a nivel episodio ────────────────────────────────────────
dnaemias_sintom <- dnaemias %>%
  mutate(
    sintomatico = as.integer(!is.na(clinica) & clinica != "VA"),
    gbghm_f     = factor(gbghmix, levels = c(0,1), labels = c("Único","Mixto")),
    gbm_f       = factor(gbmix,   levels = c(0,1), labels = c("Único","Mixto")),
    ghm_f       = factor(ghmix,   levels = c(0,1), labels = c("Único","Mixto"))
  ) %>%
  left_join(
    pacientes %>% select(NHC, serDR_pac = serDR),
    by = "NHC"
  ) %>%
  mutate(serDR_f = factor(serDR_pac))

cat(sprintf("\nTotal episodios: %d\n", nrow(dnaemias_sintom)))
cat(sprintf("Sintomáticos: %d\n",     sum(dnaemias_sintom$sintomatico)))

# ── 5. Modelo mixto por episodio ──────────────────────────────────────────────
res_mixto <- bind_rows(
  calc_or_mixto("gbm_f",   dnaemias_sintom, "gB mixtos"),
  calc_or_mixto("ghm_f",   dnaemias_sintom, "gH mixtos"),
  calc_or_mixto("gbghm_f", dnaemias_sintom, "gB+gH mixtos")
)

cat("\n── Tabla OR mixto resumen ──\n")
print(res_mixto)

# ── 6. Frecuencias genotípicas por episodio ───────────────────────────────────
res_freq <- bind_rows(
  freq_gt("gbmix",   "gB mixtos",    dnaemias_sintom),
  freq_gt("ghmix",   "gH mixtos",    dnaemias_sintom),
  freq_gt("gbghmix", "gB+gH mixtos", dnaemias_sintom)
)

cat("\n── Frecuencias genotípicas ──\n")
print(res_freq)


# ── Combinar resultados ───────────────────────────────────────────────────────
tabla_sintom <- res_or %>%
  left_join(res_mixto, by = "genotipo") %>%
  rename(
    "Genotipos mixtos"          = genotipo,
    "OR cruda\n(IC 95%)"        = or_cruda,
    "p"                         = p_cruda,
    "OR ajustada*\n(IC 95%)"    = or_ajustada,
    "p "                        = p_ajustada,
    "OR por episodio**\n(IC 95%)"= or_mixto,
    "p  "                       = p_mixto
  )

ft_sintom <- flextable(tabla_sintom) %>%
  width(j = "Genotipos mixtos",           width = 1.6) %>%
  width(j = "OR cruda\n(IC 95%)",         width = 2.0) %>%
  width(j = "p",                          width = 0.8) %>%
  width(j = "OR ajustada*\n(IC 95%)",     width = 2.0) %>%
  width(j = "p ",                         width = 0.8) %>%
  width(j = "OR por episodio**\n(IC 95%)",width = 2.0) %>%
  width(j = "p  ",                        width = 0.8) %>%
  align(align = "center", part = "all") %>%
  align(j = "Genotipos mixtos", align = "left", part = "all") %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header") %>%
  bg(bg = "white", part = "all") %>%
  bg(i = seq(2, nrow(tabla_sintom), 2), bg = "#F7F7F7", part = "body") %>%
  border_remove() %>%
  hline_top(part = "header",
            border = fp_border(color = "black", width = 1.2)) %>%
  hline(part = "header",
        border = fp_border(color = "black", width = 0.8)) %>%
  hline_bottom(part = "body",
               border = fp_border(color = "black", width = 1.2)) %>%
  padding(padding.top = 3, padding.bottom = 3,
          padding.left = 4, padding.right = 4, part = "all") %>%
  add_footer_lines(values = c(
    "*Ajustada por perfil serológico D/R.",
    "**Modelo de regresión logística mixta con efecto aleatorio por paciente."
  )) %>%
  fontsize(size = 8, part = "footer") %>%
  color(color = "grey30", part = "footer") %>%
  set_caption(
    caption = "Tabla R12. Asociación entre genotipos mixtos gB/gH y enfermedad sintomática por CMV."
  ) %>%
  set_table_properties(layout = "autofit")

ft_sintom

guardar_tabla(ft_sintom, "Tabla_19")

