# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.7.06. Figura 35. Tasa de episodios de DNAemia en pacientes con genotipos únicos y mixtos
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 35. Tasa de episodios de DNAemia en pacientes con genotipos únicos y mixtos ----
# Incluye los modelos de regresión binomial negativa (IRR) citados en el texto.
# Resolver conflictos
select <- dplyr::select
if (exists("offset")) rm(offset)

# ── 1. Preparar datos IRR ─────────────────────────────────────────────────────
df_irr <- pacientes %>%
  dplyr::mutate(
    Ftx_dt       = as.Date(Ftx),
    Ffinseg_dt   = as.Date(Ffinseg),
    persona_anos = as.numeric(Ffinseg_dt - Ftx_dt) / 365.25,
    risk_f       = factor(risk, levels = c("MEDIO", "ALTO")),
    gbghm_grupo  = factor(gbghm, levels = c(0,1),
                          labels = c("Genotipos únicos",
                                     "Genotipos mixtos"))
  ) %>%
  dplyr::filter(persona_anos > 0)

# ── 2. Estadísticos descriptivos ──────────────────────────────────────────────
cat("── Total episodios CMV ──\n")
cat(sprintf("n = %d\n", sum(pacientes$nepi_tot, na.rm = TRUE)))

cat("\n── Tasas por grupo ──\n")
df_irr %>%
  dplyr::group_by(gbghm_grupo) %>%
  dplyr::summarise(
    n               = n(),
    total_episodios = sum(nepi_tot, na.rm = TRUE),
    total_pa        = round(sum(persona_anos, na.rm = TRUE), 1),
    tasa_pa         = round(sum(nepi_tot) / sum(persona_anos), 3),
    media_ep        = round(mean(nepi_tot, na.rm = TRUE), 2),
    mediana_ep      = median(nepi_tot, na.rm = TRUE),
    .groups         = "drop"
  ) %>% print()

# ── 3. Sobredispersión ────────────────────────────────────────────────────────
cat("\n── Sobredispersión ──\n")
df_irr %>%
  dplyr::group_by(gbghm_grupo) %>%
  dplyr::summarise(
    media    = round(mean(nepi_tot), 2),
    varianza = round(var(nepi_tot), 2),
    ratio    = round(var(nepi_tot) / mean(nepi_tot), 2),
    .groups  = "drop"
  ) %>% print()

m_pois <- glm(nepi_tot ~ gbghm_grupo + offset(log(persona_anos)),
              family = poisson, data = df_irr)
cat(sprintf("Poisson deviance/df = %.2f\n",
            m_pois$deviance / m_pois$df.residual))

# ── 4. Modelos binomial negativa ──────────────────────────────────────────────
cat("\n── IRR crudo (binomial negativa) ──\n")
m_nb <- glm.nb(nepi_tot ~ gbghm_grupo + offset(log(persona_anos)),
               data = df_irr)
print(round(exp(cbind(IRR = coef(m_nb), confint(m_nb))), 3))
cat("\np-valor gbghm:\n")
print(round(summary(m_nb)$coefficients, 3))

cat("\n── IRR ajustado (binomial negativa) ──\n")
m_nb_adj <- glm.nb(nepi_tot ~ gbghm_grupo + edad_tx + risk_f +
                     offset(log(persona_anos)), data = df_irr)
print(round(exp(cbind(IRR = coef(m_nb_adj), confint(m_nb_adj))), 3))
cat("\nCoeficientes con p-valores:\n")
print(round(summary(m_nb_adj)$coefficients, 3))

select <- dplyr::select

# ── 1. Preparar datos base (a nivel paciente) ────────────────────────────────
df_irr3 <- pacientes %>%
  dplyr::mutate(
    Ftx_dt       = as.Date(Ftx),
    Ffinseg_dt   = as.Date(Ffinseg),
    persona_anos = as.numeric(Ffinseg_dt - Ftx_dt) / 365.25,
    risk_f       = factor(risk, levels = c("MEDIO", "ALTO"))
  ) %>%
  dplyr::filter(persona_anos > 0)

# ── 2. Funcion: IRR ajustado para una variable de grupo ──────────────────────
calc_irr <- function(data, var_grupo) {
  data$grupo <- data[[var_grupo]]
  m <- glm.nb(nepi_tot ~ grupo + edad_tx + risk_f + offset(log(persona_anos)),
              data = data)
  irr <- round(exp(coef(m)[2]), 2)
  ci  <- round(exp(confint.default(m)[2, ]), 2)   # Wald (rapido y estable)
  p   <- summary(m)$coefficients[2, 4]
  list(irr = irr, ci_lo = ci[1], ci_hi = ci[2], p = p)
}

# ── 3. Calcular IRR para las 3 glucoproteinas ────────────────────────────────
vars   <- c(gbm = "gbm", ghm = "ghm", gbghm = "gbghm")
labs_g <- c(gbm = "gB", ghm = "gH", gbghm = "gB/gH")

res_irr <- lapply(names(vars), function(v) {
  d <- df_irr3 %>%
    dplyr::mutate(grupo = factor(.data[[v]], levels = c(0, 1),
                                 labels = c("Únicos", "Mixtos"))) %>%
    dplyr::filter(!is.na(grupo))
  r <- calc_irr(d %>% dplyr::mutate(!!v := grupo), v)
  data.frame(
    glucoproteina = labs_g[[v]],
    irr = r$irr, ci_lo = r$ci_lo, ci_hi = r$ci_hi, p = r$p
  )
}) %>% bind_rows() %>%
  dplyr::mutate(
    glucoproteina = factor(glucoproteina, levels = c("gB", "gH", "gB/gH")),
    label_irr = sprintf("IRR = %.2f (%.2f\u2013%.2f)\np = %s",
                        irr, ci_lo, ci_hi,
                        ifelse(p < 0.001, "<0,001",
                               sub("\\.", ",", sprintf("%.3f", p))))
  )

cat("── IRR ajustados ──\n"); print(res_irr)

# ── 4. Datos en formato largo: paciente x glucoproteina ──────────────────────
df_long3 <- df_irr3 %>%
  dplyr::select(NHC, nepi_tot, persona_anos, gbm, ghm, gbghm) %>%
  tidyr::pivot_longer(c(gbm, ghm, gbghm),
                      names_to = "gen", values_to = "mixto") %>%
  dplyr::filter(!is.na(mixto)) %>%
  dplyr::mutate(
    glucoproteina = factor(gen, levels = c("gbm", "ghm", "gbghm"),
                           labels = c("gB", "gH", "gB/gH")),
    grupo     = factor(ifelse(mixto == 1, "Mixtos", "Únicos"),
                       levels = c("Únicos", "Mixtos")),
    cat_color = paste(glucoproteina,
                      ifelse(mixto == 1, "Mixto", "Unico"), sep = "_")
  )

# ── 5. Tasas (ep/paciente-año) por panel y grupo ─────────────────────────────
tasas3 <- df_long3 %>%
  dplyr::group_by(glucoproteina, grupo) %>%
  dplyr::summarise(
    tasa = round(sum(nepi_tot) / sum(persona_anos), 2),
    .groups = "drop"
  )

# ── 6. n por panel y grupo (para etiquetas del eje) ──────────────────────────
n_grupo <- df_long3 %>%
  dplyr::count(glucoproteina, grupo)

# ── 7. Paletas (gB azul | gH rosa | gB/gH verde-rojo) ────────────────────────
cols_fill <- c(
  "gB_Unico"    = "#AED6F1", "gB_Mixto"    = "#4682B4",
  "gH_Unico"    = "pink",    "gH_Mixto"    = "hotpink3",
  "gB/gH_Unico" = "darkseagreen3", "gB/gH_Mixto" = "#c0392b"
)
cols_jit <- c(
  "gB_Unico"    = "#2980B9",        "gB_Mixto"    = "#0D2B45",
  "gH_Unico"    = "palevioletred3", "gH_Mixto"    = "palevioletred4",
  "gB/gH_Unico" = "#698B69",        "gB/gH_Mixto" = "#922b21"
)

# ── 8. Figura ────────────────────────────────────────────────────────────────
figura_R35_triple <- ggplot(df_long3,
                            aes(x = grupo, y = nepi_tot)) +
  geom_boxplot(aes(fill = cat_color),
               width = 0.45, outlier.shape = NA, alpha = 0.8,
               color = "black", linewidth = 0.6) +
  geom_jitter(aes(color = cat_color), width = 0.15, size = 1.8,
              alpha = 0.6, shape = 16, show.legend = FALSE) +
  # Etiqueta de tasa sobre cada caja
  geom_label(data = tasas3,
             aes(x = grupo, y = 10,
                 label = paste0(gsub("\\.", ",", tasa))),
             fill = "white", color = "black", size = 4,
             label.size = 0.25, fontface = "bold") +
  # Barra de significacion + IRR por panel
  
  geom_text(data = res_irr,
            aes(x = 1.5, y = 13, label = label_irr),
            size = 4, color = "black", fontface = "italic",
            inherit.aes = FALSE) +
  facet_wrap(~ glucoproteina) +
  scale_fill_manual(values = cols_fill, guide = "none") +
  scale_color_manual(values = cols_jit, guide = "none") +
  scale_y_continuous(
    name   = "Número de episodios de DNAemia",
    breaks = seq(0, 12, 2),
    limits = c(0, 15),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    panel.spacing.x    = unit(14, "pt"),
    axis.title.y       = element_text(size = 13, color = "black",
                                      margin = margin(r = 10)),
    axis.text.x        = element_text(size = 12, color = "black"),
    axis.text.y        = element_text(size = 12, color = "black"),
    axis.line          = element_line(color = "black", linewidth = 0.4),
    axis.ticks         = element_line(color = "black", linewidth = 0.3),
    strip.background   = element_blank(),
    strip.text         = element_text(size = 14, color = "black",
                                      face = "bold",
                                      margin = margin(t = 2, b = 6)),
    panel.grid         = element_blank(),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    plot.margin        = margin(15, 16, 10, 12)
  )

print(figura_R35_triple)

# ════════════════════════════════════════════════════════════════════════════
#  DATOS NUMERICOS — Figura R35 triple (episodios CMV por gB, gH, gB/gH)
# ════════════════════════════════════════════════════════════════════════════

# ── 1. n, tasa, media, mediana (RIC) por panel y grupo ───────────────────────
cat("\n========== DESCRIPTIVO POR GRUPO ==========\n")
tabla_desc <- df_long3 %>%
  dplyr::group_by(glucoproteina, grupo) %>%
  dplyr::summarise(
    n          = dplyr::n(),
    episodios  = sum(nepi_tot),
    persona_a  = round(sum(persona_anos), 1),
    tasa_pa    = round(sum(nepi_tot) / sum(persona_anos), 2),
    media      = round(mean(nepi_tot), 2),
    mediana    = median(nepi_tot),
    q1         = quantile(nepi_tot, 0.25),
    q3         = quantile(nepi_tot, 0.75),
    .groups    = "drop"
  ) %>%
  dplyr::mutate(
    `mediana (RIC)` = sprintf("%g (%g\u2013%g)", mediana, q1, q3)
  )
print(tabla_desc, n = Inf)

# ── 2. IRR ajustados con IC y p (formato tesis) ──────────────────────────────
cat("\n========== IRR AJUSTADOS (NB: edad + risk) ==========\n")
res_irr %>%
  dplyr::mutate(
    `IRR (IC95%)` = sprintf("%.2f (%.2f\u2013%.2f)", irr, ci_lo, ci_hi),
    `p`           = ifelse(p < 0.001, "<0,001",
                           sub("\\.", ",", sprintf("%.3f", p)))
  ) %>%
  dplyr::select(glucoproteina, `IRR (IC95%)`, `p`) %>%
  print()

# ── 3. p crudas de comparacion de distribucion (Wilcoxon) por panel ──────────
#     (util si quieres contrastar tasa vs distribucion de episodios)
cat("\n========== WILCOXON nepi_tot ~ grupo (crudo) ==========\n")
df_long3 %>%
  dplyr::group_by(glucoproteina) %>%
  dplyr::summarise(
    p_raw = wilcox.test(nepi_tot ~ grupo, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_fmt = ifelse(p_raw < 0.001, "<0,001",
                   sub("\\.", ",", sprintf("%.3f", p_raw)))
  ) %>%
  print()

guardar_fig(figura_R35_triple, "Figura_35", ancho = 18, alto = 10)

