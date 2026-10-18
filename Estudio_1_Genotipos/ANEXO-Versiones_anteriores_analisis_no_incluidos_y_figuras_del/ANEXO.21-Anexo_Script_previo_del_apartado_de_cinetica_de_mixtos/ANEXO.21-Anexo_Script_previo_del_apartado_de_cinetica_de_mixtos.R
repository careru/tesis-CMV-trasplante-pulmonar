# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.21. Anexo: Script previo del apartado de cinética de mixtos (Tablas R9-R11, Figuras R20-R23)
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Script previo del apartado de cinética de mixtos (Tablas R9-R11, Figuras R20-R23) ----
#  APARTADO 4.5.3 — IMPACTO DE LOS GENOTIPOS MIXTOS EN LA CINÉTICA VIRAL
#
#  Tesis CMV en RTP — Carmen, HU12O
#
#  Contenido del script:
#    1. Preparación de datos
#    2. Tabla R9   — Cinética viral: únicos vs mixtos (3 glucoproteínas)
#    3. Figura R20 — Boxplots cinética estratificada por serostatus
#    4. Tabla R10  — Cinética estratificada por serostatus (p-valores)
#    5. Análisis primer episodio por paciente
#    6. Figura R21 — Tiempo TX → primer episodio
#    7. Figura R22 — Univariante mixtas
#    8. Tabla R11  — Multivariante mixtas (backward + sensibilidad mixta)
#    9. Figura R23 — Influencia de gB3 en episodios mixtos
#
# ════════════════════════════════════════════════════════════════════════════


# ════════════════════════════════════════════════════════════════════════════
#  1. PREPARACIÓN DE DATOS
# ════════════════════════════════════════════════════════════════════════════

dnaemias <- dnaemias %>%
  mutate(
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini)
  )


# ════════════════════════════════════════════════════════════════════════════
#  2. TABLA R9 — Cinética viral: únicos vs mixtos (3 glucoproteínas)
# ════════════════════════════════════════════════════════════════════════════

fmt_med <- function(x) {
  n   <- sum(!is.na(x))
  med <- median(x, na.rm = TRUE)
  q1  <- quantile(x, 0.25, na.rm = TRUE)
  q3  <- quantile(x, 0.75, na.rm = TRUE)
  sprintf("%s (%s\u2013%s)\n(n=%d)",
          formatC(med, format="f", digits=0),
          formatC(q1,  format="f", digits=0),
          formatC(q3,  format="f", digits=0), n)
}

fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0,001*")
  txt <- sprintf("%.3f", p)
  txt <- sub("\\.", ",", txt)
  if (p < 0.05) paste0(txt, "*") else txt
}

wilcox_p <- function(data, var, grupo) {
  d <- data %>% filter(!is.na(.data[[var]]), !is.na(.data[[grupo]]))
  if (nrow(d) < 4) return(NA_real_)
  tryCatch(
    wilcox.test(as.formula(paste(var, "~", grupo)), data = d)$p.value,
    error = function(e) NA_real_
  )
}

bloque <- function(data, var_mix) {
  d0 <- data %>% filter(.data[[var_mix]] == 0)
  d1 <- data %>% filter(.data[[var_mix]] == 1)
  
  vars  <- c("cvini", "cvpico", "aclaramiento_dias", "AUCCV")
  labs  <- c("CV inicial (UI/mL)", "CV pico (UI/mL)",
             "Aclaramiento (días)", "AUC viral")
  
  tibble(
    Parámetro        = labs,
    `Genotipo único` = sapply(vars, function(v) fmt_med(d0[[v]])),
    `Genotipo mixto` = sapply(vars, function(v) fmt_med(d1[[v]])),
    p                = sapply(vars, function(v) fmt_p(wilcox_p(data, v, var_mix)))
  )
}

cabecera <- function(txt) tibble(
  Parámetro        = txt,
  `Genotipo único` = NA_character_,
  `Genotipo mixto` = NA_character_,
  p                = NA_character_
)

tabla_R9 <- bind_rows(
  cabecera("Glucoproteína gB"),
  bloque(dnaemias, "gbmix"),
  cabecera("Glucoproteína gH"),
  bloque(dnaemias, "ghmix"),
  cabecera("Glucoproteínas gB y gH"),
  bloque(dnaemias, "gbghmix")
)

idx_cab_R9 <- which(tabla_R9$Parámetro %in%
                      c("Glucoproteína gB", "Glucoproteína gH",
                        "Glucoproteínas gB y gH"))

ft_R9 <- flextable(tabla_R9) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  bold(i = idx_cab_R9, j = 1, part = "body") %>%
  bg(i = idx_cab_R9, bg = "#F2F2F2", part = "body") %>%
  align(j = "Parámetro", align = "left", part = "all") %>%
  align(j = c("Genotipo único", "Genotipo mixto", "p"),
        align = "center", part = "all") %>%
  padding(i = setdiff(seq_len(nrow(tabla_R9)), idx_cab_R9),
          j = 1, padding.left = 14, part = "body") %>%
  autofit() %>%
  add_footer_lines(values = c(
    "Los datos se presentan como mediana (RIC).",
    "CV: carga viral; AUC: área bajo la curva.",
    "Comparaciones mediante la prueba U de Mann-Whitney.",
    "*p < 0,05."
  )) %>%
  fontsize(part = "footer", size = 8) %>%
  align(part = "footer", align = "left")

ft_R9


# ════════════════════════════════════════════════════════════════════════════
#  3. FIGURA R20 — Boxplots cinética por SEROSTATUS
# ════════════════════════════════════════════════════════════════════════════
##  Versión corregida: sustituido recode() por case_when()

id_patient <- "NHC"
cv_y_max   <- 7

vars_to_pivot <- c("diasdur_vir", "AUCCV", "cvpico", "cvini")

df_long <- dnaemias %>%
  mutate(
    gB   = ifelse(gbmix   == 1, "gB mixto",    "gB único"),
    gH   = ifelse(ghmix   == 1, "gH mixto",    "gH único"),
    gBGH = ifelse(gbghmix == 1, "gB+gH mixto", "gB+gH único")
  ) %>%
  pivot_longer(cols = c(gB, gH, gBGH),
               names_to  = "tipo_gen",
               values_to = "genotipo") %>%
  select(serDR, tipo_gen, genotipo,
         !!sym(id_patient), all_of(vars_to_pivot)) %>%
  pivot_longer(cols      = all_of(vars_to_pivot),
               names_to  = "variable",
               values_to = "valor_raw") %>%
  mutate(
    valor = case_when(
      variable %in% c("cvpico", "cvini") & !is.na(valor_raw) ~
        log10(as.numeric(valor_raw) + 1),
      TRUE ~ as.numeric(valor_raw)
    ),
    # ── CORREGIDO: recode() reemplazado por case_when() ─────────────────────
    variable = factor(case_when(
      variable == "diasdur_vir" ~ "Aclaramiento (días)",
      variable == "AUCCV"       ~ "AUC (UI/mL\u00b7día)",
      variable == "cvpico"      ~ "log10 CV pico (UI/mL)",
      variable == "cvini"       ~ "log10 CV inicial (UI/mL)"
    ), levels = c(
      "log10 CV inicial (UI/mL)",
      "log10 CV pico (UI/mL)",
      "Aclaramiento (días)",
      "AUC (UI/mL\u00b7día)"
    )),
    serDR    = factor(serDR, levels = c("D-/R+", "D+/R+", "D+/R-")),
    genotipo = factor(genotipo, levels = c(
      "gB único",    "gB mixto",
      "gH único",    "gH mixto",
      "gB+gH único", "gB+gH mixto"
    ))
  )

x_map <- c(
  "gB único"    = 1,   "gB mixto"    = 2,
  "gH único"    = 3.3, "gH mixto"    = 4.3,
  "gB+gH único" = 5.6, "gB+gH mixto" = 6.6
)

df_plot <- df_long %>%
  mutate(x_pos = x_map[as.character(genotipo)])

caps_dyn <- tibble(
  variable = c("log10 CV inicial (UI/mL)", "log10 CV pico (UI/mL)",
               "Aclaramiento (días)", "AUC (UI/mL\u00b7día)"),
  cap_vis  = c(cv_y_max, cv_y_max, 100, 700)
)

df_plot_vis <- df_plot %>%
  left_join(caps_dyn, by = "variable") %>%
  filter(is.na(valor) | valor <= cap_vis)

comparaciones <- df_plot %>%
  group_by(variable, serDR, tipo_gen) %>%
  wilcox_test(valor ~ genotipo) %>%
  add_significance("p") %>%
  ungroup() %>%
  filter(p.signif != "ns") %>%
  mutate(
    group1 = case_when(
      tipo_gen == "gB"   ~ "gB único",
      tipo_gen == "gH"   ~ "gH único",
      tipo_gen == "gBGH" ~ "gB+gH único"
    ),
    group2 = case_when(
      tipo_gen == "gB"   ~ "gB mixto",
      tipo_gen == "gH"   ~ "gH mixto",
      tipo_gen == "gBGH" ~ "gB+gH mixto"
    ),
    xmin = x_map[group1],
    xmax = x_map[group2]
  )

base_positions <- df_plot_vis %>%
  group_by(variable, serDR) %>%
  summarise(maxval  = max(valor, na.rm = TRUE),
            cap_vis = first(cap_vis), .groups = "drop") %>%
  mutate(base = pmin(maxval * 1.12, cap_vis * 0.90))

comp_pos <- comparaciones %>%
  left_join(base_positions, by = c("variable", "serDR")) %>%
  group_by(variable, serDR) %>%
  arrange(p) %>%
  mutate(idx = row_number(),
         y   = base + (idx - 1) * 0.05 * cap_vis) %>%
  ungroup()

counts_serDR <- dnaemias %>%
  group_by(serDR) %>%
  summarise(
    episodios = n(),
    pacientes = n_distinct(.data[[id_patient]]),
    .groups   = "drop"
  ) %>%
  mutate(label = paste0(
    "<b>", serDR, "</b><br>",
    "<span style='font-size:10pt'>",
    pacientes, " pac. / ", episodios, " ep.</span>"
  ))

serdr_labels        <- counts_serDR$label
names(serdr_labels) <- counts_serDR$serDR

colores <- c(
  "gB único"    = "#63B8FF",
  "gB mixto"    = "dodgerblue3",
  "gH único"    = "#DDA0DD",
  "gH mixto"    = "#CD69C9",
  "gB+gH único" = "#3CB371",
  "gB+gH mixto" = "#F08080"
)

figura_R20_estratificado <- ggplot(df_plot_vis,
                                   aes(x = x_pos, y = valor, fill = genotipo)) +
  geom_boxplot(
    width         = 0.55,
    color         = "grey25",
    linewidth     = 0.4,
    outlier.shape = 21,
    outlier.size  = 1.4,
    outlier.fill  = "white",
    outlier.color = "grey50",
    outlier.alpha = 0.6
  ) +
  geom_vline(xintercept = c(2.65, 4.95),
             color = "grey80", linewidth = 0.35, linetype = "dashed") +
  facet_grid(variable ~ serDR, scales = "free_y",
             labeller = labeller(serDR = serdr_labels)) +
  geom_text(data = comp_pos,
            aes(x = (xmin + xmax) / 2, y = y, label = p.signif),
            inherit.aes = FALSE, size = 3.8, vjust = 0, color = "grey20") +
  geom_segment(data = comp_pos,
               aes(x    = xmin, xend = xmax,
                   y    = y - 0.012 * cap_vis,
                   yend = y - 0.012 * cap_vis),
               inherit.aes = FALSE,
               linewidth = 0.45, color = "grey30") +
  scale_fill_manual(values = colores) +
  scale_x_continuous(
    breaks = unname(x_map),
    labels = names(x_map),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position    = "none",
    panel.border       = element_rect(color = "grey70", fill = NA,
                                      linewidth = 0.4),
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.line          = element_blank(),
    axis.ticks         = element_line(color = "grey60", linewidth = 0.3),
    axis.text.x        = element_text(size = 10, color = "grey25",
                                      angle = 45, hjust = 1),
    axis.text.y        = element_text(size = 12, color = "grey30"),
    strip.background   = element_rect(fill = "grey95", color = "grey70",
                                      linewidth = 0.4),
    strip.text.x       = element_markdown(size = 13, face = "bold",
                                          color = "grey20"),
    strip.text.y       = element_text(size = 12, face = "bold",
                                      color = "grey20"),
    panel.spacing.x    = unit(6,  "pt"),
    panel.spacing.y    = unit(10, "pt"),
    plot.margin        = margin(8, 10, 6, 6)
  )

ggsave("~/Desktop/Figuras_tesis/figura_R20_cinetica_estratificado.pdf",
       plot   = figura_R20_estratificado,
       width  = 17, height = 22, units  = "cm", device = "pdf")

figura_R20_estratificado


# ════════════════════════════════════════════════════════════════════════════
#  4. TABLA R10 — Cinética estratificada por serostatus (p-valores)
# ════════════════════════════════════════════════════════════════════════════

vars <- c("cvini", "cvpico", "diasdur_vir", "AUCCV")
var_labels <- c(
  cvini       = "CV inicial",
  cvpico      = "CV pico",
  diasdur_vir = "Aclaramiento",
  AUCCV       = "AUC viral"
)

calcular_p <- function(data_in, mixvar, mix_label) {
  map_dfr(vars, function(v) {
    d0 <- data_in %>% filter(.data[[mixvar]] == 0) %>% pull(.data[[v]])
    d1 <- data_in %>% filter(.data[[mixvar]] == 1) %>% pull(.data[[v]])
    p_global <- tryCatch(
      wilcox.test(d0, d1, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
    grupos_res <- map_dfr(c("D-/R+","D+/R+","D+/R-"), function(grupo) {
      dg  <- data_in %>% filter(serDR == grupo)
      d0g <- dg %>% filter(.data[[mixvar]] == 0) %>% pull(.data[[v]])
      d1g <- dg %>% filter(.data[[mixvar]] == 1) %>% pull(.data[[v]])
      p <- if (length(na.omit(d0g)) > 0 & length(na.omit(d1g)) > 0)
        tryCatch(wilcox.test(d0g, d1g, exact = FALSE)$p.value,
                 error = function(e) NA_real_)
      else NA_real_
      tibble(serDR = grupo, p = p)
    }) %>%
      pivot_wider(names_from = serDR, values_from = p, names_prefix = "p_")
    bind_cols(
      tibble(Glucoproteína = mix_label, Variable = var_labels[v], p_global),
      grupos_res
    )
  })
}

tabla_R10_raw <- bind_rows(
  calcular_p(dnaemias, "gbmix",   "gB"),
  calcular_p(dnaemias, "ghmix",   "gH"),
  calcular_p(dnaemias, "gbghmix", "gB+gH")
)

tabla_R10_fmt <- tabla_R10_raw %>%
  mutate(across(where(is.numeric), ~ sapply(., fmt_p)))

# (Construcción de la flextable como en tu código original)
# — omitida por brevedad, usar tu versión con tabla_R10_fmt —


# ════════════════════════════════════════════════════════════════════════════
#  5. ANÁLISIS PRIMER EPISODIO POR PACIENTE
# ════════════════════════════════════════════════════════════════════════════

##  Tiempo TX → primer episodio (estratificado por riesgo)
##  Análisis y test omitido por brevedad: usa tu código original
##  (sección "M" en tu script)


# ════════════════════════════════════════════════════════════════════════════
#  6. FIGURA R21 — Tiempo TX → primer episodio
# ════════════════════════════════════════════════════════════════════════════

##  Usa tu código original de la sección [M]
##  Solo cambiar el nombre de salida:
##  ggsave("~/Desktop/Figuras_tesis/figura_R21_tiempo_primer_episodio.png", ...)


# ════════════════════════════════════════════════════════════════════════════
#  7. FIGURA R22 — Univariante mixtas (forest plot tipo tabla con cajas)
# ════════════════════════════════════════════════════════════════════════════
##  Reescrito siguiendo el estilo de la Figura R19 de gB3
##  (cajas por dominio, columna Estimador, forest plot, columna p)

# 7.1 Preparar datos (con sexo y edad del donante) ───────────────────────────
df_uni_mix <- dnaemias %>%
  left_join(pacientes %>%
              select(NHC, edad_tx, tipotx, ebase, SexF, edad_don),
            by = "NHC") %>%
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
    gbghmix_f    = factor(gbghmix,    levels = c(0,1), labels = c("Único","Mixto")),
    ingr_f       = factor(ingr,       levels = c(0,1), labels = c("No","Sí")),
    ing_cmv_f    = factor(ing_cmv,    levels = c(0,1), labels = c("No","Sí")),
    sintom_f     = factor(sintomatico, levels = c(0,1), labels = c("No","Sí")),
    profilaxis_f = factor(antiviral,  levels = c(0,1), labels = c("No","Sí")),
    tto_f        = factor(tto,        levels = c(0,1), labels = c("No","Sí")),
    tipotx_f     = factor(tipotx,     levels = c("UNI","BI")),
    ebase_f      = factor(ebase),
    sexo_f       = factor(SexF, levels = c(0,1), labels = c("Hombre","Mujer"))
  ) %>%
  filter(!is.na(aclaramiento_dias))

cat(sprintf("\n── N episodios mixtas: %d │ N pacientes: %d ──\n",
            nrow(df_uni_mix), n_distinct(df_uni_mix$NHC)))

# 7.2 Univariante ────────────────────────────────────────────────────────────
cor_vars <- list(
  "Edad receptor (años)"     = "edad_tx",
  "Edad donante (años)"      = "edad_don",
  "CV inicial (UI/mL)"       = "cvini",
  "CV pico (UI/mL)"          = "cvpico",
  "Tiempo desde TXP (meses)" = "Tdesdetx_meses"
)

res_cor <- map_dfr(names(cor_vars), function(lab) {
  v <- cor_vars[[lab]]
  x <- as.numeric(df_uni_mix[[v]])
  if (v %in% c("cvini","cvpico")) x <- log10(x + 1)
  idx <- !is.na(x) & !is.na(df_uni_mix$aclaramiento_dias)
  test <- cor.test(x[idx], df_uni_mix$aclaramiento_dias[idx],
                   method = "spearman", exact = FALSE)
  tibble(Variable = lab, p_val = test$p.value,
         label_est = sprintf("\u03C1 = %+.3f", test$estimate))
})

cat_vars <- list(
  "Genotipos mixtos gB/gH"            = "gbghmix_f",
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
  g  <- df_uni_mix[[v]]
  ac <- df_uni_mix$aclaramiento_dias
  d0 <- ac[g == levels(g)[1] & !is.na(g)]
  d1 <- ac[g == levels(g)[2] & !is.na(g)]
  test <- wilcox.test(d1, d0, exact = FALSE)
  dif  <- median(d1, na.rm = TRUE) - median(d0, na.rm = TRUE)
  tibble(Variable = lab, p_val = test$p.value,
         label_est = sprintf("\u0394 = %+.0f d\u00edas", dif))
})

kw_eb <- kruskal.test(aclaramiento_dias ~ ebase_f, data = df_uni_mix)
res_eb <- tibble(Variable = "Enfermedad de base",
                 p_val = kw_eb$p.value,
                 label_est = "(no calculable)")

# 7.3 Ordenar y categorizar (con dominios) ───────────────────────────────────
orden_vars <- c(
  "Edad receptor (años)", "Edad donante (años)", "Sexo (Mujer vs Hombre)",
  "Tipo de trasplante (BI vs UNI)", "Enfermedad de base", "D+/R- (vs R+)",
  "Tiempo desde TXP (meses)", "CV inicial (UI/mL)", "CV pico (UI/mL)",
  "Genotipos mixtos gB/gH",
  "Enfermedad sintomática", "Profilaxis al inicio del episodio",
  "Tratamiento antiviral", "Hospitalización cualquier causa",
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
  "Genotipos mixtos gB/gH"            = "Cinéticas",
  "Enfermedad sintomática"            = "Clínicas",
  "Profilaxis al inicio del episodio" = "Clínicas",
  "Tratamiento antiviral"             = "Clínicas",
  "Hospitalización cualquier causa"   = "Clínicas",
  "Hospitalización por CMV"           = "Clínicas"
)

res_all_mix <- bind_rows(res_cor, res_cat, res_eb) %>%
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

cat("\n── Univariante mixtas: resultados ordenados por p ──\n")
print(res_all_mix %>% arrange(p_val) %>%
        select(Variable, label_est, p_fmt), n = Inf)

# 7.4 Figura R22 (mismo estilo que R19) ──────────────────────────────────────
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
                    levels = levels(res_all_mix$Variable))
)

figura_R22_univariante_mixtas <- ggplot(res_all_mix,
                                        aes(x = p_val, y = Variable)) +
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

figura_R22_univariante_mixtas

ggsave("~/Desktop/Figuras_tesis/figura_R22_univariante_mixtas.png",
       figura_R22_univariante_mixtas,
       width = 9, height = 6, dpi = 320, bg = "white")


# ════════════════════════════════════════════════════════════════════════════
#  8. TABLA R11 — Multivariante mixtas (backward + análisis de sensibilidad)
# ════════════════════════════════════════════════════════════════════════════

# 8.1 df con variables candidatas (p < 0,20) ─────────────────────────────────
df_multi_mix <- df_uni_mix %>%
  select(NHC, aclaramiento_dias,
         gbghmix_f, risk_f, cvini, cvpico, Tdesdetx_meses,
         sintom_f, profilaxis_f, ing_cmv_f, ingr_f) %>%
  filter(complete.cases(.))

cat(sprintf("\n── N multivariante mixtas: %d episodios / %d pacientes ──\n",
            nrow(df_multi_mix), n_distinct(df_multi_mix$NHC)))

# 8.2 Modelo completo ─────────────────────────────────────────────────────────
modelo_completo_mix <- lm(
  aclaramiento_dias ~ gbghmix_f + risk_f + log10(cvini + 1) + log10(cvpico + 1) +
    Tdesdetx_meses + sintom_f + profilaxis_f + ing_cmv_f + ingr_f,
  data = df_multi_mix
)

cat("\n── Modelo completo ──\n")
print(summary(modelo_completo_mix))

# 8.3 Backward ────────────────────────────────────────────────────────────────
modelo_backward_mix <- step(modelo_completo_mix, direction = "backward",
                            trace = FALSE)

cat("\n══════════════════════════════════════════════\n")
cat("  MODELO FINAL MIXTAS (post-backward)\n")
cat("══════════════════════════════════════════════\n")
print(summary(modelo_backward_mix))

cat("\n── IC 95% ──\n")
print(confint(modelo_backward_mix))

cat("\n── VIF ──\n")
print(car::vif(modelo_backward_mix))


# 8.4 Análisis de sensibilidad: modelo lineal mixto ──────────────────────────
##  Solo si el modelo final del backward es similar al de gB3
##  (mismas variables: gbghmix + cvini + cvpico + profilaxis + ing_cmv)

modelo_mixto_mix <- lmer(
  aclaramiento_dias ~ gbghmix_f + log10(cvini + 1) + log10(cvpico + 1) +
    profilaxis_f + ing_cmv_f + (1 | NHC),
  data = df_multi_mix
)

cat("\n══════════════════════════════════════════════\n")
cat("  MODELO MIXTO MIXTAS: paciente como efecto aleatorio\n")
cat("══════════════════════════════════════════════\n")
print(summary(modelo_mixto_mix))

cat("\n── IC 95% (Wald) ──\n")
print(confint(modelo_mixto_mix, method = "Wald"))

vc_mix  <- as.data.frame(VarCorr(modelo_mixto_mix))
icc_mix <- vc_mix$vcov[1] / sum(vc_mix$vcov)
cat(sprintf("\n── ICC (mixtas): %.3f ──\n", icc_mix))

# 8.5 Comparación lineal vs mixto ────────────────────────────────────────────
coef_lin_mix <- coef(modelo_backward_mix)["gbghmix_fMixto"]
ci_lin_mix   <- confint(modelo_backward_mix)["gbghmix_fMixto", ]
p_lin_mix    <- summary(modelo_backward_mix)$coefficients["gbghmix_fMixto", "Pr(>|t|)"]

coef_mix_mix <- fixef(modelo_mixto_mix)["gbghmix_fMixto"]
ci_mix_mix   <- confint(modelo_mixto_mix, method = "Wald")["gbghmix_fMixto", ]
p_mix_mix    <- summary(modelo_mixto_mix)$coefficients["gbghmix_fMixto", "Pr(>|t|)"]

cat("\n── COMPARACIÓN LINEAL vs MIXTO (efecto mixtas) ──\n")
cat(sprintf("  Lineal:  %+.1f días (IC 95%%: %+.1f a %+.1f), p = %.4f\n",
            coef_lin_mix, ci_lin_mix[1], ci_lin_mix[2], p_lin_mix))
cat(sprintf("  Mixto:   %+.1f días (IC 95%%: %+.1f a %+.1f), p = %.4f\n",
            coef_mix_mix, ci_mix_mix[1], ci_mix_mix[2], p_mix_mix))


# ════════════════════════════════════════════════════════════════════════════
#  9. FIGURA R23 — Influencia de gB3 en episodios mixtos
# ════════════════════════════════════════════════════════════════════════════

##  Reusar el código de la sección [I] y [J] de tu script original
##  Solo cambiar el nombre de salida a:
##  ggsave("~/Desktop/Figuras_tesis/figura_R23_gB3_en_mixtos.png", ...)

# ────────────────────────────────────────────────────────────────────────────
#  TABLA R10 — Flextable: p-valores del análisis estratificado por serostatus
# ────────────────────────────────────────────────────────────────────────────

# Filas cabecera por glucoproteína
cab_R10 <- function(txt) tibble(
  Glucoproteína = txt,
  Variable      = NA_character_,
  p_global      = NA_character_,
  `p_D-/R+`     = NA_character_,
  `p_D+/R+`     = NA_character_,
  `p_D+/R-`     = NA_character_
)

tabla_R10_con_cab <- bind_rows(
  cab_R10("Glucoproteína gB"),
  tabla_R10_fmt %>% filter(Glucoproteína == "gB"),
  cab_R10("Glucoproteína gH"),
  tabla_R10_fmt %>% filter(Glucoproteína == "gH"),
  cab_R10("Glucoproteínas gB+gH"),
  tabla_R10_fmt %>% filter(Glucoproteína == "gB+gH")
)

# Índices de cabeceras
idx_cab_R10 <- which(is.na(tabla_R10_con_cab$Variable))

# Construir flextable
ft_R10 <- tabla_R10_con_cab %>%
  select(-Glucoproteína) %>%
  flextable() %>%
  set_header_labels(
    Variable  = "Variable",
    p_global  = "Global",
    `p_D-/R+` = "D-/R+",
    `p_D+/R+` = "D+/R+",
    `p_D+/R-` = "D+/R-"
  ) %>%
  add_header_row(
    values    = c("", "p (único vs mixto)"),
    colwidths = c(1, 4)
  ) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  bold(i = idx_cab_R10, j = 1, part = "body") %>%
  bg(i = idx_cab_R10, bg = "#F2F2F2", part = "body") %>%
  padding(i = setdiff(seq_len(nrow(tabla_R10_con_cab)), idx_cab_R10),
          j = 1, padding.left = 14, part = "body") %>%
  align(j = "Variable", align = "left", part = "all") %>%
  align(j = c("p_global", "p_D-/R+", "p_D+/R+", "p_D+/R-"),
        align = "center", part = "all") %>%
  autofit() %>%
  add_footer_lines(values = c(
    "Comparaciones mediante la prueba U de Mann-Whitney (único vs mixto).",
    "Los p-valores se muestran para el análisis global y estratificado por perfil serológico D/R.",
    "*p < 0,05."
  )) %>%
  fontsize(part = "footer", size = 8) %>%
  align(part = "footer", align = "left")

ft_R10

# Guardar en Word
read_docx() %>%
  body_add_flextable(ft_R10) %>%
  print(target = "~/Desktop/Figuras_tesis/Tabla_R10.docx")
# ════════════════════════════════════════════════════════════════════════════
#  FIN DEL SCRIPT
# ════════════════════════════════════════════════════════════════════════════

