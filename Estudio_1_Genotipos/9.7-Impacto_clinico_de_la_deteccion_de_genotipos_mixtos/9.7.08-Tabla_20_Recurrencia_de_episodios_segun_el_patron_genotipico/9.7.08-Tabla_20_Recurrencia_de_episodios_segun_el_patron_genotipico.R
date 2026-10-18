# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.7.08. Tabla 20. Recurrencia de episodios según el patrón genotípico
# Extracto del script completo (ver README.md del estudio)

## ---- Tabla 20. Recurrencia de episodios según el patrón genotípico ----


# ── 1. Identificar primer episodio por paciente ─────────────────────────────
# Definición: el episodio con Fini más cercano a Ftx (mínimo Fini por NHC)

primer_episodio <- dnaemias %>%
  mutate(
    Fini = as.Date(Fini),
    Ffin = as.Date(Ffin)
  ) %>%
  group_by(NHC) %>%
  arrange(Fini, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    Ftx              = as.Date(Ftx),
    dias_tx_a_inicio = as.numeric(Fini - Ftx)
  )

# ── 2. Excluir pacientes sin Ffin documentado ────────────────────────────────
# (4 fallecidos antes de aclaramiento + 3 trasladados = 7 excluidos)
primer_ep_evaluables <- primer_episodio %>%
  filter(!is.na(Ffin))

cat(sprintf("── Pacientes evaluables: %d ──\n", nrow(primer_ep_evaluables)))

# ── 3. Identificar segundo episodio por paciente ────────────────────────────
segundo_episodio <- dnaemias %>%
  mutate(
    Fini = as.Date(Fini),
    Ffin = as.Date(Ffin)
  ) %>%
  inner_join(
    primer_ep_evaluables %>% select(NHC, Ffin_primer = Ffin),
    by = "NHC"
  ) %>%
  filter(Fini > Ffin_primer) %>%
  group_by(NHC) %>%
  arrange(Fini, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(NHC,
         Fini_segundo    = Fini,
         Ffin_segundo    = Ffin,
         gbghmix_segundo = gbghmix,
         gbmix_segundo   = gbmix,
         ghmix_segundo   = ghmix,
         cvini_segundo   = cvini,
         clinica_segundo = clinica)

# ── 4. Construir base df_recur ──────────────────────────────────────────────
df_recur <- primer_ep_evaluables %>%
  select(NHC,
         Ftx,
         Fini_primer    = Fini,
         Ffin_primer    = Ffin,
         gbghmix_primer = gbghmix,
         gbmix_primer   = gbmix,
         ghmix_primer   = ghmix,
         cvini_primer   = cvini,
         cvpico_primer  = cvpico,
         dias_tx_a_inicio) %>%
  left_join(segundo_episodio, by = "NHC") %>%
  left_join(
    pacientes %>%
      mutate(
        Ffinseg = as.Date(Ffinseg),
        Fexitus = as.Date(Fexitus)
      ) %>%
      select(NHC, Ffinseg, exitus, Fexitus, serDR),
    by = "NHC"
  )

# ── 5. Calcular tiempos y eventos para análisis ──────────────────────────────
df_recur <- df_recur %>%
  mutate(
    fecha_fin_seg = pmin(Ffinseg, Fexitus, na.rm = TRUE),
    dias_a_recur  = as.numeric(Fini_segundo - Ffin_primer),
    dias_seguim   = as.numeric(fecha_fin_seg - Ffin_primer),
    
    # Recurrencia a 6 meses
    recur_180 = case_when(
      !is.na(dias_a_recur) & dias_a_recur <= 180 ~ 1L,
      TRUE                                       ~ 0L
    ),
    tiempo_km = case_when(
      recur_180 == 1     ~ dias_a_recur,
      dias_seguim < 180  ~ dias_seguim,
      TRUE               ~ 180
    ),
    evento_km = recur_180,
    
    # Recurrencia a 12 meses
    recur_365 = case_when(
      !is.na(dias_a_recur) & dias_a_recur <= 365 ~ 1L,
      TRUE                                       ~ 0L
    ),
    tiempo_km_365 = case_when(
      recur_365 == 1     ~ dias_a_recur,
      dias_seguim < 365  ~ dias_seguim,
      TRUE               ~ 365
    ),
    evento_km_365 = recur_365,
    
    # Recurrencia en seguimiento completo
    recur_global = as.integer(!is.na(Fini_segundo)),
    tiempo_completo = case_when(
      recur_global == 1 ~ dias_a_recur,
      TRUE              ~ dias_seguim
    ),
    
    # Factores de exposición
    exp_grupo = factor(gbghmix_primer, levels = c(0, 1),
                       labels = c("Único", "Mixto")),
    gb_grupo  = factor(gbmix_primer,   levels = c(0, 1),
                       labels = c("gB único", "gB mixto")),
    gh_grupo  = factor(ghmix_primer,   levels = c(0, 1),
                       labels = c("gH único", "gH mixto")),
    serDR_f   = factor(serDR)
  )

# ── 5b. Excepción: aclaramiento posterior al cierre del seguimiento ──────────
# El paciente "NHC_anonimizado" se aclaró (17/01/2023) después del cierre administrativo
# del seguimiento (31/12/2022). Se incluye con tiempo 0 (censurado en el
# aclaramiento) en lugar de un tiempo negativo.
df_recur <- df_recur %>%
  mutate(across(c(tiempo_km, tiempo_km_365, tiempo_completo), ~ pmax(.x, 0)))

cat("Tiempos negativos tras la corrección:", sum(df_recur$tiempo_km_365 < 0), "\n")

# Objetos Surv
surv_180      <- Surv(time = df_recur$tiempo_km,     event = df_recur$evento_km)
surv_365      <- Surv(time = df_recur$tiempo_km_365, event = df_recur$evento_km_365)
surv_completo <- Surv(time = df_recur$tiempo_completo, event = df_recur$recur_global)

# ── 6 meses ─────────────────────────────────────────────────────────────────
cox_gbgh_6m     <- coxph(surv_180 ~ exp_grupo, data = df_recur)
cox_gb_6m       <- coxph(surv_180 ~ gb_grupo,  data = df_recur)
cox_gh_6m       <- coxph(surv_180 ~ gh_grupo,  data = df_recur)
cox_gbgh_6m_adj <- coxph(surv_180 ~ exp_grupo + serDR_f + cvpico_primer,
                         data = df_recur)
cox_gb_6m_adj   <- coxph(surv_180 ~ gb_grupo  + serDR_f + cvpico_primer,
                         data = df_recur)
cox_gh_6m_adj   <- coxph(surv_180 ~ gh_grupo  + serDR_f + cvpico_primer,
                         data = df_recur)

# ── 12 meses ────────────────────────────────────────────────────────────────
cox_gbgh_12m     <- coxph(surv_365 ~ exp_grupo, data = df_recur)
cox_gb_12m       <- coxph(surv_365 ~ gb_grupo,  data = df_recur)
cox_gh_12m       <- coxph(surv_365 ~ gh_grupo,  data = df_recur)
cox_gbgh_12m_adj <- coxph(surv_365 ~ exp_grupo + serDR_f + cvpico_primer,
                          data = df_recur)
cox_gb_12m_adj   <- coxph(surv_365 ~ gb_grupo  + serDR_f + cvpico_primer,
                          data = df_recur)
cox_gh_12m_adj   <- coxph(surv_365 ~ gh_grupo  + serDR_f + cvpico_primer,
                          data = df_recur)

# ── Seguimiento completo ────────────────────────────────────────────────────
cox_gbgh_sc     <- coxph(surv_completo ~ exp_grupo, data = df_recur)
cox_gb_sc       <- coxph(surv_completo ~ gb_grupo,  data = df_recur)
cox_gh_sc       <- coxph(surv_completo ~ gh_grupo,  data = df_recur)
cox_gbgh_sc_adj <- coxph(surv_completo ~ exp_grupo + serDR_f + cvpico_primer,
                         data = df_recur)
cox_gb_sc_adj   <- coxph(surv_completo ~ gb_grupo  + serDR_f + cvpico_primer,
                         data = df_recur)
cox_gh_sc_adj   <- coxph(surv_completo ~ gh_grupo  + serDR_f + cvpico_primer,
                         data = df_recur)

# ── Funciones auxiliares ────────────────────────────────────────────────────
extraer_hr_fmt <- function(cox) {
  hr <- exp(coef(cox)[1])
  ci <- exp(confint(cox)[1, ])
  p  <- summary(cox)$coefficients[1, 5]
  list(
    hr_ci = sprintf("%.2f (%.2f\u2013%.2f)", hr, ci[1], ci[2]),
    p     = ifelse(p < 0.001, "<0,001",
                   sub("\\.", ",", sprintf("%.3f", p)))
  )
}

fmt_n_pct <- function(n, N) {
  sub("\\.", ",", sprintf("%d (%.1f%%)", n, 100 * n / N))
}

# ── Construir filas ─────────────────────────────────────────────────────────
filas_R13 <- bind_rows(
  
  # Sección 1: 6 meses
  tibble(variable = "RECURRENCIA A 6 MESES (n = 33)",
         unicos = "", mixtos = "", hr_crudo = "", p_crudo = "",
         hr_adj = "", p_adj = ""),
  
  tibble(variable = "gB+gH",
         unicos   = fmt_n_pct(21, 99),
         mixtos   = fmt_n_pct(12, 40),
         hr_crudo = extraer_hr_fmt(cox_gbgh_6m)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gbgh_6m)$p,
         hr_adj   = extraer_hr_fmt(cox_gbgh_6m_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gbgh_6m_adj)$p),
  
  tibble(variable = "gB",
         unicos   = fmt_n_pct(25, 108),
         mixtos   = fmt_n_pct(8, 31),
         hr_crudo = extraer_hr_fmt(cox_gb_6m)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gb_6m)$p,
         hr_adj   = extraer_hr_fmt(cox_gb_6m_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gb_6m_adj)$p),
  
  tibble(variable = "gH",
         unicos   = fmt_n_pct(24, 112),
         mixtos   = fmt_n_pct(9, 27),
         hr_crudo = extraer_hr_fmt(cox_gh_6m)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gh_6m)$p,
         hr_adj   = extraer_hr_fmt(cox_gh_6m_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gh_6m_adj)$p),
  
  # Sección 2: 12 meses
  tibble(variable = "RECURRENCIA A 12 MESES (n = 42)",
         unicos = "", mixtos = "", hr_crudo = "", p_crudo = "",
         hr_adj = "", p_adj = ""),
  
  tibble(variable = "gB+gH",
         unicos   = fmt_n_pct(25, 99),
         mixtos   = fmt_n_pct(17, 40),
         hr_crudo = extraer_hr_fmt(cox_gbgh_12m)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gbgh_12m)$p,
         hr_adj   = extraer_hr_fmt(cox_gbgh_12m_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gbgh_12m_adj)$p),
  
  tibble(variable = "gB",
         unicos   = fmt_n_pct(30, 108),
         mixtos   = fmt_n_pct(12, 31),
         hr_crudo = extraer_hr_fmt(cox_gb_12m)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gb_12m)$p,
         hr_adj   = extraer_hr_fmt(cox_gb_12m_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gb_12m_adj)$p),
  
  tibble(variable = "gH",
         unicos   = fmt_n_pct(29, 112),
         mixtos   = fmt_n_pct(13, 27),
         hr_crudo = extraer_hr_fmt(cox_gh_12m)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gh_12m)$p,
         hr_adj   = extraer_hr_fmt(cox_gh_12m_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gh_12m_adj)$p),
  
  # Sección 3: Seguimiento completo
  tibble(variable = "RECURRENCIA EN SEGUIMIENTO COMPLETO (n = 56)",
         unicos = "", mixtos = "", hr_crudo = "", p_crudo = "",
         hr_adj = "", p_adj = ""),
  
  tibble(variable = "gB+gH",
         unicos   = fmt_n_pct(35, 99),
         mixtos   = fmt_n_pct(21, 40),
         hr_crudo = extraer_hr_fmt(cox_gbgh_sc)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gbgh_sc)$p,
         hr_adj   = extraer_hr_fmt(cox_gbgh_sc_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gbgh_sc_adj)$p),
  
  tibble(variable = "gB",
         unicos   = fmt_n_pct(42, 108),
         mixtos   = fmt_n_pct(14, 31),
         hr_crudo = extraer_hr_fmt(cox_gb_sc)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gb_sc)$p,
         hr_adj   = extraer_hr_fmt(cox_gb_sc_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gb_sc_adj)$p),
  
  tibble(variable = "gH",
         unicos   = fmt_n_pct(40, 112),
         mixtos   = fmt_n_pct(16, 27),
         hr_crudo = extraer_hr_fmt(cox_gh_sc)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gh_sc)$p,
         hr_adj   = extraer_hr_fmt(cox_gh_sc_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gh_sc_adj)$p)
)

colnames(filas_R13) <- c("Genotipos",
                         "Únicos\nn (%)",
                         "Mixtos\nn (%)",
                         "HR cruda\n(IC 95%)",
                         "p",
                         "HR ajustada*\n(IC 95%)",
                         "p ")

# ── Filas especiales por índice ─────────────────────────────────────────────
filas_seccion_R13 <- c(1, 5, 9)

es_p_sig <- function(x) {
  num <- suppressWarnings(as.numeric(sub(",", ".", x)))
  !is.na(num) & num < 0.05
}

filas_p_sig_crudo <- which(es_p_sig(filas_R13[[5]]))
filas_p_sig_adj   <- which(es_p_sig(filas_R13[[7]]))

# ── Construir flextable ─────────────────────────────────────────────────────
ft_R13 <- flextable(filas_R13) %>%
  width(j = 1, width = 2.0) %>%
  width(j = 2, width = 1.4) %>%
  width(j = 3, width = 1.4) %>%
  width(j = 4, width = 1.8) %>%
  width(j = 5, width = 0.9) %>%
  width(j = 6, width = 1.8) %>%
  width(j = 7, width = 0.9) %>%
  align(align = "left", part = "all") %>%
  align(j = c(2, 3, 4, 5, 6, 7), align = "center", part = "all") %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header") %>%
  bg(bg = "white", part = "all") %>%
  bold(i = filas_seccion_R13, part = "body") %>%
  bg(i = filas_seccion_R13, bg = "#F0F0F0", part = "body") %>%
  border_remove() %>%
  hline_top(part = "header",
            border = fp_border(color = "black", width = 1.2)) %>%
  hline(part = "header",
        border = fp_border(color = "black", width = 0.8)) %>%
  hline_bottom(part = "body",
               border = fp_border(color = "black", width = 1.2)) %>%
  hline(i = c(4, 8),
        border = fp_border(color = "grey70", width = 0.4),
        part = "body") %>%
  padding(padding.top = 2, padding.bottom = 2,
          padding.left = 4, padding.right = 4, part = "all") %>%
  add_footer_lines(values = c(
    "*Modelo de regresión de Cox ajustado por perfil serológico D/R y carga viral pico del primer episodio.",
    "Las HR comparan el grupo mixto frente al grupo único, considerado como referencia.",
    "Valores de p en rojo indican diferencia estadísticamente significativa (p < 0,05)."
  )) %>%
  fontsize(size = 8, part = "footer") %>%
  color(color = "grey30", part = "footer") %>%
  set_caption(
    caption = "Tabla R13. Recurrencia de DNAemia por CMV según el patrón genotípico del primer episodio en tres ventanas temporales."
  )

if (length(filas_p_sig_crudo) > 0) {
  ft_R13 <- ft_R13 %>%
    color(i = filas_p_sig_crudo, j = 5, color = "firebrick3", part = "body") %>%
    bold( i = filas_p_sig_crudo, j = 5, part = "body")
}
if (length(filas_p_sig_adj) > 0) {
  ft_R13 <- ft_R13 %>%
    color(i = filas_p_sig_adj, j = 7, color = "firebrick3", part = "body") %>%
    bold( i = filas_p_sig_adj, j = 7, part = "body")
}

ft_R13 <- ft_R13 %>%
  set_table_properties(layout = "autofit")

ft_R13

guardar_tabla(ft_R13, "Tabla_20")

