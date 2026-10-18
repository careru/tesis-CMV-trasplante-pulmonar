# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.1. Variantes asociadas a resistencia antiviral (MRA)
# Extracto del script completo (ver README.md del estudio)

# ---- 10.2.1. Variantes asociadas a resistencia antiviral (MRA) ----

## ---- Tabla 24. Características de las MRA detectadas en la paciente RTP-2 ----
nhc_rtp2 <- "NHC_RTP-2"

# ── Datos fijos (CHARMD) ──────────────────────────────────────────────────────
charmd <- tibble(
  Protein    = c("UL97", "UL97", "UL97", "UL54"),
  Mutation   = c("A594V", "T409M", "H411Y", "T503I"),
  Resistencia_fenotipica = c(
    "Ganciclovir",
    "Maribavir",
    "Maribavir",
    "Ganciclovir y Cidofovir"
  )
)

# ── Calcular metricas desde all_data_txp ──────────────────────────────────────
mra_raw <- all_data_txp %>%
  filter(class == "MRA", Prevalence >= 5)

mra_stats <- mra_raw %>%
  group_by(Protein, Mutation) %>%
  summarise(
    n_muestras   = n_distinct(Npet),
    prev_media   = round(mean(Prevalence, na.rm = TRUE), 1),
    prev_de      = round(sd(Prevalence,   na.rm = TRUE), 1),
    cv_mediana   = round(median(as.numeric(CV), na.rm = TRUE), 0),
    cv_p25       = round(quantile(as.numeric(CV), 0.25, na.rm = TRUE), 0),
    cv_p75       = round(quantile(as.numeric(CV), 0.75, na.rm = TRUE), 0),
    cob_mediana  = round(median(cobertura_total, na.rm = TRUE), 0),
    cob_p25      = round(quantile(cobertura_total, 0.25, na.rm = TRUE), 0),
    cob_p75      = round(quantile(cobertura_total, 0.75, na.rm = TRUE), 0),
    cob_min      = round(min(cobertura_total, na.rm = TRUE), 0),
    q_med        = round(median(q_score, na.rm = TRUE), 0),
    q_p25        = round(quantile(q_score, 0.25, na.rm = TRUE), 0),
    q_p75        = round(quantile(q_score, 0.75, na.rm = TRUE), 0),
    n_balance_ok = sum(balance_ok, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(charmd, by = c("Protein", "Mutation"))

# ── Umbral minimo de deteccion ────────────────────────────────────────────────
umbral_min_mra <- map_dfr(c(5, 10, 15, 20), function(thr) {
  all_data_txp %>%
    filter(class == "MRA", Prevalence >= thr) %>%
    distinct(Protein, Mutation) %>%
    mutate(umbral = thr)
}) %>%
  group_by(Protein, Mutation) %>%
  summarise(umbral_min = paste0("\u2265", min(umbral), "%"), .groups = "drop")

mra_stats <- mra_stats %>%
  left_join(umbral_min_mra, by = c("Protein", "Mutation"))

# ── Total muestras validas por gen en RTP-2 (calculado) ───────────────────────
n_muestras_gen <- all_data_nhc %>%
  filter(NHC == nhc_rtp2, Protein %in% c("UL97", "UL54")) %>%
  group_by(Protein) %>%
  summarise(n_total = n_distinct(Npet), .groups = "drop")

cat("── Denominadores de la Tabla 25 ──\n")
print(n_muestras_gen)

n_ul97 <- n_muestras_gen$n_total[n_muestras_gen$Protein == "UL97"]
n_ul54 <- n_muestras_gen$n_total[n_muestras_gen$Protein == "UL54"]

# ── Helper de formato numerico sin espacios sobrantes ─────────────────────────
fnum <- function(x) format(x, big.mark = ".", trim = TRUE)

mra_stats <- mra_stats %>%
  left_join(n_muestras_gen, by = "Protein") %>%
  mutate(
    pct   = round(n_muestras / n_total * 100, 1),
    col_n = paste0(n_muestras, " (",
                   format(pct, decimal.mark = ",", trim = TRUE), "%)"),
    col_prev = if_else(
      is.na(prev_de),
      paste0(format(prev_media, decimal.mark = ",", trim = TRUE), " (n = 1)"),
      paste0(format(prev_media, decimal.mark = ",", trim = TRUE),
             " \u00b1 ", format(prev_de, decimal.mark = ",", trim = TRUE))
    ),
    col_cv = if_else(
      is.na(cv_mediana),
      "\u2014",
      paste0(fnum(cv_mediana), " (", fnum(cv_p25), "\u2013", fnum(cv_p75), ")")
    ),
    col_cob = if_else(
      n_muestras == 1,
      paste0(fnum(cob_mediana), "\u00d7"),
      paste0(fnum(cob_mediana), "\u00d7 (", fnum(cob_p25), "\u2013",
             fnum(cob_p75), "; m\u00edn. ", fnum(cob_min), "\u00d7)")
    ),
    col_q = if_else(
      n_muestras == 1,
      as.character(q_med),
      paste0(q_med, " (", q_p25, "\u2013", q_p75, ")")
    ),
    col_balance = paste0(n_balance_ok, "/", n_muestras)
  )

# ── Orden explicito ───────────────────────────────────────────────────────────
orden_mra <- tibble(
  Protein  = c("UL97", "UL97", "UL97", "UL54"),
  Mutation = c("A594V", "T409M", "H411Y", "T503I")
)

filas_mra <- orden_mra %>%
  left_join(mra_stats, by = c("Protein", "Mutation")) %>%
  mutate(
    Gen = case_when(
      Mutation == "A594V" ~ paste0("UL97\n(n = ", n_ul97, " muestras)"),
      Mutation == "T503I" ~ paste0("UL54\n(n = ", n_ul54, " muestras)"),
      TRUE ~ ""
    )
  ) %>%
  dplyr::select(
    Gen, Mutation, col_n, col_prev, col_cv,
    col_cob, col_q, col_balance, umbral_min, Resistencia_fenotipica
  )

names(filas_mra) <- c(
  "Gen", "MRA", "Muestras, n (%)",
  "Prevalencia (media \u00b1 DE, %)",
  "CV en muestras con MRA\nmediana (RIC), UI/mL",
  "Cobertura\n(mediana, RIC; m\u00edn.)",
  "Q-score\n(mediana, RIC)",
  "Balance FW/RV\n(n muestras)",
  "Umbral m\u00ednimo de\ndetecci\u00f3n",
  "Resistencia fenot\u00edpica*"
)

# ── Bordes ────────────────────────────────────────────────────────────────────
b_azul <- officer::fp_border(color = "#2E5FA3", width = 1.5)
b_sep  <- officer::fp_border(color = "#2E5FA3", width = 1.2)
b_gris <- officer::fp_border(color = "#CCCCCC", width = 0.5)
b_head <- officer::fp_border(color = "#FFFFFF", width = 0.8)

# ── Flextable ─────────────────────────────────────────────────────────────────
ft_mra <- flextable::flextable(filas_mra) %>%
  flextable::merge_at(i = 1:3, j = 1, part = "body") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = 1:2, align = "left", part = "body") %>%
  flextable::valign(valign = "center", part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1:2, part = "body") %>%
  flextable::bg(part = "header", bg = "#2E5FA3") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::border_inner_v(border = b_head, part = "header") %>%
  flextable::bg(i = 1:3, bg = "#EEF4FB", part = "body") %>%
  flextable::bg(i = 4,   bg = "#FFFFFF", part = "body") %>%
  flextable::border_outer(border = b_azul, part = "all") %>%
  flextable::border_inner_h(border = b_gris, part = "body") %>%
  flextable::border_inner_v(border = b_gris, part = "body") %>%
  flextable::hline(i = 3, border = b_sep, part = "body") %>%
  flextable::width(j = 1,  width = 2.4, unit = "cm") %>%
  flextable::width(j = 2,  width = 1.8, unit = "cm") %>%
  flextable::width(j = 3,  width = 2.0, unit = "cm") %>%
  flextable::width(j = 4,  width = 2.5, unit = "cm") %>%
  flextable::width(j = 5,  width = 2.8, unit = "cm") %>%
  flextable::width(j = 6,  width = 3.0, unit = "cm") %>%
  flextable::width(j = 7,  width = 1.8, unit = "cm") %>%
  flextable::width(j = 8,  width = 2.1, unit = "cm") %>%
  flextable::width(j = 9,  width = 2.0, unit = "cm") %>%
  flextable::width(j = 10, width = 3.0, unit = "cm") %>%
  flextable::padding(padding = 5, part = "all") %>%
  flextable::add_footer_lines(
    "*Resistencia fenot\u00edpica seg\u00fan la base de datos CHARMD (Tilloy et al., 2024). CV: carga viral; DE: desviaci\u00f3n est\u00e1ndar; RIC: rango intercuart\u00edlico. El umbral m\u00ednimo de detecci\u00f3n indica el umbral de frecuencia de variante m\u00e1s bajo al que la variante es detectable. Cobertura: n\u00famero total de lecturas forward + reverse que cubren la posici\u00f3n, expresado como mediana, rango intercuart\u00edlico y valor m\u00ednimo. Q-score: puntuaci\u00f3n de calidad de llamada. Balance FW/RV: n\u00famero de muestras con soporte en ambas cadenas sobre el total de muestras con la variante."
  ) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft_mra

guardar_tabla(ft_mra, "Tabla_24", apaisada = TRUE)

