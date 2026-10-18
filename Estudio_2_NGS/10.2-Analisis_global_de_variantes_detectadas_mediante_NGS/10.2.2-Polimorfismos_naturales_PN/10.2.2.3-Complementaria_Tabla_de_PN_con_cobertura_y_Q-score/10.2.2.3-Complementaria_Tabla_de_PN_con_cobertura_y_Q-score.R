# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.2.3. Complementaria. Tabla de PN con cobertura y Q-score
# Extracto del script completo (ver README.md del estudio)

## ---- Complementaria. Tabla de PN con cobertura y Q-score ----
# ── Datos ─────────────────────────────────────────────────────────────────────
pol_stats <- all_data_txp %>%
  filter(class == "POL", Prevalence >= 5) %>%
  group_by(Protein, Mutation) %>%
  summarise(
    n_muestras   = n_distinct(Npet),
    n_pacientes  = n_distinct(NHC),
    prev_media   = round(mean(Prevalence, na.rm = TRUE), 1),
    prev_de      = round(sd(Prevalence, na.rm = TRUE), 1),
    cob_mediana  = round(median(cobertura_total, na.rm = TRUE), 0),
    cob_p25      = round(quantile(cobertura_total, 0.25, na.rm = TRUE), 0),
    cob_p75      = round(quantile(cobertura_total, 0.75, na.rm = TRUE), 0),
    cob_min      = round(min(cobertura_total, na.rm = TRUE), 0),
    q_med        = round(median(q_score, na.rm = TRUE), 0),
    q_p25        = round(quantile(q_score, 0.25, na.rm = TRUE), 0),
    q_p75        = round(quantile(q_score, 0.75, na.rm = TRUE), 0),
    .groups = "drop"
  )

# ── Umbral máximo de detección ────────────────────────────────────────────────
umbral_max_pol <- map_dfr(c(5, 10, 15, 20), function(thr) {
  all_data_txp %>%
    filter(class == "POL", Prevalence >= thr) %>%
    distinct(Protein, Mutation) %>%
    mutate(umbral = thr)
}) %>%
  group_by(Protein, Mutation) %>%
  summarise(
    umbral_max = paste0("\u2265", max(umbral), "%"),
    .groups = "drop"
  )

# ── Total muestras y pacientes por gen ────────────────────────────────────────
n_total_gen <- tibble(
  Protein         = c("UL97", "UL54", "UL56"),
  n_tot_muestras  = c(26, 24, 11),
  n_tot_pacientes = c(11, 11, 3)
)

# ── Unir todo ─────────────────────────────────────────────────────────────────
pol_stats <- pol_stats %>%
  left_join(umbral_max_pol, by = c("Protein", "Mutation")) %>%
  left_join(n_total_gen,    by = "Protein") %>%
  mutate(
    pct_muestras  = round(n_muestras  / n_tot_muestras  * 100, 1),
    pct_pacientes = round(n_pacientes / n_tot_pacientes * 100, 1),
    col_prev = if_else(
      is.na(prev_de),
      paste0(format(prev_media, decimal.mark = ","), " (n=1)"),
      paste0(format(prev_media, decimal.mark = ","),
             " \u00b1 ",
             format(prev_de, decimal.mark = ","))
    ),
    col_muestras  = paste0(n_muestras, " (",
                           format(pct_muestras,  decimal.mark = ","), "%)"),
    col_pacientes = paste0(n_pacientes, " (",
                           format(pct_pacientes, decimal.mark = ","), "%)"),
    col_cobertura = paste0(
      format(cob_mediana, big.mark = "."), "\u00d7",
      " (", format(cob_p25, big.mark = "."),
      "\u2013", format(cob_p75, big.mark = "."), "\u00d7;",
      " m\u00edn. ", format(cob_min, big.mark = "."), "\u00d7)"
    ),
    col_qscore = paste0(
      q_med, " (", q_p25, "\u2013", q_p75, ")"
    )
  )

# ── Orden ─────────────────────────────────────────────────────────────────────
pol_stats <- pol_stats %>%
  mutate(Protein = factor(Protein, levels = c("UL97", "UL54", "UL56"))) %>%
  arrange(Protein, desc(n_muestras))

# ── Dataframe para flextable ──────────────────────────────────────────────────
pol_ft <- pol_stats %>%
  group_by(Protein) %>%
  mutate(
    Gen = if_else(row_number() == 1,
                  paste0(as.character(Protein),
                         "\n(n = ", n_tot_muestras, " muestras)"),
                  "")
  ) %>%
  ungroup() %>%
  dplyr::select(Gen, Mutation, col_muestras, col_pacientes,
                col_prev, umbral_max, col_cobertura, col_qscore)

names(pol_ft) <- c(
  "Gen", "PN",
  "Muestras, n (%)",
  "Pacientes, n (%)",
  "Prevalencia\n(media \u00b1 DE, %)",
  "Umbral m\u00e1x.\ndetecci\u00f3n",
  "Cobertura\n(mediana, RIC; m\u00edn.)",
  "Q-Score\n(mediana, RIC)"
)

# ── Filas por bloque ──────────────────────────────────────────────────────────
n_ul97 <- sum(pol_stats$Protein == "UL97")
n_ul54 <- sum(pol_stats$Protein == "UL54")
n_ul56 <- sum(pol_stats$Protein == "UL56")

filas_ul97 <- 1:n_ul97
filas_ul54 <- (n_ul97 + 1):(n_ul97 + n_ul54)
filas_ul56 <- (n_ul97 + n_ul54 + 1):(n_ul97 + n_ul54 + n_ul56)
filas_sep  <- c(n_ul97, n_ul97 + n_ul54)

# ── Bordes ────────────────────────────────────────────────────────────────────
b_azul <- officer::fp_border(color = "#2E5FA3", width = 1.5)
b_sep  <- officer::fp_border(color = "#2E5FA3", width = 1.2)
b_gris <- officer::fp_border(color = "#CCCCCC", width = 0.5)
b_head <- officer::fp_border(color = "#FFFFFF", width = 0.8)

# ── Flextable ─────────────────────────────────────────────────────────────────
ft_pol <- flextable::flextable(pol_ft) %>%
  flextable::merge_at(i = filas_ul97, j = 1, part = "body") %>%
  flextable::merge_at(i = filas_ul54, j = 1, part = "body") %>%
  flextable::merge_at(i = filas_ul56, j = 1, part = "body") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = 1:2, align = "left", part = "body") %>%
  flextable::valign(valign = "center", part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1:2, part = "body") %>%
  flextable::bg(part = "header", bg = "#2E5FA3") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::border_inner_v(border = b_head, part = "header") %>%
  flextable::bg(i = filas_ul97, bg = "#EEF4FB", part = "body") %>%
  flextable::bg(i = filas_ul54, bg = "#FFFFFF", part = "body") %>%
  flextable::bg(i = filas_ul56, bg = "#EEF4FB", part = "body") %>%
  flextable::border_outer(border = b_azul, part = "all") %>%
  flextable::border_inner_h(border = b_gris, part = "body") %>%
  flextable::border_inner_v(border = b_gris, part = "body") %>%
  flextable::hline(i = filas_sep, border = b_sep, part = "body") %>%
  flextable::width(j = 1, width = 3.0, unit = "cm") %>%
  flextable::width(j = 2, width = 2.0, unit = "cm") %>%
  flextable::width(j = 3, width = 2.5, unit = "cm") %>%
  flextable::width(j = 4, width = 2.5, unit = "cm") %>%
  flextable::width(j = 5, width = 3.0, unit = "cm") %>%
  flextable::width(j = 6, width = 2.0, unit = "cm") %>%
  flextable::width(j = 7, width = 3.5, unit = "cm") %>%
  flextable::width(j = 8, width = 2.5, unit = "cm") %>%
  flextable::padding(padding = 6, part = "all") %>%
  flextable::add_footer_lines(paste0(
    "PN: polimorfismo natural. DE: desviaci\u00f3n est\u00e1ndar. RIC: rango intercuart\u00edlico. ",
    "Los porcentajes de muestras y pacientes se calculan sobre el total con datos v\u00e1lidos por gen ",
    "(UL97: 26 muestras, 11 pacientes; UL54: 24 muestras, 11 pacientes; UL56: 11 muestras, 3 pacientes). ",
    "El umbral m\u00e1ximo de detecci\u00f3n indica el umbral de frecuencia de variante m\u00e1s alto al que la variante es detectable. ",
    "Cobertura: suma de lecturas forward + reverse en la posici\u00f3n (mediana y RIC entre todas las muestras con la variante; m\u00ednimo observado). ",
    "Q-Score: mediana y RIC de la puntuaci\u00f3n de calidad de llamada de variante."
  )) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft_pol

guardar_tabla(ft_pol, "Complementaria_PN", apaisada = TRUE)

