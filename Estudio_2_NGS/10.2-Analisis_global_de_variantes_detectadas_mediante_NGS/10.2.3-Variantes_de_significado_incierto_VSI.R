# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.3. Variantes de significado incierto (VSI)
# Extracto del script completo (ver README.md del estudio)

# ---- 10.2.3. Variantes de significado incierto (VSI) ----

## ---- Datos y tabla de VSI (vacía: D605E se reclasificó como PN) ----
all_data_txp %>%
  filter(class == "VSI", Prevalence >= 5) %>%
  group_by(Protein, Mutation) %>%
  summarise(
    n_muestras  = n_distinct(Npet),
    n_pacientes = n_distinct(NHC),
    NHC_list    = paste(sort(unique(NHC)), collapse = ", "),
    prev_media  = round(mean(Prevalence, na.rm = TRUE), 1),
    prev_de     = round(sd(Prevalence, na.rm = TRUE), 1),
    prev_min    = round(min(Prevalence, na.rm = TRUE), 1),
    prev_max    = round(max(Prevalence, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(Protein, desc(n_muestras)) %>%
  print(n = Inf)


all_data_txp %>%
  filter(class == "VSI", Prevalence >= 5) %>%
  select(Npet, Freg, Protein, Mutation, Prevalence, CV) %>%
  arrange(Protein)

# ── Datos VSI ─────────────────────────────────────────────────────────────────
vsi_stats <- all_data_txp %>%
  filter(class == "VSI", Prevalence >= 5) %>%
  group_by(Protein, Mutation) %>%
  summarise(
    n_muestras  = n_distinct(Npet),
    n_pacientes = n_distinct(NHC),
    cob_mediana  = round(median(cobertura_total, na.rm = TRUE), 0),
    cob_p25      = round(quantile(cobertura_total, 0.25, na.rm = TRUE), 0),
    cob_p75      = round(quantile(cobertura_total, 0.75, na.rm = TRUE), 0),
    cob_min      = round(min(cobertura_total, na.rm = TRUE), 0),
    q_med        = round(median(q_score, na.rm = TRUE), 0),
    q_p25        = round(quantile(q_score, 0.25, na.rm = TRUE), 0),
    q_p75        = round(quantile(q_score, 0.75, na.rm = TRUE), 0),
    n_balance_ok = sum(balance_ok, na.rm = TRUE),
    prev_media  = round(mean(Prevalence, na.rm = TRUE), 1),
    prev_de     = round(sd(Prevalence, na.rm = TRUE), 1),
    cv_mediana  = round(median(as.numeric(CV), na.rm = TRUE), 0),
    cv_p25      = round(quantile(as.numeric(CV), 0.25, na.rm = TRUE), 0),
    cv_p75      = round(quantile(as.numeric(CV), 0.75, na.rm = TRUE), 0),
    .groups = "drop"
  )

# ── Umbral máximo de detección ────────────────────────────────────────────────
umbral_max_vsi <- map_dfr(c(5, 10, 15, 20), function(thr) {
  all_data_txp %>%
    filter(class == "VSI", Prevalence >= thr) %>%
    distinct(Protein, Mutation) %>%
    mutate(umbral = thr)
}) %>%
  group_by(Protein, Mutation) %>%
  summarise(
    umbral_max = paste0("\u2265", max(umbral), "%"),
    .groups = "drop"
  )

# ── Total muestras por gen (denominador RTP-2) ────────────────────────────────
n_muestras_gen <- tibble(
  Protein = c("UL97", "UL54"),
  n_total = c(14, 12)
)

# ── Unir todo ─────────────────────────────────────────────────────────────────
vsi_stats <- vsi_stats %>%
  left_join(umbral_max_vsi, by = c("Protein", "Mutation")) %>%
  left_join(n_muestras_gen, by = "Protein") %>%
  mutate(
    pct      = round(n_muestras / n_total * 100, 1),
    col_n    = paste0(n_muestras, " (",
                      format(pct, decimal.mark = ","), "%)"),
    col_prev = if_else(
      is.na(prev_de),
      paste0(format(prev_media, decimal.mark = ","), " (n=1)"),
      paste0(format(prev_media, decimal.mark = ","),
             " \u00b1 ",
             format(prev_de, decimal.mark = ","))
    ),
    col_cv = if_else(
      is.na(cv_mediana),
      "\u2014",
      paste0(format(cv_mediana, big.mark = "."),
             " (", format(cv_p25, big.mark = "."),
             "\u2013", format(cv_p75, big.mark = "."), ")")
    )
  )

# ── Orden: UL97 primero, luego UL54 ──────────────────────────────────────────
vsi_stats <- vsi_stats %>%
  mutate(Protein = factor(Protein, levels = c("UL97", "UL54"))) %>%
  arrange(Protein)

# ── Dataframe para flextable ──────────────────────────────────────────────────
vsi_ft <- vsi_stats %>%
  group_by(Protein) %>%
  mutate(
    col_cob = if_else(
      n_muestras == 1,
      paste0(format(cob_mediana, big.mark = "."), "×"),
      paste0(
        format(cob_mediana, big.mark = "."), "× (",
        format(cob_p25, big.mark = "."), "–",
        format(cob_p75, big.mark = "."), "; mín. ",
        format(cob_min, big.mark = "."), "×)"
      )
    ),
    col_q = if_else(
      n_muestras == 1,
      as.character(q_med),
      paste0(q_med, " (", q_p25, "–", q_p75, ")")
    ),
    col_balance = paste0(n_balance_ok, "/", n_muestras),
    Gen = if_else(row_number() == 1,
                  paste0(as.character(Protein),
                         "\n(n = ", n_total, " muestras)"),
                  "")
  ) %>%
  ungroup() %>%
  dplyr::select(
    Gen, Mutation, col_n, col_prev, col_cv,
    col_cob, col_q, col_balance, umbral_max
  )
names(vsi_ft) <- c(
  "Gen", "VSI",
  "Muestras, n (%)",
  "Prevalencia (media ± DE, %)",
  "CV en muestra con VSI\n(UI/mL)",
  "Cobertura\n(mediana, RIC; mín.)",
  "Q-score\n(mediana, RIC)",
  "Balance FW/RV\n(n muestras)",
  "Umbral máximo de\ndetección"
)

# ── Bordes ────────────────────────────────────────────────────────────────────
b_azul <- officer::fp_border(color = "#2E5FA3", width = 1.5)
b_sep  <- officer::fp_border(color = "#2E5FA3", width = 1.2)
b_gris <- officer::fp_border(color = "#CCCCCC", width = 0.5)
b_head <- officer::fp_border(color = "#FFFFFF", width = 0.8)

n_ul97 <- sum(vsi_stats$Protein == "UL97")
n_ul54 <- sum(vsi_stats$Protein == "UL54")
n_total_filas <- nrow(vsi_ft)

filas_ul97 <- if (n_ul97 > 0) 1:n_ul97 else integer(0)
filas_ul54 <- if (n_ul54 > 0) (n_ul97 + 1):n_total_filas else integer(0)

ft_vsi <- flextable::flextable(vsi_ft) %>%
  { if (length(filas_ul97) > 1) flextable::merge_at(., i = filas_ul97, j = 1, part = "body") else . } %>%
  { if (length(filas_ul54) > 1) flextable::merge_at(., i = filas_ul54, j = 1, part = "body") else . } %>%
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
  { if (length(filas_ul97) > 0) flextable::bg(., i = filas_ul97, bg = "#EEF4FB", part = "body") else . } %>%
  { if (length(filas_ul54) > 0) flextable::bg(., i = filas_ul54, bg = "#FFFFFF",  part = "body") else . } %>%
  flextable::border_outer(border = b_azul, part = "all") %>%
  flextable::border_inner_h(border = b_gris, part = "body") %>%
  flextable::border_inner_v(border = b_gris, part = "body") %>%
  { if (n_ul97 > 0 & n_ul54 > 0) flextable::hline(., i = n_ul97, border = b_sep, part = "body") else . } %>%
  flextable::width(j = 1, width = 2.6, unit = "cm") %>%
  flextable::width(j = 2, width = 1.8, unit = "cm") %>%
  flextable::width(j = 3, width = 2.2, unit = "cm") %>%
  flextable::width(j = 4, width = 2.7, unit = "cm") %>%
  flextable::width(j = 5, width = 2.6, unit = "cm") %>%
  flextable::width(j = 6, width = 3.0, unit = "cm") %>%
  flextable::width(j = 7, width = 2.0, unit = "cm") %>%
  flextable::width(j = 8, width = 2.2, unit = "cm") %>%
  flextable::width(j = 9, width = 2.1, unit = "cm") %>%
  flextable::padding(padding = 6, part = "all") %>%
  flextable::add_footer_lines(
    "VSI: variante de significado incierto. DE: desviaci\u00f3n est\u00e1ndar. Los porcentajes se calculan sobre el total de muestras v\u00e1lidas del paciente RTP-2 por gen (UL97: 14 muestras; UL54: 12 muestras). El umbral m\u00e1ximo de detecci\u00f3n indica el umbral de prevalencia al\u00e9lica m\u00e1s alto al que la variante es detectable."
  ) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft_vsi


