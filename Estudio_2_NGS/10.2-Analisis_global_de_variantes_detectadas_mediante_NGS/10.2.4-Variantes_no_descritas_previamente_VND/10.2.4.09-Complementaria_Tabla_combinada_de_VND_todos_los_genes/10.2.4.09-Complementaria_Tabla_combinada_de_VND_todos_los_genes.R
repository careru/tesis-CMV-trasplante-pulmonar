# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.4.09. Complementaria. Tabla combinada de VND (todos los genes)
# Extracto del script completo (ver README.md del estudio)

## ---- Complementaria. Tabla combinada de VND (todos los genes) ----
# ── Anotación de región de interés (prioridad: hotspot > dominio) ─────────────
annot_ul97 <- function(pos) dplyr::case_when(
  pos >= 335 & pos <= 370  ~ "Hotspot GCV/MBV-R (aa 335\u2013370)",
  pos >= 397 & pos <= 480  ~ "Hotspot MBV-R (aa 397\u2013480)",
  pos >= 519 & pos <= 521  ~ "Hotspot GCV-R (aa 519\u2013521)",
  pos >= 590 & pos <= 607  ~ "Hotspot GCV-R (aa 590\u2013607)",
  pos >= 230 & pos <= 620  ~ "Dominio quinasa",
  TRUE                     ~ ""
)
annot_ul54 <- function(pos) dplyr::case_when(
  pos >= 300 & pos <= 325  ~ "Hotspot GCV/CDV (aa 300\u2013325)",
  pos >= 408 & pos <= 420  ~ "Hotspot GCV/CDV (aa 408\u2013420)",
  pos >= 498 & pos <= 545  ~ "Hotspot GCV/CDV (aa 498\u2013545)",
  pos >= 985 & pos <= 990  ~ "Hotspot GCV/CDV (aa 985\u2013990)",
  pos >= 598 & pos <= 600  ~ "Hotspot GCV/CDV/FOS (aa 598\u2013600)",
  pos >= 802 & pos <= 825  ~ "Hotspot GCV/CDV/FOS (aa 802\u2013825)",
  pos >= 700 & pos <= 715  ~ "Hotspot FOS/GCV (aa 700\u2013715)",
  pos >= 750 & pos <= 760  ~ "Hotspot FOS (aa 750\u2013760)",
  pos >= 781 & pos <= 790  ~ "Hotspot FOS (aa 781\u2013790)",
  pos >= 945 & pos <= 955  ~ "Hotspot FOS (aa 945\u2013955)",
  pos >= 290 & pos <= 550  ~ "Dominio exonucleasa",
  pos >= 696 & pos <= 990  ~ "Dominio polimerasa",
  pos >= 1150 & pos <= 1243 ~ "Regi\u00f3n NLS / uni\u00f3n ppUL44",
  TRUE                     ~ ""
)
annot_ul56 <- function(pos) dplyr::case_when(
  pos >= 230 & pos <= 370  ~ "Hotspot LMV-R (aa 230\u2013370)",
  TRUE                     ~ ""
)

# ── Formateadores (coma decimal, punto de millar) ─────────────────────────────
pc  <- function(x) formatC(x, format = "f", digits = 1, decimal.mark = ",")
mil <- function(x) formatC(x, format = "f", digits = 0, big.mark = ".", decimal.mark = ",")

# ── Denominadores (de tu objeto n_total_gen, bloque [3b]) ─────────────────────
n_total_gen <- tibble(
  Protein         = c("UL97", "UL54", "UL56"),
  n_tot_muestras  = c(26, 24, 11),
  n_tot_pacientes = c(11, 11, 3)
)

# ── Construcción de la tabla ──────────────────────────────────────────────────
tabla_vnd <- all_data_txp %>%
  filter(class == "NEW", Prevalence >= 5) %>%
  mutate(Position = as.numeric(Position)) %>%
  group_by(Protein, Mutation, Position) %>%
  summarise(
    n_muestras  = n_distinct(Npet),
    n_pacientes = n_distinct(NHC),
    prev_media  = round(mean(Prevalence, na.rm = TRUE), 1),
    prev_de     = round(sd(Prevalence,   na.rm = TRUE), 1),
    cob_mediana = round(median(cobertura_total, na.rm = TRUE), 0),
    cob_p25     = round(quantile(cobertura_total, 0.25, na.rm = TRUE), 0),
    cob_p75     = round(quantile(cobertura_total, 0.75, na.rm = TRUE), 0),
    .groups = "drop"
  ) %>%
  left_join(n_total_gen, by = "Protein") %>%
  mutate(
    pct_muestras  = round(n_muestras  / n_tot_muestras  * 100, 1),
    pct_pacientes = round(n_pacientes / n_tot_pacientes * 100, 1),
    col_muestras  = paste0(n_muestras,  " (", pc(pct_muestras),  "%)"),
    col_pacientes = paste0(n_pacientes, " (", pc(pct_pacientes), "%)"),
    col_prev = if_else(
      is.na(prev_de),
      paste0(pc(prev_media), " (n=1)"),
      paste0(pc(prev_media), " \u00b1 ", pc(prev_de))
    ),
    col_cob = if_else(
      n_muestras == 1,
      paste0(mil(cob_mediana), "\u00d7"),
      paste0(mil(cob_mediana), "\u00d7 (", mil(cob_p25), "\u2013", mil(cob_p75), ")")
    ),
    dom_hot = case_when(
      Protein == "UL97" ~ annot_ul97(Position),
      Protein == "UL54" ~ annot_ul54(Position),
      Protein == "UL56" ~ annot_ul56(Position),
      TRUE              ~ ""
    ),
    Protein = factor(Protein, levels = c("UL97", "UL54", "UL56"))
  ) %>%
  arrange(Protein, Position)

# ── Diagnóstico (verifica que son 101) ────────────────────────────────────────
cat("Total VND (todas):", nrow(tabla_vnd), "\n")
tabla_vnd %>% count(Protein) %>% print()
cat("Con Dominio/Hotspot asignado:", sum(tabla_vnd$dom_hot != ""), "\n")

# ── Dataframe para flextable ──────────────────────────────────────────────────
vnd_ft <- tabla_vnd %>%
  transmute(
    Gen = as.character(Protein), Variante = Mutation,
    col_muestras, col_pacientes, col_prev, col_cob, dom_hot
  )
names(vnd_ft) <- c(
  "Gen", "Variante", "Muestras, n (%)", "Pacientes, n (%)",
  "Frecuencia al\u00e9lica\n(media \u00b1 DE)",
  "Cobertura de posici\u00f3n\nMediana (RIC)",
  "Dominio / Hotspot"
)

# ── Bloques por gen (para sombreado y separadores) ───────────────────────────
n_ul97 <- sum(tabla_vnd$Protein == "UL97")
n_ul54 <- sum(tabla_vnd$Protein == "UL54")
n_ul56 <- sum(tabla_vnd$Protein == "UL56")
filas_ul97 <- seq_len(n_ul97)
filas_ul54 <- n_ul97 + seq_len(n_ul54)
filas_ul56 <- n_ul97 + n_ul54 + seq_len(n_ul56)
filas_sep  <- c(n_ul97, n_ul97 + n_ul54)

b_azul <- officer::fp_border(color = "#2E5FA3", width = 1.5)
b_sep  <- officer::fp_border(color = "#2E5FA3", width = 1.2)
b_gris <- officer::fp_border(color = "#CCCCCC", width = 0.5)
b_head <- officer::fp_border(color = "#FFFFFF", width = 0.8)

# ── Flextable ─────────────────────────────────────────────────────────────────
ft_vnd <- flextable::flextable(vnd_ft) %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = c(1, 2, 7), align = "left", part = "body") %>%
  flextable::valign(valign = "center", part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 9, part = "all") %>%
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
  flextable::width(j = 1, width = 1.6, unit = "cm") %>%
  flextable::width(j = 2, width = 1.9, unit = "cm") %>%
  flextable::width(j = 3, width = 2.3, unit = "cm") %>%
  flextable::width(j = 4, width = 2.3, unit = "cm") %>%
  flextable::width(j = 5, width = 2.7, unit = "cm") %>%
  flextable::width(j = 6, width = 3.3, unit = "cm") %>%
  flextable::width(j = 7, width = 4.6, unit = "cm") %>%
  flextable::padding(padding = 4, part = "all") %>%
  flextable::add_footer_lines(paste0(
    "VND: variante no descrita en CHARMD (Tilloy et al., 2024), umbral \u22655%. ",
    "Porcentajes sobre muestras con secuenciaci\u00f3n v\u00e1lida por gen (UL97: 26; UL54: 24; UL56: 11) ",
    "y sobre pacientes con dato por gen (UL97: 11; UL54: 11; UL56: 3). ",
    "Frecuencia al\u00e9lica: media \u00b1 DE; en variantes con una sola muestra se indica el valor \u00fanico (n=1). ",
    "Cobertura: lecturas FW+RV en la posici\u00f3n, mediana (RIC). ",
    "Dominio / Hotspot: detalle solo en variantes localizadas en una regi\u00f3n de inter\u00e9s; en blanco si quedan fuera. ",
    "El asterisco (*) indica cod\u00f3n de parada."
  )) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft_vnd

guardar_tabla(ft_vnd, "Complementaria_VND_todas", apaisada = TRUE)

