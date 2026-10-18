# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.4.06. Complementaria. VND de UL54 en regiones funcionales
# Extracto del script completo (ver README.md del estudio)

## ---- Complementaria. VND de UL54 en regiones funcionales ----
# ── Funciones auxiliares UL54 (definir ANTES del bloque [5f]) ─────────────────

get_hotspot <- function(pos) {
  dplyr::case_when(
    (pos >= 300 & pos <= 325) |
      (pos >= 408 & pos <= 420) |
      (pos >= 498 & pos <= 545) |
      (pos >= 985 & pos <= 990)  ~ "GCV/CDV",
    (pos >= 598 & pos <= 600) |
      (pos >= 802 & pos <= 825)  ~ "GCV/CDV/FOS",
    (pos >= 700 & pos <= 715)  ~ "FOS > GCV",
    (pos >= 750 & pos <= 760) |
      (pos >= 781 & pos <= 790) |
      (pos >= 945 & pos <= 955)  ~ "FOS",
    TRUE                       ~ "No"
  )
}

get_dominio <- function(pos) {
  dplyr::case_when(
    pos >= 290  & pos <= 550  ~ "Exonucleasa",
    pos >= 696  & pos <= 990  ~ "Polimerasa",
    pos >= 1150 & pos <= 1243 ~ "NLS + unión ppUL44",
    TRUE                      ~ "Fuera de dominio"
  )
  
}

ul54_datos <- all_data_txp %>%
  filter(class == "NEW", Protein == "UL54", Prevalence >= 5) %>%
  mutate(Position = as.numeric(Position)) %>%
  group_by(Mutation, Position) %>%
  summarise(
    Paciente     = paste(sort(unique(RTP)), collapse = ", "),
    N_muestras   = n_distinct(Npet),
    N_pacientes  = n_distinct(NHC),
    Prev_min     = round(min(Prevalence), 1),
    Prev_max     = round(max(Prevalence), 1),
    cob_media    = round(mean(cobertura_total, na.rm = TRUE), 0),
    cob_min      = round(min(cobertura_total,  na.rm = TRUE), 0),
    q_med        = round(median(q_score, na.rm = TRUE), 0),
    n_balance_ok = sum(balance_ok, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Prevalencia = ifelse(
      Prev_min == Prev_max,
      paste0(format(Prev_max, decimal.mark = ","), "%"),
      paste0(format(Prev_min, decimal.mark = ","), "–",
             format(Prev_max, decimal.mark = ","), "%")
    ),
    Cobertura = paste0(format(cob_media, big.mark = "."), "×",
                       " (mín. ", format(cob_min, big.mark = "."), "×)"),
    QScore    = as.character(q_med),
    Balance   = paste0(n_balance_ok, "/", N_muestras),
    Hotspot   = get_hotspot(Position)
  ) %>%
  filter(
    Hotspot != "No" |
      (get_dominio(Position) %in% c("Exonucleasa", "Polimerasa", "NLS + unión ppUL44") &
         Prev_max >= 10)
  ) %>%
  arrange(Position) %>%
  select(Mutation, Paciente, N_muestras, Prevalencia,
         Cobertura, QScore, Balance, Hotspot)

names(ul54_datos) <- c(
  "Mutación", "Paciente", "N\nmuestras", "Prevalencia",
  "Cobertura\n(media, mín.)", "Q-Score\n(mediana)",
  "Balance FW/RV\n(n muestras)", "Hotspot"
)

filas_hotspot <- which(ul54_datos$Hotspot != "No")
filas_no      <- which(ul54_datos$Hotspot == "No")
filas_sanger  <- which(
  as.numeric(sub("%.*", "", sub(",", ".", ul54_datos$Prevalencia))) >= 20 |
    grepl("^[2-9][0-9]", ul54_datos$Prevalencia)
)

ft_ul54 <- flextable::flextable(ul54_datos) %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = c(1, 2, 8), align = "left", part = "body") %>%
  flextable::valign(valign = "center", part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1, part = "body") %>%
  flextable::bg(i = filas_hotspot, bg = "#FEF9E7", part = "body") %>%
  flextable::bg(i = filas_no,      bg = "#F0F4F8", part = "body") %>%
  flextable::bold(i = filas_sanger, part = "body") %>%
  flextable::bg(i = filas_sanger,   bg = "#FDECEA", part = "body") %>%
  flextable::bg(part = "header", bg = "#2E5FA3") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::border_inner_v(border = b_head, part = "header") %>%
  flextable::border_outer(border = b_azul, part = "all") %>%
  flextable::border_inner_h(border = b_gris, part = "body") %>%
  flextable::border_inner_v(border = b_gris, part = "body") %>%
  flextable::width(j = 1, width = 2.0, unit = "cm") %>%
  flextable::width(j = 2, width = 3.8, unit = "cm") %>%
  flextable::width(j = 3, width = 1.4, unit = "cm") %>%
  flextable::width(j = 4, width = 2.2, unit = "cm") %>%
  flextable::width(j = 5, width = 2.8, unit = "cm") %>%
  flextable::width(j = 6, width = 1.8, unit = "cm") %>%
  flextable::width(j = 7, width = 2.2, unit = "cm") %>%
  flextable::width(j = 8, width = 3.8, unit = "cm") %>%
  flextable::padding(padding = 5, part = "all") %>%
  flextable::add_footer_lines(paste0(
    "VND en UL54 seleccionadas al umbral \u22655% (en hotspot o en dominio funcional con prevalencia \u226515% o >3 pacientes). ",
    "Cobertura: lecturas totales FW+RV (media y mínimo entre muestras con la variante). ",
    "Q-Score: mediana de la puntuación de calidad de llamada. ",
    "Balance FW/RV: muestras con soporte en ambas cadenas. ",
    "Fondo amarillo: en hotspot de resistencia. Fondo azul: fuera de hotspot, incluida por dominio/prevalencia. ",
    "Fondo rojo: prevalencia \u226520% (equivalente Sanger). ",
    "Hotspot según regiones documentadas para GCV/CDV/FOS en UL54."
  )) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft_ul54

guardar_tabla(ft_ul54, "Complementaria_VND_UL54", apaisada = TRUE)

