# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.4.08. Complementaria. VND de UL56
# Extracto del script completo (ver README.md del estudio)

## ---- Complementaria. VND de UL56 ----
ul56_dom_datos <- all_data_txp %>%
  filter(class == "NEW", Protein == "UL56", Prevalence >= 5) %>%
  mutate(Position = as.numeric(Position)) %>%
  group_by(Mutation, Position) %>%
  summarise(
    Paciente     = paste(sort(unique(RTP)), collapse = ", "),
    N_muestras   = n_distinct(Npet),
    Prev_min     = round(min(Prevalence), 1),
    Prev_max     = round(max(Prevalence), 1),
    cob_mediana = round(median(cobertura_total, na.rm = TRUE), 0),
    cob_p25     = round(quantile(cobertura_total, 0.25, na.rm = TRUE), 0),
    cob_p75     = round(quantile(cobertura_total, 0.75, na.rm = TRUE), 0),
    cob_min     = round(min(cobertura_total, na.rm = TRUE), 0),
    q_med       = round(median(q_score, na.rm = TRUE), 0),
    q_p25       = round(quantile(q_score, 0.25, na.rm = TRUE), 0),
    q_p75       = round(quantile(q_score, 0.75, na.rm = TRUE), 0),
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
    Cobertura = ifelse(
      N_muestras == 1,
      paste0(format(cob_mediana, big.mark = "."), "×"),
      paste0(
        format(cob_mediana, big.mark = "."), "× (",
        format(cob_p25, big.mark = "."), "–",
        format(cob_p75, big.mark = "."), "×; mín. ",
        format(cob_min, big.mark = "."), "×)"
      )
    ),
    QScore = ifelse(
      N_muestras == 1,
      as.character(q_med),
      paste0(q_med, " (", q_p25, "–", q_p75, ")")
    ),
    Balance = paste0(n_balance_ok, "/", N_muestras)
  ) %>%
  arrange(Position)

hotspot_map_ul56 <- tibble::tribble(
  ~Mutation, ~Hotspot,
  "Y281H",   "LMV-R (aa 230–370)",
  "N330S",   "LMV-R (aa 230–370)",
  "D480G",   "No",
  "I508V",   "No",
  "R513C",   "No",
  "Y595C",   "No",
  "N599S",   "No",
  "L612P",   "No",
  "M697I",   "No"
)

ul56_dom_datos <- ul56_dom_datos %>%
  left_join(hotspot_map_ul56, by = "Mutation") %>%
  mutate(Hotspot = replace_na(Hotspot, "No")) %>%
  select(Mutation, Paciente, N_muestras, Prevalencia,
         Cobertura, QScore, Balance, Hotspot)

names(ul56_dom_datos) <- c(
  "Mutación", "Paciente", "N\nmuestras", "Prevalencia",
  "Cobertura\n(mediana, RIC; mín.)", "Q-Score\n(mediana, RIC)",
  "Balance FW/RV\n(n muestras)", "Hotspot"
)

filas_hotspot <- which(ul56_dom_datos$Hotspot != "No")
filas_no      <- which(ul56_dom_datos$Hotspot == "No")
fila_r513c    <- which(ul56_dom_datos[[1]] == "R513C")

ft_ul56_dom <- flextable::flextable(ul56_dom_datos) %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = c(1, 2, 8), align = "left", part = "body") %>%
  flextable::valign(valign = "center", part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1, part = "body") %>%
  flextable::bold(i = fila_r513c, part = "body") %>%
  flextable::bg(i = fila_r513c, bg = "#FDECEA", part = "body") %>%
  flextable::bg(i = filas_hotspot[!filas_hotspot %in% fila_r513c],
                bg = "#FEF9E7", part = "body") %>%
  flextable::bg(i = filas_no[!filas_no %in% fila_r513c],
                bg = "#F0F4F8", part = "body") %>%
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
  flextable::width(j = 8, width = 3.5, unit = "cm") %>%
  flextable::padding(padding = 5, part = "all") %>%
  flextable::add_footer_lines(paste0(
    "VND en UL56 al umbral \u22655%. ",
    "Cobertura: lecturas totales FW+RV en la posición (media y mínimo entre muestras con la variante). ",
    "Q-Score: mediana de la puntuación de calidad de llamada. ",
    "Balance FW/RV: muestras con soporte en ambas cadenas. ",
    "R513C: única variante con prevalencia \u226520% (equivalente Sanger). ",
    "Hotspot LMV-R según CHARMD (Tilloy et al., 2024)."
  )) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft_ul56_dom

guardar_tabla(ft_ul56_dom, "Complementaria_VND_UL56", apaisada = TRUE)

