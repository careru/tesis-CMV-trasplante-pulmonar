# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.4.04. Complementaria. VND del dominio quinasa de UL97
# Extracto del script completo (ver README.md del estudio)

## ---- Complementaria. VND del dominio quinasa de UL97 ----
ul97_quinasa_datos <- all_data_txp %>%
  filter(class == "NEW", Protein == "UL97", Prevalence >= 5) %>%
  mutate(Position = as.numeric(Position)) %>%
  filter(Position >= 230, Position <= 620) %>%
  group_by(Mutation, Position) %>%
  summarise(
    Paciente       = paste(sort(unique(RTP)), collapse = ", "),
    N_muestras     = n_distinct(Npet),
    Prev_min       = round(min(Prevalence), 1),
    Prev_max       = round(max(Prevalence), 1),
    # ── Métricas de calidad agregadas ────────────────────────────────────
    cob_mediana = round(median(cobertura_total, na.rm = TRUE), 0),
    cob_p25     = round(quantile(cobertura_total, 0.25, na.rm = TRUE), 0),
    cob_p75     = round(quantile(cobertura_total, 0.75, na.rm = TRUE), 0),
    cob_min     = round(min(cobertura_total, na.rm = TRUE), 0),
    q_med = round(median(q_score, na.rm = TRUE), 0),
    q_p25 = round(quantile(q_score, 0.25, na.rm = TRUE), 0),
    q_p75 = round(quantile(q_score, 0.75, na.rm = TRUE), 0),
    n_balance_ok   = sum(balance_ok, na.rm = TRUE),   # muestras con FW y RV > 0
    .groups = "drop"
  ) %>%
  mutate(
    Prevalencia = ifelse(
      Prev_min == Prev_max,
      paste0(format(Prev_max, decimal.mark = ","), "%"),
      paste0(format(Prev_min, decimal.mark = ","), "–",
             format(Prev_max, decimal.mark = ","), "%")
    ),
    Cobertura = paste0(
      format(cob_mediana, big.mark = "."), "× (",
      format(cob_p25, big.mark = "."), "–",
      format(cob_p75, big.mark = "."), "×; mín. ",
      format(cob_min, big.mark = "."), "×)"
    ),
    QScore = paste0(
      q_med, " (",
      q_p25, "–",
      q_p75, ")"
    ),
    Balance = paste0(n_balance_ok, "/", N_muestras)
  ) %>%
  arrange(Position)

hotspot_map <- tibble::tribble(
  ~Mutation, ~Hotspot,
  "K316R",   "No",
  "G340C",   "GCV/MBV-R (aa 335–370)",
  "R412Q",   "MBV-R (aa 397–480)",
  "W424*",   "MBV-R (aa 397–480)",
  "D481N",   "No",
  "A497P",   "No",
  "R516C",   "No",
  "D574G",   "No",
  "S604P",   "GCV-R (aa 590–607)",
  "L608P",   "No"
)

ul97_quinasa_datos <- ul97_quinasa_datos %>%
  left_join(hotspot_map, by = "Mutation") %>%
  select(Mutation, Paciente, N_muestras, Prevalencia,
         Cobertura, QScore, Balance, Hotspot)

names(ul97_quinasa_datos) <- c(
  "Mutación", "Paciente", "N\nmuestras", "Prevalencia",
  "Cobertura\n(mediana, RIC; mín.)", "Q-Score\n(mediana)",
  "Balance FW/RV\n(n muestras)", "Hotspot"
)

filas_hotspot <- which(ul97_quinasa_datos$Hotspot != "No")
filas_no      <- which(ul97_quinasa_datos$Hotspot == "No")
fila_w424     <- which(ul97_quinasa_datos[[1]] == "W424*")

b_azul <- officer::fp_border(color = "#2E5FA3", width = 1.5)
b_gris <- officer::fp_border(color = "#CCCCCC", width = 0.5)
b_head <- officer::fp_border(color = "#FFFFFF", width = 0.8)

ft_ul97_quinasa <- flextable::flextable(ul97_quinasa_datos) %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = c(1, 2, 8), align = "left", part = "body") %>%
  flextable::valign(valign = "center", part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1, part = "body") %>%
  flextable::bold(i = fila_w424, part = "body") %>%
  flextable::bg(i = fila_w424,   bg = "#FDECEA", part = "body") %>%
  flextable::bg(i = filas_hotspot[!filas_hotspot %in% fila_w424],
                bg = "#FEF9E7", part = "body") %>%
  flextable::bg(i = filas_no, bg = "#F0F4F8", part = "body") %>%
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
    "VND en dominio quinasa de UL97 (aa 230–620) al umbral \u22655%. ",
    "Cobertura: suma de lecturas forward + reverse en la posición (media entre muestras con la variante; mínimo observado). ",
    "Q-Score: puntuación de calidad de llamada de variante (mediana). ",
    "Balance FW/RV: número de muestras con lecturas en ambas cadenas sobre total de muestras con la variante. ",
    "W424*: codón de parada prematuro."
  )) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft_ul97_quinasa

guardar_tabla(ft_ul97_quinasa, "Complementaria_VND_UL97", apaisada = TRUE)

