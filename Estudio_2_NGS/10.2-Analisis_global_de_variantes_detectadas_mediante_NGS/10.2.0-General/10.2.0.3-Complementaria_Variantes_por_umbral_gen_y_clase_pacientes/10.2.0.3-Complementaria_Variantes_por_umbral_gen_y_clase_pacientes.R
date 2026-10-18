# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.0.3. Complementaria. Variantes por umbral, gen y clase (pacientes no RTP)
# Extracto del script completo (ver README.md del estudio)

## ---- Complementaria. Variantes por umbral, gen y clase (pacientes no RTP) ----
# Complemento: todos los pacientes con datos que NO son TXP
all_data_notxp <- all_data_nhc %>%
  filter(!NHC %in% pacientes_txp)

cat("── Subcohorte NO TXP ───────────────────────────\n")
cat("Pacientes NO TXP con datos:", n_distinct(all_data_notxp$NHC), "\n")
cat("Npet únicos NO TXP:        ", n_distinct(all_data_notxp$Npet), "\n")
cat("Variantes totales NO TXP:  ", nrow(all_data_notxp), "\n")

# ── Construir dataframe ───────────────────────────────────────────────────────
filas_notxp <- map_dfr(genes, function(gene) {
  total_pac <- n_distinct(filter(all_data_notxp, Protein == gene)$NHC)
  map_dfr(thresholds, function(thr) {
    tibble(
      Gen    = if (thr == thresholds[1]) paste0(gene, "  (n = ", total_pac, " pac.)") else "",
      Umbral = paste0("\u2265 ", thr, "%"),
      MRA    = calcular_celda(all_data_notxp, gene, "MRA", thr),
      POL    = calcular_celda(all_data_notxp, gene, "POL", thr),
      VSI    = calcular_celda(all_data_notxp, gene, "VSI", thr),
      Nuevas = calcular_celda(all_data_notxp, gene, "NEW", thr)
    )
  })
})

filas_mra_notxp <- which(filas_notxp$MRA != "0 / 0")

ft_notxp <- flextable::flextable(filas_notxp) %>%
  flextable::set_header_labels(
    Gen = "Gen", Umbral = "Umbral",
    MRA = "MRA", POL = "POL", VSI = "VSI", Nuevas = "Nuevas"
  ) %>%
  flextable::add_header_row(
    values    = c("", "", "Variantes \u00fanicas / Pacientes afectados", "", "", ""),
    colwidths = c(1, 1, 1, 1, 1, 1),
    top       = TRUE
  ) %>%
  flextable::merge_at(i = 1, j = 3:6, part = "header") %>%
  flextable::merge_at(i = 1:4,  j = 1, part = "body") %>%
  flextable::merge_at(i = 5:8,  j = 1, part = "body") %>%
  flextable::merge_at(i = 9:12, j = 1, part = "body") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = 1:2, align = "left", part = "body") %>%
  flextable::valign(valign = "center", part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1, part = "body") %>%
  flextable::bg(part = "header", bg = "#2E5FA3") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::border_inner_v(border = b_head, part = "header") %>%
  flextable::bg(i = 1:4,  bg = "#EEF4FB", part = "body") %>%
  flextable::bg(i = 5:8,  bg = "#FFFFFF", part = "body") %>%
  flextable::bg(i = 9:12, bg = "#EEF4FB", part = "body") %>%
  flextable::bg(j = 2, bg = "#F5F5F5", part = "body") %>%
  flextable::color(j = 2, color = "#555555", part = "body") %>%
  flextable::italic(j = 2, part = "body") %>%
  flextable::bold(i = filas_mra_notxp, j = 3, part = "body") %>%
  flextable::color(i = filas_mra_notxp, j = 3, color = "#B22222", part = "body") %>%
  flextable::border_outer(border = b_azul, part = "all") %>%
  flextable::border_inner_h(border = b_gris, part = "body") %>%
  flextable::border_inner_v(border = b_gris, part = "body") %>%
  flextable::hline(i = c(4, 8), border = b_sep, part = "body") %>%
  flextable::hline(i = 1, border = b_sep, part = "header") %>%
  flextable::width(j = 1,   width = 4.0, unit = "cm") %>%
  flextable::width(j = 2,   width = 1.8, unit = "cm") %>%
  flextable::width(j = 3:6, width = 2.8, unit = "cm") %>%
  flextable::padding(padding = 6, part = "all") %>%
  flextable::add_footer_lines(
    "Valores expresados como: variantes \u00fanicas detectadas / pacientes con al menos una variante de esa clase al umbral indicado. MRA: mutaci\u00f3n de resistencia a antivirales; POL: polimorfismo natural; VSI: variante de significado incierto; Nuevas: variantes no clasificadas en CHARMD (Tilloy et al., 2024). An\u00e1lisis restringido a pacientes no pertenecientes a la subcohorte de trasplante pulmonar."
  ) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft_notxp

guardar_tabla(ft_notxp, "Complementaria_variantes_umbral_noRTP")

