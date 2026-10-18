# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.0.2. Complementaria. Variantes por umbral, gen y clase (RTP)
# Extracto del script completo (ver README.md del estudio)

## ---- Complementaria. Variantes por umbral, gen y clase (RTP) ----
thresholds <- c(5, 10, 15, 20)
genes      <- c("UL97", "UL54", "UL56")

calcular_celda <- function(data, gene, clase, thr) {
  d <- data %>%
    filter(Protein == gene, class == clase, Prevalence >= thr)
  n_var <- n_distinct(d$Mutation)
  n_pac <- n_distinct(d$NHC)
  if (n_var == 0) return("0 / 0")
  paste0(n_var, " / ", n_pac)
}

# ── Construir dataframe ───────────────────────────────────────────────────────
filas <- map_dfr(genes, function(gene) {
  total_pac <- n_distinct(filter(all_data_txp, Protein == gene)$NHC)
  map_dfr(thresholds, function(thr) {
    tibble(
      Gen    = if (thr == thresholds[1]) paste0(gene, "  (n = ", total_pac, " pac.)") else "",
      Umbral = paste0("\u2265 ", thr, "%"),
      MRA    = calcular_celda(all_data_txp, gene, "MRA",   thr),
      POL    = calcular_celda(all_data_txp, gene, "POL",   thr),
      VSI    = calcular_celda(all_data_txp, gene, "VSI",   thr),
      Nuevas = calcular_celda(all_data_txp, gene, "NEW",   thr)
    )
  })
})

# ── Identificar filas con MRA > 0 para colorear ───────────────────────────────
filas_mra <- which(filas$MRA != "0 / 0")

# ── Bordes ────────────────────────────────────────────────────────────────────
b_azul <- officer::fp_border(color = "#2E5FA3", width = 1.5)
b_sep  <- officer::fp_border(color = "#2E5FA3", width = 1.2)
b_gris <- officer::fp_border(color = "#CCCCCC", width = 0.5)
b_head <- officer::fp_border(color = "#FFFFFF", width = 0.8)

# ── Flextable ─────────────────────────────────────────────────────────────────
ft <- flextable::flextable(filas) %>%
  
  # Cabeceras
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
  
  # Fusión Gen
  flextable::merge_at(i = 1:4,  j = 1, part = "body") %>%
  flextable::merge_at(i = 5:8,  j = 1, part = "body") %>%
  flextable::merge_at(i = 9:12, j = 1, part = "body") %>%
  
  # Alineación
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = 1:2, align = "left", part = "body") %>%
  flextable::valign(valign = "center", part = "all") %>%
  
  # Tipografía base
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1, part = "body") %>%
  
  # Cabecera: fondo azul oscuro, texto blanco
  flextable::bg(part = "header", bg = "#2E5FA3") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::border_inner_v(border = b_head, part = "header") %>%
  
  # Cuerpo: sombreado alternado por bloque gen
  flextable::bg(i = 1:4,  bg = "#EEF4FB", part = "body") %>%
  flextable::bg(i = 5:8,  bg = "#FFFFFF", part = "body") %>%
  flextable::bg(i = 9:12, bg = "#EEF4FB", part = "body") %>%
  
  # Columna Umbral: fondo gris neutro para guía visual
  flextable::bg(j = 2, bg = "#F5F5F5", part = "body") %>%
  flextable::color(j = 2, color = "#555555", part = "body") %>%
  flextable::italic(j = 2, part = "body") %>%
  
  # MRA con color si > 0
  flextable::bold(i = filas_mra, j = 3, part = "body") %>%
  flextable::color(i = filas_mra, j = 3, color = "#B22222", part = "body") %>%
  
  # Bordes
  flextable::border_outer(border = b_azul, part = "all") %>%
  flextable::border_inner_h(border = b_gris, part = "body") %>%
  flextable::border_inner_v(border = b_gris, part = "body") %>%
  flextable::hline(i = c(4, 8), border = b_sep, part = "body") %>%
  flextable::hline(i = 1, border = b_sep, part = "header") %>%
  
  # Anchos
  flextable::width(j = 1,   width = 4.0, unit = "cm") %>%
  flextable::width(j = 2,   width = 1.8, unit = "cm") %>%
  flextable::width(j = 3:6, width = 2.8, unit = "cm") %>%
  
  # Padding
  flextable::padding(padding = 6, part = "all") %>%
  
  # Nota al pie
  flextable::add_footer_lines(
    "Valores expresados como: variantes \u00fanicas detectadas / pacientes con al menos una variante de esa clase al umbral indicado. MRA: mutaci\u00f3n de resistencia a antivirales; POL: polimorfismo natural; VSI: variante de significado incierto; Nuevas: variantes no clasificadas en CHARMD (Tilloy et al., 2024). An\u00e1lisis restringido a la subcohorte de trasplante pulmonar (n\u00famero de pacientes con datos por gen indicado entre par\u00e9ntesis)."
  ) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft

guardar_tabla(ft, "Complementaria_variantes_umbral_RTP")

