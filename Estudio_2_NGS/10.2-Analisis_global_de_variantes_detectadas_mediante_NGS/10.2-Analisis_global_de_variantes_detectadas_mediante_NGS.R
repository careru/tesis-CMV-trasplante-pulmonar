# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2. Análisis global de variantes detectadas mediante NGS
# Extracto del script completo (ver README.md del estudio)

# ---- 10.2. Análisis global de variantes detectadas mediante NGS ----

## ---- Muestras analizadas por gen (26 muestras; UL97 26, UL54 24, UL56 11) ----
npet_txp_total <- all_data_nhc %>%
  filter(NHC %in% pacientes_txp) %>%
  distinct(Npet) %>%
  nrow()
cat("Total peticiones únicas (Npet) en subcohorte TXP:", npet_txp_total, "\n\n")

cat("── Npet con resultado válido por gen ────────────────────────────────────\n")
all_data_nhc %>%
  filter(NHC %in% pacientes_txp) %>%
  group_by(Protein) %>%
  summarise(
    n_npet      = n_distinct(Npet),
    pct         = round(n_distinct(Npet) / npet_txp_total * 100, 1),
    n_pacientes = n_distinct(NHC),
    .groups     = "drop"
  ) %>%
  filter(Protein %in% c("UL97", "UL54", "UL56")) %>%
  arrange(Protein) %>%
  print()

npet_ul97 <- all_data_nhc %>%
  filter(NHC %in% pacientes_txp, Protein == "UL97") %>% distinct(Npet)
npet_ul54 <- all_data_nhc %>%
  filter(NHC %in% pacientes_txp, Protein == "UL54") %>% distinct(Npet)
npet_ul56 <- all_data_nhc %>%
  filter(NHC %in% pacientes_txp, Protein == "UL56") %>% distinct(Npet)

cat("\n── Npet en UL97 pero no en UL54 (cobertura insuficiente u omitidas) ─────\n")
npet_solo_ul97 <- setdiff(npet_ul97$Npet, npet_ul54$Npet)
cat("N =", length(npet_solo_ul97), "\n")
cat(npet_solo_ul97, sep = "\n")

cat("\n── Npet en UL54 pero no en UL97 ─────────────────────────────────────────\n")
npet_solo_ul54 <- setdiff(npet_ul54$Npet, npet_ul97$Npet)
cat("N =", length(npet_solo_ul54), "\n")
cat(npet_solo_ul54, sep = "\n")

cat("\n── Npet con UL56 analizado ───────────────────────────────────────────────\n")
cat("N =", nrow(npet_ul56), "\n")
cat("Porcentaje sobre total TXP:",
    round(nrow(npet_ul56) / npet_txp_total * 100, 1), "%\n")

cat("\n── Desglose por Npet: genes analizados (TRUE/FALSE) ─────────────────────\n")
all_data_nhc %>%
  filter(NHC %in% pacientes_txp) %>%
  distinct(Npet, NHC) %>%
  left_join(npet_ul97 %>% mutate(UL97 = TRUE), by = "Npet") %>%
  left_join(npet_ul54 %>% mutate(UL54 = TRUE), by = "Npet") %>%
  left_join(npet_ul56 %>% mutate(UL56 = TRUE), by = "Npet") %>%
  mutate(across(c(UL97, UL54, UL56), ~replace_na(., FALSE))) %>%
  left_join(
    all_data_nhc %>% filter(NHC %in% pacientes_txp) %>% distinct(Npet, RTP),
    by = "Npet"
  ) %>%
  arrange(RTP, Npet) %>%
  print(n = Inf)

# ── Diagnóstico calidad: resumen de cobertura y Q-Score en subcohorte TXP ────
cat("\n── Resumen de métricas de calidad (subcohorte TXP, Prevalence ≥5%) ──────\n")
all_data_txp %>%
  filter(Protein %in% c("UL97", "UL54", "UL56")) %>%
  group_by(Protein) %>%
  summarise(
    cob_media   = round(mean(cobertura_total, na.rm = TRUE), 0),
    cob_mediana = round(median(cobertura_total, na.rm = TRUE), 0),
    cob_min     = round(min(cobertura_total, na.rm = TRUE), 0),
    q_mediana   = round(median(q_score, na.rm = TRUE), 1),
    pct_balance = round(mean(balance_ok, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) %>%
  rename(
    "Gen"             = Protein,
    "Cob. media (×)"  = cob_media,
    "Cob. mediana (×)"= cob_mediana,
    "Cob. mínima (×)" = cob_min,
    "Q-Score mediana" = q_mediana,
    "% balance FW/RV" = pct_balance
  ) %>%
  print()

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

## ---- Complementaria. Cobertura de secuenciación por paciente y gen ----
# ── Pacientes y número de muestras por paciente ───────────────────────────────
pacientes_cob <- tibble::tribble(
  ~RTP,      ~n,
  "RTP-1",    1,
  "RTP-2",   15,
  "RTP-3",    1,
  "RTP-4",    1,
  "RTP-5",    1,
  "RTP-6",    1,
  "RTP-7",    1,
  "RTP-8",    1,
  "RTP-9",    1,
  "RTP-10",   1,
  "RTP-11",   2
)
genes <- c("UL97", "UL54", "UL56")

# ── Estado por muestra (válido / cobertura insuficiente / no analizado) ───────
estado_muestra <- function(rtp, gen, n) {
  if (gen == "UL97") return(rep("valido", n))
  if (gen == "UL54") {
    v <- rep("valido", n)
    if (rtp == "RTP-2") v[c(6, 14)] <- "insuf"   # muestras 6 y 14 (días +361, +479)
    return(v)
  }
  # UL56: analizado solo en RTP-2 (8 de 15), RTP-7 (1) y RTP-11 (2)
  if (rtp == "RTP-2")  return(c(rep("valido", 8), rep("no", 7)))
  if (rtp == "RTP-7")  return(rep("valido", n))
  if (rtp == "RTP-11") return(rep("valido", n))
  rep("no", n)
}

# ── Expandir a una fila por muestra ───────────────────────────────────────────
datos <- pacientes_cob %>%
  tidyr::crossing(Gen = factor(genes, levels = genes)) %>%
  mutate(estado = pmap(list(RTP, as.character(Gen), n), estado_muestra)) %>%
  unnest(estado) %>%
  group_by(RTP, Gen) %>%
  mutate(idx = row_number()) %>%
  ungroup()

# ── Maquetación de la cuadrícula (5 casillas por fila dentro de cada celda) ───
ncol_cell <- 5
gap_pac   <- 0.7
gene_base <- c(UL97 = 0, UL54 = 7, UL56 = 14)   # x de inicio de cada gen

pac_layout <- pacientes_cob %>%
  mutate(
    RTP   = factor(RTP, levels = pacientes_cob$RTP),
    n_sub = ceiling(n / ncol_cell)
  ) %>%
  arrange(RTP) %>%
  mutate(
    y_top    = c(0, -head(cumsum(n_sub + gap_pac), -1)),
    y_centro = y_top - (n_sub - 1) / 2
  )

datos <- datos %>%
  left_join(pac_layout %>% select(RTP, y_top), by = "RTP") %>%
  mutate(
    col    = (idx - 1) %% ncol_cell,
    sub    = (idx - 1) %/% ncol_cell,
    x      = gene_base[as.character(Gen)] + col,
    y      = y_top - sub,
    estado = factor(estado, levels = c("valido", "insuf", "no"))
  )

# ── Etiquetas y resumen (resumen calculado desde los datos) ──────────────────
lab_genes <- tibble(Gen = names(gene_base),
                    x = gene_base + (ncol_cell - 1) / 2,
                    y = max(datos$y) + 1.4)

resumen <- datos %>%
  group_by(Gen) %>%
  summarise(validos = sum(estado == "valido"), total = n(), .groups = "drop") %>%
  mutate(
    pct     = validos / total * 100,
    pct_txt = ifelse(abs(pct - 100) < 1e-9, "100",
                     formatC(pct, format = "f", digits = 1, decimal.mark = ",")),
    label   = paste0(validos, "/", total, " \u00b7 ", pct_txt, "%"),
    x       = gene_base[as.character(Gen)] + (ncol_cell - 1) / 2,
    y       = min(datos$y) - 1.6
  )

# ── Colores por estado ────────────────────────────────────────────────────────
col_fill <- c(valido = "#CDEBC5", insuf = "#F6CDA0", no = "#E6E6E6")
col_line <- c(valido = "#6FB16A", insuf = "#DC8A3C", no = "#C2C2C2")

# ── Figura ────────────────────────────────────────────────────────────────────
p_cobertura <- ggplot(datos, aes(x, y)) +
  geom_tile(aes(fill = estado, color = estado),
            width = 0.82, height = 0.82, linewidth = 0.5) +
  geom_text(data = pac_layout,
            aes(x = -4.6, y = y_centro, label = RTP),
            hjust = 0, fontface = "bold", size = 3.6,
            color = "grey20", inherit.aes = FALSE) +
  geom_text(data = pac_layout,
            aes(x = -0.5, y = y_centro, label = paste0("n=", n)),
            hjust = 1, size = 3.0, color = "grey50", inherit.aes = FALSE) +
  geom_text(data = lab_genes,
            aes(x = x, y = y, label = Gen),
            fontface = "bold", size = 4.3, color = "grey20", inherit.aes = FALSE) +
  geom_text(data = resumen,
            aes(x = x, y = y, label = label),
            fontface = "bold", size = 3.4, color = "grey25", inherit.aes = FALSE) +
  scale_fill_manual(
    values = col_fill,
    labels = c("Resultado válido", "Cobertura insuficiente", "No analizado"),
    name   = NULL
  ) +
  scale_color_manual(values = col_line, guide = "none") +
  scale_x_continuous(limits = c(-5, 19), expand = c(0, 0)) +
  coord_equal(clip = "off") +
  theme_void(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size = 11),
    legend.key.size = unit(0.5, "cm"),
    plot.margin     = margin(10, 15, 10, 15)
  )

p_cobertura

# [guardado antiguo] ggsave(
# [guardado antiguo]   filename = "~/Desktop/Figuras_tesis/figura_cobertura_NGS_paciente_gen.png",
# [guardado antiguo]   plot     = p_cobertura,
# [guardado antiguo]   width    = 20, height = 22, units = "cm", dpi = 600, bg = "white"
# [guardado antiguo] )

guardar_fig(p_cobertura, "Complementaria_cobertura_paciente_gen", ancho = 20, alto = 22, dpi = 600)


