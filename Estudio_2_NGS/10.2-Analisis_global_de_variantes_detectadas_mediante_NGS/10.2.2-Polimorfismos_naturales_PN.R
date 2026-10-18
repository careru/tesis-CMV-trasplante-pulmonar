# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.2. Polimorfismos naturales (PN)
# Extracto del script completo (ver README.md del estudio)

# ---- 10.2.2. Polimorfismos naturales (PN) ----

## ---- Datos de PN ----
all_data_txp %>%
  filter(class == "POL", Prevalence >= 5) %>%
  group_by(Protein) %>%
  summarise(
    n_variantes = n_distinct(Mutation),
    n_pacientes = n_distinct(NHC),
    .groups = "drop"
  )

all_data_txp %>%
  filter(class == "POL", Prevalence >= 5) %>%
  group_by(Protein, Mutation) %>%
  summarise(
    NHC_list = paste(sort(unique(NHC)), collapse = ", "),
    n_pac    = n_distinct(NHC),
    .groups  = "drop"
  ) %>%
  arrange(Protein, desc(n_pac)) %>%
  print(n = Inf)

## ---- Anexo I. Polimorfismos naturales identificados mediante NGS ----
# El bloque guarda el Word apaisado (Anexo_I.docx) en la carpeta de salidas.
## Requiere all_data_txp (con class y columnas DeepChek). Ejecutar antes el bloque del
## Anexo II para disponer de fmt_n y fmt_1, o dejar las dos funciones de abajo.


fmt_n <- function(x) format(round(x), big.mark = ".", decimal.mark = ",", trim = TRUE)
fmt_1 <- function(x) sub("\\.", ",", sprintf("%.1f", x))

denom <- all_data_txp %>% group_by(Protein) %>%
  summarise(N_m = n_distinct(Npet), N_p = n_distinct(RTP), .groups = "drop")
etq_gen <- setNames(sprintf("p%s (%d muestras, %d pacientes)", denom$Protein, denom$N_m, denom$N_p), denom$Protein)

pn <- all_data_txp %>%
  filter(class == "POL", Prevalence >= 5) %>%
  mutate(mut_tot = mut_fw + mut_rv, balance = pmin(mut_fw, mut_rv) / mut_tot) %>%
  group_by(Protein, Mutation) %>%
  summarise(
    n_m = n_distinct(Npet), n_p = n_distinct(RTP),
    pacientes = paste(unique(RTP[order(as.numeric(sub("RTP-", "", RTP)))]), collapse = ", "),
    prev_mean = mean(Prevalence), prev_sd = sd(Prevalence),
    cob_med = median(cobertura_total), cob_q1 = quantile(cobertura_total, .25), cob_q3 = quantile(cobertura_total, .75),
    lect_min = min(mut_tot), fw_rv_min = paste0(mut_fw[which.min(mut_tot)], "/", mut_rv[which.min(mut_tot)]),
    bal_min = min(balance), q_min = min(q_score), .groups = "drop") %>%
  left_join(denom, by = "Protein") %>%
  mutate(
    Posicion   = as.numeric(gsub("[^0-9]", "", Mutation)),
    Gen        = factor(etq_gen[Protein], levels = etq_gen[c("UL97", "UL54", "UL56")]),
    Muestras   = sprintf("%d (%s)", n_m, fmt_1(100 * n_m / N_m)),
    Pacientes  = sprintf("%d (%s)", n_p, fmt_1(100 * n_p / N_p)),
    Frecuencia = ifelse(n_m == 1, fmt_1(prev_mean), sprintf("%s \u00b1 %s", fmt_1(prev_mean), fmt_1(prev_sd))),
    Cobertura  = ifelse(n_m == 1, paste0(fmt_n(cob_med), " \u00d7"),
                        sprintf("%s \u00d7 (%s\u2013%s)", fmt_n(cob_med), fmt_n(cob_q1), fmt_n(cob_q3))),
    Lecturas   = sprintf("%d (%s)", lect_min, fw_rv_min),
    Hebra      = fmt_1(100 * bal_min),
    Q          = as.character(q_min)
  ) %>%
  arrange(Gen, Posicion) %>%
  select(Gen, Mutation, Muestras, Pacientes, pacientes, Frecuencia, Cobertura, Lecturas, Hebra, Q)

cat("PN por gen:\n"); print(table(pn$Gen))

# ── Flextable con fila combinada por gen ─────────────────────────────────────
pn_grp <- as_grouped_data(pn, groups = "Gen")

ft_pn <- as_flextable(pn_grp, hide_grouplabel = TRUE) %>%
  set_header_labels(
    Mutation = "Variante", Muestras = "Muestras, n (%)", Pacientes = "Pacientes, n (%)",
    pacientes = "Pacientes (Tabla X)",
    Frecuencia = "Frecuencia de variante (%), media \u00b1 DE",
    Cobertura = "Cobertura por posici\u00f3n, mediana (RIC)",
    Lecturas = "Lecturas mutadas, m\u00ednimo (fw/rv)",
    Hebra = "Hebra minoritaria, m\u00ednimo (%)",
    Q = "Q-score, m\u00ednimo") %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "body") %>%
  bold(part = "header") %>%
  bold(i = ~ !is.na(Gen), bold = TRUE) %>%                       # filas de gen en negrita
  bg(i = ~ !is.na(Gen), bg = "#1F3B73") %>%                      # fondo azul oscuro en las filas de gen
  color(i = ~ !is.na(Gen), color = "white") %>%
  bg(part = "header", bg = "#1F3B73") %>% color(part = "header", color = "white") %>%
  bold(j = "Mutation", bold = TRUE) %>%
  flextable::font(fontname = "Source Sans Pro", part = "all") %>%
  fontsize(size = 7.5, part = "all") %>%
  border_inner(border = fp_border(color = "grey70", width = 0.5)) %>%
  border_outer(border = fp_border(color = "grey40", width = 1)) %>%
  add_footer_lines(paste0(
    "PN: polimorfismos naturales (variantes recogidas en CHARMD, versi\u00f3n actualizada a 20 de abril de 2026, sin efecto ",
    "sobre la susceptibilidad antiviral; incluye D605E, reclasificada como PN). Muestras y pacientes sobre el total con ",
    "resultado v\u00e1lido para cada gen. Frecuencia y cobertura resumen las muestras en que se detect\u00f3 la variante. ",
    "Lecturas mutadas: n\u00famero m\u00ednimo de lecturas que soportan la variante, con su reparto entre hebras directa (fw) e ",
    "inversa (rv). Hebra minoritaria: proporci\u00f3n m\u00ednima de lecturas mutadas en la hebra menos representada. ",
    "El filtro por Q-score del sistema experto estaba desactivado.")) %>%
  fontsize(size = 6.5, part = "footer") %>%
  set_table_properties(layout = "autofit", width = 1)

ft_pn

# ── Exportar a Word en horizontal (carpeta Descargas) ─────────────────────────
doc <- read_docx() %>%
  body_add_par("Anexo. Polimorfismos naturales identificados mediante NGS", style = "heading 2") %>%
  body_add_flextable(ft_pn) %>%
  body_end_section_landscape()
print(doc, target = file.path(dir_salida, "Anexo_I.docx"))

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

## ---- Complementaria. Cobertura de los PN ----
# ── Datos ─────────────────────────────────────────────────────────────────────
pol_cob <- all_data_txp %>%
  filter(class == "POL", Prevalence >= 5) %>%
  mutate(
    Position = as.integer(str_extract(Mutation, "\\d+")),
    Protein  = factor(Protein, levels = c("UL97", "UL54", "UL56"))
  ) %>%
  group_by(Protein, Mutation, Position) %>%
  summarise(
    mediana  = median(cobertura_total, na.rm = TRUE),
    p25      = quantile(cobertura_total, 0.25, na.rm = TRUE),
    p75      = quantile(cobertura_total, 0.75, na.rm = TRUE),
    n        = n_distinct(Npet),
    prev_max = max(Prevalence, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  arrange(Protein, Position) %>%
  mutate(
    Mutation  = factor(Mutation, levels = unique(Mutation)),
    unica     = n == 1,
    baja_prev = prev_max < 20
  )

# ── Colores por gen ───────────────────────────────────────────────────────────
col_gen <- c("UL97" = "#4472C4", "UL54" = "#C8880A", "UL56" = "#C0392B")

# ── Shapes: rombo (18) si prev_max < 20%, círculo (16/21) si no ──────────────
# Combinamos unica y baja_prev en un indicador de shape
pol_cob <- pol_cob %>%
  mutate(
    shape_tipo = case_when(
      baja_prev & unica  ~ "baja_unica",
      baja_prev & !unica ~ "baja_multi",
      !baja_prev & unica ~ "norm_unica",
      TRUE               ~ "norm_multi"
    )
  )

shape_map <- c(
  "norm_multi" = 16,   # círculo relleno
  "norm_unica" = 21,   # círculo hueco
  "baja_multi" = 18,   # rombo relleno
  "baja_unica" = 23    # rombo hueco
)

# ── Plot ──────────────────────────────────────────────────────────────────────
p_cob_pol <- ggplot(pol_cob,
                    aes(x = Mutation, y = mediana, color = Protein)) +
  
  # RIC (solo si n > 1)
  geom_errorbar(
    data = pol_cob %>% filter(!unica),
    aes(ymin = p25, ymax = p75),
    width = 0.35, linewidth = 0.6, alpha = 0.7
  ) +
  
  # Puntos con shape según tipo
  geom_point(
    aes(shape = shape_tipo,
        fill  = ifelse(unica | shape_tipo == "baja_unica", "white", Protein)),
    size = 3.2, stroke = 0.8
  ) +
  scale_shape_manual(
    values = shape_map,
    labels = c(
      "norm_multi" = "Frecuencia variante \u226520%",
      "norm_unica" = "Frecuencia variante \u226520%, muestra \u00fanica",
      "baja_multi" = "Frecuencia variante <20%",
      "baja_unica" = "Frecuencia variante <20%, muestra \u00fanica"
    ),
    name = NULL
  ) +
  scale_fill_manual(values = c(col_gen, "white" = "white"), guide = "none") +
  scale_color_manual(values = col_gen, guide = "none") +
  
  scale_y_continuous(
    labels  = function(x) format(x, big.mark = ".", decimal.mark = ","),
    expand  = expansion(mult = c(0.05, 0.12))
  ) +
  
  facet_wrap(~ Protein, scales = "free_x", nrow = 1) +
  
  labs(
    x       = NULL,
    y       = "Cobertura de posici\u00f3n (x)",
    caption = paste0(
    )
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1,
                                    size = 12, face = "bold"),
    axis.text.y      = element_text(size = 10),
    strip.text       = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "#EEF4FB", color = NA),
    panel.spacing    = unit(1.2, "cm"),
    legend.position  = "bottom",
    legend.text      = element_text(size = 10),
    plot.caption     = element_text(size = 10, color = "grey40",
                                    hjust = 0, margin = margin(t = 8)),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin      = margin(t = 10, r = 15, b = 10, l = 10)
  )

p_cob_pol

# [guardado antiguo] ggsave(
# [guardado antiguo]   "~/Desktop/Figuras_tesis/figura_cobertura_PN.png",
# [guardado antiguo]   p_cob_pol,
# [guardado antiguo]   width = 32, height = 14, units = "cm", dpi = 400,
# [guardado antiguo]   bg = "white"
# [guardado antiguo] )

guardar_fig(p_cob_pol, "Complementaria_cobertura_PN", ancho = 32, alto = 14, dpi = 400)


