# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.4.10. Anexo II. VND con región funcional y evidencia técnica
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo II. VND con región funcional y evidencia técnica ----
# El bloque guarda el Word apaisado (Anexo_II.docx) en la carpeta de salidas.
## Requiere all_data_txp (con class y columnas DeepChek: cobertura_total, mut_fw, mut_rv, q_score).
## Una fila por VND; frecuencia y cobertura resumen todas las muestras en que se detectó.


fmt_n <- function(x) format(round(x), big.mark = ".", decimal.mark = ",", trim = TRUE)
fmt_1 <- function(x) sub("\\.", ",", sprintf("%.1f", x))

# ── Regiones funcionales y hotspots (editar límites según Material y métodos) ──
## Se asigna la primera fila que contenga la posición: los hotspots van antes que los dominios.
regiones_funcionales <- tribble(
  ~Protein, ~ini,  ~fin,  ~etiqueta,
  # UL97
  "UL97",   335,   370,   "GCV/MBV-R (aa 335\u2013370)",
  "UL97",   397,   480,   "MBV-R (aa 397\u2013480)",
  "UL97",   590,   607,   "GCV-R (aa 590\u2013607)",
  "UL97",   335,   665,   "GCV-R (aa 335\u2013665)",
  "UL97",   337,   651,   "pUL97-KD",
  # UL54
  "UL54",   300,   325,   "GCV/CDV-R (aa 300\u2013325)",
  "UL54",   498,   545,   "GCV/CDV-R (aa 498\u2013545)",
  "UL54",   700,   715,   "FOS/GCV-R (aa 700\u2013715)",
  "UL54",   750,   760,   "FOS-R (aa 750\u2013760)",
  "UL54",   802,   825,   "GCV/CDV/FOS-R (aa 802\u2013825)",
  "UL54",   985,   990,   "GCV/CDV-R (aa 985\u2013990)",
  "UL54",   300,   545,   "pUL54-Exo",
  "UL54",   696,   988,   "pUL54-Pol",
  "UL54",  1153,  1242,   "pUL54-Cter/NLS",
  # UL56
  "UL56",   230,   370,   "LMV-R (aa 230\u2013370)"
)
asignar_region <- function(prot, pos) {
  r <- regiones_funcionales %>% filter(Protein == prot, pos >= ini, pos <= fin)
  if (nrow(r) == 0) "" else r$etiqueta[1]
}

# ── Denominadores por gen ────────────────────────────────────────────────────
denom <- all_data_txp %>% group_by(Protein) %>%
  summarise(N_m = n_distinct(Npet), N_p = n_distinct(RTP), .groups = "drop")

# ── Tabla ────────────────────────────────────────────────────────────────────
anexo2 <- all_data_txp %>%
  filter(class == "NEW", Prevalence >= 5) %>%
  mutate(mut_tot = mut_fw + mut_rv,
         balance = pmin(mut_fw, mut_rv) / mut_tot) %>%
  group_by(Protein, Mutation) %>%
  summarise(
    n_m       = n_distinct(Npet),
    n_p       = n_distinct(RTP),
    pacientes = paste(unique(RTP[order(as.numeric(sub("RTP-", "", RTP)))]), collapse = ", "),
    prev_mean = mean(Prevalence), prev_sd = sd(Prevalence),
    cob_med   = median(cobertura_total),
    cob_q1    = quantile(cobertura_total, .25), cob_q3 = quantile(cobertura_total, .75),
    lect_min  = min(mut_tot),
    fw_rv_min = paste0(mut_fw[which.min(mut_tot)], "/", mut_rv[which.min(mut_tot)]),
    bal_min   = min(balance),
    q_min     = min(q_score),
    .groups = "drop") %>%
  left_join(denom, by = "Protein") %>%
  mutate(
    Posicion  = as.numeric(gsub("[^0-9]", "", Mutation)),
    Region    = mapply(asignar_region, Protein, Posicion),
    Gen       = factor(Protein, levels = c("UL97", "UL54", "UL56"),
                       labels = c(sprintf("pUL97 (%d muestras, %d pacientes)", denom$N_m[denom$Protein == "UL97"], denom$N_p[denom$Protein == "UL97"]),
                                  sprintf("pUL54 (%d muestras, %d pacientes)", denom$N_m[denom$Protein == "UL54"], denom$N_p[denom$Protein == "UL54"]),
                                  sprintf("pUL56 (%d muestras, %d pacientes)", denom$N_m[denom$Protein == "UL56"], denom$N_p[denom$Protein == "UL56"]))),
    Muestras   = sprintf("%d (%s)", n_m, fmt_1(100 * n_m / N_m)),
    Pacientes  = sprintf("%d (%s)", n_p, fmt_1(100 * n_p / N_p)),
    Frecuencia = ifelse(n_m == 1, fmt_1(prev_mean),
                        sprintf("%s \u00b1 %s", fmt_1(prev_mean), fmt_1(prev_sd))),
    Cobertura  = ifelse(n_m == 1, paste0(fmt_n(cob_med), " \u00d7"),
                        sprintf("%s \u00d7 (%s\u2013%s)", fmt_n(cob_med), fmt_n(cob_q1), fmt_n(cob_q3))),
    Lecturas   = sprintf("%d (%s)", lect_min, fw_rv_min),
    Hebra      = fmt_1(100 * bal_min),
    Q          = as.character(q_min)
  ) %>%
  arrange(Gen, Posicion) %>%
  select(Gen, Mutation, Muestras, Pacientes, pacientes, Frecuencia, Cobertura, Lecturas, Hebra, Q, Region)

cat("VND por gen:\n"); print(table(anexo2$Gen))

# ── Flextable con fila combinada por gen ─────────────────────────────────────
anexo2_grp <- as_grouped_data(anexo2, groups = "Gen")

ft_anexo2 <- as_flextable(anexo2_grp, hide_grouplabel = TRUE) %>%
  set_header_labels(
    Mutation = "Variante",
    Muestras = "Muestras, n (%)", Pacientes = "Pacientes, n (%)",
    pacientes = "Pacientes (Tabla X)",
    Frecuencia = "Frecuencia de variante (%), media \u00b1 DE",
    Cobertura = "Cobertura por posici\u00f3n, mediana (RIC)",
    Lecturas = "Lecturas mutadas, m\u00ednimo (fw/rv)",
    Hebra = "Hebra minoritaria, m\u00ednimo (%)",
    Q = "Q-score, m\u00ednimo",
    Region = "Regi\u00f3n funcional / hotspot de resistencia"
  ) %>%
  align(align = "center", part = "all") %>%
  align(j = c(1, 10), align = "left", part = "body") %>%
  bold(part = "header") %>%
  bold(i = ~ !is.na(Gen), bold = TRUE) %>%
  bg(i = ~ !is.na(Gen), bg = "#1F3B73") %>%
  color(i = ~ !is.na(Gen), color = "white") %>%
  bg(part = "header", bg = "#1F3B73") %>% color(part = "header", color = "white") %>%
  bold(j = "Mutation", bold = TRUE) %>%
  flextable::font(fontname = "Source Sans Pro", part = "all") %>%
  fontsize(size = 7.5, part = "all") %>%
  border_inner(border = fp_border(color = "grey70", width = 0.5)) %>%
  border_outer(border = fp_border(color = "grey40", width = 1)) %>%
  add_footer_lines(paste0(
    "VND: variantes no descritas previamente (no recogidas en CHARMD, versi\u00f3n actualizada a 20 de abril de 2026). ",
    "Muestras y pacientes sobre el total con resultado v\u00e1lido para cada gen. Frecuencia y cobertura resumen las muestras en ",
    "que se detect\u00f3 la variante (valor \u00fanico cuando se detect\u00f3 en una sola). Lecturas mutadas: n\u00famero m\u00ednimo de lecturas que ",
    "soportan la variante en cualquiera de las muestras, con su reparto entre hebras directa (fw) e inversa (rv). Hebra minoritaria: ",
    "proporci\u00f3n m\u00ednima de lecturas mutadas en la hebra menos representada. El filtro por Q-score del sistema experto estaba ",
    "desactivado. pUL97-KD: dominio quinasa de pUL97; Exo: dominio exonucleasa; Pol: dominio polimerasa; Cter/NLS: extremo ",
    "C-terminal y se\u00f1al de localizaci\u00f3n nuclear; -R: regi\u00f3n asociada a resistencia; GCV: ganciclovir; CDV: cidofovir; ",
    "FOS: foscarnet; MBV: maribavir; LMV: letermovir. * Cod\u00f3n de terminaci\u00f3n prematuro.")) %>%
  fontsize(size = 6.5, part = "footer") %>%
  set_table_properties(layout = "autofit", width = 1)

ft_anexo2

# ── Exportar a Word en horizontal ────────────────────────────────────────────
doc <- read_docx() %>%
  body_add_par("Anexo II. Variantes no descritas previamente identificadas mediante NGS", style = "heading 2") %>%
  body_add_flextable(ft_anexo2) %>%
  body_end_section_landscape()
print(doc, target = file.path(dir_salida, "Anexo_II.docx"))


