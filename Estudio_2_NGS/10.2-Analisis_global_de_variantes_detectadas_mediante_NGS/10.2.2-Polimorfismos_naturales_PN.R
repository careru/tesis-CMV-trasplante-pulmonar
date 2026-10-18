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

