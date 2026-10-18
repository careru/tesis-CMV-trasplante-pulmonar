# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.4. Variantes no descritas previamente (VND)
# Extracto del script completo (ver README.md del estudio)

# ---- 10.2.4. Variantes no descritas previamente (VND) ----

## ---- Datos de VND ----
# ── Resumen actualizado de variantes NEW ──────────────────────────────────────
new_vars <- all_data_txp %>%
  filter(class == "NEW", Prevalence >= 5) %>%
  mutate(Position = as.numeric(Position)) %>%
  group_by(Protein, Mutation, Position) %>%
  summarise(
    n_muestras  = n_distinct(Npet),
    n_pacientes = n_distinct(NHC),
    prev_media  = round(mean(Prevalence, na.rm = TRUE), 1),
    prev_max    = round(max(Prevalence, na.rm = TRUE), 1),
    .groups = "drop"
  )

cat("Total NEW por gen:\n")
new_vars %>% count(Protein) %>% print()

cat("\nNEW por gen y umbral:\n")
map_dfr(c(5, 10, 15, 20), function(thr) {
  all_data_txp %>%
    filter(class == "NEW", Prevalence >= thr) %>%
    group_by(Protein) %>%
    summarise(n_variantes = n_distinct(Mutation),
              n_pacientes = n_distinct(NHC), .groups = "drop") %>%
    mutate(umbral = paste0("≥", thr, "%"))
}) %>%
  arrange(Protein, umbral) %>%
  print(n = Inf)


new_vars_dom <- new_vars %>%
  mutate(
    dominio = case_when(
      Protein == "UL97" & Position >= 335 & Position <= 370 ~ "P-loop",
      Protein == "UL97" & Position >= 409 & Position <= 411 ~ "Gatekeeper",
      Protein == "UL97" & Position >= 440 & Position <= 480 ~ "Loop catalítico",
      Protein == "UL97" & Position >= 520 & Position <= 607 ~ "Región canónica GCV",
      Protein == "UL54" & Position >= 408 & Position <= 412 ~ "Región IV",
      Protein == "UL54" & Position >= 501 & Position <= 545 ~ "Delta-C",
      Protein == "UL54" & Position >= 700 & Position <= 715 ~ "Región II",
      Protein == "UL54" & Position >= 725 & Position <= 964 ~ "Unión ATP/dNTP",
      Protein == "UL54" & Position >= 802 & Position <= 821 ~ "Región III",
      Protein == "UL54" & Position >= 910 & Position <= 913 ~ "Región I",
      Protein == "UL54" & Position >= 962 & Position <= 987 ~ "Región V",
      Protein == "UL54" & Position >= 1153 & Position <= 1159 ~ "NLS-A",
      Protein == "UL54" & Position >= 1213 & Position <= 1242 ~ "UL44 binding",
      Protein == "UL56" & Position == 25 ~ "N-terminal",
      Protein == "UL56" & Position >= 202 & Position <= 220 ~ "LAGLIDADG",
      Protein == "UL56" & Position >= 231 & Position <= 369 ~ "Región LMV",
      Protein == "UL56" & Position >= 514 & Position <= 522 ~ "Interacción UL89",
      Protein == "UL56" & Position >= 671 & Position <= 680 ~ "Interacción UL89",
      Protein == "UL56" & Position >= 757 & Position <= 764 ~ "Interacción UL89",
      TRUE ~ "Fuera de dominio funcional"
    ),
    en_dominio = dominio != "Fuera de dominio funcional"
  ) %>%
  arrange(Protein, desc(en_dominio), desc(n_pacientes), desc(prev_max))

cat("── En dominio funcional ─────────────────────────\n")
new_vars_dom %>% filter(en_dominio) %>% print(n = Inf)

cat("\n── Fuera de dominio pero relevantes ─────────────\n")
new_vars_dom %>%
  filter(!en_dominio, n_pacientes >= 2 | prev_max >= 15) %>%
  print(n = Inf)

## ---- Complementaria. VND en dominios funcionales (todos los genes) ----
# ── Datos: variantes nuevas en dominio funcional ──────────────────────────────
new_dom_resumen <- new_vars_dom %>%
  filter(en_dominio) %>%
  mutate(
    rtp = case_when(
      Protein == "UL54" & Mutation %in% c(
        "D711Y","D854G","E530D","E846*","F702L","F763L","F817L",
        "G967C","G971C","I833V","I840V","M1219I","M959V","P704S",
        "P712H","P744H","P797H","R1157S","R876L","R901L","S913G",
        "G340C","Y281H") ~ "RTP-2",
      Protein == "UL54" & Mutation %in% c("A714V","L736M","R1157M","T760A") ~ "RTP-1",
      Protein == "UL54" & Mutation %in% c("A835T","V772A") ~ "RTP-10",
      Protein == "UL54" & Mutation == "Q511R" ~ "RTP-3",
      Protein == "UL56" & Mutation == "N330S" ~ "RTP-11",
      Protein == "UL56" & Mutation == "Y281H" ~ "RTP-2",
      Protein == "UL97" & Mutation == "D574G" ~ "RTP-1",
      Protein == "UL97" & Mutation == "S604P" ~ "RTP-10",
      Protein == "UL97" & Mutation == "G340C" ~ "RTP-2",
      TRUE ~ "otro"
    )
  ) %>%
  arrange(Protein, dominio, desc(prev_max)) %>%
  mutate(
    Protein = factor(Protein, levels = c("UL54", "UL56", "UL97")),
    col_prev = paste0(format(prev_max, decimal.mark = ","), "%")
  ) %>%
  select(Protein, Mutation, dominio, prev_max, rtp)

# ── Construir dataframe para flextable ────────────────────────────────────────
ft_new_dom <- new_dom_resumen %>%
  group_by(Protein) %>%
  mutate(
    Gen = if_else(row_number() == 1, as.character(Protein), "")
  ) %>%
  ungroup() %>%
  select(Gen, Mutation, dominio, prev_max, rtp)

names(ft_new_dom) <- c(
  "Gen", "Variante",
  "Dominio funcional",
  "Prevalencia máx. (%)",
  "Paciente"
)

# ── Filas por bloque gen ──────────────────────────────────────────────────────
n_ul54 <- sum(new_dom_resumen$Protein == "UL54")
n_ul56 <- sum(new_dom_resumen$Protein == "UL56")
n_ul97 <- sum(new_dom_resumen$Protein == "UL97")

filas_ul54 <- 1:n_ul54
filas_ul56 <- (n_ul54 + 1):(n_ul54 + n_ul56)
filas_ul97 <- (n_ul54 + n_ul56 + 1):(n_ul54 + n_ul56 + n_ul97)
filas_sep  <- c(n_ul54, n_ul54 + n_ul56)

# ── Bordes ────────────────────────────────────────────────────────────────────
b_azul <- officer::fp_border(color = "#2E5FA3", width = 1.5)
b_sep  <- officer::fp_border(color = "#2E5FA3", width = 1.2)
b_gris <- officer::fp_border(color = "#CCCCCC", width = 0.5)
b_head <- officer::fp_border(color = "#FFFFFF", width = 0.8)

# ── Flextable ─────────────────────────────────────────────────────────────────
ft_new <- flextable::flextable(ft_new_dom) %>%
  flextable::merge_at(i = filas_ul54, j = 1, part = "body") %>%
  flextable::merge_at(i = filas_ul56, j = 1, part = "body") %>%
  flextable::merge_at(i = filas_ul97, j = 1, part = "body") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = 1:3, align = "left", part = "body") %>%
  flextable::valign(valign = "center", part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1:2, part = "body") %>%
  flextable::bg(part = "header", bg = "#2E5FA3") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::border_inner_v(border = b_head, part = "header") %>%
  flextable::bg(i = filas_ul54, bg = "#EEF4FB", part = "body") %>%
  flextable::bg(i = filas_ul56, bg = "#FFFFFF", part = "body") %>%
  flextable::bg(i = filas_ul97, bg = "#EEF4FB", part = "body") %>%
  flextable::border_outer(border = b_azul, part = "all") %>%
  flextable::border_inner_h(border = b_gris, part = "body") %>%
  flextable::border_inner_v(border = b_gris, part = "body") %>%
  flextable::hline(i = filas_sep, border = b_sep, part = "body") %>%
  flextable::width(j = 1, width = 2.0, unit = "cm") %>%
  flextable::width(j = 2, width = 2.5, unit = "cm") %>%
  flextable::width(j = 3, width = 4.5, unit = "cm") %>%
  flextable::width(j = 4, width = 2.5, unit = "cm") %>%
  flextable::width(j = 5, width = 2.0, unit = "cm") %>%
  flextable::padding(padding = 5, part = "all") %>%
  flextable::add_footer_lines(
    "Tabla restringida a variantes no clasificadas en CHARMD (Tilloy et al., 2024) detectadas en al menos un dominio funcional conocido, al umbral de prevalencia \u2265 5%. Prevalencia m\u00e1xima expresada como el valor m\u00e1ximo observado entre todas las muestras en que se detect\u00f3 la variante."
  ) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft_new

guardar_tabla(ft_new, "Complementaria_VND_dominios")

## ---- Figura 45. VND en UL97 ----
dominios_ul97 <- data.frame(
  start = c(230),
  end   = c(620),
  label = c("Dominio quinasa"),
  fill  = c("lightskyblue")
)

hotspots_ul97 <- data.frame(
  start = c(335, 397, 519, 590),
  end   = c(370, 480, 521, 607),
  label = c("", "", "", ""),
  fill  = c("#B03060", "#8E44AD", "#E67E22", "#E67E22")
)

ul97_umbrales <- all_data_txp %>%
  filter(class == "NEW", Protein == "UL97") %>%
  group_by(Mutation, Position) %>%
  summarise(
    prev_max   = max(Prevalence, na.rm = TRUE),
    umbral_max = case_when(
      max(Prevalence) >= 20 ~ 20,
      max(Prevalence) >= 15 ~ 15,
      max(Prevalence) >= 10 ~ 10,
      TRUE                  ~ 5
    ),
    .groups = "drop"
  ) %>%
  mutate(Position = as.numeric(Position))

ul97 <- new_vars_dom %>%
  filter(Protein == "UL97") %>%
  left_join(ul97_umbrales %>% select(Mutation, umbral_max), by = "Mutation")

nivel <- c("5" = 0.6, "10" = 1.1, "15" = 1.6, "20" = 2.1)
paso  <- 0.4

df <- ul97 %>%
  mutate(
    y_nivel    = nivel[as.character(umbral_max)],
    en_quinasa = Position >= 230 & Position <= 620,
    y_point    = ifelse(en_quinasa, -y_nivel, y_nivel),
    color      = ifelse(en_quinasa, "#C0392B", "#2980B9"),
    label_mut  = case_when(
      en_quinasa                        ~ Mutation,
      prev_max >= 10 | n_pacientes >= 4 ~ Mutation,
      TRUE                              ~ NA_character_
    )
  )

df_lab_azul <- df %>%
  filter(!is.na(label_mut), !en_quinasa) %>%
  arrange(Position) %>%
  mutate(prox_group = cumsum(c(TRUE, diff(Position) > 30))) %>%
  group_by(prox_group) %>%
  mutate(
    within_rank = row_number(),
    y_lab       = 3.0 + (within_rank - 1) * 0.60
  ) %>%
  ungroup()

df_lab_rojo <- df %>%
  filter(!is.na(label_mut), en_quinasa) %>%
  arrange(Position) %>%
  mutate(prox_group = cumsum(c(TRUE, diff(Position) > 30))) %>%
  group_by(prox_group) %>%
  mutate(
    within_rank = row_number(),
    y_lab       = -(3.0 + (within_rank - 1) * 0.60)
  ) %>%
  ungroup()

ymax_plot <- max(
  max(df_lab_azul$y_lab, na.rm = TRUE),
  max(abs(df_lab_rojo$y_lab), na.rm = TRUE)
) + 0.5

niveles_ref <- data.frame(
  y     = c(0.6, 1.1, 1.6, 2.1, -0.6, -1.1, -1.6, -2.1),
  label = rep(c("\u22655%", "\u226510%", "\u226515%", "\u226520%"), 2),
  lado  = rep(c("arriba", "abajo"), each = 4)
)

p_ul97 <- ggplot() +
  
  # Líneas de referencia — empiezan después de las etiquetas, acaban en aa 707
  geom_segment(data = niveles_ref,
               aes(x = -42, xend = 707, y = y, yend = y),
               color = "grey85", linewidth = 0.3, linetype = "dashed") +
  
  # Etiquetas de prevalencia
  geom_text(data = niveles_ref %>% filter(lado == "arriba"),
            aes(x = -55, y = y, label = label),
            hjust = 1, size = 4.0, color = "grey50", fontface = "bold") +
  geom_text(data = niveles_ref %>% filter(lado == "abajo"),
            aes(x = -55, y = y, label = label),
            hjust = 1, size = 4.0, color = "grey50", fontface = "bold") +
  
  # Dominios
  geom_rect(data = dominios_ul97,
            aes(xmin = start, xmax = end, ymin = 0.0, ymax = 0.40, fill = fill),
            alpha = 0.9, color = "white", linewidth = 0.25) +
  geom_text(data = dominios_ul97 %>% mutate(mid = (start + end) / 2),
            aes(x = mid, y = 0.20, label = label),
            size = 4.0, color = "black", fontface = "bold",
            vjust = 0.5, hjust = 0.5) +
  
  # Hotspots
  geom_rect(data = hotspots_ul97,
            aes(xmin = start, xmax = end, ymin = -0.40, ymax = 0.0, fill = fill),
            alpha = 0.9, color = "white", linewidth = 0.25) +
  geom_text(data = hotspots_ul97 %>% mutate(mid = (start + end) / 2),
            aes(x = mid, y = -0.20, label = label),
            size = 3, color = "white", fontface = "bold",
            vjust = 0.5, hjust = 0.5) +
  
  scale_fill_identity() +
  
  # Backbone aa 1–707 con tapas
  geom_segment(aes(x = 1,   xend = 707, y = 0,     yend = 0),
               color = "grey10", linewidth = 1.0) +
  geom_segment(aes(x = 1,   xend = 1,   y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  geom_segment(aes(x = 707, xend = 707, y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  
  # Etiquetas extremos backbone
  annotate("text", x = 1,   y = -0.28, label = "1",   size = 4.0,
           color = "grey10", vjust = 1) +
  annotate("text", x = 707, y = -0.28, label = "707", size = 4.0,
           color = "grey10", vjust = 1) +
  
  # Stems
  geom_segment(data = df,
               aes(x = Position, xend = Position,
                   y = ifelse(en_quinasa, -0.42, 0.42),
                   yend = y_point, color = color),
               linewidth = 0.45, alpha = 0.8) +
  
  # Puntos
  geom_point(data = df,
             aes(x = Position, y = y_point, color = color, size = prev_max),
             alpha = 0.88) +
  scale_color_identity() +
  scale_size_continuous(range = c(1.2, 4.5), guide = "none") +
  
  # Guías azules
  geom_segment(data = df_lab_azul,
               aes(x = Position, xend = Position,
                   y = y_point + 0.08, yend = y_lab - 0.08),
               color = "#2980B9", linewidth = 0.3, alpha = 0.4,
               linetype = "dotted") +
  
  # Guías rojas
  geom_segment(data = df_lab_rojo,
               aes(x = Position, xend = Position,
                   y = y_point - 0.08, yend = y_lab + 0.08),
               color = "#C0392B", linewidth = 0.3, alpha = 0.4,
               linetype = "dotted") +
  
  # Etiquetas azules
  geom_text(data = df_lab_azul,
            aes(x = Position, y = y_lab, label = label_mut),
            color = "#2980B9", size = 4.0, fontface = "bold",
            hjust = 0.5, vjust = 0) +
  
  # Etiquetas rojas
  geom_text(data = df_lab_rojo,
            aes(x = Position, y = y_lab, label = label_mut),
            color = "#C0392B", size = 4.0, fontface = "bold",
            hjust = 0.5, vjust = 1) +
  
  scale_x_continuous(
    limits = c(-100, 730),
    breaks = NULL,
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(limits = c(-ymax_plot, ymax_plot)) +
  labs(title = NULL, x = "", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(t = 8, r = 12, b = 8, l = 35)
  )

leyenda_hotspots <- ggplot(data.frame(
  x     = c(1.5, 3.7, 6),
  y     = c(1, 1, 1),
  label = c("GCV/MBV-R (335-370)",
            "MBV-R (397-480)",
            "GCV-R cod\u00f3n 520 + (590-607)"),
  color = c("#B03060", "#8E44AD", "#E67E22")
)) +
  geom_point(aes(x, y, color = color), size = 4, shape = 15) +
  geom_text(aes(x + 0.15, y, label = label),
            hjust = 0, size = 4.0, color = "grey20") +
  scale_color_identity() +
  scale_x_continuous(limits = c(0.8, 9.5)) +
  scale_y_continuous(limits = c(0.5, 1.5)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(t = 0, b = 6, l = 35)
  )

figura_ul97 <- p_ul97 / leyenda_hotspots +
  plot_layout(heights = c(10, 1))

figura_ul97
# [guardado antiguo] ggsave(
# [guardado antiguo]   "~/Desktop/Figuras_tesis/figura_lollipop_UL97.png",
# [guardado antiguo]   figura_ul97,
# [guardado antiguo]   width = 30, height = 15, units = "cm", dpi = 400,
# [guardado antiguo]   bg = "white"
# [guardado antiguo] )

guardar_fig(figura_ul97, "Figura_45", ancho = 30, alto = 15, dpi = 400)

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

## ---- Figura 46. VND en UL54 ----
# ── Dominios UL54 ─────────────────────────────────────────────────────────────
dominios_ul54 <- data.frame(
  start = c(290, 696, 1150),
  end   = c(550, 990, 1243),
  label = c("Exonucleasa", "Polimerasa", "NLS + unión ppUL44"),
  fill  = c("#82E0AA", "#7FB3D3", "#E8A0D0")
)

# ── Hotspots UL54 ─────────────────────────────────────────────────────────────
hotspots_ul54 <- data.frame(
  start = c(300, 408, 498, 985,   598, 802,   700,   750, 781, 945),
  end   = c(325, 420, 545, 990,   600, 825,   715,   760, 790, 955),
  fill  = c(rep("#E67E22", 4),
            rep("#8E44AD", 2),
            "#E74C3C",
            rep("#C0392B", 3))
)

# ── Datos UL54 con umbral_max ─────────────────────────────────────────────────
ul54_umbrales <- all_data_txp %>%
  filter(class == "NEW", Protein == "UL54") %>%
  group_by(Mutation, Position) %>%
  summarise(
    prev_max   = max(Prevalence, na.rm = TRUE),
    umbral_max = case_when(
      max(Prevalence) >= 20 ~ 20,
      max(Prevalence) >= 15 ~ 15,
      max(Prevalence) >= 10 ~ 10,
      TRUE                  ~ 5
    ),
    .groups = "drop"
  ) %>%
  mutate(Position = as.numeric(Position))

ul54 <- new_vars_dom %>%
  filter(Protein == "UL54") %>%
  left_join(ul54_umbrales %>% select(Mutation, umbral_max), by = "Mutation")

# ── Regiones de resistencia ───────────────────────────────────────────────────
en_region_resistencia <- function(pos) {
  (pos >= 300 & pos <= 325) |
    (pos >= 408 & pos <= 420) |
    (pos >= 498 & pos <= 545) |
    (pos >= 598 & pos <= 600) |
    (pos >= 700 & pos <= 715) |
    (pos >= 750 & pos <= 760) |
    (pos >= 781 & pos <= 790) |
    (pos >= 802 & pos <= 825) |
    (pos >= 945 & pos <= 955) |
    (pos >= 985 & pos <= 990)
}

# ── Niveles fijos ─────────────────────────────────────────────────────────────
nivel <- c("5" = 0.6, "10" = 1.1, "15" = 1.6, "20" = 2.1)

df <- ul54 %>%
  mutate(
    en_resist = en_region_resistencia(as.numeric(Position)),
    y_nivel   = nivel[as.character(umbral_max)],
    y_point   = ifelse(en_resist, -y_nivel, y_nivel),
    color     = ifelse(en_resist, "#C0392B", "#2980B9"),
    label_mut = case_when(
      en_resist                         ~ Mutation,
      prev_max >= 10 | n_pacientes >= 4 ~ Mutation,
      TRUE                              ~ NA_character_
    )
  )

# ── Etiquetas con agrupación por proximidad ───────────────────────────────────
df_lab_azul <- df %>%
  filter(!is.na(label_mut), !en_resist) %>%
  arrange(Position) %>%
  mutate(prox_group = cumsum(c(TRUE, diff(Position) > 40))) %>%
  group_by(prox_group) %>%
  mutate(
    within_rank = row_number(),
    y_lab       = 3.0 + (within_rank - 1) * 0.60
  ) %>%
  ungroup()

df_lab_rojo <- df %>%
  filter(!is.na(label_mut), en_resist) %>%
  arrange(Position) %>%
  mutate(
    escalon = case_when(
      Position < 817  ~ (row_number() - 1) %% 3,
      TRUE            ~ (row_number() - 1 - sum(Position < 817)) %% 3
    ),
    y_lab = -(3.0 + escalon * 0.60)
  )%>%
  ungroup()

ymax_plot <- max(df_lab_azul$y_lab, na.rm = TRUE) + 0.5
ymin_plot <- min(df_lab_rojo$y_lab, na.rm = TRUE) - 0.5

# ── Niveles de referencia ─────────────────────────────────────────────────────
niveles_ref <- data.frame(
  y     = c(0.6, 1.1, 1.6, 2.1, -0.6, -1.1, -1.6, -2.1),
  label = rep(c("\u22655%", "\u226510%", "\u226515%", "\u226520%"), 2),
  lado  = rep(c("arriba", "abajo"), each = 4)
)

# ── Plot UL54 ─────────────────────────────────────────────────────────────────
p_ul54 <- ggplot() +
  
  # Líneas de referencia — de -55 a aa 1243
  geom_segment(data = niveles_ref,
               aes(x = -55, xend = 1243, y = y, yend = y),
               color = "grey85", linewidth = 0.3, linetype = "dashed") +
  
  # Etiquetas de prevalencia
  geom_text(data = niveles_ref %>% filter(lado == "arriba"),
            aes(x = -70, y = y, label = label),
            hjust = 1, size = 4.0, color = "grey50", fontface = "bold") +
  geom_text(data = niveles_ref %>% filter(lado == "abajo"),
            aes(x = -70, y = y, label = label),
            hjust = 1, size = 4.0, color = "grey50", fontface = "bold") +
  
  # Dominios
  geom_rect(data = dominios_ul54,
            aes(xmin = start, xmax = end, ymin = 0.0, ymax = 0.40, fill = fill),
            alpha = 0.85, color = "white", linewidth = 0.25) +
  geom_text(data = dominios_ul54 %>% mutate(mid = (start + end) / 2),
            aes(x = mid, y = 0.23, label = label),
            size = 4.0, color = "grey15", fontface = "bold",
            vjust = 0.5, hjust = 0.5) +
  
  # Hotspots
  geom_rect(data = hotspots_ul54,
            aes(xmin = start, xmax = end, ymin = -0.40, ymax = 0.0, fill = fill),
            alpha = 0.9, color = "white", linewidth = 0.25) +
  
  scale_fill_identity() +
  
  # Backbone aa 1–1243 con tapas
  geom_segment(aes(x = 1,    xend = 1243, y = 0,     yend = 0),
               color = "grey10", linewidth = 1.0) +
  geom_segment(aes(x = 1,    xend = 1,    y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  geom_segment(aes(x = 1243, xend = 1243, y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  
  # Etiquetas extremos backbone
  annotate("text", x = 1,    y = -0.28, label = "1",    size = 4.0,
           color = "grey10", vjust = 1) +
  annotate("text", x = 1243, y = -0.28, label = "1243", size = 4.0,
           color = "grey10", vjust = 1) +
  
  # Stems
  geom_segment(data = df,
               aes(x = Position, xend = Position,
                   y = ifelse(en_resist, -0.42, 0.42),
                   yend = y_point, color = color),
               linewidth = 0.45, alpha = 0.8) +
  
  # Puntos
  geom_point(data = df,
             aes(x = Position, y = y_point, color = color, size = prev_max),
             alpha = 0.88) +
  scale_color_identity() +
  scale_size_continuous(range = c(1.2, 4.5), guide = "none") +
  
  # Guías azules
  geom_segment(data = df_lab_azul,
               aes(x = Position, xend = Position,
                   y = y_point + 0.08, yend = y_lab - 0.08),
               color = "#2980B9", linewidth = 0.3, alpha = 0.4,
               linetype = "dotted") +
  
  # Guías rojas
  geom_segment(data = df_lab_rojo,
               aes(x = Position, xend = Position,
                   y = y_point - 0.08, yend = y_lab + 0.08),
               color = "#C0392B", linewidth = 0.3, alpha = 0.4,
               linetype = "dotted") +
  
  # Etiquetas azules
  geom_text(data = df_lab_azul,
            aes(x = Position, y = y_lab, label = label_mut),
            color = "#2980B9", size = 4.0, fontface = "bold",
            hjust = 0.5, vjust = 0) +
  
  # Etiquetas rojas
  geom_text(data = df_lab_rojo,
            aes(x = Position, y = y_lab, label = label_mut),
            color = "#C0392B", size = 4.0, fontface = "bold",
            hjust = 0.5, vjust = 1) +
  
  scale_x_continuous(
    limits = c(-130, 1380),
    breaks = NULL,
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(limits = c(ymin_plot, ymax_plot)) +
  labs(title = NULL, x = "", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(t = 8, r = 12, b = 8, l = 50)
  )

# ── Leyenda hotspots ──────────────────────────────────────────────────────────
leyenda_hotspots_ul54 <- ggplot(data.frame(
  x     = c(1.0, 4.2, 7.4, 10.6),
  y     = 1,
  label = c("GCV/CDV", "GCV/CDV/FOS", "FOS > GCV", "FOS"),
  color = c("#E67E22", "#8E44AD", "#E74C3C", "#C0392B")
)) +
  geom_point(aes(x, y, color = color), size = 4, shape = 15) +
  geom_text(aes(x + 0.18, y, label = label),
            hjust = 0, size = 4.0, color = "grey20") +
  scale_color_identity() +
  scale_x_continuous(limits = c(0.5, 13)) +
  scale_y_continuous(limits = c(0.5, 1.5)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(t = 0, b = 6, l = 50)
  )

# ── Composición final ─────────────────────────────────────────────────────────
figura_ul54 <- p_ul54 / leyenda_hotspots_ul54 +
  plot_layout(heights = c(10, 1)) +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  )

figura_ul54

# [guardado antiguo] ggsave(
# [guardado antiguo]   "~/Desktop/Figuras_tesis/figura_lollipop_UL54.png",
# [guardado antiguo]   figura_ul54,
# [guardado antiguo]   width = 34, height = 18, units = "cm", dpi = 400
# [guardado antiguo] )

guardar_fig(figura_ul54, "Figura_46", ancho = 34, alto = 18, dpi = 400)

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

## ---- Figura 47. VND en UL56 ----
dominios_ul56 <- data.frame(
  start = c(4,   117, 280, 514, 671, 757),
  end   = c(20,  130, 310, 522, 680, 764),
  label = c("NLS", "endonucleasa", "Dedo de zinc", "Interacción pUL89", "Estab. 1", "Estab. 2"),
  fill  = c("cadetblue3", "cadetblue3", "cadetblue3", "#F0A500", "#F0A500", "#F0A500")
)

hotspots_ul56 <- data.frame(
  start = c(20,  230, 392),
  end   = c(30,  370, 398),
  label = c("",  "", ""),
  fill  = c("#B03060", "#B03060", "#B03060")
)

ul56_umbrales <- all_data_txp %>%
  filter(class == "NEW", Protein == "UL56") %>%
  group_by(Mutation, Position) %>%
  summarise(
    prev_max   = max(Prevalence, na.rm = TRUE),
    umbral_max = case_when(
      max(Prevalence) >= 20 ~ 20,
      max(Prevalence) >= 15 ~ 15,
      max(Prevalence) >= 10 ~ 10,
      TRUE                  ~ 5
    ),
    .groups = "drop"
  ) %>%
  mutate(Position = as.numeric(Position))

ul56 <- new_vars_dom %>%
  filter(Protein == "UL56") %>%
  left_join(ul56_umbrales %>% select(Mutation, umbral_max), by = "Mutation")

nivel <- c("5" = 0.6, "10" = 1.1, "15" = 1.6, "20" = 2.1)
paso  <- 0.4

df <- ul56 %>%
  mutate(
    y_nivel   = nivel[as.character(umbral_max)],
    en_lmv    = Position >= 231 & Position <= 412,
    y_point   = ifelse(en_lmv, -y_nivel, y_nivel),
    color     = ifelse(en_lmv, "#C0392B", "#2980B9"),
    label_mut = Mutation
  )

df_lab_azul_base <- df %>%
  filter(!is.na(label_mut), !en_lmv) %>%
  arrange(Position) %>%
  mutate(row_idx = row_number())

df_lab_azul <- df_lab_azul_base %>%
  mutate(
    y_lab_raw = 2.5 + (row_idx - 1) * paso,
    y_lab_raw = case_when(
      Position <= 30 ~ y_lab_raw + 0.8,
      TRUE           ~ y_lab_raw
    ),
    y_lab = case_when(
      Position >= 595 ~ y_lab_raw - 1.6,
      TRUE            ~ y_lab_raw
    )
  )

df_lab_rojo <- df %>%
  filter(!is.na(label_mut), en_lmv) %>%
  arrange(Position) %>%
  mutate(y_lab = -(0.9 + (row_number() - 1) * paso))

ymax_plot <- max(df_lab_azul$y_lab, na.rm = TRUE) + 0.5

niveles_ref <- data.frame(
  y     = c(0.6, 1.1, 1.6, 2.1, -0.6, -1.1, -1.6, -2.1),
  label = rep(c("\u22655%", "\u226510%", "\u226515%", "\u226520%"), 2),
  lado  = rep(c("arriba", "abajo"), each = 4)
)

p_ul56 <- ggplot() +
  geom_segment(data = niveles_ref,
               aes(x = -42, xend = 850, y = y, yend = y),
               color = "grey85", linewidth = 0.3, linetype = "dashed") +
  # Etiquetas prevalencia — más cerca del aa 1 (x = -55, límite -75)
  geom_text(data = niveles_ref %>% filter(lado == "arriba"),
            aes(x = -55, y = y, label = label),
            hjust = 1, size = 3.5, color = "grey50", fontface = "bold") +
  geom_text(data = niveles_ref %>% filter(lado == "abajo"),
            aes(x = -55, y = y, label = label),
            hjust = 1, size = 3.5, color = "grey50", fontface = "bold") +
  # Dominios
  geom_rect(data = dominios_ul56,
            aes(xmin = start, xmax = end, ymin = 0.0, ymax = 0.40, fill = fill),
            alpha = 0.9, color = "white", linewidth = 0.25) +
  geom_text(data = dominios_ul56 %>% mutate(mid = (start + end) / 2),
            aes(x = mid, y = 0.20, label = label),
            size = 3.5, color = "black", fontface = "bold",
            vjust = 0.5, hjust = 0.5, lineheight = 0.85) +
  # Hotspots
  geom_rect(data = hotspots_ul56,
            aes(xmin = start, xmax = end, ymin = -0.40, ymax = 0.0, fill = fill),
            alpha = 0.9, color = "white", linewidth = 0.25) +
  scale_fill_identity() +
  # Backbone aa 1–850
  geom_segment(aes(x = 1, xend = 850, y = 0, yend = 0),
               color = "grey10", linewidth = 1.0) +
  # Tapas extremos (palitos)
  geom_segment(aes(x = 1,   xend = 1,   y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  geom_segment(aes(x = 850, xend = 850, y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  # Stems
  geom_segment(data = df,
               aes(x = Position, xend = Position,
                   y = ifelse(en_lmv, -0.42, 0.42),
                   yend = y_point, color = color),
               linewidth = 0.45, alpha = 0.8) +
  # Puntos
  geom_point(data = df,
             aes(x = Position, y = y_point, color = color, size = prev_max),
             alpha = 0.88) +
  scale_color_identity() +
  scale_size_continuous(range = c(1.2, 4.5), guide = "none") +
  # Guías azules
  geom_segment(data = df_lab_azul,
               aes(x = Position, xend = Position,
                   y = y_point + 0.08, yend = y_lab - 0.08),
               color = "#2980B9", linewidth = 0.3, alpha = 0.4,
               linetype = "dotted") +
  # Guías rojas
  geom_segment(data = df_lab_rojo,
               aes(x = Position, xend = Position,
                   y = y_point - 0.08, yend = y_lab + 0.08),
               color = "#C0392B", linewidth = 0.3, alpha = 0.4,
               linetype = "dotted") +
  # Etiquetas azules
  geom_text(data = df_lab_azul,
            aes(x = Position, y = y_lab, label = label_mut),
            color = "#2980B9", size = 3.5, fontface = "bold",
            hjust = 0.5, vjust = 0) +
  # Etiquetas rojas
  geom_text(data = df_lab_rojo,
            aes(x = Position, y = y_lab, label = label_mut),
            color = "#C0392B", size = 3.5, fontface = "bold",
            hjust = 0.5, vjust = 1) +
  annotate("text", x = 1,   y = -0.28, label = "1",   size = 3.5, color = "grey10", vjust = 1) +
  annotate("text", x = 850, y = -0.28, label = "850", size = 3.5, color = "grey10", vjust = 1) +
  scale_x_continuous(
    limits = c(-100, 870),
    breaks = NULL,
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(limits = c(-2.5, ymax_plot)) +
  # Sin título de panel — el nombre "UL56" va en el pie de figura
  labs(title = NULL, x = "", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(t = 8, r = 12, b = 8, l = 35)
  )

leyenda_hotspots_ul56 <- ggplot(data.frame(
  x     = 1.5,
  y     = 1,
  label = "LMV-R",
  color = "#B03060"
)) +
  geom_point(aes(x, y, color = color), size = 4, shape = 15) +
  geom_text(aes(x + 0.15, y, label = label),
            hjust = 0, size = 3.5, color = "grey20") +
  scale_color_identity() +
  scale_x_continuous(limits = c(0.8, 4)) +
  scale_y_continuous(limits = c(0.5, 1.5)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(t = 0, b = 6, l = 35)
  )

figura_ul56 <- p_ul56 / leyenda_hotspots_ul56 +
  plot_layout(heights = c(10, 1))

figura_ul56

# [guardado antiguo] ggsave(
# [guardado antiguo]   "~/Desktop/Figuras_tesis/figura_lollipop_UL56.png",
# [guardado antiguo]   figura_ul56,
# [guardado antiguo]   width = 28, height = 15, units = "cm", dpi = 400,
# [guardado antiguo]   bg = "white"
# [guardado antiguo] )

guardar_fig(figura_ul56, "Figura_47", ancho = 28, alto = 15, dpi = 400)

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


