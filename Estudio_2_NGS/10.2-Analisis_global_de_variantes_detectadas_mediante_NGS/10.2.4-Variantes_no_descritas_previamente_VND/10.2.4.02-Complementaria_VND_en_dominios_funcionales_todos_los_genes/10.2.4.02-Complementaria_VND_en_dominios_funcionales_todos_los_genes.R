# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.4.02. Complementaria. VND en dominios funcionales (todos los genes)
# Extracto del script completo (ver README.md del estudio)

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

