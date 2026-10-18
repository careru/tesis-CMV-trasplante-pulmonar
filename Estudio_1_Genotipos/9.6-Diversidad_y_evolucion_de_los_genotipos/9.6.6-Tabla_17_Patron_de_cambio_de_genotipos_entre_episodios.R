# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.6.6. Tabla 17. Patrón de cambio de genotipos entre episodios
# Extracto del script completo (ver README.md del estudio)

## ---- Tabla 17. Patrón de cambio de genotipos entre episodios ----
# ── 8. Construir tabla ────────────────────────────────────────────────────────
tabla_inter <- estab_inter %>%
  left_join(cambios_inter, by = "NHC") %>%
  left_join(pacientes %>% select(NHC, serDR), by = "NHC") %>%
  mutate(
    estables = pmap_chr(
      list(gb1, gb2, gb3, gb4, gh1, gh2),
      function(gb1, gb2, gb3, gb4, gh1, gh2) {
        gt <- c(
          if (gb1 == "estable") "gB1",
          if (gb2 == "estable") "gB2",
          if (gb3 == "estable") "gB3",
          if (gb4 == "estable") "gB4",
          if (gh1 == "estable") "gH1",
          if (gh2 == "estable") "gH2"
        )
        if (length(gt) == 0) "—" else paste(gt, collapse = "/")
      }
    ),
    variables = pmap_chr(
      list(gb1, gb2, gb3, gb4, gh1, gh2),
      function(gb1, gb2, gb3, gb4, gh1, gh2) {
        gt <- c(
          if (gb1 == "variable") "gB1",
          if (gb2 == "variable") "gB2",
          if (gb3 == "variable") "gB3",
          if (gb4 == "variable") "gB4",
          if (gh1 == "variable") "gH1",
          if (gh2 == "variable") "gH2"
        )
        if (length(gt) == 0) "—" else paste(gt, collapse = "/")
      }
    )
  ) %>%
  arrange(desc(n_ep), desc(n_cambios)) %>%
  mutate(pac_num = row_number()) %>%
  select(pac_num, serDR, n_ep, estables, variables, n_cambios)

# ── 9. Flextable ──────────────────────────────────────────────────────────────
tabla_inter_ft <- tabla_inter %>%
  rename(
    "Paciente"              = pac_num,
    "Serología D/R"         = serDR,
    "N episodios"           = n_ep,
    "Patrón dominante"      = estables,
    "Genotipos fluctuantes" = variables,
    "N cambios"             = n_cambios
  )

fila_sin_estable <- which(tabla_inter_ft$`Patrón dominante` == "—")

ft_inter <- flextable(tabla_inter_ft) %>%
  width(j = "Paciente",              width = 0.7) %>%
  width(j = "Serología D/R",         width = 1.2) %>%
  width(j = "N episodios",           width = 0.9) %>%
  width(j = "Patrón dominante",      width = 2.0) %>%
  width(j = "Genotipos fluctuantes", width = 2.2) %>%
  width(j = "N cambios",             width = 0.9) %>%
  align(align = "center", part = "all") %>%
  align(j = c("Patrón dominante", "Genotipos fluctuantes"),
        align = "left", part = "body") %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header") %>%
  bg(bg = "white", part = "all") %>%
  bold(i = fila_sin_estable, part = "body") %>%
  color(i = fila_sin_estable, color = "firebrick3", part = "body") %>%
  border_remove() %>%
  hline_top(part = "header",
            border = fp_border(color = "black", width = 1.2)) %>%
  hline(part = "header",
        border = fp_border(color = "black", width = 0.8)) %>%
  hline_bottom(part = "body",
               border = fp_border(color = "black", width = 1.2)) %>%
  padding(padding.top = 2, padding.bottom = 2,
          padding.left = 4, padding.right = 4, part = "all") %>%
  add_footer_lines(values = c(
    "Patrón dominante: genotipos presentes en todos los episodios del paciente.",
    "Genotipos fluctuantes: genotipos con aparición o desaparición entre episodios.",
    "El paciente en rojo corresponde al único caso sin ningún genotipo estable entre episodios.",
    "N cambios: número de transiciones de patrón genotípico entre episodios consecutivos."
  )) %>%
  fontsize(size = 8, part = "footer") %>%
  color(color = "grey30", part = "footer") %>%
  set_caption(
    caption = "Tabla R9. Patrón de cambio de genotipos gB/gH en pacientes con infecciones múltiples que presentaron variabilidad interepisodio."
  ) %>%
  set_table_properties(layout = "autofit")

ft_inter

guardar_tabla(ft_inter, "Tabla_17")

