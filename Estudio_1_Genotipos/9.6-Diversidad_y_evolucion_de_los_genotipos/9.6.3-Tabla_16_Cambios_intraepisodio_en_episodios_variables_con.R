# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.6.3. Tabla 16. Cambios intraepisodio en episodios variables con genotipos mixtos gB/gH
# Extracto del script completo (ver README.md del estudio)

## ---- Tabla 16. Cambios intraepisodio en episodios variables con genotipos mixtos gB/gH ----
muestras_mv <- muestras_ep %>%
  semi_join(mixtos_var, by = c("NHC", "Fini")) %>%
  mutate(patron = paste(GB1, GB2, GB3, GB4, GH1, GH2, sep = "-")) %>%
  arrange(NHC, Fini, Freg)

# Número de cambios por episodio
cambios_ep <- muestras_mv %>%
  group_by(NHC, Fini) %>%
  arrange(Freg, .by_group = TRUE) %>%
  mutate(cambio = patron != lag(patron)) %>%
  summarise(
    n_muestras = n(),
    n_cambios  = sum(cambio, na.rm = TRUE),
    .groups    = "drop"
  )

cat("── Número de cambios por episodio ──\n")
cat(sprintf("Rango: %d – %d\n", min(cambios_ep$n_cambios), max(cambios_ep$n_cambios)))
cat(sprintf("Mediana: %.0f (RIQ %.0f–%.0f)\n",
            median(cambios_ep$n_cambios),
            quantile(cambios_ep$n_cambios, 0.25),
            quantile(cambios_ep$n_cambios, 0.75)))

# Estabilidad por genotipo
genotipos <- c("GB1","GB2","GB3","GB4","GH1","GH2")
nombres   <- c("gB1","gB2","gB3","gB4","gH1","gH2")

estabilidad_gt <- muestras_mv %>%
  group_by(NHC, Fini) %>%
  summarise(
    across(all_of(genotipos), ~ {
      presente <- sum(. == 1, na.rm = TRUE)
      n_tot    <- n()
      if (presente == 0) "ausente"
      else if (presente == n_tot) "estable"
      else "variable"
    }),
    .groups = "drop"
  )

cat("\n── Frecuencia genotipos estables ──\n")
for (i in seq_along(genotipos)) {
  g <- genotipos[i]
  n_est <- sum(estabilidad_gt[[g]] == "estable")
  cat(sprintf("%s estable: %d/35 (%.1f%%)\n", nombres[i], n_est, 100*n_est/35))
}

cat("\n── Frecuencia genotipos variables ──\n")
for (i in seq_along(genotipos)) {
  g <- genotipos[i]
  n_var <- sum(estabilidad_gt[[g]] == "variable")
  cat(sprintf("%s variable: %d/35 (%.1f%%)\n", nombres[i], n_var, 100*n_var/35))
}

estabilidad_gt <- estabilidad_gt %>%
  mutate(tiene_estable = rowSums(across(all_of(genotipos), ~ . == "estable")) > 0)

cat(sprintf("\n── Episodios con ≥1 genotipo estable: %d/35 (%.1f%%)\n",
            sum(estabilidad_gt$tiene_estable),
            100 * mean(estabilidad_gt$tiene_estable)))

cat(sprintf("── Episodios con cambio completo: %d\n",
            sum(!estabilidad_gt$tiene_estable)))

tabla_mv <- estabilidad_gt %>%
  select(-any_of("n_muestras")) %>%
  left_join(
    cambios_ep %>% select(NHC, Fini, n_cambios, n_muestras),
    by = c("NHC", "Fini")
  ) %>%
  left_join(pacientes %>% select(NHC, serDR), by = "NHC") %>%
  mutate(
    estables = pmap_chr(
      list(GB1, GB2, GB3, GB4, GH1, GH2),
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
      list(GB1, GB2, GB3, GB4, GH1, GH2),
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
  arrange(NHC, Fini) %>%
  mutate(
    pac_num   = as.integer(factor(NHC, levels = unique(NHC))),
    ep_num    = row_number(),
    pac_label = paste0(pac_num, " ", serDR)
  ) %>%
  select(ep_num, pac_label, n_muestras, estables, variables, n_cambios)

tabla_ft <- tabla_mv %>%
  rename(
    "Episodio"              = ep_num,
    "Paciente"              = pac_label,
    "N muestras"            = n_muestras,
    "Patrón dominante"      = estables,
    "Genotipos fluctuantes" = variables,
    "N cambios"             = n_cambios
  )

ft <- flextable(tabla_ft) %>%
  width(j = "Episodio",              width = 0.7) %>%
  width(j = "Paciente",              width = 1.6) %>%
  width(j = "N muestras",            width = 0.9) %>%
  width(j = "Patrón dominante",      width = 2.0) %>%
  width(j = "Genotipos fluctuantes", width = 2.0) %>%
  width(j = "N cambios",             width = 0.9) %>%
  align(align = "center", part = "all") %>%
  align(j = c("Patrón dominante", "Genotipos fluctuantes"),
        align = "left", part = "body") %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header") %>%
  bg(bg = "white", part = "all") %>%
  bold(i = 1, part = "body") %>%
  color(i = 1, color = "firebrick3", part = "body") %>%
  border_remove() %>%
  hline_top(part = "header",
            border = fp_border(color = "black", width = 1.2)) %>%
  hline(part = "header",
        border = fp_border(color = "black", width = 0.8)) %>%
  hline_bottom(part = "body",
               border = fp_border(color = "black", width = 1.2)) %>%
  hline(
    i = which(diff(as.integer(factor(tabla_ft$Paciente))) != 0),
    border = fp_border(color = "grey50", width = 0.5),
    part = "body"
  ) %>%
  padding(padding.top = 2, padding.bottom = 2,
          padding.left = 4, padding.right = 4, part = "all") %>%
  add_footer_lines(values = c(
    "Patrón dominante: genotipos presentes en todas las muestras del episodio.",
    "Genotipos fluctuantes: genotipos con aparición o desaparición entre muestras consecutivas.",
    "El episodio en rojo corresponde al único caso con cambio completo del patrón dominante.",
    "N cambios: número de transiciones de patrón genotípico entre muestras consecutivas."
  )) %>%
  fontsize(size = 8, part = "footer") %>%
  color(color = "grey30", part = "footer") %>%
  set_caption(caption = "Tabla R8. Episodios con variabilidad intraepisodio en genotipos mixtos gB/gH.") %>%
  set_table_properties(layout = "autofit")

ft

guardar_tabla(ft, "Tabla_16")

