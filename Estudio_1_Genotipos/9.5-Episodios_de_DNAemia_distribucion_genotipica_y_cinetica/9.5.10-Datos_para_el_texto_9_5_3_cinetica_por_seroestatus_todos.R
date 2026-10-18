# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.5.10. Datos para el texto (9.5.3): cinética por seroestatus, todos los episodios y primer episodio
# Extracto del script completo (ver README.md del estudio)

## ---- Datos para el texto (9.5.3): cinética por seroestatus, todos los episodios y primer episodio ----
# ── 1. Datos ──────────────────────────────────────────────────────────────────
df <- dnaemias

vars <- c("cvini", "cvpico", "diasdur_vir", "AUCCV")
var_labels <- c(
  cvini       = "CV inicial",
  cvpico      = "CV pico",
  diasdur_vir = "Aclaramiento",
  AUCCV       = "AUC viral"
)

# ── 2. Función de cálculo ─────────────────────────────────────────────────────
calcular_p <- function(mixvar, mix_label) {
  map_dfr(vars, function(v) {
    d0 <- df %>% filter(.data[[mixvar]] == 0) %>% pull(.data[[v]])
    d1 <- df %>% filter(.data[[mixvar]] == 1) %>% pull(.data[[v]])
    p_global <- tryCatch(
      wilcox.test(d0, d1, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
    grupos_res <- map_dfr(c("D-/R+","D+/R+","D+/R-"), function(grupo) {
      dg <- df %>% filter(serDR == grupo)
      d0g <- dg %>% filter(.data[[mixvar]] == 0) %>% pull(.data[[v]])
      d1g <- dg %>% filter(.data[[mixvar]] == 1) %>% pull(.data[[v]])
      p <- if (length(na.omit(d0g)) > 0 & length(na.omit(d1g)) > 0)
        tryCatch(wilcox.test(d0g, d1g, exact = FALSE)$p.value,
                 error = function(e) NA_real_)
      else NA_real_
      tibble(serDR = grupo, p = p)
    }) %>%
      pivot_wider(names_from = serDR, values_from = p,
                  names_prefix = "p_")
    bind_cols(
      tibble(Glucoproteína = mix_label, Variable = var_labels[v], p_global),
      grupos_res
    )
  })
}

# ── 3. Calcular y combinar ────────────────────────────────────────────────────
tabla_raw <- bind_rows(
  calcular_p("gbmix",   "gB"),
  calcular_p("ghmix",   "gH"),
  calcular_p("gbghmix", "gB+gH")
)

# ── 4. Formatear p-valores ────────────────────────────────────────────────────
fmt_p <- function(p) {
  if (is.na(p))    return("")
  if (p < 0.001)   return("<0,001")
  txt <- sub("\\.", ",", sprintf("%.3f", p))
  if (p < 0.05) paste0(txt, "*") else txt
}

tabla_fmt <- tabla_raw %>%
  mutate(across(where(is.numeric), ~ sapply(., fmt_p)))

# ── 5. Añadir filas cabecera por glucoproteína ────────────────────────────────
cab <- function(txt) tibble(
  Glucoproteína = txt,
  Variable      = NA_character_,
  p_global      = NA_character_,
  `p_D-/R+`     = NA_character_,
  `p_D+/R+`     = NA_character_,
  `p_D+/R-`     = NA_character_
)

tabla_con_cab <- bind_rows(
  cab("Glucoproteína gB"),
  tabla_fmt %>% filter(Glucoproteína == "gB"),
  cab("Glucoproteína gH"),
  tabla_fmt %>% filter(Glucoproteína == "gH"),
  cab("Glucoproteínas gB+gH"),
  tabla_fmt %>% filter(Glucoproteína == "gB+gH")
)

# Recalcular índices sobre la tabla con cabeceras
idx_cab <- which(is.na(tabla_con_cab$Variable))

# Recalcular índices sig sobre tabla_raw reordenada
tabla_raw_ord <- bind_rows(
  tabla_raw %>% filter(Glucoproteína == "gB"),
  tabla_raw %>% filter(Glucoproteína == "gH"),
  tabla_raw %>% filter(Glucoproteína == "gB+gH")
)

# Posiciones en tabla_con_cab (desplazadas por las filas cabecera)
offset <- function(idx_raw, idx_cab_vec) {
  sapply(idx_raw, function(i) {
    i + sum(idx_cab_vec <= (i + sum(idx_cab_vec < i)))
  })
}

idx_sig_global <- which(!is.na(tabla_raw_ord$p_global) & tabla_raw_ord$p_global < 0.05)
idx_sig_DmR    <- which(!is.na(tabla_raw_ord$`p_D-/R+`) & tabla_raw_ord$`p_D-/R+` < 0.05)
idx_sig_DpR    <- which(!is.na(tabla_raw_ord$`p_D+/R+`) & tabla_raw_ord$`p_D+/R+` < 0.05)
idx_sig_DpRm   <- which(!is.na(tabla_raw_ord$`p_D+/R-`) & tabla_raw_ord$`p_D+/R-` < 0.05)

# Convertir índices raw a índices en tabla con cabeceras
to_table_idx <- function(raw_idx) {
  # Cada bloque de 4 filas tiene 1 cabecera antes
  raw_idx + ceiling(raw_idx / 4)
}

ft <- tabla_con_cab %>%
  select(-Glucoproteína) %>%
  flextable() %>%
  set_header_labels(
    Variable  = "Variable",
    p_global  = "Global",
    `p_D-/R+` = "D-/R+",
    `p_D+/R+` = "D+/R+",
    `p_D+/R-` = "D+/R-"
  ) %>%
  add_header_row(
    values     = c("", "p (único vs mixto)"),
    colwidths  = c(1, 4)
  ) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  bold(i = idx_cab, j = 1, part = "body") %>%
  bg(i = idx_cab, bg = "#F2F2F2", part = "body") %>%
  padding(i = setdiff(seq_len(nrow(tabla_con_cab)), idx_cab),
          j = 1, padding.left = 14, part = "body") %>%
  bold(i = to_table_idx(idx_sig_global), j = "p_global",  part = "body") %>%
  bold(i = to_table_idx(idx_sig_DmR),    j = "p_D-/R+",  part = "body") %>%
  bold(i = to_table_idx(idx_sig_DpR),    j = "p_D+/R+",  part = "body") %>%
  bold(i = to_table_idx(idx_sig_DpRm),   j = "p_D+/R-",  part = "body") %>%
  align(j = "Variable", align = "left",   part = "all") %>%
  align(j = c("p_global","p_D-/R+","p_D+/R+","p_D+/R-"),
        align = "center", part = "all") %>%
  autofit() %>%
  add_footer_lines(values = c(
    "Comparaciones mediante la prueba U de Mann-Whitney (único vs mixto).",
    "Los p-valores se muestran para el análisis global y estratificado por perfil serológico D/R.",
    "*p < 0,05."
  )) %>%
  fontsize(part = "footer", size = 8) %>%
  align(part = "footer", align = "left")

ft

# ── 1. Datos — solo primer episodio por paciente ──────────────────────────────
df <- dnaemias %>% filter(primer_mil == 1)

vars <- c("cvini", "cvpico", "diasdur_vir", "AUCCV")
var_labels <- c(
  cvini       = "CV inicial",
  cvpico      = "CV pico",
  diasdur_vir = "Aclaramiento",
  AUCCV       = "AUC viral"
)

# ── 2. Función de cálculo ─────────────────────────────────────────────────────
calcular_p <- function(mixvar, mix_label) {
  map_dfr(vars, function(v) {
    d0 <- df %>% filter(.data[[mixvar]] == 0) %>% pull(.data[[v]])
    d1 <- df %>% filter(.data[[mixvar]] == 1) %>% pull(.data[[v]])
    p_global <- tryCatch(
      wilcox.test(d0, d1, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
    grupos_res <- map_dfr(c("D-/R+","D+/R+","D+/R-"), function(grupo) {
      dg <- df %>% filter(serDR == grupo)
      d0g <- dg %>% filter(.data[[mixvar]] == 0) %>% pull(.data[[v]])
      d1g <- dg %>% filter(.data[[mixvar]] == 1) %>% pull(.data[[v]])
      p <- if (length(na.omit(d0g)) > 0 & length(na.omit(d1g)) > 0)
        tryCatch(wilcox.test(d0g, d1g, exact = FALSE)$p.value,
                 error = function(e) NA_real_)
      else NA_real_
      tibble(serDR = grupo, p = p)
    }) %>%
      pivot_wider(names_from = serDR, values_from = p,
                  names_prefix = "p_")
    bind_cols(
      tibble(Glucoproteína = mix_label, Variable = var_labels[v], p_global),
      grupos_res
    )
  })
}

# ── 3. Calcular y combinar ────────────────────────────────────────────────────
tabla_raw <- bind_rows(
  calcular_p("gbmix",   "gB"),
  calcular_p("ghmix",   "gH"),
  calcular_p("gbghmix", "gB+gH")
)

# ── 4. Formatear p-valores ────────────────────────────────────────────────────
fmt_p <- function(p) {
  if (is.na(p))    return("")
  if (p < 0.001)   return("<0,001")
  txt <- sub("\\.", ",", sprintf("%.3f", p))
  if (p < 0.05) paste0(txt, "*") else txt
}

tabla_fmt <- tabla_raw %>%
  mutate(across(where(is.numeric), ~ sapply(., fmt_p)))

# ── 5. Añadir filas cabecera por glucoproteína ────────────────────────────────
cab <- function(txt) tibble(
  Glucoproteína = txt,
  Variable      = NA_character_,
  p_global      = NA_character_,
  `p_D-/R+`     = NA_character_,
  `p_D+/R+`     = NA_character_,
  `p_D+/R-`     = NA_character_
)

tabla_con_cab <- bind_rows(
  cab("Glucoproteína gB"),
  tabla_fmt %>% filter(Glucoproteína == "gB"),
  cab("Glucoproteína gH"),
  tabla_fmt %>% filter(Glucoproteína == "gH"),
  cab("Glucoproteínas gB+gH"),
  tabla_fmt %>% filter(Glucoproteína == "gB+gH")
)

# Recalcular índices sobre la tabla con cabeceras
idx_cab <- which(is.na(tabla_con_cab$Variable))

# Recalcular índices sig sobre tabla_raw reordenada
tabla_raw_ord <- bind_rows(
  tabla_raw %>% filter(Glucoproteína == "gB"),
  tabla_raw %>% filter(Glucoproteína == "gH"),
  tabla_raw %>% filter(Glucoproteína == "gB+gH")
)

# Posiciones en tabla_con_cab (desplazadas por las filas cabecera)
offset <- function(idx_raw, idx_cab_vec) {
  sapply(idx_raw, function(i) {
    i + sum(idx_cab_vec <= (i + sum(idx_cab_vec < i)))
  })
}

idx_sig_global <- which(!is.na(tabla_raw_ord$p_global) & tabla_raw_ord$p_global < 0.05)
idx_sig_DmR    <- which(!is.na(tabla_raw_ord$`p_D-/R+`) & tabla_raw_ord$`p_D-/R+` < 0.05)
idx_sig_DpR    <- which(!is.na(tabla_raw_ord$`p_D+/R+`) & tabla_raw_ord$`p_D+/R+` < 0.05)
idx_sig_DpRm   <- which(!is.na(tabla_raw_ord$`p_D+/R-`) & tabla_raw_ord$`p_D+/R-` < 0.05)

# Convertir índices raw a índices en tabla con cabeceras
to_table_idx <- function(raw_idx) {
  # Cada bloque de 4 filas tiene 1 cabecera antes
  raw_idx + ceiling(raw_idx / 4)
}

ft <- tabla_con_cab %>%
  select(-Glucoproteína) %>%
  flextable() %>%
  set_header_labels(
    Variable  = "Variable",
    p_global  = "Global",
    `p_D-/R+` = "D-/R+",
    `p_D+/R+` = "D+/R+",
    `p_D+/R-` = "D+/R-"
  ) %>%
  add_header_row(
    values     = c("", "PRIMER EPISODIO p (único vs mixto)"),
    colwidths  = c(1, 4)
  ) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  bold(i = idx_cab, j = 1, part = "body") %>%
  bg(i = idx_cab, bg = "#F2F2F2", part = "body") %>%
  padding(i = setdiff(seq_len(nrow(tabla_con_cab)), idx_cab),
          j = 1, padding.left = 14, part = "body") %>%
  bold(i = to_table_idx(idx_sig_global), j = "p_global",  part = "body") %>%
  bold(i = to_table_idx(idx_sig_DmR),    j = "p_D-/R+",  part = "body") %>%
  bold(i = to_table_idx(idx_sig_DpR),    j = "p_D+/R+",  part = "body") %>%
  bold(i = to_table_idx(idx_sig_DpRm),   j = "p_D+/R-",  part = "body") %>%
  align(j = "Variable", align = "left",   part = "all") %>%
  align(j = c("p_global","p_D-/R+","p_D+/R+","p_D+/R-"),
        align = "center", part = "all") %>%
  autofit() %>%
  add_footer_lines(values = c(
    "Comparaciones mediante la prueba U de Mann-Whitney (único vs mixto).",
    "Los p-valores se muestran para el análisis global y estratificado por perfil serológico D/R.",
    "*p < 0,05."
  )) %>%
  fontsize(part = "footer", size = 8) %>%
  align(part = "footer", align = "left")

ft

