# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.5.06. Tabla 12. Cinética viral según genotipo individual (genotipos únicos de gB y gH)
# Extracto del script completo (ver README.md del estudio)

## ---- Tabla 12. Cinética viral según genotipo individual (genotipos únicos de gB y gH) ----
# ── 1. Función genérica ───────────────────────────────────────────────────────
resumen_bloque <- function(data, geno_var, geno_levels, geno_labels,
                           vars, var_labels, tipo_comparacion) {
  map_dfr(seq_along(vars), function(i) {
    v   <- vars[i]
    lab <- var_labels[i]
    
    celdas <- map_chr(geno_levels, function(g) {
      x <- data %>% filter(.data[[geno_var]] == g) %>%
        pull(.data[[v]]) %>% na.omit()
      if (length(x) == 0) return("–")
      sprintf("%s (%s–%s)\nn = %d",
              formatC(median(x),          format = "f", digits = 0),
              formatC(quantile(x, 0.25),  format = "f", digits = 0),
              formatC(quantile(x, 0.75),  format = "f", digits = 0),
              length(x))
    })
    
    grupos <- map(geno_levels, function(g)
      data %>% filter(.data[[geno_var]] == g) %>%
        pull(.data[[v]]) %>% na.omit())
    
    p_val <- tryCatch({
      if (tipo_comparacion == "kruskal") kruskal.test(grupos)$p.value
      else wilcox.test(grupos[[1]], grupos[[2]], exact = FALSE)$p.value
    }, error = function(e) NA_real_)
    
    p_fmt <- if (is.na(p_val)) "" else
      if (p_val < 0.001) "<0,001" else
        sub("\\.", ",", sprintf("%.3f", p_val))
    
    row        <- as_tibble_row(setNames(celdas, geno_labels))
    row$Parámetro <- lab
    row$p         <- p_fmt
    row$p_raw     <- p_val
    row
  })
}

# ── 2. Preparar datos ─────────────────────────────────────────────────────────
dnaemias <- dnaemias %>%
  mutate(
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini)
  )

vars     <- c("cvini", "cvpico", "aclaramiento_dias", "AUCCV")
var_labs <- c("CV inicial (UI/mL)", "CV pico (UI/mL)",
              "Aclaramiento (días)", "AUC viral (log10(UI/mL)·día)")

# ── 3. Bloque gB ──────────────────────────────────────────────────────────────
df_gb <- dnaemias %>%
  filter(gbsum == 1) %>%
  mutate(gb_tipo = case_when(
    gb1 == 1 ~ "gB1", gb2 == 1 ~ "gB2",
    gb3 == 1 ~ "gB3", gb4 == 1 ~ "gB4"
  ))

bloque_gb <- resumen_bloque(
  data             = df_gb,
  geno_var         = "gb_tipo",
  geno_levels      = c("gB1","gB2","gB3","gB4"),
  geno_labels      = c("gB1","gB2","gB3","gB4"),
  vars             = vars,
  var_labels       = var_labs,
  tipo_comparacion = "kruskal"
)

# ── 4. Bloque gH ──────────────────────────────────────────────────────────────
df_gh <- dnaemias %>%
  filter(ghsum == 1) %>%
  mutate(gh_tipo = case_when(
    gh1 == 1 ~ "gH1", gh2 == 1 ~ "gH2"
  ))

bloque_gh_raw <- resumen_bloque(
  data             = df_gh,
  geno_var         = "gh_tipo",
  geno_levels      = c("gH1","gH2"),
  geno_labels      = c("gH1","gH2"),
  vars             = vars,
  var_labels       = var_labs,
  tipo_comparacion = "wilcox"
)

# Añadir columnas gB3/gB4 vacías y renombrar gH1/gH2 → gB1/gB2
# (el flextable usará los nombres de columna, las cabeceras se personalizan)
bloque_gh <- bloque_gh_raw %>%
  rename(gB1 = gH1, gB2 = gH2) %>%
  mutate(gB3 = NA_character_, gB4 = NA_character_) %>%
  select(Parámetro, gB1, gB2, gB3, gB4, p, p_raw)

# ── 5. Filas cabecera ─────────────────────────────────────────────────────────
cab <- function(txt) tibble(
  Parámetro = txt,
  gB1 = NA_character_, gB2 = NA_character_,
  gB3 = NA_character_, gB4 = NA_character_,
  p   = NA_character_,  p_raw = NA_real_
)

# ── 6. Tabla completa ─────────────────────────────────────────────────────────
tabla_final <- bind_rows(
  cab("Glucoproteína gB"),
  bloque_gb,
  cab("Glucoproteína gH"),
  bloque_gh
)

# ── 7. Índices para formato ───────────────────────────────────────────────────
idx_cab <- which(tabla_final$Parámetro %in%
                   c("Glucoproteína gB", "Glucoproteína gH"))
idx_sig <- which(!is.na(tabla_final$p_raw) & tabla_final$p_raw < 0.05)
idx_gh  <- which(tabla_final$Parámetro == "Glucoproteína gH") +
  seq_len(length(vars))   # filas de datos de gH

# ── 8. Flextable ──────────────────────────────────────────────────────────────
ft <- tabla_final %>%
  select(-p_raw) %>%
  flextable() %>%
  # cabecera con nombres correctos por sección
  set_header_labels(
    Parámetro = "Parámetro",
    gB1 = "gB1 / gH1",
    gB2 = "gB2 / gH2",
    gB3 = "gB3",
    gB4 = "gB4",
    p   = "p"
  ) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  # cabeceras de sección
  bold(i = idx_cab, j = 1, part = "body") %>%
  bg(i = idx_cab,   bg = "#F2F2F2", part = "body") %>%
  # p significativas en negrita
  bold(i = idx_sig, j = "p", part = "body") %>%
  # alineación
  align(j = "Parámetro", align = "left",   part = "all") %>%
  align(j = c("gB1","gB2","gB3","gB4","p"),
        align = "center", part = "all") %>%
  # sangría en filas de datos
  padding(i = setdiff(seq_len(nrow(tabla_final)), idx_cab),
          j = 1, padding.left = 14, part = "body") %>%
  # gris para celdas vacías (gB3/gB4 en filas de gH)
  bg(i = idx_gh, j = c("gB3","gB4"), bg = "#F8F8F8", part = "body") %>%
  autofit() %>%
  add_footer_lines(values = c(
    "Datos presentados como mediana (RIC). CV: carga viral; AUC: área bajo la curva viral.",
    "Solo se incluyen episodios con genotipo único de cada glucoproteína.",
    "Glucoproteína gB: comparación entre gB1, gB2, gB3 y gB4 (Kruskal-Wallis).",
    "Glucoproteína gH: comparación entre gH1 y gH2 (U de Mann-Whitney).",
    "Las columnas gB3 y gB4 no aplican a gH (–). Los valores de p < 0,05 en negrita."
  )) %>%
  fontsize(part = "footer", size = 8) %>%
  align(part = "footer", align = "left")

ft

guardar_tabla(ft, "Tabla_12")

