# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.5.09. Tabla 14. Cinética viral según presencia de genotipos mixtos de gB, gH y gB/gH
# Extracto del script completo (ver README.md del estudio)

## ---- Tabla 14. Cinética viral según presencia de genotipos mixtos de gB, gH y gB/gH ----
# ── 1. Calcular aclaramiento si no existe ─────────────────────────────────────
dnaemias <- dnaemias %>%
  mutate(
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini)
  )

# ── 2. Funciones auxiliares ───────────────────────────────────────────────────
fmt_med <- function(x) {
  n   <- sum(!is.na(x))
  med <- median(x, na.rm = TRUE)
  q1  <- quantile(x, 0.25, na.rm = TRUE)
  q3  <- quantile(x, 0.75, na.rm = TRUE)
  sprintf("%s (%s\u2013%s)\n(n=%d)",
          formatC(med, format="f", digits=0),
          formatC(q1,  format="f", digits=0),
          formatC(q3,  format="f", digits=0), n)
}

fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0,001*")
  txt <- sprintf("%.3f", p)
  txt <- sub("\\.", ",", txt)
  if (p < 0.05) paste0(txt, "*") else txt
}

wilcox_p <- function(data, var, grupo) {
  d <- data %>% filter(!is.na(.data[[var]]), !is.na(.data[[grupo]]))
  if (nrow(d) < 4) return(NA_real_)
  tryCatch(
    wilcox.test(as.formula(paste(var, "~", grupo)), data = d)$p.value,
    error = function(e) NA_real_
  )
}

# ── 3. Función que construye 4 filas para un grupo genotípico ─────────────────
bloque <- function(data, var_mix) {
  d0 <- data %>% filter(.data[[var_mix]] == 0)
  d1 <- data %>% filter(.data[[var_mix]] == 1)
  
  vars  <- c("cvini", "cvpico", "aclaramiento_dias", "AUCCV")
  labs  <- c("CV inicial (UI/mL)", "CV pico (UI/mL)",
             "Aclaramiento (días)", "AUC viral")
  
  tibble(
    Parámetro = labs,
    `Genotipo único` = sapply(vars, function(v) fmt_med(d0[[v]])),
    `Genotipo mixto`  = sapply(vars, function(v) fmt_med(d1[[v]])),
    p = sapply(vars, function(v) fmt_p(wilcox_p(data, v, var_mix)))
  )
}

# ── 4. Construir tabla completa ───────────────────────────────────────────────
cabecera <- function(txt) tibble(
  Parámetro        = txt,
  `Genotipo único` = NA_character_,
  `Genotipo mixto` = NA_character_,
  p                = NA_character_
)

tabla_cin <- bind_rows(
  cabecera("Glucoproteína gB"),
  bloque(dnaemias, "gbmix"),
  cabecera("Glucoproteína gH"),
  bloque(dnaemias, "ghmix"),
  cabecera("Glucoproteínas gB y gH"),
  bloque(dnaemias, "gbghmix")
)

# ── 5. Flextable ──────────────────────────────────────────────────────────────
# Índices de filas cabecera
idx_cab <- which(tabla_cin$Parámetro %in%
                   c("Glucoproteína gB", "Glucoproteína gH",
                     "Glucoproteínas gB y gH"))

ft_cin <- flextable(tabla_cin) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  bold(i = idx_cab, j = 1, part = "body") %>%
  bg(i = idx_cab, bg = "#F2F2F2", part = "body") %>%
  align(j = "Parámetro", align = "left",   part = "all") %>%
  align(j = c("Genotipo único", "Genotipo mixto", "p"),
        align = "center", part = "all") %>%
  padding(i = setdiff(seq_len(nrow(tabla_cin)), idx_cab),
          j = 1, padding.left = 14, part = "body") %>%
  autofit() %>%
  add_footer_lines(values = c(
    "Los datos se presentan como mediana (RIC).",
    "CV: carga viral; AUC: área bajo la curva.",
    "Comparaciones mediante la prueba U de Mann-Whitney.",
    "*p < 0,05."
  )) %>%
  fontsize(part = "footer", size = 8) %>%
  align(part = "footer", align = "left")

ft_cin

guardar_tabla(ft_cin, "Tabla_14")

