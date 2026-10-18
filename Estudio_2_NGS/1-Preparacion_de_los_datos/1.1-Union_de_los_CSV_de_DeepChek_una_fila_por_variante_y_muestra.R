# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 1.1. Unión de los CSV de DeepChek (una fila por variante y muestra)
# Extracto del script completo (ver README.md del estudio)

## ---- Unión de los CSV de DeepChek (una fila por variante y muestra) ----
# ── Configuración ─────────────────────────────────────────────────────────────
path        <- dir_aa
output_file <- file.path(dir_salida, "1_out_mutaciones.xlsx")
extra_mra   <- c("L545F")

# ── Carga y fusión de archivos ────────────────────────────────────────────────
files <- list.files(path, pattern = "\\.csv$", full.names = TRUE)
cat("Archivos encontrados:", length(files), "\n")

all_data <- map(files, ~ {
  df <- read_delim(.x, delim = ";", col_types = cols(.default = "c")) %>%
    rename(Npet = 1) %>%
    mutate(Npet = trimws(Npet))
  
  names(df) <- str_replace(
    names(df),
    "Is mutation interest \\(CMV variant database .*\\)",
    "Resistencia"
  )
  names(df) <- str_replace(names(df), "Passed Expert System", "Filtro")
  
  return(df)
}) %>%
  bind_rows()

cat("Muestras únicas (pre-filtro):", n_distinct(all_data$Npet), "\n")

# ── Selección de columnas ─────────────────────────────────────────────────────
all_data <- all_data %>%
  dplyr::select(
    Npet, Protein, Position, Mutation, Prevalence, Resistencia, Filtro,
    matches("Mutation count FW|Mutation\\.count\\.FW"),
    matches("Mutation count RV|Mutation\\.count\\.RV"),
    matches("Position coverage FW|Position\\.coverage\\.FW"),
    matches("Position coverage RV|Position\\.coverage\\.RV"),
    matches("^Q-Score$|^Q\\.Score$")
  )

# ── Renombrar a nombres cortos ────────────────────────────────────────────────
old_names <- names(all_data)
names(all_data)[str_detect(old_names, "Mutation.*count.*FW|Mutation\\.count\\.FW")]       <- "mut_fw"
names(all_data)[str_detect(old_names, "Mutation.*count.*RV|Mutation\\.count\\.RV")]       <- "mut_rv"
names(all_data)[str_detect(old_names, "Position.*coverage.*FW|Position\\.coverage\\.FW")] <- "cov_fw"
names(all_data)[str_detect(old_names, "Position.*coverage.*RV|Position\\.coverage\\.RV")] <- "cov_rv"
names(all_data)[str_detect(old_names, "Q.Score|Q-Score")]                                 <- "q_score"

# ── Convertir a numérico y calcular métricas derivadas ───────────────────────
all_data <- all_data %>%
  mutate(
    across(c(mut_fw, mut_rv, cov_fw, cov_rv, q_score), as.numeric),
    cobertura_total = cov_fw + cov_rv,
    balance_ok      = mut_fw > 0 & mut_rv > 0
  )

cat("Columnas de calidad disponibles:",
    paste(c("mut_fw","mut_rv","cov_fw","cov_rv","q_score",
            "cobertura_total","balance_ok"), collapse = ", "), "\n")

# ── Filtros ───────────────────────────────────────────────────────────────────
n_dash <- sum(all_data$Mutation == "-", na.rm = TRUE)
cat("Filas con Mutation == '-' eliminadas:", n_dash, "\n")

all_data <- all_data %>%
  filter(Mutation != "-") %>%
  filter(Filtro == "Yes")

cat("Muestras únicas (post-filtro Expert System):", n_distinct(all_data$Npet), "\n")
cat("Variantes totales (post-filtro Expert System):", nrow(all_data), "\n")

# ── Corrección de Resistencia ─────────────────────────────────────────────────
all_data <- all_data %>%
  mutate(Resistencia = if_else(str_trim(Mutation) %in% extra_mra, "Yes", Resistencia))
cat("Variantes marcadas como MRA adicional:",
    sum(str_trim(all_data$Mutation) %in% extra_mra), "\n")

nas <- sum(is.na(all_data$Resistencia))
cat("NAs en Resistencia:", nas, "\n")
if (nas > 0) {
  all_data <- all_data %>% mutate(Resistencia = replace_na(Resistencia, "No"))
  cat("NAs corregidos a 'No'. Comprobación:", sum(is.na(all_data$Resistencia)), "\n")
}

# ── Exportar ──────────────────────────────────────────────────────────────────
write.xlsx(all_data, file = output_file)
cat("Archivo exportado:", output_file, "\n")

