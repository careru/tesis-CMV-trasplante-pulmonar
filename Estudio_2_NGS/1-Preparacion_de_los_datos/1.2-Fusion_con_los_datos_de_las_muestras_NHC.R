# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 1.2. Fusión con los datos de las muestras (NHC)
# Extracto del script completo (ver README.md del estudio)

## ---- Fusión con los datos de las muestras (NHC) ----
# ── Exclusiones ───────────────────────────────────────────────────────────────
npet_excluidos <- unique(c(
  "102564799", "102344108", "102498654", "102500899", "213553115",
  "213645499", "334497189", "666063012", "102578305", "102473733",
  "102523258", "610121903", "610104837", "888328133", "102527866",
  "213797155", "666062877", "102775723", "102554490", "666075028",
  "111636931", "111631934", "213640305", "210744884", "210744746",
  "210744882", "210707909", "210734985", "210740947", "210740946",
  "666047983", "666050147", "666072527", "130384461", "210740916",
  "210743105", "120395862", "102856101", "102379827", "610099233",
  "102432523", "102344508", "610126260"
))
cat("Exclusiones únicas definidas:", length(npet_excluidos), "\n")

# ── Importar datos clínicos ───────────────────────────────────────────────────
pacientes <- read_dta(f_muestras)
setDT(pacientes)
pacientes <- pacientes %>% mutate(Npet = as.character(Npet))

# ── Aplicar exclusiones antes del join ───────────────────────────────────────
n_antes        <- n_distinct(all_data$Npet)
all_data_clean <- all_data %>% filter(!Npet %in% npet_excluidos)
n_despues      <- n_distinct(all_data_clean$Npet)
cat("Muestras antes de exclusiones:", n_antes,             "\n")
cat("Muestras excluidas:            ", n_antes - n_despues, "\n")
cat("Muestras tras exclusiones:     ", n_despues,           "\n")

# ── Fusión con datos clínicos ─────────────────────────────────────────────────
all_data_nhc <- all_data_clean %>%
  left_join(pacientes %>% dplyr::select(Npet, NHC, CV, Freg), by = "Npet") %>%
  filter(!is.na(NHC))

# ── Diagnóstico ───────────────────────────────────────────────────────────────
cat("\n── Pacientes (allsamples.dta) ──────────────────\n")
cat("Npet únicos:", n_distinct(pacientes$Npet), "\n")
cat("NHC únicos: ", n_distinct(pacientes$NHC),  "\n")
cat("\n── all_data_nhc (tras exclusiones + join) ───────\n")
cat("Npet únicos:", n_distinct(all_data_nhc$Npet), "\n")
cat("NHC únicos: ", n_distinct(all_data_nhc$NHC),  "\n")

# ── Comprobación: no deben quedar muestras sin match ─────────────────────────
sin_match_restantes <- setdiff(unique(all_data_clean$Npet), unique(pacientes$Npet))
if (length(sin_match_restantes) > 0) {
  cat("\n⚠️  Muestras aún sin match en allsamples.dta:",
      paste(sin_match_restantes, collapse = ", "), "\n")
} else {
  cat("\nTodas las muestras restantes tienen match en allsamples.dta.\n")
}

