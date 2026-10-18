# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 1.4. Selección de los receptores de trasplante pulmonar (RTP)
# Extracto del script completo (ver README.md del estudio)

## ---- Selección de los receptores de trasplante pulmonar (RTP) ----
pacientes_txp <- c(
  "NHC_RTP-2", "NHC_RTP-11", "NHC_RTP-1", "NHC_RTP-9", "NHC_RTP-4",
  "NHC_RTP-10", "NHC_RTP-6", "NHC_RTP-5", "NHC_RTP-3", "NHC_RTP-8", "NHC_RTP-7"
)

all_data_txp <- all_data_nhc %>%
  filter(NHC %in% pacientes_txp)

cat("── Subcohorte TXP ──────────────────────────────\n")
cat("Pacientes TXP definidos:", length(pacientes_txp),         "\n")
cat("Pacientes TXP con datos:", n_distinct(all_data_txp$NHC),  "\n")
cat("Npet únicos TXP:        ", n_distinct(all_data_txp$Npet), "\n")
cat("Variantes totales TXP:  ", nrow(all_data_txp),            "\n")

cat("\nDistribución por gen (TXP):\n")
all_data_txp %>%
  count(Protein, name = "n_variantes") %>%
  arrange(Protein) %>%
  print()

sin_datos <- setdiff(pacientes_txp, unique(all_data_txp$NHC))
if (length(sin_datos) > 0) {
  cat("\n⚠️  Pacientes TXP sin variantes tras filtro:",
      paste(sin_datos, collapse = ", "), "\n")
} else {
  cat("\nTodos los pacientes TXP tienen datos tras el filtro.\n")
}

