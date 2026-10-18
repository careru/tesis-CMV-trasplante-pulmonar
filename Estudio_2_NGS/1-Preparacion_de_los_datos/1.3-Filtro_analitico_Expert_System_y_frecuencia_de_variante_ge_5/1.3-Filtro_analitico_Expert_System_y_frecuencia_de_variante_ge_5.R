# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 1.3. Filtro analítico: Expert System y frecuencia de variante ≥ 5%
# Extracto del script completo (ver README.md del estudio)

## ---- Filtro analítico: Expert System y frecuencia de variante ≥ 5% ----
all_data_nhc <- all_data_nhc %>%
  mutate(Prevalence = as.numeric(Prevalence)) %>%
  filter(Filtro == "Yes", Prevalence >= 5, cobertura_total >= 100)

cat("── Tras filtro Prevalence ≥5% ──────────────────\n")
cat("Variantes totales:", nrow(all_data_nhc),          "\n")
cat("Npet únicos:      ", n_distinct(all_data_nhc$Npet), "\n")
cat("NHC únicos:       ", n_distinct(all_data_nhc$NHC),  "\n")

cat("\nDistribución por gen:\n")
all_data_nhc %>%
  count(Protein, name = "n_variantes") %>%
  arrange(Protein) %>%
  print()

