# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.2.1. Datos de PN
# Extracto del script completo (ver README.md del estudio)

## ---- Datos de PN ----
all_data_txp %>%
  filter(class == "POL", Prevalence >= 5) %>%
  group_by(Protein) %>%
  summarise(
    n_variantes = n_distinct(Mutation),
    n_pacientes = n_distinct(NHC),
    .groups = "drop"
  )

all_data_txp %>%
  filter(class == "POL", Prevalence >= 5) %>%
  group_by(Protein, Mutation) %>%
  summarise(
    NHC_list = paste(sort(unique(NHC)), collapse = ", "),
    n_pac    = n_distinct(NHC),
    .groups  = "drop"
  ) %>%
  arrange(Protein, desc(n_pac)) %>%
  print(n = Inf)

