# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.2.3. Selección de las muestras con genotipado completo (n = 483)
# Extracto del script completo (ver README.md del estudio)

## ---- Selección de las muestras con genotipado completo (n = 483) ----
muestras_sel <- muestras %>%
  mutate(
    gb = rowSums(across(c(GB1, GB2, GB3, GB4)), na.rm = TRUE),
    gh = rowSums(across(c(GH1, GH2)), na.rm = TRUE),
    gt_cat = case_when(
      gb >= 1 & gh >= 1 ~ "Genotipado completo",
      xor(gb >= 1, gh >= 1) ~ "Genotipado parcial",
      TRUE ~ "Genotipado negativo"
    )
  ) %>%
  filter(gt_cat == "Genotipado completo")

#Número de pacientes
n_distinct(muestras_sel$NHC)

# Número de muestras por paciente mínimo y maximo
muestras_por_paciente <- muestras_sel %>%
  count(NHC, name = "n_muestras")

min(muestras_por_paciente$n_muestras)
max(muestras_por_paciente$n_muestras)

# Mediana + IQR
quantile(muestras_por_paciente$n_muestras,
         probs = c(0.25, 0.5, 0.75),
         na.rm = TRUE)


