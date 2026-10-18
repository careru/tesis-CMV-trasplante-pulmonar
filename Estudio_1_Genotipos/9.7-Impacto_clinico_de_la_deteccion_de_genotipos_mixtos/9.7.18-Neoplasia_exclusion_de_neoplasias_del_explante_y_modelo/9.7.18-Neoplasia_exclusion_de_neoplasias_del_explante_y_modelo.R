# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.7.18. Neoplasia: exclusión de neoplasias del explante y modelo dependiente del tiempo
# Extracto del script completo (ver README.md del estudio)

## ---- Neoplasia: exclusión de neoplasias del explante y modelo dependiente del tiempo ----
explante <- !is.na(pacientes$F_onco) & as.Date(pacientes$F_onco) == as.Date(pacientes$Ftx)
pacientes <- pacientes %>%
  mutate(onco       = ifelse(explante, 0L, onco),
         F_onco     = if_else(explante, as.Date(NA), as.Date(F_onco)),
         onco_tipo  = ifelse(explante, NA, onco_tipo),
         onco_tipo2 = ifelse(explante, NA, onco_tipo2),
         onco_det   = ifelse(explante, NA, onco_det))

res_onco <- list(
  gB   = analisis_td("gbm",   "onco", "F_onco", "Neoplasia global | gB mixtos"),
  gH   = analisis_td("ghm",   "onco", "F_onco", "Neoplasia global | gH mixtos"),
  gBgH = analisis_td("gbghm", "onco", "F_onco", "Neoplasia global | gB/gH mixtos")
)

## Descomentar para ver los nueve modelos de cada clasificación
# analisis_td("gbm",   "onco", "F_onco", "Neoplasia global | gB mixtos", solo_principal = FALSE)
# analisis_td("ghm",   "onco", "F_onco", "Neoplasia global | gH mixtos", solo_principal = FALSE)
# analisis_td("gbghm", "onco", "F_onco", "Neoplasia global | gB/gH mixtos", solo_principal = FALSE)

