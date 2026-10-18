# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.7.20. Tumor sólido no cutáneo: modelo dependiente del tiempo
# Extracto del script completo (ver README.md del estudio)

## ---- Tumor sólido no cutáneo: modelo dependiente del tiempo ----
res_solido <- list(
  gB   = analisis_td("gbm",   "solido", "F_onco", "Tumor sólido no cutáneo | gB mixtos"),
  gH   = analisis_td("ghm",   "solido", "F_onco", "Tumor sólido no cutáneo | gH mixtos"),
  gBgH = analisis_td("gbghm", "solido", "F_onco", "Tumor sólido no cutáneo | gB/gH mixtos")
)

## Descomentar para ver los nueve modelos de cada clasificación
# analisis_td("gbm",   "solido", "F_onco", "Tumor sólido no cutáneo | gB mixtos", solo_principal = FALSE)
# analisis_td("ghm",   "solido", "F_onco", "Tumor sólido no cutáneo | gH mixtos", solo_principal = FALSE)
# analisis_td("gbghm", "solido", "F_onco", "Tumor sólido no cutáneo | gB/gH mixtos", solo_principal = FALSE)

