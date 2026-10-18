# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.7.12. Modelo dependiente del tiempo: infecciones oportunistas (para la tabla resumen)
# Extracto del script completo (ver README.md del estudio)

## ---- Modelo dependiente del tiempo: infecciones oportunistas (para la tabla resumen) ----
res_inf_opo <- list(
  gB   = analisis_td("gbm",   "inf_opo", "F_opo", "Infecciones oportunistas no CMV | gB mixtos"),
  gH   = analisis_td("ghm",   "inf_opo", "F_opo", "Infecciones oportunistas no CMV | gH mixtos"),
  gBgH = analisis_td("gbghm", "inf_opo", "F_opo", "Infecciones oportunistas no CMV | gB/gH mixtos")
)

## Descomentar para ver los nueve modelos de cada clasificación
# analisis_td("gbm",   "inf_opo", "F_opo", "Infecciones oportunistas no CMV | gB mixtos", solo_principal = FALSE)
# analisis_td("ghm",   "inf_opo", "F_opo", "Infecciones oportunistas no CMV | gH mixtos", solo_principal = FALSE)
# analisis_td("gbghm", "inf_opo", "F_opo", "Infecciones oportunistas no CMV | gB/gH mixtos", solo_principal = FALSE)

