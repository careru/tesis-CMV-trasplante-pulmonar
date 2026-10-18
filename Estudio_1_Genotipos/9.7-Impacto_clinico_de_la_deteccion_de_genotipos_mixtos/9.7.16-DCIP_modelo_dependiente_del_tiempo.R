# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.7.16. DCIP: modelo dependiente del tiempo
# Extracto del script completo (ver README.md del estudio)

## ---- DCIP: modelo dependiente del tiempo ----
res_DCIP <- list(
  gB   = analisis_td("gbm",   "DCIP", "F_DCIP", "Disfunción crónica del injerto pulmonar | gB mixtos"),
  gH   = analisis_td("ghm",   "DCIP", "F_DCIP", "Disfunción crónica del injerto pulmonar | gH mixtos"),
  gBgH = analisis_td("gbghm", "DCIP", "F_DCIP", "Disfunción crónica del injerto pulmonar | gB/gH mixtos")
)

## Descomentar para ver los nueve modelos de cada clasificación
# analisis_td("gbm",   "DCIP", "F_DCIP", "Disfunción crónica del injerto pulmonar | gB mixtos", solo_principal = FALSE)
# analisis_td("ghm",   "DCIP", "F_DCIP", "Disfunción crónica del injerto pulmonar | gH mixtos", solo_principal = FALSE)
# analisis_td("gbghm", "DCIP", "F_DCIP", "Disfunción crónica del injerto pulmonar | gB/gH mixtos", solo_principal = FALSE)

