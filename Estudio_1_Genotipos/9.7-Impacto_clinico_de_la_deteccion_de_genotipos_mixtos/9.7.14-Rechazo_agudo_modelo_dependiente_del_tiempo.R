# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.7.14. Rechazo agudo: modelo dependiente del tiempo
# Extracto del script completo (ver README.md del estudio)

## ---- Rechazo agudo: modelo dependiente del tiempo ----
## 21 de 26 rechazos ocurren antes del primer genotipado: quedan 5 eventos evaluables.
## El modelo se ejecuta para documentarlo, pero no es interpretable y no se dibuja KM.
res_RA <- list(
  gB   = analisis_td("gbm",   "RA", "F_RA1", "Rechazo agudo | gB mixtos"),
  gH   = analisis_td("ghm",   "RA", "F_RA1", "Rechazo agudo | gH mixtos"),
  gBgH = analisis_td("gbghm", "RA", "F_RA1", "Rechazo agudo | gB/gH mixtos")
)

## Descomentar para ver los nueve modelos de cada clasificación
# analisis_td("gbm",   "RA", "F_RA1", "Rechazo agudo | gB mixtos", solo_principal = FALSE)
# analisis_td("ghm",   "RA", "F_RA1", "Rechazo agudo | gH mixtos", solo_principal = FALSE)
# analisis_td("gbghm", "RA", "F_RA1", "Rechazo agudo | gB/gH mixtos", solo_principal = FALSE)

