# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.7.22. Mortalidad: modelo dependiente del tiempo
# Extracto del script completo (ver README.md del estudio)

## ---- Mortalidad: modelo dependiente del tiempo ----
## Ajuste homogéneo con el resto de desenlaces: edad al trasplante + perfil serológico (risk_f)
res_exitus <- list(
  gB   = analisis_td("gbm",   "exitus", "Fexitus", "Mortalidad | gB mixtos"),
  gH   = analisis_td("ghm",   "exitus", "Fexitus", "Mortalidad | gH mixtos"),
  gBgH = analisis_td("gbghm", "exitus", "Fexitus", "Mortalidad | gB/gH mixtos")
)

## Descomentar para ver los nueve modelos de cada clasificación
# analisis_td("gbm",   "exitus", "Fexitus", "Mortalidad | gB mixtos", solo_principal = FALSE)
# analisis_td("ghm",   "exitus", "Fexitus", "Mortalidad | gH mixtos", solo_principal = FALSE)
# analisis_td("gbghm", "exitus", "Fexitus", "Mortalidad | gB/gH mixtos", solo_principal = FALSE)

