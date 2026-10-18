# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 0.1. Rutas y archivos de datos
# Extracto del script completo (ver README.md del estudio)

## ---- Rutas y archivos de datos ----
dir_datos  <- "~/Documents/CMV/ANALISIS"     # carpeta con los tres archivos de datos
dir_salida <- file.path(dir_datos, "salidas_tesis")

f_bd        <- "_BD CMV TXP CLINICA 22-08-24.xlsx"     # base clínica (pacientes, dnaemias...)
f_genotipos <- "Resultados_genotiposMIL_COMPLETO.xls"  # resultados de genotipado (618 muestras)
f_cargas    <- "DTA_CV_TXP.dta"                        # cargas virales seriadas

ejecutar_anexo <- FALSE   # TRUE para ejecutar también los bloques del anexo

setwd(dir_datos)
dir.create(dir_salida, showWarnings = FALSE)

