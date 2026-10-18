# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 0.1. Rutas y archivos de datos
# Extracto del script completo (ver README.md del estudio)

## ---- Rutas y archivos de datos ----
dir_datos  <- "~/Documents/CMV/ESTUDIO VARIANTES NGS"   # carpeta de trabajo
dir_aa     <- file.path(dir_datos, "NGS_CMV/AA")       # CSV de DeepChek (uno por muestra)
dir_salida <- file.path(dir_datos, "salidas_tesis")

f_muestras <- file.path(dir_datos, "allsamples.dta")             # muestras con NHC
f_clinica  <- file.path(dir_datos, "CMV_NGS_pacientes_DEF.xlsx")  # datos clínicos (Tabla 23)

ejecutar_anexo <- FALSE   # TRUE para ejecutar también los bloques del anexo

dir.create(dir_salida, showWarnings = FALSE)

