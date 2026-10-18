# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.15. Anexo: Importación de dnaemias SIN la AUC recalculada (no usar)
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Importación de dnaemias SIN la AUC recalculada (no usar) ----
dnaemias <- read_excel(f_bd, 
                       sheet = "dnaemias_sel")


episodios_por_serostatus <- dnaemias %>%
  filter(!is.na(serDR)) %>%     # por si hay algún NA
  count(serDR) %>%
  rename(
    seroestatus = serDR,
    n_episodios = n
  )

episodios_por_serostatus

