# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.5.01. Importación de los episodios (n = 234) con la AUC recalculada
# Extracto del script completo (ver README.md del estudio)

## ---- Importación de los episodios (n = 234) con la AUC recalculada ----
# Única importación de dnaemias del script: incluye incorporar_AUC().
# (En el script antiguo, el bloque [4.1][a] repetido volvía a leerla sin la AUC recalculada.)
dnaemias <- read_excel(f_bd, 
                       sheet = "dnaemias_sel") %>%
  incorporar_AUC()

episodios_por_serostatus <- dnaemias %>%
  filter(!is.na(serDR)) %>%
  count(serDR) %>%
  rename(seroestatus = serDR, n_episodios = n)
episodios_por_serostatus

