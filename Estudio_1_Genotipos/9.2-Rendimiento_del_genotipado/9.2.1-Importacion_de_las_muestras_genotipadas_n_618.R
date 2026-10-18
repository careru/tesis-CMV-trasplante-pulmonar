# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.2.1. Importación de las muestras genotipadas (n = 618)
# Extracto del script completo (ver README.md del estudio)

## ---- Importación de las muestras genotipadas (n = 618) ----
muestras <- read_excel(f_genotipos,
                       col_types = c(
                         rep("numeric", 2), "date", "text", "numeric", "text",
                         rep("numeric", 14), "text", rep("skip", 3)))

length(unique(muestras$Nreg)) #692 muestras de 177 pacientes
length(unique(muestras$NHC))

muestras <- muestras %>% #quito 74 que no se encontraron / volumen insuf. (NE=.)
  filter(NE != ".")

head(muestras)
length(unique(muestras$Nreg)) #618 muestras de 169 pacientes
length(unique(muestras$NHC))

#Borrar entradas vacias  
muestras <- muestras %>% 
  filter(!is.na(Freg))

muestras %>%
  mutate(
    gb_sum = rowSums(across(c(GB1, GB2, GB3, GB4)), na.rm = TRUE),
    gh_sum = rowSums(across(c(GH1, GH2)), na.rm = TRUE),
    categoria = case_when(
      gb_sum > 0 & gh_sum > 0 ~ "Genotipado completo",
      gb_sum > 0 & gh_sum == 0 ~ "Parcial - Solo gB",
      gh_sum > 0 & gb_sum == 0 ~ "Parcial - solo gH",
      gb_sum == 0 & gh_sum == 0 ~ "Negativo gB y gH" )) %>%
  count(categoria) %>%
  mutate(porcentaje = round(100 * n / sum(n), 1))

