# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.4.03. Resúmenes de frecuencias genotípicas por muestra (datos para las Figuras 19-21)
# Extracto del script completo (ver README.md del estudio)

## ---- Resúmenes de frecuencias genotípicas por muestra (datos para las Figuras 19-21) ----
df_gb <- muestras_sel %>%
  mutate(
    gb_sum = rowSums(across(GB1:GB4, ~replace_na(., 0))),  # suma de GB1-GB4
    combo  = apply(select(., GB1:GB4) == 1, 1, function(r)
      paste(colnames(select(., GB1:GB4))[which(r)], collapse = "+"))
  ) %>%
  filter(gb_sum %in% 1:3)

tabla_resumen <- df_gb %>%
  count(gb_sum, name = "n") %>%
  mutate(
    tipo = factor(paste(gb_sum, "genotipo(s)"),
                  levels = c("1 genotipo(s)","2 genotipo(s)","3 genotipo(s)")),
    porcentaje = round(100 * n / sum(n), 1)
  ) %>%
  arrange(tipo) %>%
  select(tipo, n, porcentaje)

tabla_resumen

desglose <- df_gb %>%
  count(gb_sum, combo, name = "n") %>%
  group_by(gb_sum) %>%
  mutate(porcentaje = round(100 * n / sum(n), 1)) %>%
  arrange(gb_sum, desc(n)) %>%
  ungroup()

desglose_1 <- desglose %>% filter(gb_sum == 1) %>% select(combo, n, porcentaje)
desglose_2 <- desglose %>% filter(gb_sum == 2) %>% select(combo, n, porcentaje)
desglose_3 <- desglose %>% filter(gb_sum == 3) %>% select(combo, n, porcentaje)

desglose_1  # combinaciones cuando hay 1 genotipo
desglose_2  # combinaciones cuando hay 2 genotipos
desglose_3  # combinaciones cuando hay 3 genotipos

df_gh <- muestras_sel %>%
  mutate(
    gh_sum = rowSums(across(GH1:GH2, ~replace_na(., 0))),  # suma GH1–GH2
    combo  = apply(select(., GH1:GH2) == 1, 1, function(r)
      paste(colnames(select(., GH1:GH2))[which(r)], collapse = "+"))
  ) %>%
  filter(gh_sum %in% 1:2)

tabla_resumen_gh <- df_gh %>%
  count(gh_sum, name = "n") %>%
  mutate(
    tipo = factor(paste(gh_sum, "genotipo(s)"),
                  levels = c("1 genotipo(s)", "2 genotipo(s)")),
    porcentaje = round(100 * n / sum(n), 1)
  ) %>%
  arrange(tipo) %>%
  select(tipo, n, porcentaje)

tabla_resumen_gh

desglose_gh <- df_gh %>%
  count(gh_sum, combo, name = "n") %>%
  group_by(gh_sum) %>%
  mutate(porcentaje = round(100 * n / sum(n), 1)) %>%
  arrange(gh_sum, desc(n)) %>%
  ungroup()

desglose_1_gh <- desglose_gh %>%
  filter(gh_sum == 1) %>%
  select(combo, n, porcentaje)

desglose_2_gh <- desglose_gh %>%
  filter(gh_sum == 2) %>%
  select(combo, n, porcentaje)

desglose_1_gh   # combinaciones cuando hay 1 genotipo gH
desglose_2_gh   # combinaciones cuando hay 2 genotipos gH

# Partimos de `muestras_sel` (solo muestras con genotipado completo)

df_gbgh <- muestras_sel %>%
  mutate(
    # nº de genotipos presentes
    gb_sum = rowSums(across(GB1:GB4, ~ tidyr::replace_na(., 0))),
    gh_sum = rowSums(across(GH1:GH2, ~ tidyr::replace_na(., 0))),
    
    # combinaciones gB (formato gB1, gB1 + gB3, etc.)
    gB_combo = apply(select(., GB1:GB4) == 1, 1, function(r) {
      if (!any(r)) return(NA_character_)
      paste0("gB", which(r), collapse = " + ")
    }),
    
    # combinaciones gH (formato gH1, gH1 + gH2)
    gH_combo = apply(select(., GH1:GH2) == 1, 1, function(r) {
      if (!any(r)) return(NA_character_)
      paste0("gH", which(r), collapse = " + ")
    }),
    
    # tipo único/mixto por proteína
    gB_tipo = if_else(gb_sum == 1, "gBunico", "gBmixto"),
    gH_tipo = if_else(gh_sum == 1, "gHunico", "gHmixto"),
    
    # grupo combinado
    grupo = paste(gB_tipo, gH_tipo, sep = "_"),
    
    # combinación conjunta gB/gH
    combo_total = paste(gB_combo, gH_combo, sep = " / ")
  ) %>%
  # por seguridad, nos quedamos solo con filas con gB y gH definidos
  filter(!is.na(gB_combo), !is.na(gH_combo))

## 1) Tabla resumen de grupos (gBunico_gHunico, etc.)

niveles_grupo <- c(
  "gBunico_gHunico",
  "gBmixto_gHunico",
  "gBunico_gHmixto",
  "gBmixto_gHmixto"
)

tabla_grupos <- df_gbgh %>%
  count(grupo, name = "n") %>%
  mutate(
    grupo = factor(grupo, levels = niveles_grupo),
    porcentaje = round(100 * n / sum(n), 1)
  ) %>%
  arrange(grupo)

tabla_grupos

## Subdesglose gBunico_gHunico
desglose_gBunico_gHunico <- df_gbgh %>%
  filter(grupo == "gBunico_gHunico") %>%
  count(gB_combo, gH_combo, combo_total, name = "n") %>%
  mutate(
    porcentaje = round(100 * n / sum(n), 1)
  ) %>%
  arrange(desc(n))

desglose_gBunico_gHunico

## Subdesglose gBmixto_gHunico
desglose_gBmixto_gHunico <- df_gbgh %>%
  filter(grupo == "gBmixto_gHunico") %>%
  count(gB_combo, gH_combo, combo_total, name = "n") %>%
  mutate(
    porcentaje = round(100 * n / sum(n), 1)
  ) %>%
  arrange(desc(n))

desglose_gBmixto_gHunico

## Subdesglose gBunico_gHmixto
desglose_gBunico_gHmixto <- df_gbgh %>%
  filter(grupo == "gBunico_gHmixto") %>%
  count(gB_combo, gH_combo, combo_total, name = "n") %>%
  mutate(
    porcentaje = round(100 * n / sum(n), 1)
  ) %>%
  arrange(desc(n))

desglose_gBunico_gHmixto

## 5) Subdesglose gBmixto_gHmixto
desglose_gBmixto_gHmixto <- df_gbgh %>%
  filter(grupo == "gBmixto_gHmixto") %>%
  count(gB_combo, gH_combo, combo_total, name = "n") %>%
  mutate(
    porcentaje = round(100 * n / sum(n), 1)
  ) %>%
  arrange(desc(n))

desglose_gBmixto_gHmixto

