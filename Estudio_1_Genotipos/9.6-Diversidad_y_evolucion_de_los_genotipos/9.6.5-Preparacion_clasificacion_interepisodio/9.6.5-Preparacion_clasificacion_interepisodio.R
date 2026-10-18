# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.6.5. Preparación: clasificación interepisodio
# Extracto del script completo (ver README.md del estudio)

## ---- Preparación: clasificación interepisodio ----
# ── 1. Pacientes con >1 episodio genotipado ───────────────────────────────────
ep_por_paciente <- dnaemias %>%
  filter(!is.na(gbghmix)) %>%
  mutate(Fini = as.Date(Fini)) %>%
  group_by(NHC) %>%
  summarise(n_episodios = n(), .groups = "drop") %>%
  filter(n_episodios > 1)

cat(sprintf("Pacientes con >1 episodio: %d\n", nrow(ep_por_paciente)))
cat(sprintf("Total episodios: %d\n", sum(ep_por_paciente$n_episodios)))

cat("\n── Distribución por número de episodios ──\n")
ep_por_paciente %>%
  count(n_episodios) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()

# ── 2. Episodios de pacientes con >1 episodio ─────────────────────────────────
ep_multi <- dnaemias %>%
  filter(!is.na(gbghmix)) %>%
  mutate(Fini = as.Date(Fini)) %>%
  semi_join(ep_por_paciente, by = "NHC")

cat("\n── Clasificación único/mixto ──\n")
ep_multi %>%
  mutate(clase = ifelse(gbghmix == 0, "Único", "Mixto")) %>%
  count(clase) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()

# ── 3. Clasificar pacientes según comportamiento interepisodio ────────────────
patron_paciente <- ep_multi %>%
  group_by(NHC) %>%
  summarise(
    n_episodios = n(),
    n_unicos    = sum(gbghmix == 0),
    n_mixtos    = sum(gbghmix == 1),
    .groups     = "drop"
  ) %>%
  mutate(
    patron = case_when(
      n_unicos > 0 & n_mixtos == 0 ~ "Siempre único",
      n_unicos == 0 & n_mixtos > 0 ~ "Siempre mixto",
      TRUE                          ~ "Mixtos y únicos"
    )
  )

cat("\n── Clasificación pacientes ──\n")
patron_paciente %>%
  count(patron) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()

cat("\n── Episodios por grupo ──\n")
patron_paciente %>%
  group_by(patron) %>%
  summarise(n_episodios = sum(n_episodios), .groups = "drop") %>%
  mutate(pct = 100 * n_episodios / sum(n_episodios)) %>%
  print()

# ── 4. Perfil serológico por grupo ────────────────────────────────────────────
cat("\n── Perfil serológico por grupo ──\n")
patron_paciente %>%
  left_join(pacientes %>% select(NHC, serDR), by = "NHC") %>%
  group_by(patron, serDR) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(patron) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()

# ── 5. Tiempo entre episodios por grupo ───────────────────────────────────────
cat("\n── Tiempo entre episodios (Ffin → Fini siguiente, días) ──\n")
ep_multi %>%
  mutate(Ffin = as.Date(Ffin)) %>%
  arrange(NHC, Fini) %>%
  group_by(NHC) %>%
  mutate(tiempo_entre = as.numeric(Fini - lag(Ffin))) %>%
  filter(!is.na(tiempo_entre), tiempo_entre >= 0) %>%
  left_join(patron_paciente %>% select(NHC, patron), by = "NHC") %>%
  group_by(patron) %>%
  summarise(
    mediana = median(tiempo_entre, na.rm = TRUE),
    q1      = quantile(tiempo_entre, 0.25, na.rm = TRUE),
    q3      = quantile(tiempo_entre, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

# ── 6. Análisis genotipos estables/variables interepisodio ────────────────────
pacientes_cambio <- patron_paciente %>%
  filter(patron == "Mixtos y únicos") %>%
  pull(NHC)

ep_cambio <- ep_multi %>%
  filter(NHC %in% pacientes_cambio) %>%
  arrange(NHC, Fini)

genotipos <- c("gb1","gb2","gb3","gb4","gh1","gh2")
nombres   <- c("gB1","gB2","gB3","gB4","gH1","gH2")

estab_inter <- ep_cambio %>%
  group_by(NHC) %>%
  summarise(
    n_ep = n(),
    across(all_of(genotipos), ~ {
      presente <- sum(. == 1, na.rm = TRUE)
      n_tot    <- n()
      if (presente == 0) "ausente"
      else if (presente == n_tot) "estable"
      else "variable"
    }),
    .groups = "drop"
  ) %>%
  mutate(
    tiene_estable = rowSums(across(all_of(genotipos), ~ . == "estable")) > 0
  )

cat(sprintf("\n── Pacientes con ≥1 genotipo estable: %d/21 (%.1f%%)\n",
            sum(estab_inter$tiene_estable),
            100 * mean(estab_inter$tiene_estable)))

cat(sprintf("── Pacientes sin ningún genotipo estable: %d\n",
            sum(!estab_inter$tiene_estable)))

cat("\n── Frecuencia genotipos estables ──\n")
for (i in seq_along(genotipos)) {
  n_est <- sum(estab_inter[[genotipos[i]]] == "estable")
  cat(sprintf("%s estable: %d/21 (%.1f%%)\n", nombres[i], n_est, 100*n_est/21))
}

cat("\n── Frecuencia genotipos variables ──\n")
for (i in seq_along(genotipos)) {
  n_var <- sum(estab_inter[[genotipos[i]]] == "variable")
  cat(sprintf("%s variable: %d/21 (%.1f%%)\n", nombres[i], n_var, 100*n_var/21))
}

# ── 7. Número de cambios entre episodios consecutivos ─────────────────────────
cambios_inter <- ep_cambio %>%
  arrange(NHC, Fini) %>%
  group_by(NHC) %>%
  mutate(
    patron_ep = paste(gb1, gb2, gb3, gb4, gh1, gh2, sep = "-"),
    cambio    = patron_ep != lag(patron_ep)
  ) %>%
  summarise(
    n_cambios = sum(cambio, na.rm = TRUE),
    .groups   = "drop"
  )

