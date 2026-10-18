# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.6.1. Preparación: muestras genotipadas por episodio
# Extracto del script completo (ver README.md del estudio)

## ---- Preparación: muestras genotipadas por episodio ----
# ── 1. Cargar muestras ────────────────────────────────────────────────────────
muestras <- read_excel(f_bd,
                       sheet = "genotipado",
                       col_types = c("numeric", "numeric", "date", "numeric",
                                     "numeric", "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "text"))

# ── 1b. Filtrar solo muestras con genotipado completo ─────────────────────────
muestras <- muestras %>%
  mutate(
    gb = rowSums(across(c(GB1, GB2, GB3, GB4)), na.rm = TRUE),
    gh = rowSums(across(c(GH1, GH2)), na.rm = TRUE)
  ) %>%
  filter(gb >= 1 & gh >= 1)

cat(sprintf("Muestras con genotipado completo: %d\n", nrow(muestras)))
cat(sprintf("Pacientes: %d\n", n_distinct(muestras$NHC)))

# ── 2. Asignar cada muestra a su episodio ─────────────────────────────────────
muestras_ep <- muestras %>%
  inner_join(
    dnaemias %>%
      select(NHC, Fini, Ffin, gbghmix) %>%
      mutate(Fini = as.Date(Fini), Ffin = as.Date(Ffin)),
    by = "NHC",
    relationship = "many-to-many"
  ) %>%
  filter(
    Freg >= Fini &
      (is.na(Ffin) | Freg <= Ffin)
  )

cat(sprintf("Muestras asignadas: %d\n", nrow(muestras_ep)))
cat(sprintf("Episodios representados: %d\n",
            n_distinct(muestras_ep %>% select(NHC, Fini))))

# ── 3. Episodios con ≥2 muestras ──────────────────────────────────────────────
ep_n_muestras <- muestras_ep %>%
  group_by(NHC, Fini, Ffin, gbghmix) %>%
  summarise(n_muestras = n(), .groups = "drop")

n_total_ep <- 234
n_intra_ep <- ep_n_muestras %>% filter(n_muestras >= 2) %>% nrow()

cat(sprintf("Episodios con ≥2 muestras: %d (%.1f%%)\n",
            n_intra_ep, 100 * n_intra_ep / n_total_ep))

ep_intra <- ep_n_muestras %>% filter(n_muestras >= 2)

cat(sprintf("Mediana muestras por episodio: %.0f (RIQ %.0f-%.0f)\n",
            median(ep_intra$n_muestras),
            quantile(ep_intra$n_muestras, 0.25),
            quantile(ep_intra$n_muestras, 0.75)))

# ── 4. Clasificar episodios únicos vs mixtos ──────────────────────────────────
cat("\n── Clasificación episodios con ≥2 muestras ──\n")
ep_intra %>%
  mutate(clase = ifelse(gbghmix == 0, "Único", "Mixto")) %>%
  count(clase) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()

# ── 5. Variabilidad intraepisodio ─────────────────────────────────────────────
muestras_intra <- muestras_ep %>%
  semi_join(ep_intra, by = c("NHC", "Fini")) %>%
  mutate(patron = paste(GB1, GB2, GB3, GB4, GH1, GH2, sep = "-"))

variabilidad_ep <- muestras_intra %>%
  group_by(NHC, Fini, gbghmix) %>%
  summarise(
    n_muestras = n(),
    n_patrones = n_distinct(patron),
    .groups    = "drop"
  ) %>%
  mutate(
    clase    = ifelse(gbghmix == 0, "Único", "Mixto"),
    variable = ifelse(n_patrones > 1, "Variable", "Estable")
  )

cat("\n── Variabilidad por clase ──\n")
variabilidad_ep %>%
  group_by(clase, variable) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(clase) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()

# ── 6. Episodios mixtos variables ─────────────────────────────────────────────
mixtos_var <- variabilidad_ep %>%
  filter(clase == "Mixto", variable == "Variable")

cat(sprintf("\n── Episodios mixtos variables: %d episodios, %d pacientes ──\n",
            nrow(mixtos_var), n_distinct(mixtos_var$NHC)))

mixtos_var %>%
  left_join(pacientes %>% select(NHC, serDR), by = "NHC") %>%
  distinct(NHC, serDR) %>%
  count(serDR) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()

