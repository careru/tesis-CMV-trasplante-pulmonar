# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2. Análisis global de variantes detectadas mediante NGS
# Extracto del script completo (ver README.md del estudio)

# ---- 10.2. Análisis global de variantes detectadas mediante NGS ----

## ---- Muestras analizadas por gen (26 muestras; UL97 26, UL54 24, UL56 11) ----
npet_txp_total <- all_data_nhc %>%
  filter(NHC %in% pacientes_txp) %>%
  distinct(Npet) %>%
  nrow()
cat("Total peticiones únicas (Npet) en subcohorte TXP:", npet_txp_total, "\n\n")

cat("── Npet con resultado válido por gen ────────────────────────────────────\n")
all_data_nhc %>%
  filter(NHC %in% pacientes_txp) %>%
  group_by(Protein) %>%
  summarise(
    n_npet      = n_distinct(Npet),
    pct         = round(n_distinct(Npet) / npet_txp_total * 100, 1),
    n_pacientes = n_distinct(NHC),
    .groups     = "drop"
  ) %>%
  filter(Protein %in% c("UL97", "UL54", "UL56")) %>%
  arrange(Protein) %>%
  print()

npet_ul97 <- all_data_nhc %>%
  filter(NHC %in% pacientes_txp, Protein == "UL97") %>% distinct(Npet)
npet_ul54 <- all_data_nhc %>%
  filter(NHC %in% pacientes_txp, Protein == "UL54") %>% distinct(Npet)
npet_ul56 <- all_data_nhc %>%
  filter(NHC %in% pacientes_txp, Protein == "UL56") %>% distinct(Npet)

cat("\n── Npet en UL97 pero no en UL54 (cobertura insuficiente u omitidas) ─────\n")
npet_solo_ul97 <- setdiff(npet_ul97$Npet, npet_ul54$Npet)
cat("N =", length(npet_solo_ul97), "\n")
cat(npet_solo_ul97, sep = "\n")

cat("\n── Npet en UL54 pero no en UL97 ─────────────────────────────────────────\n")
npet_solo_ul54 <- setdiff(npet_ul54$Npet, npet_ul97$Npet)
cat("N =", length(npet_solo_ul54), "\n")
cat(npet_solo_ul54, sep = "\n")

cat("\n── Npet con UL56 analizado ───────────────────────────────────────────────\n")
cat("N =", nrow(npet_ul56), "\n")
cat("Porcentaje sobre total TXP:",
    round(nrow(npet_ul56) / npet_txp_total * 100, 1), "%\n")

cat("\n── Desglose por Npet: genes analizados (TRUE/FALSE) ─────────────────────\n")
all_data_nhc %>%
  filter(NHC %in% pacientes_txp) %>%
  distinct(Npet, NHC) %>%
  left_join(npet_ul97 %>% mutate(UL97 = TRUE), by = "Npet") %>%
  left_join(npet_ul54 %>% mutate(UL54 = TRUE), by = "Npet") %>%
  left_join(npet_ul56 %>% mutate(UL56 = TRUE), by = "Npet") %>%
  mutate(across(c(UL97, UL54, UL56), ~replace_na(., FALSE))) %>%
  left_join(
    all_data_nhc %>% filter(NHC %in% pacientes_txp) %>% distinct(Npet, RTP),
    by = "Npet"
  ) %>%
  arrange(RTP, Npet) %>%
  print(n = Inf)

# ── Diagnóstico calidad: resumen de cobertura y Q-Score en subcohorte TXP ────
cat("\n── Resumen de métricas de calidad (subcohorte TXP, Prevalence ≥5%) ──────\n")
all_data_txp %>%
  filter(Protein %in% c("UL97", "UL54", "UL56")) %>%
  group_by(Protein) %>%
  summarise(
    cob_media   = round(mean(cobertura_total, na.rm = TRUE), 0),
    cob_mediana = round(median(cobertura_total, na.rm = TRUE), 0),
    cob_min     = round(min(cobertura_total, na.rm = TRUE), 0),
    q_mediana   = round(median(q_score, na.rm = TRUE), 1),
    pct_balance = round(mean(balance_ok, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) %>%
  rename(
    "Gen"             = Protein,
    "Cob. media (×)"  = cob_media,
    "Cob. mediana (×)"= cob_mediana,
    "Cob. mínima (×)" = cob_min,
    "Q-Score mediana" = q_mediana,
    "% balance FW/RV" = pct_balance
  ) %>%
  print()

