# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.4.01. Datos de VND
# Extracto del script completo (ver README.md del estudio)

## ---- Datos de VND ----
# ── Resumen actualizado de variantes NEW ──────────────────────────────────────
new_vars <- all_data_txp %>%
  filter(class == "NEW", Prevalence >= 5) %>%
  mutate(Position = as.numeric(Position)) %>%
  group_by(Protein, Mutation, Position) %>%
  summarise(
    n_muestras  = n_distinct(Npet),
    n_pacientes = n_distinct(NHC),
    prev_media  = round(mean(Prevalence, na.rm = TRUE), 1),
    prev_max    = round(max(Prevalence, na.rm = TRUE), 1),
    .groups = "drop"
  )

cat("Total NEW por gen:\n")
new_vars %>% count(Protein) %>% print()

cat("\nNEW por gen y umbral:\n")
map_dfr(c(5, 10, 15, 20), function(thr) {
  all_data_txp %>%
    filter(class == "NEW", Prevalence >= thr) %>%
    group_by(Protein) %>%
    summarise(n_variantes = n_distinct(Mutation),
              n_pacientes = n_distinct(NHC), .groups = "drop") %>%
    mutate(umbral = paste0("≥", thr, "%"))
}) %>%
  arrange(Protein, umbral) %>%
  print(n = Inf)


new_vars_dom <- new_vars %>%
  mutate(
    dominio = case_when(
      Protein == "UL97" & Position >= 335 & Position <= 370 ~ "P-loop",
      Protein == "UL97" & Position >= 409 & Position <= 411 ~ "Gatekeeper",
      Protein == "UL97" & Position >= 440 & Position <= 480 ~ "Loop catalítico",
      Protein == "UL97" & Position >= 520 & Position <= 607 ~ "Región canónica GCV",
      Protein == "UL54" & Position >= 408 & Position <= 412 ~ "Región IV",
      Protein == "UL54" & Position >= 501 & Position <= 545 ~ "Delta-C",
      Protein == "UL54" & Position >= 700 & Position <= 715 ~ "Región II",
      Protein == "UL54" & Position >= 725 & Position <= 964 ~ "Unión ATP/dNTP",
      Protein == "UL54" & Position >= 802 & Position <= 821 ~ "Región III",
      Protein == "UL54" & Position >= 910 & Position <= 913 ~ "Región I",
      Protein == "UL54" & Position >= 962 & Position <= 987 ~ "Región V",
      Protein == "UL54" & Position >= 1153 & Position <= 1159 ~ "NLS-A",
      Protein == "UL54" & Position >= 1213 & Position <= 1242 ~ "UL44 binding",
      Protein == "UL56" & Position == 25 ~ "N-terminal",
      Protein == "UL56" & Position >= 202 & Position <= 220 ~ "LAGLIDADG",
      Protein == "UL56" & Position >= 231 & Position <= 369 ~ "Región LMV",
      Protein == "UL56" & Position >= 514 & Position <= 522 ~ "Interacción UL89",
      Protein == "UL56" & Position >= 671 & Position <= 680 ~ "Interacción UL89",
      Protein == "UL56" & Position >= 757 & Position <= 764 ~ "Interacción UL89",
      TRUE ~ "Fuera de dominio funcional"
    ),
    en_dominio = dominio != "Fuera de dominio funcional"
  ) %>%
  arrange(Protein, desc(en_dominio), desc(n_pacientes), desc(prev_max))

cat("── En dominio funcional ─────────────────────────\n")
new_vars_dom %>% filter(en_dominio) %>% print(n = Inf)

cat("\n── Fuera de dominio pero relevantes ─────────────\n")
new_vars_dom %>%
  filter(!en_dominio, n_pacientes >= 2 | prev_max >= 15) %>%
  print(n = Inf)

