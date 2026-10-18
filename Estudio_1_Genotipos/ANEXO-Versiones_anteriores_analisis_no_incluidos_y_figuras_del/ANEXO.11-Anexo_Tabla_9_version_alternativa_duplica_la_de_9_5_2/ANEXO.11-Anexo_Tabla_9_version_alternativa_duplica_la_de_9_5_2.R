# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.11. Anexo: Tabla 9, versión alternativa (duplica la de 9.5.2)
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Tabla 9, versión alternativa (duplica la de 9.5.2) ----
## ── 1. gB a nivel de muestra ─────────────────────────────────────────────
df_gb_m <- muestras_sel %>%
  mutate(
    gb_sum   = rowSums(across(GB1:GB4, ~ replace_na(., 0))),
    gB_combo = apply(select(., GB1:GB4) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gB", which(r), collapse = " + ")
    })
  ) %>%
  filter(gb_sum >= 1, !is.na(gB_combo))

n_m_total  <- nrow(df_gb_m)
n_m_gb_uni <- sum(df_gb_m$gb_sum == 1)
n_m_gb_mix <- sum(df_gb_m$gb_sum > 1)

freq_gb_uni_m <- df_gb_m %>%
  filter(gb_sum == 1) %>%
  count(gB_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_m_gb_uni, 1))

freq_gb_mix_m <- df_gb_m %>%
  filter(gb_sum > 1) %>%
  count(gB_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_m_gb_mix, 1)) %>%
  arrange(desc(n))

## ── 2. gB a nivel de episodio ────────────────────────────────────────────
df_gb_e <- dnaemias %>%
  mutate(
    gB_combo = apply(select(., gb1:gb4) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gB", which(r), collapse = " + ")
    })
  ) %>%
  filter(!is.na(gB_combo))

n_e_total  <- nrow(df_gb_e)
n_e_gb_uni <- sum(df_gb_e$gbsum == 1)
n_e_gb_mix <- sum(df_gb_e$gbsum > 1)

freq_gb_uni_e <- df_gb_e %>%
  filter(gbsum == 1) %>%
  count(gB_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_e_gb_uni, 1))

freq_gb_mix_e <- df_gb_e %>%
  filter(gbsum > 1) %>%
  count(gB_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_e_gb_mix, 1)) %>%
  arrange(desc(n))

## ── 3. Comprobacion en consola ───────────────────────────────────────────
cat(sprintf("Muestras: %d | unicos %d (%.1f%%) | mixtos %d (%.1f%%)\n",
            n_m_total, n_m_gb_uni, 100 * n_m_gb_uni / n_m_total,
            n_m_gb_mix, 100 * n_m_gb_mix / n_m_total))
cat(sprintf("Episodios: %d | unicos %d (%.1f%%) | mixtos %d (%.1f%%)\n",
            n_e_total, n_e_gb_uni, 100 * n_e_gb_uni / n_e_total,
            n_e_gb_mix, 100 * n_e_gb_mix / n_e_total))

## ── 4. Formato con coma decimal ──────────────────────────────────────────
fmt_np <- function(n, pct) {
  sprintf("%d (%s)", n,
          formatC(pct, format = "f", digits = 1, decimal.mark = ","))
}

get_m <- function(df, key) {
  r <- df[df[[1]] == key, ]
  if (nrow(r) == 0) return("\u2013")
  fmt_np(r$n, r$pct)
}

## ── 5. Construccion de las filas ─────────────────────────────────────────
combos_gb     <- sort(union(freq_gb_uni_m$gB_combo, freq_gb_uni_e$gB_combo))
combos_gb_mix <- union(freq_gb_mix_m$gB_combo, freq_gb_mix_e$gB_combo)

filas_t9 <- bind_rows(
  tibble(Variable  = "Genotipos únicos gB, n (%)",
         Muestras  = fmt_np(n_m_gb_uni, round(100 * n_m_gb_uni / n_m_total, 1)),
         Episodios = fmt_np(n_e_gb_uni, round(100 * n_e_gb_uni / n_e_total, 1)),
         header    = TRUE),
  
  map_dfr(combos_gb, function(k) {
    tibble(Variable  = k,
           Muestras  = get_m(freq_gb_uni_m, k),
           Episodios = get_m(freq_gb_uni_e, k),
           header    = FALSE)
  }),
  
  tibble(Variable  = "Genotipos mixtos gB, n (%)",
         Muestras  = fmt_np(n_m_gb_mix, round(100 * n_m_gb_mix / n_m_total, 1)),
         Episodios = fmt_np(n_e_gb_mix, round(100 * n_e_gb_mix / n_e_total, 1)),
         header    = TRUE),
  
  map_dfr(combos_gb_mix, function(k) {
    tibble(Variable  = k,
           Muestras  = get_m(freq_gb_mix_m, k),
           Episodios = get_m(freq_gb_mix_e, k),
           header    = FALSE)
  })
)

## ── 6. Flextable ─────────────────────────────────────────────────────────
idx_headers <- which(filas_t9$header)
idx_indent  <- which(!filas_t9$header)

ft_t9 <- filas_t9 %>%
  select(Variable, Muestras, Episodios) %>%
  flextable() %>%
  theme_booktabs() %>%
  set_header_labels(
    Variable  = "Genotipo gB",
    Muestras  = paste0("Muestras\n(n = ", n_m_total, ")"),
    Episodios = paste0("Episodios\n(n = ", n_e_total, ")")
  ) %>%
  align(j = "Variable", align = "left", part = "all") %>%
  align(j = c("Muestras", "Episodios"), align = "center", part = "all") %>%
  bold(i = idx_headers, j = "Variable", part = "body") %>%
  bg(i = idx_headers, bg = "#F0F0F0", part = "body") %>%
  padding(i = idx_indent, j = "Variable", padding.left = 14, part = "body") %>%
  bold(part = "header") %>%
  bg(bg = "#F0F0F0", part = "header") %>%
  add_footer_lines(values = c(
    "Los porcentajes de los genotipos únicos gB1 a gB4 se calcularon respecto al total de genotipos únicos.",
    "Los porcentajes de las combinaciones mixtas se calcularon respecto al total de genotipos mixtos."
  )) %>%
  fontsize(part = "footer", size = 8) %>%
  align(part = "footer", align = "left") %>%
  fontsize(part = "body",   size = 10) %>%
  fontsize(part = "header", size = 10) %>%
  autofit()

ft_t9

## ── 7. Exportar a Word ───────────────────────────────────────────────────
read_docx() %>%
  body_add_par("Tabla 9. Frecuencias genotípicas de gB en muestras y en episodios.",
               style = "Normal") %>%
  body_add_flextable(ft_t9) %>%
  print(target = "~/Desktop/Figuras_tesis/Tabla_9_gB.docx")

