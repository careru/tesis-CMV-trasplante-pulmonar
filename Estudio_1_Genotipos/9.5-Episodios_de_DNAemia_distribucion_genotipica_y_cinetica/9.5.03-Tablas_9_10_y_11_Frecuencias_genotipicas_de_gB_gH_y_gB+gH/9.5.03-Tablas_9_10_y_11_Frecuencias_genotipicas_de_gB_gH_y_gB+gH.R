# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.5.03. Tablas 9, 10 y 11. Frecuencias genotípicas de gB, gH y gB+gH en muestras y episodios
# Extracto del script completo (ver README.md del estudio)

## ---- Tablas 9, 10 y 11. Frecuencias genotípicas de gB, gH y gB+gH en muestras y episodios ----
## ── gB muestras ──────────────────────────────────────────────────────────

df_gb_m <- muestras_sel %>%
  mutate(
    gb_sum  = rowSums(across(GB1:GB4, ~ replace_na(., 0))),
    gB_combo = apply(select(., GB1:GB4) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gB", which(r), collapse = " + ")
    })
  ) %>%
  filter(gb_sum >= 1, !is.na(gB_combo))

n_m_total   <- nrow(df_gb_m)                          # 483
n_m_gb_uni  <- sum(df_gb_m$gb_sum == 1)               # únicos
n_m_gb_mix  <- sum(df_gb_m$gb_sum > 1)                # mixtos

freq_gb_uni_m <- df_gb_m %>%
  filter(gb_sum == 1) %>%
  count(gB_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_m_gb_uni, 1))

freq_gb_mix_m <- df_gb_m %>%
  filter(gb_sum > 1) %>%
  count(gB_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_m_gb_mix, 1)) %>%
  arrange(desc(n))

## ── gH muestras ──────────────────────────────────────────────────────────

df_gh_m <- muestras_sel %>%
  mutate(
    gh_sum   = rowSums(across(GH1:GH2, ~ replace_na(., 0))),
    gH_combo = apply(select(., GH1:GH2) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gH", which(r), collapse = " + ")
    })
  ) %>%
  filter(gh_sum >= 1, !is.na(gH_combo))

n_m_gh_uni <- sum(df_gh_m$gh_sum == 1)
n_m_gh_mix <- sum(df_gh_m$gh_sum > 1)

freq_gh_uni_m <- df_gh_m %>%
  filter(gh_sum == 1) %>%
  count(gH_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_m_gh_uni, 1))

freq_gh_mix_m <- df_gh_m %>%
  filter(gh_sum > 1) %>%
  count(gH_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_m_gh_mix, 1))

## ── gB+gH muestras ───────────────────────────────────────────────────────

df_gbgh_m <- muestras_sel %>%
  mutate(
    gb_sum   = rowSums(across(GB1:GB4, ~ replace_na(., 0))),
    gh_sum   = rowSums(across(GH1:GH2, ~ replace_na(., 0))),
    gB_combo = apply(select(., GB1:GB4) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gB", which(r), collapse = " + ")
    }),
    gH_combo = apply(select(., GH1:GH2) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gH", which(r), collapse = " + ")
    }),
    combo_total = paste(gB_combo, gH_combo, sep = " / "),
    mixto = gb_sum > 1 | gh_sum > 1
  ) %>%
  filter(!is.na(gB_combo), !is.na(gH_combo))

n_m_gbgh_uni <- sum(!df_gbgh_m$mixto)
n_m_gbgh_mix <- sum(df_gbgh_m$mixto)

freq_gbgh_uni_m <- df_gbgh_m %>%
  filter(!mixto) %>%
  count(combo_total, name = "n") %>%
  mutate(pct = round(100 * n / n_m_gbgh_uni, 1)) %>%
  arrange(desc(n))

freq_gbgh_mix_m <- df_gbgh_m %>%
  filter(mixto) %>%
  count(combo_total, name = "n") %>%
  mutate(pct = round(100 * n / n_m_gbgh_mix, 1)) %>%
  arrange(desc(n))


## ── gB episodios ─────────────────────────────────────────────────────────

df_gb_e <- dnaemias %>%
  mutate(
    gB_combo = apply(select(., gb1:gb4) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gB", which(r), collapse = " + ")
    })
  ) %>%
  filter(!is.na(gB_combo))

n_e_total  <- nrow(df_gb_e)                           # 234
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

## ── gH episodios ─────────────────────────────────────────────────────────

df_gh_e <- dnaemias %>%
  mutate(
    gH_combo = apply(select(., gh1:gh2) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gH", which(r), collapse = " + ")
    })
  ) %>%
  filter(!is.na(gH_combo))

n_e_gh_uni <- sum(df_gh_e$ghsum == 1)
n_e_gh_mix <- sum(df_gh_e$ghsum > 1)

freq_gh_uni_e <- df_gh_e %>%
  filter(ghsum == 1) %>%
  count(gH_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_e_gh_uni, 1))

freq_gh_mix_e <- df_gh_e %>%
  filter(ghsum > 1) %>%
  count(gH_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_e_gh_mix, 1))

## ── gB+gH episodios ──────────────────────────────────────────────────────

df_gbgh_e <- dnaemias %>%
  mutate(
    gB_combo = apply(select(., gb1:gb4) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gB", which(r), collapse = " + ")
    }),
    gH_combo = apply(select(., gh1:gh2) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gH", which(r), collapse = " + ")
    }),
    combo_total = paste(gB_combo, gH_combo, sep = " / "),
    mixto = gbsum > 1 | ghsum > 1
  ) %>%
  filter(!is.na(gB_combo), !is.na(gH_combo))

n_e_gbgh_uni <- sum(!df_gbgh_e$mixto)
n_e_gbgh_mix <- sum(df_gbgh_e$mixto)

freq_gbgh_uni_e <- df_gbgh_e %>%
  filter(!mixto) %>%
  count(combo_total, name = "n") %>%
  mutate(pct = round(100 * n / n_e_gbgh_uni, 1)) %>%
  arrange(desc(n))

freq_gbgh_mix_e <- df_gbgh_e %>%
  filter(mixto) %>%
  count(combo_total, name = "n") %>%
  mutate(pct = round(100 * n / n_e_gbgh_mix, 1)) %>%
  arrange(desc(n))


fmt_np <- function(n, pct) sprintf("%d (%s)", n, format(pct, nsmall = 1))

hacer_ft <- function(df, titulo, nota) {
  
  idx_headers <- which(df$header)
  idx_indent  <- which(!df$header)
  
  ft <- df %>%
    select(Variable, Muestras, Episodios) %>%
    flextable() %>%
    theme_booktabs() %>%
    set_header_labels(
      Variable  = titulo,
      Muestras  = paste0("Muestras\n(n = ", n_m_total, ")"),
      Episodios = paste0("Episodios\n(n = ", n_e_total, ")")
    ) %>%
    align(j = "Variable",              align = "left",   part = "all") %>%
    align(j = c("Muestras","Episodios"), align = "center", part = "all") %>%
    bold(i = idx_headers, j = "Variable", part = "body") %>%
    bg(i = idx_headers, bg = "#F2F2F2", part = "body") %>%
    padding(i = idx_indent, j = "Variable", padding.left = 14, part = "body") %>%
    bold(part = "header") %>%
    bg(bg = "#F2F2F2", part = "header") %>%
    add_footer_lines(values = nota) %>%
    fontsize(part = "footer", size = 8) %>%
    align(part = "footer", align = "left") %>%
    fontsize(part = "body",   size = 10) %>%
    fontsize(part = "header", size = 10) %>%
    autofit()
  
  ft
}


combos_gb <- sort(union(freq_gb_uni_m$gB_combo, freq_gb_uni_e$gB_combo))
combos_gb_mix <- union(freq_gb_mix_m$gB_combo, freq_gb_mix_e$gB_combo)

get_m <- function(df, key) {
  r <- df[df[[1]] == key, ]
  if (nrow(r) == 0) return("\u2013")
  fmt_np(r$n, r$pct)
}

filas_r3 <- bind_rows(
  tibble(Variable = "Genotipos únicos gB, n (%)",
         Muestras  = fmt_np(n_m_gb_uni, round(100*n_m_gb_uni/n_m_total,1)),
         Episodios = fmt_np(n_e_gb_uni, round(100*n_e_gb_uni/n_e_total,1)),
         header = TRUE),
  
  purrr::map_dfr(combos_gb, function(k) {
    tibble(Variable  = k,
           Muestras  = get_m(freq_gb_uni_m, k),
           Episodios = get_m(freq_gb_uni_e, k),
           header    = FALSE)
  }),
  
  tibble(Variable  = "Genotipos mixtos gB, n (%)",
         Muestras  = fmt_np(n_m_gb_mix, round(100*n_m_gb_mix/n_m_total,1)),
         Episodios = fmt_np(n_e_gb_mix, round(100*n_e_gb_mix/n_e_total,1)),
         header    = TRUE),
  
  purrr::map_dfr(combos_gb_mix, function(k) {
    tibble(Variable  = k,
           Muestras  = get_m(rename(freq_gb_mix_m, gB_combo=gB_combo), k),
           Episodios = get_m(rename(freq_gb_mix_e, gB_combo=gB_combo), k),
           header    = FALSE)
  })
)

ft_r3 <- hacer_ft(
  filas_r3,
  "Genotipo gB",
  c("Los porcentajes de gB1\u2013gB4 se calculan respecto al total de genotipos únicos.",
    "Los porcentajes de genotipos mixtos se calculan respecto al total de genotipos mixtos.")
)

ft_r3


filas_r4 <- bind_rows(
  tibble(Variable  = "Genotipos únicos gH, n (%)",
         Muestras  = fmt_np(n_m_gh_uni, round(100*n_m_gh_uni/n_m_total,1)),
         Episodios = fmt_np(n_e_gh_uni, round(100*n_e_gh_uni/n_e_total,1)),
         header    = TRUE),
  
  purrr::map_dfr(sort(union(freq_gh_uni_m$gH_combo, freq_gh_uni_e$gH_combo)), function(k) {
    tibble(Variable  = k,
           Muestras  = get_m(freq_gh_uni_m, k),
           Episodios = get_m(freq_gh_uni_e, k),
           header    = FALSE)
  }),
  
  tibble(Variable  = "Genotipos mixtos gH, n (%)",
         Muestras  = fmt_np(n_m_gh_mix, round(100*n_m_gh_mix/n_m_total,1)),
         Episodios = fmt_np(n_e_gh_mix, round(100*n_e_gh_mix/n_e_total,1)),
         header    = TRUE),
  
  purrr::map_dfr(union(freq_gh_mix_m$gH_combo, freq_gh_mix_e$gH_combo), function(k) {
    tibble(Variable  = k,
           Muestras  = get_m(freq_gh_mix_m, k),
           Episodios = get_m(freq_gh_mix_e, k),
           header    = FALSE)
  })
)

ft_r4 <- hacer_ft(
  filas_r4,
  "Genotipo gH",
  c("Los porcentajes de gH1 y gH2 se calculan respecto al total de genotipos únicos.",
    "Los porcentajes de genotipos mixtos se calculan respecto al total de genotipos mixtos.")
)

ft_r4


combos_gbgh_uni  <- union(freq_gbgh_uni_m$combo_total,  freq_gbgh_uni_e$combo_total)
combos_gbgh_mix  <- union(freq_gbgh_mix_m$combo_total,  freq_gbgh_mix_e$combo_total)

## Ordenar por frecuencia en muestras (desc)
combos_gbgh_uni <- freq_gbgh_uni_m %>%
  full_join(freq_gbgh_uni_e, by = "combo_total", suffix = c("_m","_e")) %>%
  arrange(desc(replace_na(n_m, 0))) %>%
  pull(combo_total)

combos_gbgh_mix <- freq_gbgh_mix_m %>%
  full_join(freq_gbgh_mix_e, by = "combo_total", suffix = c("_m","_e")) %>%
  arrange(desc(replace_na(n_m, 0))) %>%
  pull(combo_total)

filas_r5 <- bind_rows(
  tibble(Variable  = "Genotipos únicos gB / gH, n (%)",
         Muestras  = fmt_np(n_m_gbgh_uni, round(100*n_m_gbgh_uni/n_m_total,1)),
         Episodios = fmt_np(n_e_gbgh_uni, round(100*n_e_gbgh_uni/n_e_total,1)),
         header    = TRUE),
  
  purrr::map_dfr(combos_gbgh_uni, function(k) {
    tibble(Variable  = k,
           Muestras  = get_m(freq_gbgh_uni_m, k),
           Episodios = get_m(freq_gbgh_uni_e, k),
           header    = FALSE)
  }),
  
  tibble(Variable  = "Genotipos mixtos gB / gH, n (%)",
         Muestras  = fmt_np(n_m_gbgh_mix, round(100*n_m_gbgh_mix/n_m_total,1)),
         Episodios = fmt_np(n_e_gbgh_mix, round(100*n_e_gbgh_mix/n_e_total,1)),
         header    = TRUE),
  
  purrr::map_dfr(combos_gbgh_mix, function(k) {
    tibble(Variable  = k,
           Muestras  = get_m(freq_gbgh_mix_m, k),
           Episodios = get_m(freq_gbgh_mix_e, k),
           header    = FALSE)
  })
)

ft_r5 <- hacer_ft(
  filas_r5,
  "Combinación gB / gH",
  c("Los porcentajes de genotipos únicos se calculan respecto al total de genotipos únicos.",
    "Los porcentajes de genotipos mixtos se calculan respecto al total de genotipos mixtos.")
)

ft_r3

ft_r4

ft_r5


  

guardar_tabla(ft_r3, "Tabla_09")

guardar_tabla(ft_r4, "Tabla_10")

guardar_tabla(ft_r5, "Tabla_11")

