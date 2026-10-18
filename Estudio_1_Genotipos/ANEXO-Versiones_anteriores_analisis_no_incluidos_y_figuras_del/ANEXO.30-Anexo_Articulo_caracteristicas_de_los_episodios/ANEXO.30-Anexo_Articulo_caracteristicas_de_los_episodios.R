# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.30. Anexo: Artículo: características de los episodios
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Artículo: características de los episodios ----
## ── Data preparation ──────────────────────────────────────────────────────
dnaemias_t <- dnaemias %>%
  mutate(
    Fini      = as.Date(Fini),
    Ffin      = as.Date(Ffin),
    FultCVpos = as.Date(FultCVpos),
    aclaramiento_dias = as.numeric(Ffin - Fini),
    no_tratado  = case_when(tto == 0 ~ 1L, tto == 1 ~ 0L, TRUE ~ NA_integer_),
    sintomatico = case_when(clinica == "VA" ~ 0L, !is.na(clinica) ~ 1L, TRUE ~ NA_integer_),
    alt_pp_bin  = if_else(!is.na(alt_pp) & alt_pp == 1, 1L, 0L)
  )

grupos <- c("Global", "D-/R+", "D+/R+", "D+/R-")

## ── Helper functions (point decimal for English) ─────────────────────────
fmt_iqr <- function(data, var, grp) {
  d <- if (grp == "Global") data else filter(data, serDR == grp)
  d <- filter(d, !is.na(.data[[var]]))
  if (nrow(d) == 0) return("\u2013")
  med <- median(d[[var]]); q1 <- quantile(d[[var]], 0.25); q3 <- quantile(d[[var]], 0.75)
  sprintf("%.1f (%.1f\u2013%.1f)", med, q1, q3)
}

fmt_npct <- function(data, var, grp) {
  d <- if (grp == "Global") data else filter(data, serDR == grp)
  d <- filter(d, !is.na(.data[[var]]))
  if (nrow(d) == 0) return("\u2013")
  n <- sum(d[[var]] == 1); tot <- nrow(d)
  if (n == 0) return("\u2013")
  sprintf("%d (%.1f)", n, 100 * n / tot)
}

fmt_p <- function(p) {
  if (is.na(p) || is.null(p)) return("")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}

base_pac <- dnaemias_t %>% select(NHC, serDR, nepi_gtc) %>% distinct()

fmt_iqr_pac <- function(var, grp) {
  d <- if (grp == "Global") base_pac else filter(base_pac, serDR == grp)
  d <- filter(d, !is.na(.data[[var]]))
  med <- median(d[[var]]); q1 <- quantile(d[[var]], 0.25); q3 <- quantile(d[[var]], 0.75)
  sprintf("%.1f (%.1f\u2013%.1f)", med, q1, q3)
}

fmt_mas1 <- function(grp) {
  d <- if (grp == "Global") dnaemias_t else filter(dnaemias_t, serDR == grp)
  n <- sum(d$PCR_gt > 1, na.rm = TRUE); tot <- nrow(d)
  sprintf("%d (%.1f)", n, 100 * n / tot)
}

fmt_altpp <- function(grp) {
  if (grp == "D-/R+") return("\u2013")
  d <- if (grp == "Global") dnaemias_t else filter(dnaemias_t, serDR == grp)
  n <- sum(d$alt_pp_bin, na.rm = TRUE); tot <- nrow(d)
  if (n == 0) return("\u2013")
  sprintf("%d (%.1f)", n, 100 * n / tot)
}

fmt_resist <- function(grp) {
  d <- if (grp == "Global") dnaemias_t else filter(dnaemias_t, serDR == grp)
  n <- sum(d$CMV_R == 1, na.rm = TRUE); tot <- nrow(d)
  if (n == 0) return("\u2013")
  sprintf("%d (%.1f)", n, 100 * n / tot)
}

fmt_notrat <- function(grp) {
  d <- if (grp == "Global") dnaemias_t else filter(dnaemias_t, serDR == grp)
  d <- filter(d, !is.na(no_tratado))
  n <- sum(d$no_tratado == 1); tot <- nrow(d)
  if (n == 0) return("\u2013")
  sprintf("%d (%.1f)", n, 100 * n / tot)
}

## ── p-values ──────────────────────────────────────────────────────────────
calc_p_cat <- function(data, var) {
  tab <- data %>%
    filter(!is.na(serDR), !is.na(.data[[var]]), serDR != "Global") %>%
    group_by(serDR) %>%
    summarise(ev = sum(.data[[var]] == 1), no = sum(.data[[var]] == 0), .groups = "drop") %>%
    column_to_rownames("serDR") %>% as.matrix()
  tryCatch(chisq.test(tab)$p.value, error = function(e) fisher.test(tab)$p.value)
}

calc_p_cont <- function(data, var) {
  d <- filter(data, !is.na(.data[[var]]), !is.na(serDR))
  tryCatch(kruskal.test(as.formula(paste(var, "~ serDR")), data = d)$p.value,
           error = function(e) NA_real_)
}

p_epi_pac  <- kruskal.test(nepi_gtc ~ serDR, data = base_pac)$p.value
p_pcr_gt   <- calc_p_cont(dnaemias_t, "PCR_gt")
p_mas1_pcr <- calc_p_cat(dnaemias_t %>% mutate(mas1 = as.integer(PCR_gt > 1)), "mas1")
p_primer   <- calc_p_cat(dnaemias_t, "primer_mil")

p_altpp <- {
  d <- filter(dnaemias_t, serDR %in% c("D+/R+", "D+/R-"))
  tab <- d %>% group_by(serDR) %>%
    summarise(ev = sum(alt_pp_bin), no = n() - sum(alt_pp_bin), .groups = "drop") %>%
    column_to_rownames("serDR") %>% as.matrix()
  tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_)
}

p_sint   <- calc_p_cat(dnaemias_t, "sintomatico")
p_ingr   <- calc_p_cat(dnaemias_t, "ingr")
p_ingcmv <- calc_p_cat(dnaemias_t, "ing_cmv")

p_notrat <- {
  d <- filter(dnaemias_t, !is.na(no_tratado), serDR %in% c("D-/R+", "D+/R+"))
  tab <- d %>% group_by(serDR) %>%
    summarise(ev = sum(no_tratado), no = n() - sum(no_tratado), .groups = "drop") %>%
    column_to_rownames("serDR") %>% as.matrix()
  tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_)
}

p_cvini  <- calc_p_cont(dnaemias_t, "cvini")
p_cvpico <- calc_p_cont(dnaemias_t, "cvpico")
p_aclar  <- calc_p_cont(filter(dnaemias_t, !is.na(aclaramiento_dias)), "aclaramiento_dias")
p_auc    <- calc_p_cont(filter(dnaemias_t, !is.na(AUCCV)), "AUCCV")

## ── Build rows ─────────────────────────────────────────────────────────────
make_row <- function(variable, vals, p_val, es_header = FALSE) {
  tibble(
    Variable = variable,
    Global   = vals[1], `D-/R+` = vals[2], `D+/R+` = vals[3], `D+/R-` = vals[4],
    p        = if (es_header) "" else fmt_p(p_val),
    header   = es_header,
    p_raw    = if (es_header) NA_real_ else p_val
  )
}

tabla_r2 <- bind_rows(
  
  make_row("Episode and sample distribution", rep("", 4), NA, es_header = TRUE),
  make_row("Episodes per patient, median (IQR)",
           sapply(grupos, fmt_iqr_pac, var = "nepi_gtc"), p_epi_pac),
  make_row("Genotyped samples per episode, median (IQR)",
           sapply(grupos, function(g) fmt_iqr(dnaemias_t, "PCR_gt", g)), p_pcr_gt),
  make_row("Episodes with >1 genotyped sample, n (%)",
           sapply(grupos, fmt_mas1), p_mas1_pcr),
  
  make_row("Clinical characteristics", rep("", 4), NA, es_header = TRUE),
  make_row("First episode with VL \u2265 1,000 IU/mL, n (%)",
           sapply(grupos, function(g) fmt_npct(dnaemias_t, "primer_mil", g)), p_primer),
  make_row("Early discontinuation of primary prophylaxis, n (%)",
           sapply(grupos, fmt_altpp), p_altpp),
  make_row("Symptomatic episodes, n (%)",
           sapply(grupos, function(g) fmt_npct(dnaemias_t, "sintomatico", g)), p_sint),
  make_row("Hospitalisation during the episode, n (%)",
           sapply(grupos, function(g) fmt_npct(dnaemias_t, "ingr", g)), p_ingr),
  make_row("CMV-attributable hospitalisation, n (%)",
           sapply(grupos, function(g) fmt_npct(dnaemias_t, "ing_cmv", g)), p_ingcmv),
  make_row("Demonstrated antiviral resistance, n (%)",
           sapply(grupos, fmt_resist), NA_real_),
  make_row("Episode not treated with antivirals, n (%)",
           sapply(grupos, fmt_notrat), p_notrat),
  
  make_row("Viral kinetics", rep("", 4), NA, es_header = TRUE),
  make_row("Initial VL (IU/mL), median (IQR)",
           sapply(grupos, function(g) fmt_iqr(dnaemias_t, "cvini", g)), p_cvini),
  make_row("Peak VL (IU/mL), median (IQR)",
           sapply(grupos, function(g) fmt_iqr(dnaemias_t, "cvpico", g)), p_cvpico),
  make_row("Time to clearance (days), median (IQR)",
           sapply(grupos, function(g) fmt_iqr(filter(dnaemias_t, !is.na(aclaramiento_dias)),
                                              "aclaramiento_dias", g)), p_aclar),
  make_row("Cumulative viral exposure (AUC), median (IQR)",
           sapply(grupos, function(g) fmt_iqr(filter(dnaemias_t, !is.na(AUCCV)), "AUCCV", g)), p_auc)
)

## ── Flextable ─────────────────────────────────────────────────────────────
idx_headers <- which(tabla_r2$header)
idx_sig     <- which(!is.na(tabla_r2$p_raw) & tabla_r2$p_raw < 0.05)

ft_r2 <- tabla_r2 %>%
  select(Variable, Global, `D-/R+`, `D+/R+`, `D+/R-`, p) %>%
  flextable() %>%
  theme_booktabs() %>%
  set_header_labels(
    Variable = "Characteristic",
    Global   = "Overall\nn = 234",
    `D-/R+`  = "D\u2212/R+\nn = 21",
    `D+/R+`  = "D+/R+\nn = 176",
    `D+/R-`  = "D+/R\u2212\nn = 37",
    p        = "p"
  ) %>%
  align(j = "Variable", align = "left", part = "all") %>%
  align(j = c("Global","D-/R+","D+/R+","D+/R-","p"), align = "center", part = "all") %>%
  bold(i = idx_headers, j = "Variable", part = "body") %>%
  bg(i = idx_headers, bg = "#F2F2F2", part = "body") %>%
  bold(i = idx_sig, j = "p", part = "body") %>%
  padding(i = setdiff(seq_len(nrow(tabla_r2)), idx_headers),
          j = "Variable", padding.left = 14, part = "body") %>%
  bold(part = "header") %>%
  bg(bg = "#F2F2F2", part = "header") %>%
  add_footer_lines(values = c(
    "Data are median (IQR) or n (%), as appropriate.",
    "VL, viral load; AUC, area under the cumulative viral exposure curve.",
    "Between-group comparisons: Kruskal-Wallis test for continuous variables; chi-square or Fisher exact test for categorical variables.",
    "\u2013, no cases recorded in that group.",
    "p < 0.05 shown in bold."
  )) %>%
  fontsize(part = "footer", size = 8) %>%
  align(part = "footer", align = "left") %>%
  fontsize(part = "body", size = 10) %>%
  fontsize(part = "header", size = 10) %>%
  autofit()

ft_r2

