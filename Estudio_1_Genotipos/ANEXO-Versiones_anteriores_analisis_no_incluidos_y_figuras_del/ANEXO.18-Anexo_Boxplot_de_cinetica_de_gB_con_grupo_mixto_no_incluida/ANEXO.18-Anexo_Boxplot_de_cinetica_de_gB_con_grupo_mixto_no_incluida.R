# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.18. Anexo: Boxplot de cinética de gB con grupo mixto (no incluida)
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Boxplot de cinética de gB con grupo mixto (no incluida) ----
# ── 1. Datos ──────────────────────────────────────────────────────────────────
df_gb <- dnaemias %>%
  mutate(
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini),
    cvini_log         = log10(cvini),
    cvpico_log        = log10(cvpico)
  ) %>%
  filter(gbsum >= 1) %>%
  mutate(
    gb_tipo = factor(case_when(
      gbsum >  1 ~ "gB mixto",
      gb1   == 1 ~ "gB1",
      gb2   == 1 ~ "gB2",
      gb3   == 1 ~ "gB3",
      gb4   == 1 ~ "gB4"
    ), levels = c("gB1","gB2","gB3","gB4","gB mixto"))
  )

# ── 2. Colores ────────────────────────────────────────────────────────────────
cols_gb <- c(
  "gB1"      = "#AED6F1",
  "gB2"      = "#5DADE2",
  "gB3"      = "#1A5276",
  "gB4"      = "#85929E",
  "gB mixto" = "#D7BDE2"
)

cols_jitter <- c(
  "gB1"      = "#2980B9",
  "gB2"      = "#1A6FA8",
  "gB3"      = "#0D2B45",
  "gB4"      = "#555F66",
  "gB mixto" = "#7D3C98"
)

col_mixto <- "#7D3C98"   # color para asteriscos/corchetes con gB mixto
col_otros <- "grey20"    # color para el resto

# ── 3. Tests Kruskal-Wallis (informativos) ────────────────────────────────────
cat("\n--- Kruskal-Wallis ---\n")
print(df_gb %>% kruskal_test(cvini_log         ~ gb_tipo))
print(df_gb %>% kruskal_test(cvpico_log        ~ gb_tipo))
print(df_gb %>% filter(aclaramiento_dias <= 200) %>%
        kruskal_test(aclaramiento_dias   ~ gb_tipo))
print(df_gb %>% filter(AUCCV <= 500) %>%
        kruskal_test(AUCCV               ~ gb_tipo))

# ── 4. Función para calcular comparaciones post-hoc ──────────────────────────
calcular_comp <- function(data, var, y_base, y_step) {
  d <- data %>%
    select(gb_tipo, valor = !!sym(var)) %>%
    filter(!is.na(valor))
  
  res <- d %>%
    dunn_test(valor ~ gb_tipo, p.adjust.method = "BH") %>%
    filter(p.adj < 0.05) %>%
    mutate(
      p.adj.signif = case_when(
        p.adj < 0.001 ~ "***",
        p.adj < 0.01  ~ "**",
        p.adj < 0.05  ~ "*",
        TRUE          ~ "ns"
      ),
      involucra_mixto = group1 == "gB mixto" | group2 == "gB mixto"
    ) %>%
    arrange(involucra_mixto, p.adj) %>%       # primero las no-mixto, luego mixto
    mutate(y.position = y_base + (row_number() - 1) * y_step) %>%
    select(group1, group2, p.adj.signif, y.position, involucra_mixto)
  
  res
}

# Calcular comparaciones para cada variable
comp_cvini  <- calcular_comp(df_gb, "cvini_log",
                             y_base =   6.8, y_step = 0.35)
comp_cvpico <- calcular_comp(df_gb, "cvpico_log",
                             y_base =   7.0, y_step = 0.35)
comp_aclar  <- calcular_comp(df_gb %>% filter(aclaramiento_dias <= 200),
                             "aclaramiento_dias",
                             y_base = 155, y_step = 13)
comp_auc    <- calcular_comp(df_gb %>% filter(AUCCV <= 500),
                             "AUCCV",
                             y_base = 380, y_step = 32)

cat("\n--- Comparaciones significativas Dunn (BH) ---\n")
print(list(cvini = comp_cvini, cvpico = comp_cvpico,
           aclar = comp_aclar,    auc    = comp_auc))

# ── 5. Función para construir cada panel ──────────────────────────────────────
make_panel <- function(data, var, ylab, comparaciones,
                       cap_vis = NULL, ylim_max = NULL) {
  
  d <- data %>%
    select(gb_tipo, valor = !!sym(var)) %>%
    filter(!is.na(valor))
  
  if (!is.null(cap_vis)) d <- d %>% filter(valor <= cap_vis)
  
  p <- ggplot(d, aes(x = gb_tipo, y = valor, fill = gb_tipo)) +
    geom_boxplot(
      width         = 0.55,
      color         = "grey25",
      linewidth     = 0.45,
      outlier.shape = 21,
      outlier.size  = 1.8,
      outlier.fill  = "white",
      outlier.color = "grey50",
      outlier.alpha = 0.7
    ) +
    geom_jitter(
      aes(color = gb_tipo),
      width = 0.12, size = 0.9, alpha = 0.35,
      show.legend = FALSE
    ) +
    scale_fill_manual(values  = cols_gb) +
    scale_color_manual(values = cols_jitter) +
    labs(x = NULL, y = ylab) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position    = "none",
      panel.border       = element_blank(),
      axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
      axis.line.y        = element_line(color = "grey60", linewidth = 0.4),
      axis.ticks         = element_line(color = "grey60", linewidth = 0.3),
      axis.text.x        = element_text(size = 10, color = "grey25"),
      axis.text.y        = element_text(size = 10, color = "grey30"),
      axis.title.y       = element_text(size = 11, color = "grey20"),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.margin        = margin(8, 12, 6, 6)
    )
  
  if (!is.null(ylim_max)) {
    p <- p + scale_y_continuous(
      name   = ylab,
      limits = c(NA, ylim_max),
      expand = expansion(mult = c(0.02, 0.05))
    )
  } else {
    p <- p + scale_y_continuous(
      name   = ylab,
      expand = expansion(mult = c(0.02, 0.20))
    )
  }
  
  # Separar comparaciones por color y añadir cada subconjunto en una capa
  if (nrow(comparaciones) > 0) {
    
    comp_otros <- comparaciones %>% filter(!involucra_mixto)
    comp_mix   <- comparaciones %>% filter( involucra_mixto)
    
    if (nrow(comp_otros) > 0) {
      p <- p + stat_pvalue_manual(
        comp_otros,
        label        = "p.adj.signif",
        xmin         = "group1",
        xmax         = "group2",
        y.position   = "y.position",
        tip.length   = 0.01,
        bracket.size = 0.45,
        size         = 4.5,
        color        = col_otros,
        inherit.aes  = FALSE
      )
    }
    
    if (nrow(comp_mix) > 0) {
      p <- p + stat_pvalue_manual(
        comp_mix,
        label        = "p.adj.signif",
        xmin         = "group1",
        xmax         = "group2",
        y.position   = "y.position",
        tip.length   = 0.01,
        bracket.size = 0.45,
        size         = 4.5,
        color        = col_mixto,
        inherit.aes  = FALSE
      )
    }
  }
  
  p
}

# ── 6. Construir los 4 paneles ────────────────────────────────────────────────
p_cvini <- make_panel(
  data          = df_gb,
  var           = "cvini_log",
  ylab          = "log10 CV inicial",
  comparaciones = comp_cvini
)

p_cvpico <- make_panel(
  data          = df_gb,
  var           = "cvpico_log",
  ylab          = "log10 CV pico",
  comparaciones = comp_cvpico
)

p_aclar <- make_panel(
  data          = df_gb,
  var           = "aclaramiento_dias",
  ylab          = "Aclaramiento (días)",
  comparaciones = comp_aclar,
  cap_vis       = 200,
  ylim_max      = 215
)

p_auc <- make_panel(
  data          = df_gb,
  var           = "AUCCV",
  ylab          = "AUC viral (UI/mL\u00b7d\u00eda)",
  comparaciones = comp_auc,
  cap_vis       = 500,
  ylim_max      = 520
)

# ── 7. Combinar 2x2 ───────────────────────────────────────────────────────────
figura_gb_cinetica <- (p_cvini | p_cvpico) / (p_aclar | p_auc)

figura_gb_cinetica

