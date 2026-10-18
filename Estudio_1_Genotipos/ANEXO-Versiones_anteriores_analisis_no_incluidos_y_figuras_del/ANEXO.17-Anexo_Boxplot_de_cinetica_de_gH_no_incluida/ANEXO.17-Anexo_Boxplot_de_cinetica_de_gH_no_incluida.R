# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.17. Anexo: Boxplot de cinética de gH (no incluida)
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Boxplot de cinética de gH (no incluida) ----
# ── 1. Datos ──────────────────────────────────────────────────────────────────
df_gh <- dnaemias %>%
  mutate(
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini),
    cvini_log         = log10(cvini),
    cvpico_log        = log10(cvpico)
  ) %>%
  filter(ghsum == 1) %>%
  mutate(
    gh_tipo = factor(case_when(
      gh1 == 1 ~ "gH1",
      gh2 == 1 ~ "gH2"
    ), levels = c("gH1","gH2"))
  )

# ── 2. Comparaciones post-hoc significativas ──────────────────────────────────
# Sólo hay una comparación posible (gH1 vs gH2). Rellena cuando tengas
# los resultados de tus tests; en caso contrario se deja vacío.
comp_vacio <- tibble(
  group1       = character(0),
  group2       = character(0),
  p.adj.signif = character(0),
  y.position   = numeric(0)
)

# Plantilla para activar la comparación cuando proceda:
# comp_gh <- tibble(
#   group1       = "gH1",
#   group2       = "gH2",
#   p.adj.signif = "*",
#   y.position   = 175
# )

# ── 3. Colores ────────────────────────────────────────────────────────────────
cols_gh <- c(
  "gH1" = "pink",
  "gH2" = "hotpink3"
)

cols_jitter <- c(
  "gH1" = "pink3",
  "gH2" = "palevioletred4"
)

# ── 4. Función para construir cada panel ──────────────────────────────────────
make_panel <- function(data, var, ylab, comparaciones,
                       cap_vis = NULL, ylim_max = NULL) {
  
  d <- data %>%
    select(gh_tipo, valor = !!sym(var)) %>%
    filter(!is.na(valor))
  
  if (!is.null(cap_vis)) d <- d %>% filter(valor <= cap_vis)
  
  p <- ggplot(d, aes(x = gh_tipo, y = valor, fill = gh_tipo)) +
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
      aes(color = gh_tipo),
      width = 0.12, size = 0.9, alpha = 0.35,
      show.legend = FALSE
    ) +
    scale_fill_manual(values  = cols_gh) +
    scale_color_manual(values = cols_jitter) +
    labs(x = NULL, y = ylab) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position    = "none",
      panel.border       = element_blank(),
      axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
      axis.line.y        = element_line(color = "grey60", linewidth = 0.4),
      axis.ticks         = element_line(color = "grey60", linewidth = 0.3),
      axis.text.x        = element_text(size = 11, color = "grey25"),
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
      expand = expansion(mult = c(0.02, 0.15))
    )
  }
  
  if (nrow(comparaciones) > 0) {
    p <- p + stat_pvalue_manual(
      comparaciones,
      label        = "p.adj.signif",
      xmin         = "group1",
      xmax         = "group2",
      y.position   = "y.position",
      tip.length   = 0.01,
      bracket.size = 0.45,
      size         = 4.5,
      inherit.aes  = FALSE
    )
  }
  
  p
}

# ── 5. Construir los 4 paneles ────────────────────────────────────────────────
p_cvini <- make_panel(
  data          = df_gh,
  var           = "cvini_log",
  ylab          = "log10 CV inicial",
  comparaciones = comp_vacio
)

p_cvpico <- make_panel(
  data          = df_gh,
  var           = "cvpico_log",
  ylab          = "log10 CV pico",
  comparaciones = comp_vacio
)

p_aclar <- make_panel(
  data          = df_gh,
  var           = "aclaramiento_dias",
  ylab          = "Aclaramiento (días)",
  comparaciones = comp_vacio,
  cap_vis       = 200,
  ylim_max      = 200
)

p_auc <- make_panel(
  data          = df_gh,
  var           = "AUCCV",
  ylab          = "AUC viral (UI/mL\u00b7d\u00eda)",
  comparaciones = comp_vacio,
  cap_vis       = 500,
  ylim_max      = 500
)

# ── 6. Combinar 2x2 ──────────────────────────────────────────────────────────
figura_gh_cinetica <- (p_cvini | p_cvpico) / (p_aclar | p_auc)

figura_gh_cinetica

