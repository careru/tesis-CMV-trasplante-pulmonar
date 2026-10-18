# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.5.07. Figura 27. Parámetros de cinética viral entre los genotipos únicos de gB
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 27. Parámetros de cinética viral entre los genotipos únicos de gB ----
# ── 1. Datos ──────────────────────────────────────────────────────────────────
df_gb <- dnaemias %>%
  mutate(
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini),
    cvini_log         = log10(cvini),
    cvpico_log        = log10(cvpico)
  ) %>%
  filter(gbsum == 1) %>%
  mutate(
    gb_tipo = factor(case_when(
      gb1 == 1 ~ "gB1", gb2 == 1 ~ "gB2",
      gb3 == 1 ~ "gB3", gb4 == 1 ~ "gB4"
    ), levels = c("gB1","gB2","gB3","gB4"))
  )

# ── 2. Comparaciones post-hoc significativas ──────────────────────────────────
comp_aclar <- tibble(
  group1       = c("gB1", "gB2"),
  group2       = c("gB3", "gB3"),
  p.adj.signif = c("**", "**"),
  y.position   = c(160, 175)
)

comp_auc <- tibble(
  group1       = c("gB1", "gB2"),
  group2       = c("gB3", "gB3"),
  p.adj.signif = c("**", "**"),
  y.position   = c(320, 355)
)

comp_vacio <- tibble(
  group1       = character(0),
  group2       = character(0),
  p.adj.signif = character(0),
  y.position   = numeric(0)
)

# ── 3. Colores: paleta de azules equilibrada (gB3 no tan oscuro) ───────────────
cols_gb <- c(
  "gB1" = "#D6EAF8",
  "gB2" = "#85C1E9",
  "gB3" = "#3498DB",
  "gB4" = "#A9CCE3"
)

cols_jitter <- c(
  "gB1" = "#2E86C1",
  "gB2" = "#2471A3",
  "gB3" = "#1B4F72",
  "gB4" = "#5499C7"
)

# ── 4. Función para construir cada panel ──────────────────────────────────────
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
    # Mediana resaltada en blanco por encima del relleno
    stat_summary(
      fun = median, geom = "crossbar",
      width = 0.55, fatten = 0,
      color = "white", linewidth = 0.7
    ) +
    geom_jitter(
      aes(color = gb_tipo),
      width = 0.12, size = 0.9, alpha = 0.35,
      show.legend = FALSE
    ) +
    scale_fill_manual(values  = cols_gb) +
    scale_color_manual(values = cols_jitter) +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = 15) +
    theme(
      legend.position    = "none",
      axis.line          = element_line(color = "black", linewidth = 0.5),
      axis.ticks         = element_line(color = "black", linewidth = 0.4),
      axis.text.x        = element_text(size = 14, color = "black"),
      axis.text.y        = element_text(size = 12, color = "black"),
      axis.title.y       = element_text(size = 14, color = "black", face = "bold"),
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
      size         = 6,
      inherit.aes  = FALSE
    )
  }
  
  p
}

# ── 5. Construir los 4 paneles ────────────────────────────────────────────────
p_cvini <- make_panel(
  data          = df_gb,
  var           = "cvini_log",
  ylab          = expression(bold(log[10]*" CV inicial (UI/mL)")),
  comparaciones = comp_vacio
)

p_cvpico <- make_panel(
  data          = df_gb,
  var           = "cvpico_log",
  ylab          = expression(bold(log[10]*" CV pico (UI/mL)")),
  comparaciones = comp_vacio
)

p_aclar <- make_panel(
  data          = df_gb,
  var           = "aclaramiento_dias",
  ylab          = "Aclaramiento (d\u00edas)",
  comparaciones = comp_aclar,
  cap_vis       = 200,
  ylim_max      = 200
)

p_auc <- make_panel(
  data          = df_gb,
  var           = "AUCCV",
  ylab          = expression(bold("AUC (log"[10]*"(UI/mL)\u00b7d\u00eda)")),
  comparaciones = comp_auc,
  cap_vis       = 400,
  ylim_max      = 400
)

# ── 6. Combinar 2x2 ──────────────────────────────────────────────────────────
figura_gb_cinetica <- (p_cvini | p_cvpico) / (p_aclar | p_auc)

figura_gb_cinetica


# ── 7. Guardar en SVG ─────────────────────────────────────────────────────────
# [guardado antiguo] ggsave(
# [guardado antiguo]   filename      = "~/Desktop/Figuras_tesis/Figura_27_cinetica_gB.svg",
# [guardado antiguo]   plot          = figura_gb_cinetica,
# [guardado antiguo]   device        = svglite::svglite,
# [guardado antiguo]   fix_text_size = FALSE,
# [guardado antiguo]   width         = 16,
# [guardado antiguo]   height        = 15,
# [guardado antiguo]   units         = "cm",
# [guardado antiguo]   bg            = "white"
# [guardado antiguo] )

# ── Preparar datos ────────────────────────────────────────────────────────────
df_gb <- dnaemias %>%
  mutate(
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini)
  ) %>%
  filter(gbsum == 1) %>%
  mutate(gb_tipo = case_when(
    gb1 == 1 ~ "gB1", gb2 == 1 ~ "gB2",
    gb3 == 1 ~ "gB3", gb4 == 1 ~ "gB4"
  ))

# ── Post-hoc Dunn — Aclaramiento ──────────────────────────────────────────────
dunn.test(
  x      = df_gb$aclaramiento_dias,
  g      = df_gb$gb_tipo,
  method = "bonferroni",
  kw     = TRUE,
  label  = TRUE
)

# ── Post-hoc Dunn — AUC ───────────────────────────────────────────────────────
dunn.test(
  x      = df_gb$AUCCV,
  g      = df_gb$gb_tipo,
  method = "bonferroni",
  kw     = TRUE,
  label  = TRUE
)

guardar_fig(figura_gb_cinetica, "Figura_27", ancho = 16, alto = 15)

