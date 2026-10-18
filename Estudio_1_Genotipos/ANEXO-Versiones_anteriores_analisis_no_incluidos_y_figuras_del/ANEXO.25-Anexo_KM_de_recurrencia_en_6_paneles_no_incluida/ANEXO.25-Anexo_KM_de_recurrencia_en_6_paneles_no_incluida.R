# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.25. Anexo: KM de recurrencia en 6 paneles (no incluida)
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: KM de recurrencia en 6 paneles (no incluida) ----
# ── Función para formatear el número con coma decimal ───────────────────────
fmt_num_coma <- function(p) {
  sub("\\.", ",", formatC(p, format = "f", digits = 3))
}

# ── Función para generar un panel KM consistente ────────────────────────────
panel_km <- function(km_fit, p_logrank, titulo, etiqueta_legend,
                     xlim_max, x_breaks, mostrar_y = TRUE,
                     mostrar_legend = FALSE,
                     pval_x = 0, pval_y = 0.10,
                     vline_at = NULL) {
  
  ggsurv <- ggsurvplot(
    km_fit,
    pval              = FALSE,
    conf.int          = TRUE,
    conf.int.alpha    = 0.15,
    conf.int.style    = "ribbon",
    risk.table        = FALSE,
    palette           = c("#3498db", "#e74c3c"),
    legend.title      = "",
    legend.labs       = etiqueta_legend,
    xlab              = "Días desde aclaramiento",
    ylab              = if (mostrar_y) "Probabilidad libre de recurrencia" else "",
    title             = titulo,
    ggtheme           = theme_classic(base_size = 11),
    font.legend       = c(9, "plain", "black"),
    font.x            = c(10, "plain", "black"),
    font.y            = c(10, "plain", "black"),
    font.tickslab     = c(9, "plain", "black"),
    break.time.by     = x_breaks,
    xlim              = c(0, xlim_max),
    ylim              = c(0, 1)
  )
  
  p <- ggsurv$plot +
    annotate("text",
             x = pval_x, y = pval_y,
             label = paste0("italic('p')*' = '*'", fmt_num_coma(p_logrank), "'"),
             parse = TRUE,
             hjust = 0, size = 3.4,
             color = "black") +
    theme(
      legend.position    = if (mostrar_legend) "bottom" else "none",
      plot.title         = element_text(size = 11, face = "bold",
                                        hjust = 0.5),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.background   = element_rect(fill = "white", color = NA),
      plot.background    = element_rect(fill = "white", color = NA),
      axis.line          = element_line(color = "grey40", linewidth = 0.4),
      axis.ticks         = element_line(color = "grey40", linewidth = 0.3),
      plot.margin        = margin(8, 8, 8, 8)
    )
  
  if (!is.null(vline_at)) {
    p <- p + geom_vline(xintercept = vline_at,
                        linetype  = "dashed",
                        color     = "grey50",
                        linewidth = 0.4)
  }
  
  p
}

# ── Preparar ajustes KM ─────────────────────────────────────────────────────
# Fila A: 12 meses
km_gbgh_12m <- survfit(surv_365 ~ exp_grupo, data = df_recur)
km_gb_12m   <- survfit(surv_365 ~ gb_grupo,  data = df_recur)
km_gh_12m   <- survfit(surv_365 ~ gh_grupo,  data = df_recur)

# Fila B: seguimiento completo
km_gbgh_glob <- survfit(surv_completo ~ exp_grupo, data = df_recur)
km_gb_glob   <- survfit(surv_completo ~ gb_grupo,  data = df_recur)
km_gh_glob   <- survfit(surv_completo ~ gh_grupo,  data = df_recur)

# ── Calcular p de log-rank ──────────────────────────────────────────────────
p_lr <- function(formula, data) {
  ld <- survdiff(formula, data = data)
  1 - pchisq(ld$chisq, df = length(ld$n) - 1)
}

p_gbgh_12m  <- p_lr(surv_365 ~ exp_grupo, df_recur)
p_gb_12m    <- p_lr(surv_365 ~ gb_grupo,  df_recur)
p_gh_12m    <- p_lr(surv_365 ~ gh_grupo,  df_recur)
p_gbgh_glob <- p_lr(surv_completo ~ exp_grupo, df_recur)
p_gb_glob   <- p_lr(surv_completo ~ gb_grupo,  df_recur)
p_gh_glob   <- p_lr(surv_completo ~ gh_grupo,  df_recur)

# ── Determinar xlim del seguimiento completo ────────────────────────────────
xlim_completo     <- max(df_recur$tiempo_completo, na.rm = TRUE)
xlim_completo_red <- ceiling(xlim_completo / 365) * 365

# ── Generar los 6 paneles ───────────────────────────────────────────────────
pval_x_factor <- 0.1
pval_x_A      <- 365 * pval_x_factor
pval_x_B      <- xlim_completo_red * pval_x_factor

# Fila A: 12 meses, línea discontinua a los 6 meses
A1 <- panel_km(km_gbgh_12m, p_gbgh_12m, "A1. gB+gH (12 meses)",
               c("gB+gH únicos", "gB+gH mixtos"), 365, 60,
               mostrar_y = TRUE,
               pval_x = pval_x_A, pval_y = 0.10,
               vline_at = 6 * 30.44)
A2 <- panel_km(km_gb_12m,   p_gb_12m,   "A2. gB (12 meses)",
               c("gB únicos", "gB mixtos"), 365, 60,
               mostrar_y = FALSE,
               pval_x = pval_x_A, pval_y = 0.10,
               vline_at = 6 * 30.44)
A3 <- panel_km(km_gh_12m,   p_gh_12m,   "A3. gH (12 meses)",
               c("gH únicos", "gH mixtos"), 365, 60,
               mostrar_y = FALSE,
               pval_x = pval_x_A, pval_y = 0.10,
               vline_at = 6 * 30.44)

# Fila B: seguimiento completo
B1 <- panel_km(km_gbgh_glob, p_gbgh_glob, "B1. gB+gH (seguimiento completo)",
               c("gB+gH únicos", "gB+gH mixtos"), xlim_completo_red, 365,
               mostrar_y = TRUE, mostrar_legend = TRUE,
               pval_x = pval_x_B, pval_y = 0.10)
B2 <- panel_km(km_gb_glob,   p_gb_glob,   "B2. gB (seguimiento completo)",
               c("gB únicos", "gB mixtos"), xlim_completo_red, 365,
               mostrar_y = FALSE, mostrar_legend = TRUE,
               pval_x = pval_x_B, pval_y = 0.10)
B3 <- panel_km(km_gh_glob,   p_gh_glob,   "B3. gH (seguimiento completo)",
               c("gH únicos", "gH mixtos"), xlim_completo_red, 365,
               mostrar_y = FALSE, mostrar_legend = TRUE,
               pval_x = pval_x_B, pval_y = 0.10)

# ── Eje X de fila A en meses redondos ───────────────────────────────────────
meses_breaks <- seq(0, 12, 2)
dias_breaks  <- meses_breaks * 30.44

A1 <- A1 + scale_x_continuous(
  breaks = dias_breaks,
  labels = meses_breaks,
  limits = c(0, 365),
  expand = c(0, 0)) +
  labs(x = "Meses desde aclaramiento")
A2 <- A2 + scale_x_continuous(
  breaks = dias_breaks,
  labels = meses_breaks,
  limits = c(0, 365),
  expand = c(0, 0)) +
  labs(x = "Meses desde aclaramiento")
A3 <- A3 + scale_x_continuous(
  breaks = dias_breaks,
  labels = meses_breaks,
  limits = c(0, 365),
  expand = c(0, 0)) +
  labs(x = "Meses desde aclaramiento")

# ── Eje X de fila B en años ─────────────────────────────────────────────────
B1 <- B1 + scale_x_continuous(
  breaks = seq(0, xlim_completo_red, 365),
  labels = seq(0, xlim_completo_red / 365, 1),
  limits = c(0, xlim_completo_red),
  expand = c(0, 0)) +
  labs(x = "Años desde aclaramiento")
B2 <- B2 + scale_x_continuous(
  breaks = seq(0, xlim_completo_red, 365),
  labels = seq(0, xlim_completo_red / 365, 1),
  limits = c(0, xlim_completo_red),
  expand = c(0, 0)) +
  labs(x = "Años desde aclaramiento")
B3 <- B3 + scale_x_continuous(
  breaks = seq(0, xlim_completo_red, 365),
  labels = seq(0, xlim_completo_red / 365, 1),
  limits = c(0, xlim_completo_red),
  expand = c(0, 0)) +
  labs(x = "Años desde aclaramiento")

# ── Componer figura 2x3 ─────────────────────────────────────────────────────
figura_R31bis <- (A1 | A2 | A3) / (B1 | B2 | B3) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  )

print(figura_R31bis)

ggsave("~/Desktop/Figuras_tesis/Figura_R31bis_recurrencia.png",
       figura_R31bis,
       width = 32, height = 16, units = "cm",
       dpi = 300, bg = "white")

