# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.28. Anexo: Artículo: Figura X3, KM de recurrencia
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Artículo: Figura X3, KM de recurrencia ----
# ── Format p-value (point decimal for English text) ─────────────────────────
fmt_num_point <- function(p) {
  formatC(p, format = "f", digits = 3)
}

# ── Function to build a consistent KM panel (B/W, line type by stratum) ─────
panel_km <- function(km_fit, p_logrank, titulo, etiqueta_legend,
                     xlim_max, x_breaks, mostrar_y = TRUE,
                     mostrar_legend = FALSE,
                     pval_x = 0, pval_y = 0.10) {
  
  ggsurv <- ggsurvplot(
    km_fit,
    pval              = FALSE,
    conf.int          = FALSE,              # no CI for clean B/W
    censor            = FALSE,              # no censoring marks
    linetype          = "strata",           # single = solid, mixed = dashed
    palette           = c("black", "black"),# both curves black
    size              = 0.6,
    legend.title      = "",
    legend.labs       = etiqueta_legend,
    xlab              = "Months since clearance",
    ylab              = if (mostrar_y) "Recurrence-free probability" else "",
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
             label = paste0("italic('p')*' = '*'", fmt_num_point(p_logrank), "'"),
             parse = TRUE,
             hjust = 0, size = 3.4,
             color = "black") +
    theme(
      legend.position    = if (mostrar_legend) "bottom" else "none",
      plot.title         = element_text(size = 11, face = "bold", hjust = 0.5),
      panel.grid         = element_blank(),           # no grid (journal style)
      panel.background   = element_rect(fill = "white", color = NA),
      plot.background    = element_rect(fill = "white", color = NA),
      axis.line          = element_line(color = "black", linewidth = 0.4),
      axis.ticks         = element_line(color = "black", linewidth = 0.3),
      plot.margin        = margin(8, 8, 8, 8)
    )
  
  p
}

# ── Fit KM models at 12 months ──────────────────────────────────────────────
km_gbgh_12m <- survfit(surv_365 ~ exp_grupo, data = df_recur)
km_gb_12m   <- survfit(surv_365 ~ gb_grupo,  data = df_recur)
km_gh_12m   <- survfit(surv_365 ~ gh_grupo,  data = df_recur)

# ── Log-rank p-values ───────────────────────────────────────────────────────
p_lr <- function(formula, data) {
  ld <- survdiff(formula, data = data)
  1 - pchisq(ld$chisq, df = length(ld$n) - 1)
}

p_gbgh_12m <- p_lr(surv_365 ~ exp_grupo, df_recur)
p_gb_12m   <- p_lr(surv_365 ~ gb_grupo,  df_recur)
p_gh_12m   <- p_lr(surv_365 ~ gh_grupo,  df_recur)

# ── Build the 3 panels (no legend; described in figure caption) ─────────────
pval_x_factor <- 0.1
pval_x_A      <- 365 * pval_x_factor

A1 <- panel_km(km_gbgh_12m, p_gbgh_12m, "A. gB/gH",
               c("Single gB/gH", "Mixed gB/gH"), 365, 60,
               mostrar_y = TRUE, mostrar_legend = FALSE,
               pval_x = pval_x_A, pval_y = 0.10)
A2 <- panel_km(km_gb_12m,   p_gb_12m,   "B. gB",
               c("Single gB", "Mixed gB"), 365, 60,
               mostrar_y = FALSE, mostrar_legend = FALSE,
               pval_x = pval_x_A, pval_y = 0.10)
A3 <- panel_km(km_gh_12m,   p_gh_12m,   "C. gH",
               c("Single gH", "Mixed gH"), 365, 60,
               mostrar_y = FALSE, mostrar_legend = FALSE,
               pval_x = pval_x_A, pval_y = 0.10)

# ── 1. Generar paneles (xlim 366 para que quepa la marca de 12) ─────────────
pval_x_factor <- 0.1
pval_x_A      <- 365 * pval_x_factor

A1 <- panel_km(km_gbgh_12m, p_gbgh_12m, "A. gB/gH",
               c("Single gB/gH", "Mixed gB/gH"), 366, 60,
               mostrar_y = TRUE, mostrar_legend = FALSE,
               pval_x = pval_x_A, pval_y = 0.10)
A2 <- panel_km(km_gb_12m,   p_gb_12m,   "B. gB",
               c("Single gB", "Mixed gB"), 366, 60,
               mostrar_y = FALSE, mostrar_legend = FALSE,
               pval_x = pval_x_A, pval_y = 0.10)
A3 <- panel_km(km_gh_12m,   p_gh_12m,   "C. gH",
               c("Single gH", "Mixed gH"), 366, 60,
               mostrar_y = FALSE, mostrar_legend = FALSE,
               pval_x = pval_x_A, pval_y = 0.10)

# ── 2. Convertir el eje X a meses (DESPUÉS de generar los paneles) ──────────
meses_breaks <- seq(0, 12, 2)
dias_breaks  <- meses_breaks * 30.44

A1 <- A1 + scale_x_continuous(
  breaks = dias_breaks, labels = meses_breaks,
  limits = c(0, 366), expand = c(0, 0)) +
  labs(x = "Months since clearance")
A2 <- A2 + scale_x_continuous(
  breaks = dias_breaks, labels = meses_breaks,
  limits = c(0, 366), expand = c(0, 0)) +
  labs(x = "Months since clearance")
A3 <- A3 + scale_x_continuous(
  breaks = dias_breaks, labels = meses_breaks,
  limits = c(0, 366), expand = c(0, 0)) +
  labs(x = "Months since clearance")


# ── 1. Paneles sin leyenda (los A1/A2/A3 que ya tienes, en meses) ───────────
fila_paneles <- plot_grid(A1, A2, A3, nrow = 1)

# ── 2. Construir una leyenda genérica común (continua = único, discont. = mixto)
df_leyenda <- data.frame(
  x = c(1, 1, 2, 2),
  y = c(1, 2, 1, 2),
  Genotype = factor(rep(c("Single genotype", "Mixed genotype"), each = 2),
                    levels = c("Single genotype", "Mixed genotype"))
)

p_leyenda <- ggplot(df_leyenda, aes(x, y, linetype = Genotype)) +
  geom_line(linewidth = 0.6, color = "black") +
  scale_linetype_manual(values = c("Single genotype" = "solid",
                                   "Mixed genotype"   = "dashed")) +
  theme_void(base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 10, color = "black"),
    legend.key.width = unit(1.2, "cm")
  )

leyenda_comun <- get_legend(p_leyenda)

# ── 3. Componer: fila de paneles arriba, leyenda común abajo ────────────────
figura_X3 <- plot_grid(
  fila_paneles,
  leyenda_comun,
  ncol        = 1,
  rel_heights = c(1, 0.10)
)

print(figura_X3)

ggsave("~/Desktop/Figuras_tesis/figure_X3_recurrence_12m.png",
       figura_X3,
       width = 30, height = 11, units = "cm",
       dpi = 300, bg = "white")

