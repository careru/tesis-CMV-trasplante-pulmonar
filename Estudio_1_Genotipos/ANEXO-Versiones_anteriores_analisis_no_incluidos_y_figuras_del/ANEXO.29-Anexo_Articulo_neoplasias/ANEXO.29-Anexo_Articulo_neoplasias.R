# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.29. Anexo: Artículo: neoplasias
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Artículo: neoplasias ----
# ── Función panel_km: estilo limpio B/N (líneas continua/discontinua) ───────
panel_km <- function(km_fit, p_logrank, titulo, etiqueta_legend,
                     xlim_max, x_breaks, mostrar_y = TRUE,
                     pval_x = 0, pval_y = 0.10) {
  
  ggsurv <- ggsurvplot(
    km_fit,
    pval              = FALSE,
    conf.int          = FALSE,              # sin IC
    censor            = FALSE,              # sin marcas de censura
    linetype          = "strata",           # único = continua, mixto = discontinua
    palette           = c("black", "black"),# ambas negras
    size              = 0.6,
    legend.title      = "",
    legend.labs       = etiqueta_legend,
    xlab              = "Years since transplant",
    ylab              = if (mostrar_y) "Malignancy-free survival" else "",
    title             = titulo,
    ggtheme           = theme_classic(base_size = 11),
    font.x            = c(10, "plain", "black"),
    font.y            = c(10, "plain", "black"),
    font.tickslab     = c(9, "plain", "black"),
    break.time.by     = x_breaks,
    xlim              = c(0, xlim_max),
    ylim              = c(0, 1)
  )
  
  ggsurv$plot +
    annotate("text",
             x = pval_x, y = pval_y,
             label = paste0("italic('p')*' = '*'",
                            formatC(p_logrank, format = "f", digits = 3), "'"),
             parse = TRUE, hjust = 0, size = 3.4, color = "black") +
    theme(
      legend.position    = "none",
      plot.title         = element_text(size = 11, face = "bold", hjust = 0.5),
      panel.grid         = element_blank(),          # sin rejilla
      panel.background   = element_rect(fill = "white", color = NA),
      plot.background    = element_rect(fill = "white", color = NA),
      axis.line          = element_line(color = "black", linewidth = 0.4),
      axis.ticks         = element_line(color = "black", linewidth = 0.3),
      plot.margin        = margin(8, 8, 8, 8)
    )
}

# ── p de log-rank usando los mismos objetos Surv y data del survfit ─────────
ld_any    <- survdiff(surv_onco ~ gbghm_grupo, data = df_onco)
p_neo_any <- 1 - pchisq(ld_any$chisq, df = length(ld_any$n) - 1)

ld_sol    <- survdiff(surv_sol ~ gbghm_grupo, data = df_onco)
p_neo_sol <- 1 - pchisq(ld_sol$chisq, df = length(ld_sol$n) - 1)

cat(sprintf("Any malignancy log-rank p = %.3f\n", p_neo_any))
cat(sprintf("Solid tumour log-rank p  = %.3f\n", p_neo_sol))

# ── Límite del eje X en años ────────────────────────────────────────────────
xmax_anios <- ceiling(max(df_mort_gb$t_analisis_anios, na.rm = TRUE))
pX1_x <- xmax_anios * 0.05

# ── Paneles (supervivencia libre de evento: descendente) ────────────────────
X1A <- panel_km(km_onco, p_neo_any, "A. Any malignancy",
                c("Single", "Mixed"), xmax_anios, 2,
                mostrar_y = TRUE,
                pval_x = pX1_x, pval_y = 0.10)
X1B <- panel_km(km_sol,  p_neo_sol, "B. Non-cutaneous solid tumour",
                c("Single", "Mixed"), xmax_anios, 2,
                mostrar_y = FALSE,
                pval_x = pX1_x, pval_y = 0.10)

# ── Etiquetas de eje X (años, marcas cada 2) ────────────────────────────────
X1A <- X1A + scale_x_continuous(breaks = seq(0, xmax_anios, 2),
                                limits = c(0, xmax_anios), expand = c(0, 0)) +
  labs(x = "Years since transplant", y = "Malignancy-free survival")
X1B <- X1B + scale_x_continuous(breaks = seq(0, xmax_anios, 2),
                                limits = c(0, xmax_anios), expand = c(0, 0)) +
  labs(x = "Years since transplant", y = NULL)

# ── Leyenda común genérica (continua = único, discontinua = mixto) ──────────
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
  theme(legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.text      = element_text(size = 10, color = "black"),
        legend.key.width = unit(1.2, "cm"))
leyenda_comun <- get_legend(p_leyenda)

# ── Componer A|B con leyenda común abajo ────────────────────────────────────
fila_paneles <- plot_grid(X1A, X1B, nrow = 1)
figura_X1 <- plot_grid(fila_paneles, leyenda_comun,
                       ncol = 1, rel_heights = c(1, 0.10))

print(figura_X1)

ggsave("~/Desktop/Figuras_tesis/figure_X1_malignancy.png",
       figura_X1,
       width = 22, height = 11, units = "cm",
       dpi = 300, bg = "white")

