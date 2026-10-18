# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.32. Anexo: Artículo: KM de tumores
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Artículo: KM de tumores ----
# Supervivencia libre de neoplasia (estilo NEJM): A. cualquier neoplasia
# B. tumor solido no cutaneo. Autocontenida.
library(dplyr); library(survival); library(survminer); library(ggplot2)

col_cmp <- c("#4A7DA6", "#E27D5F")   # azul NEJM / coral
tema_panel <- theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0),
    axis.title      = element_text(size = 11),
    axis.text       = element_text(size = 10, color = "black"),
    axis.line       = element_line(linewidth = 0.4, color = "black"),
    axis.ticks      = element_line(linewidth = 0.4, color = "black"),
    legend.position = "top",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 10)
  )
hline50 <- geom_hline(yintercept = 0.5, linetype = "dashed",
                      color = "grey60", linewidth = 0.35)
p_logrank <- function(fml, data) {
  ld <- survdiff(fml, data = data); 1 - pchisq(ld$chisq, df = length(ld$n) - 1)
}
anota_p <- function(p_val, x, y) {
  annotate("text", x = x, y = y, hjust = 0, size = 4, parse = TRUE,
           label = paste0("italic('p')*' = '*'",
                          formatC(p_val, format = "f", digits = 3), "'"))
}

df_onco <- pacientes %>%
  mutate(
    Ftx = as.Date(Ftx), F_onco = as.Date(F_onco), Ffinseg = as.Date(Ffinseg),
    fecha_fin_real   = if_else(onco == 1, F_onco, Ffinseg),
    t_analisis_anios = as.numeric(fecha_fin_real - Ftx) / 365.25,
    evento_onco      = if_else(onco == 1, 1L, 0L),
    solido         = as.integer(onco == 1 & onco_tipo == "SOL"),
    fecha_fin_sol  = if_else(solido == 1, F_onco, Ffinseg),
    t_analisis_sol = as.numeric(fecha_fin_sol - Ftx) / 365.25,
    evento_solido  = if_else(solido == 1, 1L, 0L),
    gbghm_grupo = factor(gbghm, levels = c(0, 1),
                         labels = c("Single genotype", "Mixed genotype"))
  )

km_onco <- survfit(Surv(t_analisis_anios, evento_onco)  ~ gbghm_grupo, data = df_onco)
km_sol  <- survfit(Surv(t_analisis_sol,  evento_solido) ~ gbghm_grupo, data = df_onco)
p_any <- p_logrank(Surv(t_analisis_anios, evento_onco)  ~ gbghm_grupo, df_onco)
p_sol <- p_logrank(Surv(t_analisis_sol,  evento_solido) ~ gbghm_grupo, df_onco)

xmax <- ceiling(max(df_onco$t_analisis_anios, na.rm = TRUE) / 2) * 2
breaks_x <- sort(unique(c(0.5, 1, seq(0, xmax, 2))))
labels_x <- ifelse(breaks_x %in% seq(0, xmax, 2), as.character(breaks_x), "")
escala_x <- scale_x_continuous(
  limits = c(0, xmax), breaks = breaks_x, labels = labels_x,
  expand = expansion(mult = c(0.02, 0.02))
)

panel_neo <- function(km_fit, p_val, titulo, ylab) {
  g <- ggsurvplot(
    km_fit, data = df_onco, palette = col_cmp, size = 0.9,
    conf.int = TRUE, conf.int.alpha = 0.15, conf.int.style = "ribbon",
    censor = FALSE, surv.scale = "percent", pval = FALSE,
    xlab = "Years since transplant", ylab = ylab,
    xlim = c(0, xmax), break.time.by = 2,
    legend.labs = c("Single genotype", "Mixed genotype"),
    risk.table = TRUE, risk.table.title = "No. at risk",
    risk.table.height = 0.22, risk.table.fontsize = 3.3,
    risk.table.y.text = TRUE, tables.theme = theme_cleantable(),
    ggtheme = tema_panel
  )
  g$plot <- g$plot + ggtitle(titulo) + hline50 + anota_p(p_val, 0.3, 0.10) + escala_x
  g
}

gA <- panel_neo(km_onco, p_any, "A. Any malignancy", "Malignancy-free survival (%)")
gB <- panel_neo(km_sol,  p_sol, "B. Non-cutaneous solid tumor", "Solid-tumor-free survival (%)")
arrange_ggsurvplots(list(gA, gB), ncol = 2, nrow = 1, print = TRUE)
cat(sprintf("Log-rank any = %.3f | solid = %.3f\n", p_any, p_sol))

