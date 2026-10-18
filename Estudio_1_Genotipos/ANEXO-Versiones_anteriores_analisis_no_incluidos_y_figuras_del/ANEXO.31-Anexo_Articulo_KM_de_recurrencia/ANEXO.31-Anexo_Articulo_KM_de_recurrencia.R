# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.31. Anexo: Artículo: KM de recurrencia
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Artículo: KM de recurrencia ----
# Recurrencia de DNAemia a 12 meses (estilo NEJM, 2x2). Autocontenida.
library(dplyr); library(survival); library(survminer); library(ggplot2)

col_cmp  <- c("#4A7DA6", "#E27D5F")   # azul NEJM / coral (unico vs mixto)
col_glob <- "#4D4D4D"                  # gris neutro (cohorte global)

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

df_rec <- df_recur %>% mutate(meses_365 = tiempo_km_365 / 30.44)

km_glob <- survfit(Surv(meses_365, evento_km_365) ~ 1,         data = df_rec)
km_gb   <- survfit(Surv(meses_365, evento_km_365) ~ gb_grupo,  data = df_rec)
km_gh   <- survfit(Surv(meses_365, evento_km_365) ~ gh_grupo,  data = df_rec)
km_gbgh <- survfit(Surv(meses_365, evento_km_365) ~ exp_grupo, data = df_rec)

p_gb   <- p_logrank(Surv(meses_365, evento_km_365) ~ gb_grupo,  df_rec)
p_gh   <- p_logrank(Surv(meses_365, evento_km_365) ~ gh_grupo,  df_rec)
p_gbgh <- p_logrank(Surv(meses_365, evento_km_365) ~ exp_grupo, df_rec)

xmax <- 12

panel_cmp <- function(km_fit, p_val, titulo, labs, ylab) {
  g <- ggsurvplot(
    km_fit, data = df_rec, palette = col_cmp, size = 0.9,
    conf.int = TRUE, conf.int.alpha = 0.15, conf.int.style = "ribbon",
    censor = FALSE, pval = FALSE,
    xlab = "Months since clearance", ylab = ylab,
    xlim = c(0, xmax), ylim = c(0, 1), break.time.by = 2, legend.labs = labs,
    risk.table = TRUE, risk.table.title = "No. at risk",
    risk.table.height = 0.24, risk.table.fontsize = 3.1,
    risk.table.y.text = TRUE, tables.theme = theme_cleantable(),
    ggtheme = tema_panel
  )
  g$plot <- g$plot + ggtitle(titulo) + hline50 + anota_p(p_val, xmax * 0.04, 0.08)
  g
}
panel_glob <- function(km_fit, titulo, ylab) {
  g <- ggsurvplot(
    km_fit, data = df_rec, palette = col_glob, size = 0.9,
    conf.int = TRUE, conf.int.alpha = 0.15, conf.int.style = "ribbon",
    censor = FALSE, legend = "none",
    xlab = "Months since clearance", ylab = ylab,
    xlim = c(0, xmax), ylim = c(0, 1), break.time.by = 2,
    risk.table = TRUE, risk.table.title = "No. at risk",
    risk.table.height = 0.24, risk.table.fontsize = 3.1,
    risk.table.y.text = FALSE, tables.theme = theme_cleantable(),
    ggtheme = tema_panel
  )
  g$plot <- g$plot + ggtitle(titulo) + hline50
  g
}

gA <- panel_glob(km_glob,         "A. Overall recurrence", "Recurrence-free probability")
gB <- panel_cmp (km_gbgh, p_gbgh, "B. gB/gH", c("Single gB/gH", "Mixed gB/gH"), NULL)
gC <- panel_cmp (km_gb,   p_gb,   "C. gB",    c("Single gB", "Mixed gB"), "Recurrence-free probability")
gD <- panel_cmp (km_gh,   p_gh,   "D. gH",    c("Single gH", "Mixed gH"), NULL)

arrange_ggsurvplots(list(gA, gC, gB, gD), ncol = 2, nrow = 2, print = TRUE)
cat(sprintf("Log-rank  gB = %.3f | gH = %.3f | gB/gH = %.3f\n", p_gb, p_gh, p_gbgh))

