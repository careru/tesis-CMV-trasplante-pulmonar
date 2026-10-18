# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.27. Anexo: KM de infecciones oportunistas (no incluida)
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: KM de infecciones oportunistas (no incluida) ----
p_km_opo <- ggsurvplot(
  km_opo,
  data              = df_opo,
  pval              = paste0("p = ", formatC(p_lr, format="f", digits=3)),
  pval.coord        = c(0, 0.15),
  pval.size         = 4,
  conf.int          = TRUE,
  conf.int.alpha    = 0.15,
  conf.int.style    = "ribbon",
  risk.table        = TRUE,
  risk.table.height = 0.25,
  risk.table.y.text = FALSE,
  ncensor.plot      = FALSE,
  palette           = c("#3498db", "#e74c3c"),
  legend.title      = "",
  legend.labs       = c("Genotipos únicos", "Genotipos mixtos"),
  xlab              = "Tiempo desde trasplante (años)",
  ylab              = "Supervivencia libre de IO no-CMV",
  title             = " ",
  risk.table.title  = "Pacientes en riesgo",
  ggtheme           = theme_minimal(base_size = 12),
  font.legend       = c(11, "plain", "black"),
  font.x            = c(11, "plain", "black"),
  font.y            = c(11, "plain", "black"),
  font.tickslab     = c(10, "plain", "black"),
  break.time.by     = 1,
  xlim              = c(0, max(df_opo$t_analisis_anios, na.rm = TRUE))
)

p_km_opo$plot <- p_km_opo$plot +
  theme(
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA)
  )

p_km_opo$table <- p_km_opo$table +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

print(p_km_opo)


 

