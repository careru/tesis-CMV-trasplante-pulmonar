# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.7.11. Figura 38. Infecciones oportunistas no relacionadas con CMV
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 38. Infecciones oportunistas no relacionadas con CMV ----
# ── 6. Preparar datos ─────────────────────────────────────────────────────────
df_opo <- pacientes %>%
  mutate(
    Ftx_dt           = as.Date(Ftx),
    F_opo_dt         = as.Date(F_opo),
    Ffinseg_dt       = as.Date(Ffinseg),
    fecha_fin_real   = if_else(inf_opo == 1, F_opo_dt, Ffinseg_dt),
    t_analisis_anios = as.numeric(fecha_fin_real - Ftx_dt) / 365.25,
    evento_opo       = if_else(inf_opo == 1, 1L, 0L),
    risk_f           = factor(risk, levels = c("MEDIO", "ALTO")),
    gbghm_grupo      = factor(gbghm, levels = c(0,1),
                              labels = c("Genotipos únicos",
                                         "Genotipos mixtos"))
  ) %>%
  filter(t_analisis_anios > 0)

# ── 7. Descriptivo ────────────────────────────────────────────────────────────
cat("\n── Infecciones oportunistas por grupo ──\n")
df_opo %>%
  group_by(gbghm_grupo) %>%
  summarise(
    n     = n(),
    n_opo = sum(evento_opo),
    pct   = round(n_opo / n * 100, 1),
    .groups = "drop"
  ) %>% print()

cat("\n── Total pacientes con IO ──\n")
cat(sprintf("n = %d (%.1f%%)\n",
            sum(pacientes$inf_opo == 1, na.rm = TRUE),
            100 * mean(pacientes$inf_opo == 1, na.rm = TRUE)))

cat("\n── Número de episodios por paciente ──\n")
table(pacientes$n_opo[pacientes$inf_opo == 1], useNA = "always") %>% print()

cat("\n── Etiología IO ──\n")
table(pacientes$det_opo, useNA = "always") %>%
  sort(decreasing = TRUE) %>% print()

cat("\n── Test chi-cuadrado IO vs grupo ──\n")
print(chisq.test(table(df_opo$gbghm_grupo, df_opo$evento_opo)))

# ── 8. KM + log-rank ─────────────────────────────────────────────────────────
surv_opo <- Surv(df_opo$t_analisis_anios, df_opo$evento_opo)
km_opo   <- survfit(surv_opo ~ gbghm_grupo, data = df_opo)
lr_opo   <- survdiff(surv_opo ~ gbghm_grupo, data = df_opo)
p_lr     <- 1 - pchisq(lr_opo$chisq, df = 1)
cat(sprintf("\n── Log-rank IO: p = %.3f ──\n", p_lr))

cat("\n── Supervivencia libre IO a 5 años ──\n")
print(summary(km_opo, times = 5))

# ── 9. Cox univariante y ajustado ─────────────────────────────────────────────
cox_opo_crudo <- coxph(surv_opo ~ gbghm_grupo, data = df_opo)
cox_opo_ajust <- coxph(surv_opo ~ gbghm_grupo + edad_tx + risk_f,
                       data = df_opo)

cat("\n── Cox IO crudo ──\n")
hr_oc <- exp(coef(cox_opo_crudo)[1])
ci_oc <- exp(confint(cox_opo_crudo)[1,])
p_oc  <- summary(cox_opo_crudo)$coefficients[1,5]
cat(sprintf("HR: %.2f (%.2f–%.2f)  p = %.3f\n",
            hr_oc, ci_oc[1], ci_oc[2], p_oc))

cat("\n── Cox IO ajustado ──\n")
data.frame(
  HR    = round(exp(coef(cox_opo_ajust)), 2),
  IC_lo = round(exp(confint(cox_opo_ajust))[,1], 2),
  IC_hi = round(exp(confint(cox_opo_ajust))[,2], 2),
  p     = round(summary(cox_opo_ajust)$coefficients[,5], 3)
) %>% print()

# ── IO segun genotipo mixto: gB, gH y gB/gH ──────────────────────────────────

# Funcion: descriptivo + test para una variable de mezcla
io_por_genotipo <- function(data, var_mix, etiqueta) {
  d <- data %>%
    mutate(grupo = factor(.data[[var_mix]], levels = c(0, 1),
                          labels = c("\u00danicos", "Mixtos"))) %>%
    filter(!is.na(grupo))
  
  tab <- table(d$grupo, d$evento_opo)
  
  resumen <- d %>%
    group_by(grupo) %>%
    summarise(
      n     = n(),
      n_opo = sum(evento_opo),
      pct   = round(100 * n_opo / n, 1),
      .groups = "drop"
    )
  
  # Fisher (robusto con n pequenos) y chi-cuadrado de Yates
  p_fisher <- fisher.test(tab)$p.value
  p_chisq  <- suppressWarnings(chisq.test(tab)$p.value)
  
  cat(sprintf("\n══════ %s ══════\n", etiqueta))
  print(resumen)
  cat(sprintf("p (Fisher)        = %s\n",
              sub("\\.", ",", sprintf("%.3f", p_fisher))))
  cat(sprintf("p (chi-cuadrado)  = %s\n",
              sub("\\.", ",", sprintf("%.3f", p_chisq))))
  
  invisible(resumen)
}

# Aplicar a las tres clasificaciones (usa df_opo, que ya tiene evento_opo)
io_por_genotipo(df_opo, "gbm",   "gB mixtos")
io_por_genotipo(df_opo, "ghm",   "gH mixtos")
io_por_genotipo(df_opo, "gbghm", "gB/gH mixtos")

# [Nuevo] Datos de la Figura 38: episodios por etiología a partir de det_opo.
# En el script antiguo df_counts y colores no se definían (venían de otra sesión).
# Cada paciente aporta n_opo episodios; las etiologías se leen de det_opo.
# Si det_opo menciona más etiologías que n_opo, se conservan las últimas n_opo
# (afecta a un paciente: "VEB (PTLD) + IFI pulmonar", con n_opo = 1, que cuenta como IFI).
clasificar_opo <- function(txt) {
  partes <- trimws(unlist(strsplit(txt, ",|\\+|;")))
  partes <- partes[partes != ""]
  case_when(
    grepl("VVZ", partes, ignore.case = TRUE)            ~ "VVZ",
    grepl("candidiasis", partes, ignore.case = TRUE)    ~ "Candidiasis",
    grepl("IFI", partes)                                ~ "IFI",
    grepl("VEB|PTLD|PLTD", partes)                      ~ "VEB/PTLD",
    grepl("TBC", partes)                                ~ "TBC",
    grepl("Cryptosporidium", partes)                    ~ "Cryptosporidium",
    grepl("VHH8|Kaposi", partes)                        ~ "VHH8/Sarcoma de Kaposi",
    grepl("JC", partes)                                 ~ "Virus JC",
    TRUE                                                ~ NA_character_
  )
}

episodios_opo <- pacientes %>%
  filter(inf_opo == 1, !is.na(det_opo)) %>%
  select(NHC, det_opo, n_opo) %>%
  mutate(categoria = lapply(det_opo, clasificar_opo)) %>%
  mutate(categoria = map2(categoria, n_opo, ~ tail(.x[!is.na(.x)], .y))) %>%
  tidyr::unnest(categoria)

df_counts <- episodios_opo %>%
  count(categoria) %>%
  mutate(pct = round(100 * n / sum(n), 1))

stopifnot(sum(df_counts$n) == sum(pacientes$n_opo, na.rm = TRUE))
print(df_counts %>% arrange(desc(n)))

colores <- c(
  "VVZ"                    = "#4E79A7",
  "Candidiasis"            = "#F28E2B",
  "IFI"                    = "#8CD17D",
  "VEB/PTLD"               = "#E15759",
  "TBC"                    = "#7B3F9E",
  "VHH8/Sarcoma de Kaposi" = "#2CB5A0",
  "Virus JC"               = "#8C8C8C",
  "Cryptosporidium"        = "#EDC948"
)
# ── Etiquetas personalizadas (edita aquí a tu gusto) ─────────────────────────
etiquetas_custom <- c(
  "VVZ"                    = "Herpes zóster",
  "Candidiasis"            = "Candidiasis mucocutánea",
  "IFI"                    = "Infección fúngica invasiva",
  "VEB/PTLD"               = "PTLD asociado a VEB",
  "TBC"                    = "Tuberculosis diseminada",
  "Cryptosporidium"        = "Criptosporidiosis",
  "VHH8/Sarcoma de Kaposi" = "Sarcoma de Kaposi (VHH-8)",
  "Virus JC"               = "LMP por virus JC"
)

figura_opo_barras <- ggplot(df_counts %>% arrange(n),
                            aes(x = reorder(categoria, n), y = n,
                                fill = categoria)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.2) +
  geom_text(
    aes(label = paste0("n=", n, " (",
                       gsub("\\.", ",", as.character(pct)), "%)")),
    hjust = -0.1, size = 3.8, color = "grey20"
  )  +
  scale_fill_manual(values = colores, labels = etiquetas_custom) +
  scale_x_discrete(labels = etiquetas_custom) +  # ← etiquetas eje Y
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.45))
  ) +
  coord_flip() +
  labs(x = NULL, y = "Número de episodios") +
  theme_classic(base_size = 13) +
  theme(
    legend.position    = "none",
    axis.text.y        = element_text(size = 11, color = "grey20"),
    axis.text.x        = element_text(size = 10, color = "grey40"),
    axis.title.x       = element_text(size = 11, color = "grey20"),
    axis.line          = element_line(color = "grey60", linewidth = 0.4),
    axis.ticks         = element_line(color = "grey60", linewidth = 0.3),
    plot.margin        = margin(10, 20, 10, 10)
  )

figura_opo_barras

guardar_fig(figura_opo_barras, "Figura_38", ancho = 16, alto = 9)

