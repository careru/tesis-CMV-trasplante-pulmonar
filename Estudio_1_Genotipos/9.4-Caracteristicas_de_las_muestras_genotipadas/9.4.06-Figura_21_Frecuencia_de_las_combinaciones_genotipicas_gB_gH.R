# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.4.06. Figura 21. Frecuencia de las combinaciones genotípicas gB/gH a nivel de muestra
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 21. Frecuencia de las combinaciones genotípicas gB/gH a nivel de muestra ----
## Helper: porcentaje con coma decimal
fmt_pct <- function(x) paste0(formatC(x, format = "f", digits = 1, decimal.mark = ","), "%")

## 0) Base gB+gH
df_gbgh <- muestras_sel %>%
  mutate(
    gb_sum = rowSums(across(GB1:GB4, ~ replace_na(., 0))),
    gh_sum = rowSums(across(GH1:GH2, ~ replace_na(., 0))),
    gB_combo = apply(select(., GB1:GB4) == 1, 1, function(r) {
      if (!any(r)) return(NA_character_)
      paste0("gB", which(r), collapse = " + ")
    }),
    gH_combo = apply(select(., GH1:GH2) == 1, 1, function(r) {
      if (!any(r)) return(NA_character_)
      paste0("gH", which(r), collapse = " + ")
    }),
    combo_total = paste(gB_combo, gH_combo, sep = " / "),
    # única solo si gB y gH son ambos únicos
    es_unica = (gb_sum == 1 & gh_sum == 1)
  ) %>%
  filter(!is.na(gB_combo), !is.na(gH_combo))

## 1) Combinaciones ÚNICAS (gB único / gH único)
tabla_uni <- df_gbgh %>%
  filter(es_unica) %>%
  count(combo_total, name = "n") %>%
  arrange(n) %>%
  mutate(
    combo_total = factor(combo_total, levels = combo_total),
    pct   = round(100 * n / sum(n), 1),
    etiqueta = paste0(n, " (", fmt_pct(pct), ")")
  )
n_uni <- sum(tabla_uni$n)

cols_uni <- colorRampPalette(c("#C7E9C0", "#00441B"))(nrow(tabla_uni))

## 2) Combinaciones MIXTAS (mezcla en gB, en gH o en ambas)
tabla_mix <- df_gbgh %>%
  filter(!es_unica) %>%
  count(combo_total, name = "n") %>%
  arrange(n) %>%
  mutate(
    combo_total = factor(combo_total, levels = combo_total),
    pct   = round(100 * n / sum(n), 1),
    etiqueta = paste0(n, " (", fmt_pct(pct), ")")
  )
n_mix <- sum(tabla_mix$n)

cols_mix <- colorRampPalette(c("#FCBBA1", "#67000D"))(nrow(tabla_mix))

## 3) Panel A: combinaciones únicas
p_uni <- ggplot(tabla_uni, aes(x = n, y = combo_total, fill = combo_total)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = etiqueta), hjust = -0.15, size = 4) +
  scale_fill_manual(values = cols_uni, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.5))) +
  labs(x = "", y = "",
       title = paste0("(A) Combinaciones únicas gB/gH (n = ", n_uni, ")")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0, size = 12,
                                margin = margin(b = 7, l = -30)),  
    plot.margin  = margin(t = 6, r = 30, b = 6, l = 6),
    legend.position = "none",
    axis.text    = element_text(size = 11),
    axis.line    = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks   = element_line(linewidth = 0.3, colour = "grey50")
  )

## 4) Panel B: combinaciones mixtas
p_mix <- ggplot(tabla_mix, aes(x = n, y = combo_total, fill = combo_total)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = etiqueta), hjust = -0.15, size = 4) +
  scale_fill_manual(values = cols_mix, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.50))) +
  labs(x = "Número de muestras", y = "",
       title = paste0("(B) Combinaciones mixtas gB/gH (n = ", n_mix, ")")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0, size = 12,
                                margin = margin(b = 7, l = -30)),
    plot.margin  = margin(t = 6, r = 30, b = 6, l = 6),
    legend.position = "none",
    axis.text    = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    axis.line    = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks   = element_line(linewidth = 0.3, colour = "grey50")
  )

## 5) Composición: únicas arriba, mixtas abajo (alturas según nº de barras)
panel_gbgh <- p_uni / p_mix +
  plot_layout(heights = c(nrow(tabla_uni), nrow(tabla_mix)))
panel_gbgh

guardar_fig(panel_gbgh, "Figura_21", ancho = 16, alto = 20)

