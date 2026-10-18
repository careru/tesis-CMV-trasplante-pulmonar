# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.4.05. Figura 20. Frecuencia de los genotipos de gB y gH a nivel de muestra
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 20. Frecuencia de los genotipos de gB y gH a nivel de muestra ----
## Helper: porcentaje con coma decimal
fmt_pct <- function(x) paste0(formatC(x, format = "f", digits = 1, decimal.mark = ","), "%")

## 0) Bases comunes
df_gb <- muestras_sel %>%
  mutate(
    gb_sum = rowSums(across(GB1:GB4, ~ replace_na(., 0))),
    combo  = apply(select(., GB1:GB4) == 1, 1, function(r)
      paste0("gB", which(r), collapse = " + "))
  ) %>%
  filter(gb_sum %in% 1:3)

df_gh <- muestras_sel %>%
  mutate(
    gh_sum = rowSums(across(GH1:GH2, ~ replace_na(., 0))),
    combo  = apply(select(., GH1:GH2) == 1, 1, function(r)
      paste0("gH", which(r), collapse = " + "))
  ) %>%
  filter(gh_sum %in% 1:2)

## 1) Únicos gB (gB1-gB4)
tabla_uni_gb <- df_gb %>%
  filter(gb_sum == 1) %>%
  count(combo, name = "n") %>%
  arrange(n) %>%
  mutate(
    combo = factor(combo, levels = combo),
    pct   = round(100 * n / sum(n), 1),
    etiqueta = paste0(n, " (", fmt_pct(pct), ")")
  )
n_uni_gb <- sum(tabla_uni_gb$n)

cols_uni_gb <- c("gB1" = "#C6DBEF", "gB2" = "#6BAED6",
                 "gB3" = "#3182BD", "gB4" = "#08519C")

## 2) Únicos gH (gH1-gH2)
tabla_uni_gh <- df_gh %>%
  filter(gh_sum == 1) %>%
  count(combo, name = "n") %>%
  arrange(n) %>%
  mutate(
    combo = factor(combo, levels = combo),
    pct   = round(100 * n / sum(n), 1),
    etiqueta = paste0(n, " (", fmt_pct(pct), ")")
  )
n_uni_gh <- sum(tabla_uni_gh$n)

cols_uni_gh <- c("gH1" = "#DDA0DD", "gH2" = "#CD69C9")

## 3) Mixtos gB (todas las combinaciones, incl. n=1)
tabla_mix_gb <- df_gb %>%
  filter(gb_sum %in% 2:3) %>%
  count(combo, name = "n") %>%
  arrange(n) %>%
  mutate(
    combo = factor(combo, levels = combo),
    pct   = round(100 * n / sum(n), 1),
    etiqueta = paste0(n, " (", fmt_pct(pct), ")")
  )
n_mix_gb <- sum(tabla_mix_gb$n)

cols_mix_gb <- grDevices::hcl.colors(nrow(tabla_mix_gb), "Dark 3")
names(cols_mix_gb) <- as.character(tabla_mix_gb$combo)

## 4) Panel únicos gB
p_uni_gb <- ggplot(tabla_uni_gb, aes(x = n, y = combo, fill = combo)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = etiqueta), hjust = -0.15, size = 4.2) +
  scale_fill_manual(values = cols_uni_gb, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 1.30))) +
  labs(x = "", y = "",
       title = paste0("(A) Genotipos únicos gB (n = ", n_uni_gb, ")")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.margin  = margin(t = 6, r = 30, b = 6, l = 6),
    legend.position = "none",
    axis.text    = element_text(size = 11),
    axis.line    = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks   = element_line(linewidth = 0.3, colour = "grey50")
  )

## 5) Panel únicos gH
p_uni_gh <- ggplot(tabla_uni_gh, aes(x = n, y = combo, fill = combo)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = etiqueta), hjust = -0.15, size = 4.2) +
  scale_fill_manual(values = cols_uni_gh, guide = "none") +
  scale_x_continuous(
    breaks = c(0, 150, 300),
    limits = c(0, 500),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = "", y = "",
       title = paste0("(B) Genotipos únicos gH (n = ", n_uni_gh, ")")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.margin  = margin(t = 6, r = 30, b = 6, l = 6),
    legend.position = "none",
    axis.text    = element_text(size = 11),
    axis.line    = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks   = element_line(linewidth = 0.3, colour = "grey50")
  )

## 6) Panel mixtos gB
p_mix_gb <- ggplot(tabla_mix_gb, aes(x = n, y = combo, fill = combo)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = etiqueta), hjust = -0.15, size = 4.2) +
  scale_fill_manual(values = cols_mix_gb, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.30))) +
  labs(x = "Número de muestras", y = "",
       title = paste0("(C) Genotipos mixtos gB (n = ", n_mix_gb, ")")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0, size = 13,
                                margin = margin(b = 9, l = -35)),
    plot.margin  = margin(t = 6, r = 30, b = 6, l = 6),
    legend.position = "none",
    axis.text    = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    axis.line    = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks   = element_line(linewidth = 0.3, colour = "grey50")
  )

## 7) Composición: fila superior (gB | gH) y debajo mixtos gB a todo el ancho
panel_gb <- (p_uni_gb | p_uni_gh) / p_mix_gb +
  plot_layout(heights = c(max(nrow(tabla_uni_gb), nrow(tabla_uni_gh)),
                          nrow(tabla_mix_gb)))
panel_gb

guardar_fig(panel_gb, "Figura_20", ancho = 16, alto = 18)

