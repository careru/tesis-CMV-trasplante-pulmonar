# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.4.01. Figura 17. Muestras genotipadas por año y número acumulado de receptores
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 17. Muestras genotipadas por año y número acumulado de receptores ----
pacientes <- read_excel(f_bd,
                        sheet = "pacientes_sel")

## 1) Muestras completamente genotipadas por año
df_comp <- muestras %>%
  mutate(
    anio = year(as.Date(Freg)),
    gb   = rowSums(across(c(GB1, GB2, GB3, GB4)), na.rm = TRUE),
    gh   = rowSums(across(c(GH1, GH2)), na.rm = TRUE),
    gt_cat = case_when(
      gb >= 1 & gh >= 1 ~ "Genotipado completo",
      xor(gb >= 1, gh >= 1) ~ "Genotipado parcial",
      TRUE ~ "Genotipado negativo"
    )
  ) %>%
  filter(gt_cat == "Genotipado completo") %>%
  count(anio, name = "n_muestras")

## 2) Trasplantes por año y acumulados a partir de `pacientes$Ftx`
tx_por_anio <- pacientes %>%
  mutate(anio = year(as_date(Ftx))) %>%
  count(anio, name = "n_tx") %>%
  arrange(anio) %>%
  mutate(tx_acum = cumsum(n_tx))

## 3) Unir tablas, asegurar años 2008–2022 y asignar metodología
anios_levels <- 2008:2022

metodologia_levels <- c(
  "ARGENE SANGRE (2008-2013)",
  "ALTONA SANGRE (2014-2016)",
  "VERIS PLASMA (2017-2019)",
  "ALTONA PLASMA (2020-2022)"
)

plot_df <- full_join(df_comp, tx_por_anio, by = "anio") %>%
  complete(anio = anios_levels) %>%
  arrange(anio) %>%
  mutate(
    anio_factor = factor(anio, levels = anios_levels),
    metodologia = case_when(
      anio %in% 2008:2013 ~ metodologia_levels[1],
      anio %in% 2014:2016 ~ metodologia_levels[2],
      anio %in% 2017:2019 ~ metodologia_levels[3],
      anio %in% 2020:2022 ~ metodologia_levels[4]
    ),
    metodologia = factor(metodologia, levels = metodologia_levels)
  )

## 4) Paleta: sangre en rojos, plasma en naranjas y amarillos
colores_metodologia <- c(
  "ARGENE SANGRE (2008-2013)" = "#8B3A3A",  # rojo granate
  "ALTONA SANGRE (2014-2016)" = "lightsalmon",  # rojo medio
  
  "VERIS PLASMA (2017-2019)"  = "darkorange3",  # naranja intenso
  "ALTONA PLASMA (2020-2022)" = "#FDD049"   # amarillo dorado
)

## 5) Ratio para el eje secundario
ratio <- max(plot_df$n_muestras, na.rm = TRUE) /
  max(plot_df$tx_acum, na.rm = TRUE)

## 6) Gráfico final
graph_acumulados <- ggplot(plot_df, aes(x = anio_factor)) +
  geom_col(aes(y = n_muestras, fill = metodologia),
           width = 0.8, na.rm = TRUE) +
  geom_text(aes(y = 0, label = n_muestras),
            vjust = 1.4, size = 4, fontface = "bold", na.rm = TRUE) +
  geom_line(aes(y = tx_acum * ratio, group = 1),
            linewidth = 0.8, na.rm = TRUE) +
  geom_point(aes(y = tx_acum * ratio),
             size = 2, na.rm = TRUE) +
  geom_text(aes(y = tx_acum * ratio, label = tx_acum),
            vjust = -1, size = 4, na.rm = TRUE) +
  scale_fill_manual(name = "Metodología:", values = colores_metodologia) +
  scale_y_continuous(
    name = "Número de muestras",
    expand = expansion(mult = c(0.05, 0.12)),
    sec.axis = sec_axis(~ . / ratio,
                        name = "Número de RTP (acumulado)")
  ) +
  coord_cartesian(clip = "off") +
  labs(x = "", y = "") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.line    = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks   = element_line(linewidth = 0.3, colour = "grey50")
  )
graph_acumulados

guardar_fig(graph_acumulados, "Figura_17", ancho = 18, alto = 12)

