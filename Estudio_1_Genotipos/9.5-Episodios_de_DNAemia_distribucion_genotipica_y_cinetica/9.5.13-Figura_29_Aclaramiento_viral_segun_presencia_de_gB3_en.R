# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.5.13. Figura 29. Aclaramiento viral según presencia de gB3, en genotipos únicos o mixtos de gB
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 29. Aclaramiento viral según presencia de gB3, en genotipos únicos o mixtos de gB ----
# ── 1. Crear variables de clasificación ───────────────────────────────────────
df_gb3 <- df_uni %>%
  mutate(
    # gB3 presente en el episodio (independientemente de si es único o mixto)
    gb3_presente = as.integer(gb3 == 1),
    
    # Cuatro grupos: único/mixto x gB3 sí/no
    grupo_gb3 = case_when(
      gbmix == 0 & gb3 == 0 ~ "Único sin gB3",
      gbmix == 0 & gb3 == 1 ~ "Único con gB3",
      gbmix == 1 & gb3 == 0 ~ "Mixto sin gB3",
      gbmix == 1 & gb3 == 1 ~ "Mixto con gB3",
      TRUE ~ NA_character_
    ),
    grupo_gb3 = factor(grupo_gb3,
                       levels = c("Único sin gB3", "Único con gB3",
                                  "Mixto sin gB3", "Mixto con gB3"))
  ) %>%
  filter(!is.na(aclaramiento_dias), !is.na(grupo_gb3))

# ── 2. Análisis A: entre episodios mixtos, gB3 vs no gB3 ─────────────────────
df_mixtos <- df_gb3 %>% filter(gbmix == 1)

cat("\n── Análisis A: episodios mixtos, gB3 vs no gB3 ──\n")
df_mixtos %>%
  group_by(gb3_presente) %>%
  summarise(
    n        = n(),
    mediana  = median(aclaramiento_dias, na.rm = TRUE),
    q1       = quantile(aclaramiento_dias, 0.25, na.rm = TRUE),
    q3       = quantile(aclaramiento_dias, 0.75, na.rm = TRUE)
  ) %>% print()

wilcox.test(aclaramiento_dias ~ gb3_presente, data = df_mixtos, exact = FALSE)

# ── 3. Análisis B: cuatro grupos ──────────────────────────────────────────────
cat("\n── Análisis B: cuatro grupos ──\n")
df_gb3 %>%
  group_by(grupo_gb3) %>%
  summarise(
    n        = n(),
    mediana  = median(aclaramiento_dias, na.rm = TRUE),
    q1       = quantile(aclaramiento_dias, 0.25, na.rm = TRUE),
    q3       = quantile(aclaramiento_dias, 0.75, na.rm = TRUE)
  ) %>% print()

kruskal.test(aclaramiento_dias ~ grupo_gb3, data = df_gb3)

# Post-hoc si Kruskal significativo
dunn.test(df_gb3$aclaramiento_dias, df_gb3$grupo_gb3,
          method = "bonferroni", kw = TRUE)

etiquetas_n <- c(
  "Único sin gB3" = paste0("Único sin gB3\n(n=", sum(df_gb3$grupo_gb3 == "Único sin gB3", na.rm=TRUE), ")"),
  "Único con gB3" = paste0("Único con gB3\n(n=", sum(df_gb3$grupo_gb3 == "Único con gB3", na.rm=TRUE), ")"),
  "Mixto sin gB3" = paste0("Mixto sin gB3\n(n=", sum(df_gb3$grupo_gb3 == "Mixto sin gB3", na.rm=TRUE), ")"),
  "Mixto con gB3" = paste0("Mixto con gB3\n(n=", sum(df_gb3$grupo_gb3 == "Mixto con gB3", na.rm=TRUE), ")")
)

figura_gb3 <- ggplot(
  df_gb3,
  aes(
    x = grupo_gb3,
    y = aclaramiento_dias,
    fill = ifelse(grepl("con gB3", grupo_gb3), "con gB3", "sin gB3")
  )
) +
  geom_boxplot(
    alpha         = 0.75,
    outlier.shape = NA,      # oculta outliers nativos del boxplot
    width         = 0.55,
    color         = "grey25",
    linewidth     = 0.4
  ) +
  # Puntos individuales solo hasta 150
  geom_jitter(
    data        = df_gb3 %>% filter(aclaramiento_dias <= 150),
    aes(x = grupo_gb3, y = aclaramiento_dias,
        fill = ifelse(grepl("con gB3", grupo_gb3), "con gB3", "sin gB3")),
    shape       = 21,
    size        = 1.5,
    color       = "grey40",
    alpha       = 0.6,
    width       = 0.15,
    height      = 0,
    inherit.aes = FALSE
  ) + 
  geom_vline(
    xintercept = 2.5,
    linetype   = "dashed",
    color      = "grey70",
    linewidth  = 0.4
  ) +
  # Barras de significación
  annotate("segment", x = 1,   xend = 2, y = 170, yend = 170, linewidth = 0.4) +
  annotate("segment", x = 1,   xend = 1, y = 165, yend = 170, linewidth = 0.4) +
  annotate("segment", x = 2,   xend = 2, y = 165, yend = 170, linewidth = 0.4) +
  annotate("text",    x = 1.5, y = 174,  label = "***",  size = 5) +
  annotate("segment", x = 1,   xend = 4, y = 185, yend = 185, linewidth = 0.4) +
  annotate("segment", x = 1,   xend = 1, y = 180, yend = 185, linewidth = 0.4) +
  annotate("segment", x = 4,   xend = 4, y = 180, yend = 185, linewidth = 0.4) +
  annotate("text",    x = 2.5, y = 189,  label = "****", size = 5) +
  # Etiquetas de bloque
  annotate("text", x = 1.5, y = 208, label = "Genotipo único",
           size = 4, color = "grey30", fontface = "bold") +
  annotate("text", x = 3.5, y = 208, label = "Genotipo mixto",
           size = 4, color = "grey30", fontface = "bold") +
  scale_fill_manual(
    values = c(
      "sin gB3" = "#CDC0B0",
      "con gB3" = "#4F94CD"
    )
  ) +
  scale_x_discrete(labels = etiquetas_n) +
  scale_y_continuous(
    name   = "Tiempo hasta aclaramiento viral (días)",
    breaks = seq(0, 150, by = 50)
  ) +
  coord_cartesian(ylim = c(0, 215), clip = "off") +
  labs(x = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35),
    axis.line.y        = element_line(color = "grey60", linewidth = 0.4),
    axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
    axis.text.x        = element_text(size = 11, color = "grey20"),
    axis.text.y        = element_text(size = 11, color = "grey40"),
    axis.title.y       = element_text(size = 11, color = "grey20"),
    plot.margin        = margin(8, 8, 8, 8)
  )

figura_gb3

guardar_fig(figura_gb3, "Figura_29", ancho = 14, alto = 11)


