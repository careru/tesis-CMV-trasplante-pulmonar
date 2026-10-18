# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.4.11. Datos para el texto (9.4.4): Ct en genotipos únicos frente a mixtos
# Extracto del script completo (ver README.md del estudio)

## ---- Datos para el texto (9.4.4): Ct en genotipos únicos frente a mixtos ----
# ── 1. Construir dataset largo ────────────────────────────────────────────────
df_ct_comp <- bind_rows(
  muestras_sel %>%
    filter(GB1 == 1) %>%
    transmute(genotipo = "gB1", CT = GB1CT,
              tipo = if_else(gb == 1, "Único", "Mixto")),
  muestras_sel %>%
    filter(GB2 == 1) %>%
    transmute(genotipo = "gB2", CT = GB2CT,
              tipo = if_else(gb == 1, "Único", "Mixto")),
  muestras_sel %>%
    filter(GB3 == 1) %>%
    transmute(genotipo = "gB3", CT = GB3CT,
              tipo = if_else(gb == 1, "Único", "Mixto")),
  muestras_sel %>%
    filter(GB4 == 1) %>%
    transmute(genotipo = "gB4", CT = GB4CT,
              tipo = if_else(gb == 1, "Único", "Mixto")),
  muestras_sel %>%
    filter(GH1 == 1) %>%
    transmute(genotipo = "gH1", CT = GH1CT,
              tipo = if_else(gh == 1, "Único", "Mixto")),
  muestras_sel %>%
    filter(GH2 == 1) %>%
    transmute(genotipo = "gH2", CT = GH2CT,
              tipo = if_else(gh == 1, "Único", "Mixto"))
) %>%
  filter(!is.na(CT)) %>%
  mutate(
    genotipo = factor(genotipo,
                      levels = c("gB1", "gB2", "gB3", "gB4", "gH1", "gH2")),
    tipo = factor(tipo, levels = c("Único", "Mixto"))
  )

# ── 2. Resumen numérico ───────────────────────────────────────────────────────
cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  Ct por genotipo: único vs mixto (mediana [RIC], n)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

df_ct_comp %>%
  group_by(genotipo, tipo) %>%
  summarise(
    n       = n(),
    mediana = median(CT, na.rm = TRUE),
    p25     = quantile(CT, 0.25, na.rm = TRUE),
    p75     = quantile(CT, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(resumen = sprintf("%.1f (%.1f–%.1f), n=%d", mediana, p25, p75, n)) %>%
  select(genotipo, tipo, resumen) %>%
  pivot_wider(names_from = tipo, values_from = resumen) %>%
  print()

# ── 3. Tests de Wilcoxon por genotipo ─────────────────────────────────────────
cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  Wilcoxon único vs mixto por genotipo\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

stats_ct <- df_ct_comp %>%
  group_by(genotipo) %>%
  wilcox_test(CT ~ tipo) %>%
  add_significance("p")

print(stats_ct %>% select(genotipo, group1, group2, n1, n2, p, p.signif))

# ── 4. Posiciones y para las anotaciones ──────────────────────────────────────
y_pos <- df_ct_comp %>%
  group_by(genotipo) %>%
  summarise(y.position = max(CT, na.rm = TRUE) + 1.5, .groups = "drop")

stats_plot <- left_join(stats_ct, y_pos, by = "genotipo") %>%
  mutate(
    xmin = 1,
    xmax = 2
  ) %>%
  filter(p < 0.05)   # ← solo significativos

# ── 5. Figura ─────────────────────────────────────────────────────────────────
cols_tipo   <- c("Único" = "seagreen3", "Mixto" = "#E07070")
cols_jitter <- c("Único" = "seagreen4", "Mixto" = "#A03030")

figura_ct_comp <- ggplot(df_ct_comp,
                         aes(x = tipo, y = CT, fill = tipo)) +
  geom_boxplot(
    width         = 0.55,
    color         = "grey25",
    linewidth     = 0.45,
    outlier.shape = 21,
    outlier.size  = 1.8,
    outlier.fill  = "white",
    outlier.color = "grey50",
    outlier.alpha = 0.7
  ) +
  geom_jitter(
    aes(color = tipo),
    width = 0.12, size = 0.9, alpha = 0.35,
    show.legend = FALSE
  ) +
  stat_pvalue_manual(
    stats_plot,
    label        = "p.signif",
    xmin         = "xmin",
    xmax         = "xmax",
    y.position   = "y.position",
    tip.length   = 0.01,
    bracket.size = 0.45,
    size         = 4,
    inherit.aes  = FALSE
  ) +
  facet_wrap(~ genotipo, nrow = 1, scales = "free_y") +
  scale_fill_manual(values  = cols_tipo) +
  scale_color_manual(values = cols_jitter) +
  scale_y_continuous(
    name   = "Valor Ct",
    breaks = seq(25, 42, 5),
    expand = expansion(mult = c(0.02, 0.15))
  ) +
  labs(x = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "bottom",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_blank(),
    axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
    axis.line.y        = element_line(color = "grey60", linewidth = 0.4),
    axis.ticks         = element_line(color = "grey60", linewidth = 0.3),
    axis.text.x        = element_text(size = 10, color = "grey30"),
    axis.text.y        = element_text(size = 10, color = "grey30"),
    strip.background   = element_blank(),
    strip.text         = element_text(size = 12, face = "bold", color = "grey20"),
    panel.spacing.x    = unit(16, "pt"),
    plot.margin        = margin(10, 15, 8, 8)
  )

figura_ct_comp


