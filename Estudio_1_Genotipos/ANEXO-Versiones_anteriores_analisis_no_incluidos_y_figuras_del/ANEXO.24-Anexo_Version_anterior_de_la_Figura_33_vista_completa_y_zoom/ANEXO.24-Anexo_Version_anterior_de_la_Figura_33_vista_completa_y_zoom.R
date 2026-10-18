# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.24. Anexo: Versión anterior de la Figura 33 (vista completa y zoom)
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Versión anterior de la Figura 33 (vista completa y zoom) ----
# ── 1. Preparar datos ────────────────────────────────────────────────────────
df_fig <- dnaemias %>%
  filter(primer_mil == 1) %>%
  mutate(
    Ftx              = as.Date(Ftx),
    Fini             = as.Date(Fini),
    meses_tx_inicio  = as.numeric(Fini - Ftx) / 30.44,
    anios_tx_inicio  = meses_tx_inicio / 12       # ← nueva columna
  ) %>%
  select(-any_of(c("gbm", "ghm", "gbghm"))) %>%
  left_join(pacientes %>% select(NHC, gbm, ghm, gbghm), by = "NHC") %>%
  mutate(
    serDR = factor(serDR, levels = c("D-/R+", "D+/R+", "D+/R-"))
  ) %>%
  filter(!is.na(serDR), !is.na(meses_tx_inicio))

df_long <- df_fig %>%
  select(NHC, serDR, meses_tx_inicio, anios_tx_inicio, gbm, ghm, gbghm) %>%
  pivot_longer(
    cols      = c(gbm, ghm, gbghm),
    names_to  = "glucoproteina",
    values_to = "mixto"
  ) %>%
  filter(!is.na(mixto)) %>%
  mutate(
    glucoproteina = factor(glucoproteina,
                           levels = c("gbm", "ghm", "gbghm"),
                           labels = c("gB", "gH", "gB+gH")),
    genotipo      = factor(ifelse(mixto == 1, "Mixto", "Único"),
                           levels = c("Único", "Mixto")),
    cat_color = paste(glucoproteina, genotipo, sep = "_"),
    cat_color = factor(cat_color,
                       levels = c("gB_Único",   "gB_Mixto",
                                  "gH_Único",   "gH_Mixto",
                                  "gB+gH_Único","gB+gH_Mixto"))
  )

# ── 2. N por panel: serología y n en una sola línea ─────────────────────────
counts_panel <- df_fig %>%
  group_by(serDR) %>%
  summarise(n = n_distinct(NHC), .groups = "drop")

# ── CAMBIO: serología y n en MISMA línea, no en dos líneas ────────────────
etiquetas_riesgo <- setNames(
  sprintf("%s (n = %d)", counts_panel$serDR, counts_panel$n),
  as.character(counts_panel$serDR)
)

# ── 3. Líneas de profilaxis (solo Panel B) ──────────────────────────────────
df_prof <- tibble(
  serDR = factor(c("D-/R+", "D+/R+", "D+/R-"),
                 levels = c("D-/R+", "D+/R+", "D+/R-")),
  prof_meses = c(6, 6, 12)
)

# ── 4. Tests significativos (solo en zoom 0-24m) ───────────────────────────
DODGE <- 0.7

stat_test <- df_long %>%
  group_by(serDR, glucoproteina) %>%
  rstatix::wilcox_test(meses_tx_inicio ~ genotipo) %>%
  rstatix::add_significance("p") %>%
  ungroup() %>%
  filter(p.signif != "ns") %>%
  mutate(
    x_pos      = as.numeric(glucoproteina),
    xmin       = x_pos - DODGE / 4,
    xmax       = x_pos + DODGE / 4,
    y.position = 22
  )

# ── 5. Paletas ──────────────────────────────────────────────────────────────
cols_fill <- c(
  "gB_Único"     = "#AED6F1",  "gB_Mixto"     = "#4682B4",
  "gH_Único"     = "pink",     "gH_Mixto"     = "hotpink3",
  "gB+gH_Único"  = "#4CAF7D",  "gB+gH_Mixto"  = "#E07070"
)

cols_jit <- c(
  "gB_Único"     = "#2980B9",  "gB_Mixto"     = "#0D2B45",
  "gH_Único"     = "pink3",    "gH_Mixto"     = "palevioletred4",
  "gB+gH_Único"  = "#2E7D52",  "gB+gH_Mixto"  = "#B04040"
)

# ════════════════════════════════════════════════════════════════════════════
#  PANEL A — Vista completa (eje en AÑOS, sin línea 24m)
# ════════════════════════════════════════════════════════════════════════════
fig_A_completa <- ggplot(df_long,
                         aes(x = glucoproteina, y = anios_tx_inicio,
                             fill = cat_color,
                             group = interaction(glucoproteina, genotipo))) +
  geom_boxplot(
    width          = 0.55,
    position       = position_dodge(width = DODGE),
    color          = "grey25",
    linewidth      = 0.4,
    outlier.shape  = 21,
    outlier.size   = 1.8,
    outlier.fill   = "white",
    outlier.color  = "grey50",
    outlier.stroke = 0.3,
    outlier.alpha  = 0.8
  ) +
  geom_jitter(
    aes(color = cat_color),
    position    = position_jitterdodge(jitter.width = 0.10,
                                       dodge.width  = DODGE),
    size        = 1.0,
    alpha       = 0.4,
    show.legend = FALSE
  ) +
  # ── CAMBIO: línea 24m ELIMINADA ───────────────────────────────────────────
  scale_fill_manual(values = cols_fill, guide = "none") +
  scale_color_manual(values = cols_jit, guide = "none") +
  facet_wrap(~ serDR, labeller = labeller(serDR = etiquetas_riesgo)) +
  labs(
    x = NULL,
    y = "Años desde el TX al 1er episodio",   # ← AÑOS y 1er episodio
    title = "A. Distribución completa del tiempo al primer episodio"
  ) +
  scale_y_continuous(
    breaks = c(0, 2, 4, 6, 8, 10),               # ← marcas en años
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  coord_cartesian(ylim = c(0, 11)) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "none",
    panel.border       = element_blank(),
    panel.spacing.x    = unit(15, "pt"),
    axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
    axis.line.y        = element_line(color = "grey60", linewidth = 0.4),
    axis.ticks         = element_line(color = "grey60", linewidth = 0.3),
    axis.text.x        = element_text(size = 12, color = "grey20"),
    axis.text.y        = element_text(size = 11, color = "grey30"),
    axis.title.y       = element_text(size = 12, color = "grey15",
                                      margin = margin(r = 8)),
    plot.title         = element_text(size = 13, color = "grey15",
                                      face = "bold",
                                      margin = margin(b = 8)),
    strip.background   = element_rect(fill = "grey94",
                                      color = "grey70", linewidth = 0.4),
    strip.text         = element_text(size = 12, color = "grey15",
                                      face = "bold",
                                      margin = margin(t = 5, b = 5)),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.margin        = margin(10, 12, 6, 10)
  )

# ════════════════════════════════════════════════════════════════════════════
#  PANEL B — Zoom 0-24 meses (sin la palabra "Zoom")
# ════════════════════════════════════════════════════════════════════════════
fig_B_zoom <- ggplot(df_long,
                     aes(x = glucoproteina, y = meses_tx_inicio,
                         fill = cat_color,
                         group = interaction(glucoproteina, genotipo))) +
  geom_hline(data = df_prof,
             aes(yintercept = prof_meses),
             color = "grey55", linetype = "dashed", linewidth = 0.5) +
  geom_boxplot(
    width          = 0.55,
    position       = position_dodge(width = DODGE),
    color          = "grey25",
    linewidth      = 0.4,
    outlier.shape  = 21,
    outlier.size   = 1.6,
    outlier.fill   = "white",
    outlier.color  = "grey50",
    outlier.stroke = 0.3,
    outlier.alpha  = 0.7
  ) +
  geom_jitter(
    aes(color = cat_color),
    position    = position_jitterdodge(jitter.width = 0.10,
                                       dodge.width  = DODGE),
    size        = 1.0,
    alpha       = 0.4,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = cols_fill, guide = "none") +
  scale_color_manual(values = cols_jit, guide = "none") +
  facet_wrap(~ serDR, labeller = labeller(serDR = etiquetas_riesgo)) +
  {
    if (nrow(stat_test) > 0)
      stat_pvalue_manual(
        stat_test,
        label        = "p.signif",
        xmin         = "xmin",
        xmax         = "xmax",
        y.position   = "y.position",
        tip.length   = 0.01,
        bracket.size = 0.5,
        size         = 6,
        color        = "grey20",
        inherit.aes  = FALSE
      )
  } +
  labs(
    x = NULL,
    y = "Meses desde el TX al 1er episodio",          # ← MESES y 1er episodio
    title = "B. Primeros 24 meses post-trasplante"     # ← sin "Zoom"
  ) +
  scale_y_continuous(
    breaks = seq(0, 24, 6),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  coord_cartesian(ylim = c(0, 24)) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "none",
    panel.border       = element_blank(),
    panel.spacing.x    = unit(15, "pt"),
    axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
    axis.line.y        = element_line(color = "grey60", linewidth = 0.4),
    axis.ticks         = element_line(color = "grey60", linewidth = 0.3),
    axis.text.x        = element_text(size = 12, color = "grey20"),
    axis.text.y        = element_text(size = 11, color = "grey30"),
    axis.title.y       = element_text(size = 12, color = "grey15",
                                      margin = margin(r = 8)),
    plot.title         = element_text(size = 13, color = "grey15",
                                      face = "bold",
                                      margin = margin(b = 8)),
    strip.background   = element_rect(fill = "grey94",
                                      color = "grey70", linewidth = 0.4),
    strip.text         = element_text(size = 12, color = "grey15",
                                      face = "bold",
                                      margin = margin(t = 5, b = 5)),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.margin        = margin(10, 12, 6, 10)
  )

# ── Combinar paneles ────────────────────────────────────────────────────────
figura_R20_doble <- fig_A_completa / fig_B_zoom +
  plot_layout(heights = c(1, 1))

figura_R20_doble

ggsave("~/Desktop/Figuras_tesis/figura_R20_doble_vista.png",
       figura_R20_doble,
       width = 10, height = 9, dpi = 320, bg = "white")

