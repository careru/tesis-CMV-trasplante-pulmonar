# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.4.08. Figura 23. Carga viral según genotipado completo, parcial o negativo
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 23. Carga viral según genotipado completo, parcial o negativo ----
df <- muestras %>%
  mutate(
    gb = rowSums(across(c(GB1, GB2, GB3, GB4)), na.rm = TRUE),
    gh = rowSums(across(c(GH1, GH2)), na.rm = TRUE),
    gt_cat = case_when(
      gb >= 1 & gh >= 1 ~ "Genotipado completo",
      xor(gb >= 1, gh >= 1) ~ "Genotipado parcial",
      TRUE ~ "Genotipado negativo"
    ),
    gt_cat = factor(gt_cat, levels = c(
      "Genotipado completo",
      "Genotipado parcial",
      "Genotipado negativo"
    )),
    logCV = log10(CV)
  ) %>%
  drop_na(gt_cat, logCV)

comparaciones <- list(
  c("Genotipado completo", "Genotipado parcial"),
  c("Genotipado completo", "Genotipado negativo"),
  c("Genotipado parcial", "Genotipado negativo")
)

stats <- df %>%
  pairwise_wilcox_test(logCV ~ gt_cat, comparisons = comparaciones, p.adjust.method = "none") %>%
  add_significance("p") %>%
  add_y_position(fun = "max", step.increase = 0.15)

y_top <- max(stats$y.position, na.rm = TRUE) * 1.05

p <- ggplot(df, aes(x = gt_cat, y = logCV, fill = gt_cat)) +
  geom_boxplot(
    width = 0.50, color = "grey25", linewidth = 0.45,
    outlier.shape = 21, outlier.size = 1.8,
    outlier.fill = "white", outlier.color = "grey50",
    outlier.stroke = 0.4, outlier.alpha = 0.7
  ) +
  geom_jitter(
    aes(color = gt_cat),
    width = 0.12, size = 1.0, alpha = 0.35, show.legend = FALSE
  ) +
  scale_fill_manual(values = c(
    "Genotipado completo" = "#63B8FF",
    "Genotipado parcial"  = "#FFA54F",
    "Genotipado negativo" = "#BDBDBD"
  )) +
  scale_color_manual(values = c(
    "Genotipado completo" = "#1A6EBF",
    "Genotipado parcial"  = "#CC7A1A",
    "Genotipado negativo" = "#757575"
  )) +
  labs(x = NULL, y = "Carga viral CMV Log10 (UI/mL)") +
  scale_y_continuous(
    expand = expansion(mult = c(0.04, 0.12)),
    limits = c(NA, y_top)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position  = "none",
    panel.grid.major.y = element_line(color = "white", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_blank(),
    axis.line.x        = element_line(color = "grey30", linewidth = 0.5),
    axis.line.y        = element_line(color = "grey30", linewidth = 0.5),
    axis.ticks         = element_line(color = "grey30", linewidth = 0.4),
    axis.ticks.length  = unit(3, "pt"),
    axis.text.x        = element_text(size = 11, color = "grey20", margin = margin(t = 4)),
    axis.text.y        = element_text(size = 10, color = "grey30"),
    axis.title.y       = element_text(size = 11, color = "grey20", margin = margin(r = 8)),
    plot.margin        = margin(10, 15, 8, 8)
  ) +
  coord_cartesian(clip = "off") +
  stat_pvalue_manual(
    data = stats,
    label = "p.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.012,
    bracket.size = 0.4,
    size = 4.0,
    color = "grey30",
    inherit.aes = FALSE
  )

p

guardar_fig(p, "Figura_23", ancho = 14, alto = 12)

