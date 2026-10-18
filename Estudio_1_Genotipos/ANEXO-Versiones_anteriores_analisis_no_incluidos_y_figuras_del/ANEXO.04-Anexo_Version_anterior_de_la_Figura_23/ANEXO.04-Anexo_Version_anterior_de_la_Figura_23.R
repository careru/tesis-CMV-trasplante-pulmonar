# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.04. Anexo: Versión anterior de la Figura 23
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Versión anterior de la Figura 23 ----
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
    width = 0.55, color = "black", size = 0.4,
    outlier.shape = 21, outlier.size = 1.6,
    outlier.stroke = 0, outlier.alpha = 0.5
  ) +
  scale_fill_manual(values = c(
    "Genotipado completo" = "#63B8FF",
    "Genotipado parcial"  = "#FFA54F",
    "Genotipado negativo" = "grey70"
  )) +
  labs(x = NULL, y = "Carga viral CMV log10 (UI/mL)") +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.10)),
    limits = c(NA, y_top)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position  = "none",
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "grey40", fill = NA, linewidth = 0.4),
    axis.line.x      = element_line(color = "black", linewidth = 0.4),
    axis.line.y      = element_line(color = "black", linewidth = 0.4),
    panel.spacing.x  = grid::unit(16, "pt"),
    plot.margin      = margin(5, 5, 5, 5)
  ) +
  coord_cartesian(clip = "off") +
  stat_pvalue_manual(
    data = stats,
    label = "p.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.45,
    size = 3.6,
    inherit.aes = FALSE
  )

p

# Pvalores CV en GENOTIPADO- Mann–Whitney U test
stats %>%
  select(group1, group2, p, p.signif) %>%
  mutate(
    p = formatC(p, format = "f", digits = 5, decimal.mark = ",")) %>%
  print()


# DATOS CV según genotipo completo, parcial y negativo
df %>%
  group_by(gt_cat) %>%
  summarise(
    n         = n(),
    mediana   = median(CV, na.rm = TRUE),
    p25       = quantile(CV, 0.25, na.rm = TRUE),
    p75       = quantile(CV, 0.75, na.rm = TRUE))

