# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.4.09. Figura 24. Carga viral según genotipos únicos y mixtos de gB, gH y gB+gH
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 24. Carga viral según genotipos únicos y mixtos de gB, gH y gB+gH ----
## Base (solo genotipado completo gB y gH presentes)
dfc <- muestras_sel %>%
  mutate(
    gb_sum  = rowSums(across(GB1:GB4, ~ tidyr::replace_na(., 0))),
    gh_sum  = rowSums(across(GH1:GH2, ~ tidyr::replace_na(., 0))),
    CVlog10 = log10(CV)
  ) %>%
  filter(gb_sum >= 1, gh_sum >= 1) %>%
  drop_na(CV, CVlog10)

## Secciones
gb_ind <- dfc %>%
  filter(gb_sum == 1) %>%
  pivot_longer(GB1:GB4, names_to = "categoria_raw", values_to = "v") %>%
  filter(v == 1) %>%
  transmute(
    seccion   = "gB (individuales únicos)",
    categoria = paste0("g", substr(categoria_raw, 2, nchar(categoria_raw))),
    CV, CVlog10
  )

gh_ind <- dfc %>%
  filter(gh_sum == 1) %>%
  pivot_longer(GH1:GH2, names_to = "categoria_raw", values_to = "v") %>%
  filter(v == 1) %>%
  transmute(
    seccion   = "gH (individuales únicos)",
    categoria = paste0("g", substr(categoria_raw, 2, nchar(categoria_raw))),
    CV, CVlog10
  )

gb_cat <- dfc %>%
  transmute(
    seccion   = "gB (único vs mixto)",
    categoria = ifelse(gb_sum == 1, "gB único", "gB mixto"),
    CV, CVlog10
  )

gh_cat <- dfc %>%
  transmute(
    seccion   = "gH (único vs mixto)",
    categoria = ifelse(gh_sum == 1, "gH único", "gH mixto"),
    CV, CVlog10
  )

gbgh_cat <- dfc %>%
  transmute(
    seccion   = "gB+gH (único vs mixto)",
    categoria = ifelse(gb_sum == 1 & gh_sum == 1, "gB+gH único", "gB+gH mixto"),
    CV, CVlog10
  )

## Niveles
sec_levels <- c(
  "gB (individuales únicos)",
  "gH (individuales únicos)",
  "gB (único vs mixto)",
  "gH (único vs mixto)",
  "gB+gH (único vs mixto)"
)

cat_levels <- c(
  "gB1", "gB2", "gB3", "gB4",
  "gH1", "gH2",
  "gB único", "gB mixto",
  "gH único", "gH mixto",
  "gB+gH único", "gB+gH mixto"
)

df_plot <- bind_rows(gb_ind, gh_ind, gb_cat, gh_cat, gbgh_cat) %>%
  mutate(
    seccion   = factor(seccion, levels = sec_levels),
    categoria = factor(categoria, levels = cat_levels)
  )

## Paleta refinada
cols <- c(
  "gB1"         = "#5B9BD5",
  "gB2"         = "#5B9BD5",
  "gB3"         = "#5B9BD5",
  "gB4"         = "#5B9BD5",
  "gH1"         = "#D98FD6",
  "gH2"         = "#D98FD6",
  "gB único"    = "#4CAF7D",
  "gB mixto"    = "#E07070",
  "gH único"    = "#4CAF7D",
  "gH mixto"    = "#E07070",
  "gB+gH único" = "#4CAF7D",
  "gB+gH mixto" = "#E07070"
)

## Estadísticas
max_log_by_seccion <- function(data, sec_label) {
  data %>% filter(seccion == sec_label) %>% pull(CVlog10) %>% max(na.rm = TRUE)
}

st_gb_ind <- gb_ind %>%
  pairwise_wilcox_test(CV ~ categoria, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(seccion = sec_levels[1]) %>%
  arrange(group1, group2)

max_gb_log <- max_log_by_seccion(df_plot, sec_levels[1])
st_gb_ind$y.position <- max_gb_log + seq(0.12, by = 0.12, length.out = nrow(st_gb_ind))

st_gh_ind <- gh_ind %>%
  wilcox_test(CV ~ categoria) %>%
  add_significance("p") %>%
  mutate(seccion = sec_levels[2], y.position = max_log_by_seccion(df_plot, sec_levels[2]) + 0.12)

st_gb <- gb_cat %>%
  wilcox_test(CV ~ categoria) %>%
  add_significance("p") %>%
  mutate(seccion = sec_levels[3], y.position = max_log_by_seccion(df_plot, sec_levels[3]) + 0.12)

st_gh <- gh_cat %>%
  wilcox_test(CV ~ categoria) %>%
  add_significance("p") %>%
  mutate(seccion = sec_levels[4], y.position = max_log_by_seccion(df_plot, sec_levels[4]) + 0.12)

st_gbgh <- gbgh_cat %>%
  wilcox_test(CV ~ categoria) %>%
  add_significance("p") %>%
  mutate(seccion = sec_levels[5], y.position = max_log_by_seccion(df_plot, sec_levels[5]) + 0.12)

stats_all <- bind_rows(st_gb_ind, st_gh_ind, st_gb, st_gh, st_gbgh) %>%
  mutate(
    p.signif = dplyr::coalesce(p.signif, p.adj.signif),
    seccion  = factor(seccion, levels = sec_levels)
  )

stats_sig <- stats_all %>% filter(p.signif != "ns")

## Plot
p <- ggplot(df_plot, aes(categoria, CVlog10, fill = categoria)) +
  geom_boxplot(
    width          = 0.55,
    color          = "grey20",
    linewidth      = 0.42,
    outlier.shape  = 21,
    outlier.size   = 1.7,
    outlier.fill   = "white",
    outlier.color  = "grey50",
    outlier.stroke = 0.35,
    outlier.alpha  = 0.75
  ) +
  geom_jitter(
    aes(color = categoria),
    width = 0.10, size = 0.9, alpha = 0.30, show.legend = FALSE
  ) +
  facet_grid(
    ~ seccion,
    scales = "free_x",
    space  = "free_x",
    switch = "x",
    drop   = FALSE
  ) +
  scale_fill_manual(values = cols, breaks = names(cols)) +
  scale_color_manual(
    values = c(
      "gB1" = "#1A5FA8", "gB2" = "#1A5FA8", "gB3" = "#1A5FA8", "gB4" = "#1A5FA8",
      "gH1" = "#A04EA0", "gH2" = "#A04EA0",
      "gB único" = "#2E7D52", "gB mixto" = "#B04040",
      "gH único" = "#2E7D52", "gH mixto" = "#B04040",
      "gB+gH único" = "#2E7D52", "gB+gH mixto" = "#B04040"
    )
  ) +
  labs(x = NULL, y = "Carga viral CMV Log10 (UI/mL)") +
  scale_y_continuous(
    expand = expansion(mult = c(0.03, 0.12)),
    breaks = scales::pretty_breaks(n = 5)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "none",
    # Grid: solo horizontales suaves, sin verticales
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    # Cuadrante con borde fino y fondo muy levemente tintado
    panel.border       = element_rect(color = "grey55", fill = NA, linewidth = 0.45),
    panel.background   = element_rect(fill = "#FAFAFA", color = NA),
    # Strip labels debajo
    strip.placement    = "outside",
    strip.background   = element_rect(fill = "grey90", color = "grey70", linewidth = 0.35),
    strip.text         = element_text(size = 10.5, color = "grey20", face = "bold",
                                      margin = margin(t = 4, b = 4)),
    # Ejes
    axis.line          = element_blank(),   # el borde del panel ya hace de eje
    axis.ticks         = element_line(color = "grey50", linewidth = 0.35),
    axis.ticks.length  = unit(3, "pt"),
    axis.text.x        = element_text(size = 10, color = "grey20", margin = margin(t = 3)),
    axis.text.y        = element_text(size = 10, color = "grey30"),
    axis.title.y       = element_text(size = 11, color = "grey20", margin = margin(r = 8)),
    # Separación entre paneles
    panel.spacing.x    = grid::unit(10, "pt"),
    plot.margin        = margin(10, 15, 8, 8)
  ) +
  coord_cartesian(clip = "off") +
  stat_pvalue_manual(
    stats_sig,
    label        = "p.signif",
    xmin         = "group1",
    xmax         = "group2",
    y.position   = "y.position",
    tip.length   = 0.012,
    bracket.size = 0.42,
    size         = 3.8,
    color        = "grey25",
    inherit.aes  = FALSE
  )

p

## Resumen numérico
resumen_CV_gt <- df_plot %>%
  group_by(seccion, categoria) %>%
  summarise(
    n   = dplyr::n(),
    med = median(CV, na.rm = TRUE),
    p25 = quantile(CV, .25, na.rm = TRUE),
    p75 = quantile(CV, .75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(seccion, categoria) %>%
  mutate(
    med_fmt = scales::number(med, accuracy = 1, big.mark = ".", decimal.mark = ","),
    p25_fmt = scales::number(p25, accuracy = 1, big.mark = ".", decimal.mark = ","),
    p75_fmt = scales::number(p75, accuracy = 1, big.mark = ".", decimal.mark = ","),
    resumen_texto = paste0(
      "mediana de ", med_fmt,
      " UI/mL (IQR: ", p25_fmt, "–", p75_fmt,
      "; n = ", n, ")"
    )
  )

resumen_CV_gt

## Tabla de p-valores
tabla_pvalores <- stats_all %>%
  select(seccion, group1, group2, p, p.adj, p.signif) %>%
  mutate(
    p_lab       = ifelse(is.na(p), NA_character_,
                         formatC(p, format = "f", digits = 3, decimal.mark = ",")),
    p_adj_lab   = ifelse(is.na(p.adj), NA_character_,
                         formatC(p.adj, format = "f", digits = 3, decimal.mark = ",")),
    comparacion = paste0(group1, " vs ", group2)
  ) %>%
  select(seccion, comparacion, p_lab, p_adj_lab, p.signif)

tabla_pvalores

guardar_fig(p, "Figura_24", ancho = 22, alto = 12)

