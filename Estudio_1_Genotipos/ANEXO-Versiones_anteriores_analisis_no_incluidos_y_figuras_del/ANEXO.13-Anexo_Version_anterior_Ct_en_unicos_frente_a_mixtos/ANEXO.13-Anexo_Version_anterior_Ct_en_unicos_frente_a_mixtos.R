# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.13. Anexo: Versión anterior: Ct en únicos frente a mixtos
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Versión anterior: Ct en únicos frente a mixtos ----
# 1. Crear base con número de genotipos gB y gH
df_ct_tipo <- muestras_sel %>%
  mutate(
    gb_sum = rowSums(across(GB1:GB4, ~replace_na(., 0))),
    gh_sum = rowSums(across(GH1:GH2, ~replace_na(., 0))),
    gB_tipo = if_else(gb_sum == 1, "gB único", "gB mixto"),
    gH_tipo = if_else(gh_sum == 1, "gH único", "gH mixto"),
    grupo = case_when(
      gb_sum == 1 & gh_sum == 1 ~ "gB+gH únicos",
      gb_sum >= 2 & gh_sum == 1 ~ "gB mixto / gH único",
      gb_sum == 1 & gh_sum >= 2 ~ "gB único / gH mixto",
      gb_sum >= 2 & gh_sum >= 2 ~ "gB+gH mixtos"
    )
  ) %>%
  select(GB1CT, GB2CT, GB3CT, GB4CT, GH1CT, GH2CT, gB_tipo, gH_tipo, grupo)

# 2. Crear variable "CT medio" para gB y gH
df_ct_tipo <- df_ct_tipo %>%
  rowwise() %>%
  mutate(
    CT_gB = median(c(GB1CT, GB2CT, GB3CT, GB4CT), na.rm = TRUE),
    CT_gH = median(c(GH1CT, GH2CT), na.rm = TRUE)
  ) %>%
  ungroup()

# 3. Dataset largo para los tres paneles
df_plot_ct <- bind_rows(
  df_ct_tipo %>% select(CT = CT_gB, categoria = gB_tipo) %>%
    mutate(seccion = "gB (único vs mixto)"),
  df_ct_tipo %>% select(CT = CT_gH, categoria = gH_tipo) %>%
    mutate(seccion = "gH (único vs mixto)"),
  df_ct_tipo %>%
    mutate(
      categoria = if_else(grepl("mixto", grupo), "gB+gH mixtos", "gB+gH únicos")
    ) %>%
    select(CT = CT_gB, categoria) %>%
    mutate(seccion = "gB+gH (único vs mixto)")
)

# 4. Orden de secciones
sec_levels <- c("gB (único vs mixto)", "gH (único vs mixto)", "gB+gH (único vs mixto)")

df_plot_ct <- df_plot_ct %>%
  mutate(
    seccion = factor(seccion, levels = sec_levels),
    categoria = factor(categoria,
                       levels = c("gB único", "gB mixto",
                                  "gH único", "gH mixto",
                                  "gB+gH únicos", "gB+gH mixtos"))
  )

# 5. Colores y estilo coherentes con tu figura de CV
cols <- c(
  "gB único"       = "#008B00",
  "gB mixto"       = "#CD5555",
  "gH único"       = "#008B00",
  "gH mixto"       = "#CD5555",
  "gB+gH únicos"   = "#008B00",
  "gB+gH mixtos"   = "#CD5555"
)

# 6. Comparaciones estadísticas: Wilcoxon
stats_ct <- df_plot_ct %>%
  group_by(seccion) %>%
  wilcox_test(CT ~ categoria) %>%
  add_significance("p") %>%
  mutate(
    y.position = ave(rep(max(df_plot_ct$CT, na.rm = TRUE), n()), seccion, FUN = function(y) y + 0.8)
  )

# 7. Boxplot
ggplot(df_plot_ct, aes(x = categoria, y = CT, fill = categoria)) +
  geom_boxplot(
    width = 0.55,
    color = "black",
    size = 0.4,
    outlier.shape = 21,
    outlier.size = 1.8,
    outlier.stroke = 0,
    outlier.alpha = 0.6,
    na.rm = TRUE
  ) +
  facet_grid(~ seccion, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = cols) +
  coord_cartesian(ylim = c(20, 50)) +
  labs(
    x = NULL,
    y = "Valor CT",
    title = "Comparación de valores CT según genotipos únicos y mixtos"
  ) +
  stat_pvalue_manual(
    stats_ct,
    label = "p.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.45,
    size = 3.6,
    inherit.aes = FALSE
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.4),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = grid::unit(16, "pt")
  )

