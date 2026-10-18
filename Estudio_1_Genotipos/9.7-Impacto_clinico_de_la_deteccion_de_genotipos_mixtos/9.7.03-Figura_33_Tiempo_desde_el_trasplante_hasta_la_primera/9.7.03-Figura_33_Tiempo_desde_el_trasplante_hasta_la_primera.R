# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.7.03. Figura 33. Tiempo desde el trasplante hasta la primera DNAemia genotipada
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 33. Tiempo desde el trasplante hasta la primera DNAemia genotipada ----
## Panel A: Genotipos mixtos a lo largo del seguimiento (nivel paciente)
##          -> variables gbm/ghm/gbghm desde tabla `pacientes`
## Panel B: Genotipos mixtos en el primer episodio (nivel episodio)
##          -> variables gbmix/ghmix/gbghmix desde tabla `dnaemias`
## Filtro: primer_ev = 1
## Orden de paneles: D+/R- -> D+/R+ -> D-/R+


DODGE <- 0.7

#  Paletas (compartidas)
cols_fill <- c(
  "gB_\u00danico"     = "#AED6F1",  "gB_Mixto"     = "#4682B4",
  "gH_\u00danico"     = "pink",     "gH_Mixto"     = "hotpink3",
  "gB+gH_\u00danico"  = "#4CAF7D",  "gB+gH_Mixto"  = "#E07070"
)

cols_jit <- c(
  "gB_\u00danico"     = "#2980B9",  "gB_Mixto"     = "#0D2B45",
  "gH_\u00danico"     = "pink3",    "gH_Mixto"     = "palevioletred4",
  "gB+gH_\u00danico"  = "#2E7D52",  "gB+gH_Mixto"  = "#B04040"
)

df_prof <- tibble::tibble(
  serDR = factor(c("D+/R-", "D+/R+", "D-/R+"),
                 levels = c("D+/R-", "D+/R+", "D-/R+")),
  prof_meses = c(12, 6, 6)
)

#  FUNCION CONSTRUCTORA DEL PANEL
construir_panel <- function(df_long, etiquetas_riesgo, titulo_panel) {
  
  # Tests significativos del panel
  stat_test <- df_long %>%
    dplyr::group_by(serDR, glucoproteina) %>%
    rstatix::wilcox_test(meses_tx_inicio ~ genotipo) %>%
    rstatix::add_significance("p") %>%
    dplyr::ungroup() %>%
    dplyr::filter(p.signif != "ns") %>%
    dplyr::mutate(
      x_pos      = as.numeric(glucoproteina),
      xmin       = x_pos - DODGE / 4,
      xmax       = x_pos + DODGE / 4,
      y.position = 22
    )
  
  ggplot(df_long,
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
        ggpubr::stat_pvalue_manual(
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
      y = "Meses desde TX a 1er episodio",
      title = titulo_panel
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
}

#  PANEL A - Nivel PACIENTE (gbm, ghm, gbghm desde tabla pacientes)
df_fig_A <- dnaemias %>%
  dplyr::filter(primer_ev == 1) %>%
  dplyr::mutate(
    Ftx              = as.Date(Ftx),
    Fini             = as.Date(Fini),
    meses_tx_inicio  = as.numeric(Fini - Ftx) / 30.44
  ) %>%
  dplyr::select(-any_of(c("gbm", "ghm", "gbghm"))) %>%
  dplyr::left_join(pacientes %>% dplyr::select(NHC, gbm, ghm, gbghm), by = "NHC") %>%
  dplyr::mutate(
    serDR = factor(serDR, levels = c("D+/R-", "D+/R+", "D-/R+"))
  ) %>%
  dplyr::filter(!is.na(serDR), !is.na(meses_tx_inicio))

df_long_A <- df_fig_A %>%
  dplyr::select(NHC, serDR, meses_tx_inicio, gbm, ghm, gbghm) %>%
  tidyr::pivot_longer(
    cols      = c(gbm, ghm, gbghm),
    names_to  = "glucoproteina",
    values_to = "mixto"
  ) %>%
  dplyr::filter(!is.na(mixto)) %>%
  dplyr::mutate(
    glucoproteina = factor(glucoproteina,
                           levels = c("gbm", "ghm", "gbghm"),
                           labels = c("gB", "gH", "gB+gH")),
    genotipo      = factor(ifelse(mixto == 1, "Mixto", "\u00danico"),
                           levels = c("\u00danico", "Mixto")),
    cat_color = paste(glucoproteina, genotipo, sep = "_"),
    cat_color = factor(cat_color,
                       levels = c("gB_\u00danico",   "gB_Mixto",
                                  "gH_\u00danico",   "gH_Mixto",
                                  "gB+gH_\u00danico","gB+gH_Mixto"))
  )

counts_A <- df_fig_A %>%
  dplyr::group_by(serDR) %>%
  dplyr::summarise(n = dplyr::n_distinct(NHC), .groups = "drop")

etiquetas_A <- setNames(
  sprintf("%s (n = %d)", counts_A$serDR, counts_A$n),
  as.character(counts_A$serDR)
)

panel_A <- construir_panel(
  df_long          = df_long_A,
  etiquetas_riesgo = etiquetas_A,
  titulo_panel     = "A. Genotipos únicos/mixtos a lo largo del seguimiento"
)

#  PANEL B - Nivel EPISODIO (gbmix, ghmix, gbghmix desde tabla dnaemias)
df_fig_B <- dnaemias %>%
  dplyr::filter(primer_ev == 1) %>%
  dplyr::mutate(
    Ftx              = as.Date(Ftx),
    Fini             = as.Date(Fini),
    meses_tx_inicio  = as.numeric(Fini - Ftx) / 30.44,
    serDR = factor(serDR, levels = c("D+/R-", "D+/R+", "D-/R+"))
  ) %>%
  dplyr::filter(!is.na(serDR), !is.na(meses_tx_inicio))

df_long_B <- df_fig_B %>%
  dplyr::select(NHC, serDR, meses_tx_inicio, gbmix, ghmix, gbghmix) %>%
  tidyr::pivot_longer(
    cols      = c(gbmix, ghmix, gbghmix),
    names_to  = "glucoproteina",
    values_to = "mixto"
  ) %>%
  dplyr::filter(!is.na(mixto)) %>%
  dplyr::mutate(
    glucoproteina = factor(glucoproteina,
                           levels = c("gbmix", "ghmix", "gbghmix"),
                           labels = c("gB", "gH", "gB+gH")),
    genotipo      = factor(ifelse(mixto == 1, "Mixto", "\u00danico"),
                           levels = c("\u00danico", "Mixto")),
    cat_color = paste(glucoproteina, genotipo, sep = "_"),
    cat_color = factor(cat_color,
                       levels = c("gB_\u00danico",   "gB_Mixto",
                                  "gH_\u00danico",   "gH_Mixto",
                                  "gB+gH_\u00danico","gB+gH_Mixto"))
  )

counts_B <- df_fig_B %>%
  dplyr::group_by(serDR) %>%
  dplyr::summarise(n = dplyr::n_distinct(NHC), .groups = "drop")

etiquetas_B <- setNames(
  sprintf("%s (n = %d)", counts_B$serDR, counts_B$n),
  as.character(counts_B$serDR)
)

panel_B <- construir_panel(
  df_long          = df_long_B,
  etiquetas_riesgo = etiquetas_B,
  titulo_panel     = "B. Genotipos únicos/mixtos en el primer episodio"
)

#  COMBINAR PANELES A y B
figura_R20_doble <- panel_A / panel_B +
  patchwork::plot_layout(heights = c(1, 1))

figura_R20_doble

# [guardado antiguo] ggsave("~/Desktop/Figuras_tesis/figura_R20_doble_paciente_vs_episodio.png",
# [guardado antiguo]        figura_R20_doble,
# [guardado antiguo]        width = 10, height = 9, dpi = 320, bg = "white")

guardar_fig(figura_R20_doble, "Figura_33", ancho = 25.4, alto = 22.9)

