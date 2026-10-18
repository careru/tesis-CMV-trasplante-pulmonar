# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.16. Anexo: Figura de mixtos por seroestatus en el primer episodio (no incluida)
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Figura de mixtos por seroestatus en el primer episodio (no incluida) ----
# ── df ya definido arriba con filter(primer_mil == 1) ─────────────────────────

# ── 1. Etiquetas con n= ───────────────────────────────────────────────────────
n_total <- nrow(df)
n_DR    <- sum(df$serDR == "D+/R-",  na.rm = TRUE)
n_DpR   <- sum(df$serDR == "D+/R+",  na.rm = TRUE)
n_DmR   <- sum(df$serDR == "D-/R+",  na.rm = TRUE)

etiq_x <- c(
  "Global" = paste0("Global\n(n = ", n_total, ")"),
  "D+/R-"  = paste0("D+/R-\n(n = ", n_DR,    ")"),
  "D+/R+"  = paste0("D+/R+\n(n = ", n_DpR,   ")"),
  "D-/R+"  = paste0("D-/R+\n(n = ", n_DmR,   ")")
)

# ── 2a. Por seroestatus ───────────────────────────────────────────────────────
tab_ser <- df %>%
  group_by(serDR) %>%
  summarise(
    episodios     = n(),
    `gB mixto`    = 100 * sum(gb_mixto,   na.rm = TRUE) / episodios,
    `gH mixto`    = if (all(is.na(gh_mixto)))    NA_real_ else
      100 * sum(gh_mixto,   na.rm = TRUE) / episodios,
    `gB/gH mixto` = if (all(is.na(gbgh_mixto))) NA_real_ else
      100 * sum(gbgh_mixto, na.rm = TRUE) / episodios,
    .groups = "drop"
  )

# ── 2b. Global ────────────────────────────────────────────────────────────────
tab_global <- df %>%
  summarise(
    serDR         = "Global",
    episodios     = n(),
    `gB mixto`    = 100 * sum(gb_mixto,   na.rm = TRUE) / episodios,
    `gH mixto`    = if (all(is.na(gh_mixto)))    NA_real_ else
      100 * sum(gh_mixto,   na.rm = TRUE) / episodios,
    `gB/gH mixto` = if (all(is.na(gbgh_mixto))) NA_real_ else
      100 * sum(gbgh_mixto, na.rm = TRUE) / episodios
  )

# ── 2c. Combinar y pivotar ────────────────────────────────────────────────────
tab_fig <- bind_rows(tab_global, tab_ser) %>%
  pivot_longer(
    cols      = c(`gB mixto`, `gH mixto`, `gB/gH mixto`),
    names_to  = "tipo_genotipo",
    values_to = "porcentaje"
  ) %>%
  filter(!is.na(porcentaje)) %>%
  mutate(
    serDR = factor(serDR, levels = c("Global", "D+/R-", "D+/R+", "D-/R+")),
    tipo_genotipo = factor(tipo_genotipo,
                           levels = c("gB mixto", "gH mixto", "gB/gH mixto"))
  )

# ── 3. Figura R15b (primeros episodios) ───────────────────────────────────────
figura_R15b <- ggplot(tab_fig,
                      aes(x = serDR, y = porcentaje, fill = tipo_genotipo)) +
  geom_col(
    position  = position_dodge(width = 0.75),
    width     = 0.65,
    color     = "white",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 1.5,
    color      = "grey75",
    linewidth  = 0.4,
    linetype   = "dashed"
  ) +
  geom_text(
    aes(label = paste0(formatC(porcentaje, digits = 1, format = "f",
                               decimal.mark = ","), "%")),
    position  = position_dodge(width = 0.75),
    vjust     = -0.5,
    size      = 4,
    color     = "grey20",
    na.rm     = TRUE
  ) +
  scale_fill_manual(
    values = c(
      "gB mixto"    = "#56B4E9",
      "gH mixto"    = "palevioletred1",
      "gB/gH mixto" = "firebrick3"
    ),
    name = NULL
  ) +
  scale_x_discrete(labels = etiq_x) +
  scale_y_continuous(
    limits = c(0, 50),
    expand = expansion(mult = c(0, 0.08)),
    breaks = seq(0, 50, 10),
    labels = scales::label_number(suffix = "%", decimal.mark = ",")
  ) +
  labs(x = "", y = "Frecuencia (%)") +
  theme_minimal(base_size = 13) +
  theme(
    panel.border       = element_blank(),
    axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
    axis.line.y        = element_line(color = "grey60", linewidth = 0.4),
    axis.ticks.x       = element_blank(),
    axis.ticks.y       = element_line(color = "grey60", linewidth = 0.3),
    axis.text          = element_text(color = "grey30", size = 13),
    axis.title         = element_text(color = "grey20", size = 13),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "white", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    legend.position    = "top",
    legend.text        = element_text(size = 11),
    legend.key.size    = unit(0.45, "cm"),
    legend.spacing.x   = unit(0.4, "cm"),
    plot.margin        = margin(t = 8, r = 10, b = 4, l = 4)
  )

figura_R15b


  

