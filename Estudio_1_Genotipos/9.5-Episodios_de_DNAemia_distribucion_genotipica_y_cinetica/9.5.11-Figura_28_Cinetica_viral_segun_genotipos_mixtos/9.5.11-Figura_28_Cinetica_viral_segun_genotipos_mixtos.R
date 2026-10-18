# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.5.11. Figura 28. Cinética viral según genotipos mixtos, estratificada por perfil serológico
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 28. Cinética viral según genotipos mixtos, estratificada por perfil serológico ----
# ── Variables de control ──────────────────────────────────────────────────────
id_patient <- "NHC"
cv_y_max   <- 7

## Etiquetas de las filas (orden de arriba abajo, como en la tesis)
## Se escriben como expresiones plotmath (log[10] = subíndice) y se
## interpretan con label_parsed en facet_grid()
## atop() reparte cada etiqueta en dos líneas para que quepa en la cabecera
lab_acl  <- 'atop(bold("Aclaramiento"), bold("(días)"))'
lab_auc  <- 'atop(bold("AUC"), bold("(log"[10]*"(UI/mL)·día)"))'
lab_cvi  <- 'atop(bold("log"[10]*" CV inicial"), bold("(UI/mL)"))'
lab_cvp  <- 'atop(bold("log"[10]*" CV pico"), bold("(UI/mL)"))'
orden_var <- c(lab_acl, lab_auc, lab_cvi, lab_cvp)

df_sub <- dnaemias

vars_to_pivot <- c("diasdur_vir", "AUCCV", "cvpico", "cvini")

# ── Formato largo ─────────────────────────────────────────────────────────────
## case_when() en lugar de recode(): car::recode() enmascara a dplyr::recode()
df_long <- df_sub %>%
  mutate(
    gB   = ifelse(gbmix   == 1, "gB mixto",    "gB único"),
    gH   = ifelse(ghmix   == 1, "gH mixto",    "gH único"),
    gBGH = ifelse(gbghmix == 1, "gB+gH mixto", "gB+gH único")
  ) %>%
  pivot_longer(cols = c(gB, gH, gBGH),
               names_to  = "tipo_gen",
               values_to = "genotipo") %>%
  select(serDR, tipo_gen, genotipo,
         !!sym(id_patient), all_of(vars_to_pivot)) %>%
  pivot_longer(cols      = all_of(vars_to_pivot),
               names_to  = "variable",
               values_to = "valor_raw") %>%
  mutate(
    valor = case_when(
      variable %in% c("cvpico", "cvini") & !is.na(valor_raw) ~
        log10(as.numeric(valor_raw) + 1),
      TRUE ~ as.numeric(valor_raw)
    ),
    variable = case_when(
      variable == "diasdur_vir" ~ lab_acl,
      variable == "AUCCV"       ~ lab_auc,
      variable == "cvini"       ~ lab_cvi,
      variable == "cvpico"      ~ lab_cvp
    ),
    serDR    = factor(serDR, levels = c("D-/R+", "D+/R+", "D+/R-")),
    genotipo = factor(genotipo, levels = c(
      "gB único",    "gB mixto",
      "gH único",    "gH mixto",
      "gB+gH único", "gB+gH mixto"
    ))
  )

# ── Posiciones X ──────────────────────────────────────────────────────────────
x_map <- c(
  "gB único"    = 1,   "gB mixto"    = 2,
  "gH único"    = 3.3, "gH mixto"    = 4.3,
  "gB+gH único" = 5.6, "gB+gH mixto" = 6.6
)
df_plot <- df_long %>%
  mutate(x_pos = x_map[as.character(genotipo)])

# ── Límites visuales por panel ────────────────────────────────────────────────
caps_dyn <- tibble(
  variable = orden_var,
  cap_vis  = c(100, 400, cv_y_max, cv_y_max)
)

df_plot_vis <- df_plot %>%
  left_join(caps_dyn, by = "variable") %>%
  filter(is.na(valor) | valor <= cap_vis) %>%
  mutate(variable = factor(variable, levels = orden_var))

# ── Comparaciones Wilcoxon sobre datos SIN filtrar (solo significativas) ──────
comparaciones <- df_plot %>%
  group_by(variable, serDR, tipo_gen) %>%
  wilcox_test(valor ~ genotipo) %>%
  add_significance("p") %>%
  ungroup() %>%
  filter(p.signif != "ns") %>%
  mutate(
    group1 = case_when(
      tipo_gen == "gB"   ~ "gB único",
      tipo_gen == "gH"   ~ "gH único",
      tipo_gen == "gBGH" ~ "gB+gH único"
    ),
    group2 = case_when(
      tipo_gen == "gB"   ~ "gB mixto",
      tipo_gen == "gH"   ~ "gH mixto",
      tipo_gen == "gBGH" ~ "gB+gH mixto"
    ),
    xmin = x_map[group1],
    xmax = x_map[group2]
  )

base_positions <- df_plot_vis %>%
  mutate(variable = as.character(variable)) %>%
  group_by(variable, serDR) %>%
  summarise(maxval  = max(valor,   na.rm = TRUE),
            cap_vis = first(cap_vis), .groups = "drop") %>%
  mutate(base = pmin(maxval * 1.12, cap_vis * 0.90))

comp_pos <- comparaciones %>%
  left_join(base_positions, by = c("variable", "serDR")) %>%
  group_by(variable, serDR) %>%
  arrange(p) %>%
  mutate(idx = row_number(),
         y   = base + (idx - 1) * 0.05 * cap_vis) %>%
  ungroup() %>%
  mutate(variable = factor(variable, levels = orden_var))

print(comparaciones %>% select(variable, serDR, tipo_gen, n1, n2, p, p.signif))

# ── Etiquetas serológicas con n ───────────────────────────────────────────────
counts_serDR <- df_sub %>%
  group_by(serDR) %>%
  summarise(
    episodios = n(),
    pacientes = n_distinct(.data[[id_patient]]),
    .groups   = "drop"
  ) %>%
  mutate(label = paste0(
    "<b>", serDR, "</b><br>",
    "<span style='font-size:10pt'>",
    pacientes, " pac. / ", episodios, " ep.</span>"
  ))

serdr_labels        <- counts_serDR$label
names(serdr_labels) <- counts_serDR$serDR

# ── Colores ───────────────────────────────────────────────────────────────────
colores <- c(
  "gB único"    = "#63B8FF",
  "gB mixto"    = "dodgerblue3",
  "gH único"    = "#DDA0DD",
  "gH mixto"    = "#CD69C9",
  "gB+gH único" = "#3CB371",
  "gB+gH mixto" = "#F08080"
)

# ── Figura ────────────────────────────────────────────────────────────────────
p_global <- ggplot(df_plot_vis,
                   aes(x = x_pos, y = valor, fill = genotipo)) +
  geom_boxplot(
    width         = 0.55,
    color         = "grey25",
    linewidth     = 0.4,
    outlier.shape = 21,
    outlier.size  = 1.4,
    outlier.fill  = "white",
    outlier.color = "grey50",
    outlier.alpha = 0.6
  ) +
  geom_vline(xintercept = c(2.65, 4.95),
             color = "grey80", linewidth = 0.35, linetype = "dashed") +
  facet_grid(variable ~ serDR, scales = "free_y",
             labeller = labeller(serDR    = as_labeller(serdr_labels),
                                 variable = label_parsed)) +
  geom_text(data = comp_pos,
            aes(x = (xmin + xmax) / 2, y = y, label = p.signif),
            inherit.aes = FALSE, size = 3.8, vjust = 0, color = "grey20") +
  geom_segment(data = comp_pos,
               aes(x    = xmin, xend = xmax,
                   y    = y - 0.012 * cap_vis,
                   yend = y - 0.012 * cap_vis),
               inherit.aes = FALSE,
               linewidth = 0.45, color = "grey30") +
  scale_fill_manual(values = colores) +
  scale_x_continuous(
    breaks = unname(x_map),
    labels = names(x_map),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position    = "none",
    panel.border       = element_rect(color = "grey70", fill = NA,
                                      linewidth = 0.4),
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.line          = element_blank(),
    axis.ticks         = element_line(color = "grey60", linewidth = 0.3),
    axis.text.x        = element_text(size = 10, color = "grey25",
                                      angle = 45, hjust = 1),
    axis.text.y        = element_text(size = 12, color = "grey30"),
    strip.background   = element_rect(fill = "grey95", color = "grey70",
                                      linewidth = 0.4),
    strip.text.x       = element_markdown(size = 13, face = "bold",
                                          color = "grey20"),
    strip.text.y       = element_text(size = 11, face = "bold",
                                      color = "grey20"),
    panel.spacing.x    = unit(6,  "pt"),
    panel.spacing.y    = unit(10, "pt"),
    plot.margin        = margin(8, 10, 6, 6)
  )

p_global

# ── Guardar en SVG ────────────────────────────────────────────────────────────
# [guardado antiguo] ggsave(
# [guardado antiguo]   filename      = "~/Desktop/Figuras_tesis/Figura_28_cinetica_serologia.svg",
# [guardado antiguo]   plot          = p_global,
# [guardado antiguo]   device        = svglite::svglite,
# [guardado antiguo]   fix_text_size = FALSE,
# [guardado antiguo]   width         = 16,
# [guardado antiguo]   height        = 22,
# [guardado antiguo]   units         = "cm",
# [guardado antiguo]   bg            = "white"
# [guardado antiguo] )

guardar_fig(p_global, "Figura_28", ancho = 16, alto = 22)

