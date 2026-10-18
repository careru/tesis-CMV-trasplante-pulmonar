# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.6.4. Figura 31. Distancia al percentil 75 del Ct en genotipos estables y fluctuantes
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 31. Distancia al percentil 75 del Ct en genotipos estables y fluctuantes ----
ct_cols  <- c("GB1CT","GB2CT","GB3CT","GB4CT","GH1CT","GH2CT")
gt_cols  <- c("GB1","GB2","GB3","GB4","GH1","GH2")
gt_names <- c("gB1","gB2","gB3","gB4","gH1","gH2")

gt_largo <- muestras_mv %>%
  select(NHC, Fini, Freg, all_of(gt_cols)) %>%
  pivot_longer(cols = all_of(gt_cols),
               names_to = "genotipo_var", values_to = "presente") %>%
  filter(presente == 1) %>%
  mutate(genotipo = gt_names[match(genotipo_var, gt_cols)])

ct_largo <- muestras_mv %>%
  select(NHC, Fini, Freg, all_of(ct_cols)) %>%
  pivot_longer(cols = all_of(ct_cols),
               names_to = "ct_var", values_to = "Ct") %>%
  mutate(genotipo_var = gt_cols[match(ct_var, ct_cols)])

det_largo <- gt_largo %>%
  left_join(ct_largo %>% select(NHC, Fini, Freg, genotipo_var, Ct),
            by = c("NHC", "Fini", "Freg", "genotipo_var")) %>%
  filter(!is.na(Ct))

estabilidad_largo <- estabilidad_gt %>%
  select(NHC, Fini, all_of(gt_cols)) %>%
  pivot_longer(cols = all_of(gt_cols),
               names_to = "genotipo_var", values_to = "estabilidad") %>%
  filter(estabilidad != "ausente")

det_largo <- det_largo %>%
  left_join(estabilidad_largo, by = c("NHC", "Fini", "genotipo_var")) %>%
  filter(!is.na(estabilidad))

# P75 por genotipo sobre todas las muestras
todas_ct <- muestras %>%
  select(NHC, Freg, all_of(gt_cols), all_of(ct_cols)) %>%
  pivot_longer(cols = all_of(gt_cols),
               names_to = "genotipo_var", values_to = "presente") %>%
  filter(presente == 1) %>%
  left_join(
    muestras %>%
      select(NHC, Freg, all_of(ct_cols)) %>%
      pivot_longer(cols = all_of(ct_cols),
                   names_to = "ct_var", values_to = "Ct") %>%
      mutate(genotipo_var = gt_cols[match(ct_var, ct_cols)]) %>%
      select(NHC, Freg, genotipo_var, Ct),
    by = c("NHC", "Freg", "genotipo_var")
  ) %>%
  filter(!is.na(Ct))

p75_genotipo <- todas_ct %>%
  group_by(genotipo_var) %>%
  summarise(P75 = quantile(Ct, 0.75, na.rm = TRUE), .groups = "drop")

det_largo <- det_largo %>%
  left_join(p75_genotipo, by = "genotipo_var") %>%
  mutate(dist_p75 = P75 - Ct)

cat("\n── Distancia al P75: estable vs fluctuante ──\n")
det_largo %>%
  group_by(estabilidad) %>%
  summarise(
    n             = n(),
    mediana       = median(dist_p75, na.rm = TRUE),
    q1            = quantile(dist_p75, 0.25, na.rm = TRUE),
    q3            = quantile(dist_p75, 0.75, na.rm = TRUE),
    pct_sobre_p75 = 100 * mean(dist_p75 < 0, na.rm = TRUE)
  ) %>% print()

wilcox.test(dist_p75 ~ estabilidad, data = det_largo, exact = FALSE)

etiq <- det_largo %>%
  group_by(estabilidad) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = paste0(ifelse(estabilidad == "estable",
                               "Genotipos estables", "Genotipos fluctuantes"),
                        "\n(n=", n, ")"))

det_fig <- det_largo %>%
  mutate(
    grupo = ifelse(estabilidad == "estable",
                   "Patrón estable", "Genotipos fluctuantes"),
    grupo = factor(grupo, levels = c("Patrón estable", "Genotipos fluctuantes"))
  )

etiq_x <- etiq %>%
  mutate(grupo = ifelse(estabilidad == "estable",
                        "Patrón estable", "Genotipos fluctuantes"),
         grupo = factor(grupo, levels = c("Patrón estable", "Genotipos fluctuantes")))

figura_R22 <- ggplot(det_fig, aes(x = grupo, y = dist_p75, fill = grupo)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "firebrick3", linewidth = 0.6) +
  annotate("text", x = 2.42, y = 0.8, label = "P75",
           size = 5, color = "firebrick3", fontface = "bold") +
  geom_jitter(aes(color = grupo), width = 0.15, size = 0.9,
              alpha = 0.35, show.legend = FALSE) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA, width = 0.5,
               color = "grey15", linewidth = 0.7) +
  # Mediana resaltada en blanco
  stat_summary(fun = median, geom = "crossbar",
               width = 0.5, fatten = 0,
               color = "white", linewidth = 0.7) +
  annotate("segment", x = 1, xend = 2, y = 9,   yend = 9,
           color = "black", linewidth = 0.5) +
  annotate("segment", x = 1, xend = 1, y = 8.5, yend = 9,
           color = "black", linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2, y = 8.5, yend = 9,
           color = "black", linewidth = 0.5) +
  annotate("text", x = 1.5, y = 9.7, label = "***",
           size = 7, color = "black", fontface = "bold") +
  scale_x_discrete(labels = setNames(etiq_x$label, etiq_x$grupo)) +
  scale_fill_manual(
    values = c("Patrón estable"        = "#2E7D32",
               "Genotipos fluctuantes" = "#E65100")
  ) +
  scale_color_manual(
    values = c("Patrón estable"        = "#2E7D32",
               "Genotipos fluctuantes" = "#E65100")
  ) +
  scale_y_continuous(
    name   = "Distancia al P75 del genotipo (ciclos Ct)",
    breaks = seq(-6, 10, by = 2)
  ) +
  coord_cartesian(ylim = c(-7, 11)) +
  labs(x = NULL) +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "none",
    axis.line       = element_line(color = "black", linewidth = 0.6),
    axis.ticks      = element_line(color = "black", linewidth = 0.5),
    axis.text.x     = element_text(size = 14, color = "black"),
    axis.text.y     = element_text(size = 13, color = "black"),
    axis.title.y    = element_text(size = 14, color = "black", face = "bold"),
    plot.margin     = margin(8, 8, 8, 8)
  )

figura_R22

guardar_fig(figura_R22, "Figura_31", ancho = 12, alto = 13)

