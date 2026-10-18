# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.7.07. Figura 36. AUC acumulada de exposición viral por paciente
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 36. AUC acumulada de exposición viral por paciente ----
# ── 1. Episodios: todos los de dnaemias, excluyendo los sin aclaramiento ──────
n_inicial <- nrow(dnaemias)

df_ep <- dnaemias %>%
  filter(!is.na(aclaramiento_dias), !is.na(AUCCV))

n_excluidos <- n_inicial - nrow(df_ep)

cat(sprintf("Episodios totales en dnaemias: %d\n", n_inicial))
cat(sprintf("Excluidos (sin aclaramiento/AUC): %d\n", n_excluidos))
cat(sprintf("Episodios analizados: %d\n", nrow(df_ep)))

# ── 2. AUC acumulada por paciente (suma de AUCCV de sus episodios) ────────────
auc_pac <- df_ep %>%
  group_by(NHC) %>%
  summarise(AUCa = sum(AUCCV, na.rm = TRUE), .groups = "drop") %>%
  left_join(pacientes %>% select(NHC, gbm, ghm, gbghm), by = "NHC")

cat(sprintf("Pacientes con AUC acumulada: %d\n", nrow(auc_pac)))

# ── 3. Formato largo: una fila por paciente x glucoproteina ──────────────────
df_long <- auc_pac %>%
  pivot_longer(c(gbm, ghm, gbghm),
               names_to = "gen", values_to = "mixto") %>%
  filter(!is.na(mixto)) %>%
  mutate(
    glucoproteina = factor(gen,
                           levels = c("gbm", "ghm", "gbghm"),
                           labels = c("gB", "gH", "gB/gH")),
    genotipo  = factor(ifelse(mixto == 1, "Mixtos", "\u00danicos"),
                       levels = c("\u00danicos", "Mixtos")),
    cat_color = paste(glucoproteina, ifelse(mixto == 1, "Mixto", "\u00danico"),
                      sep = "_")
  )

# ── 4. Medianas, RIC y test de Wilcoxon por glucoproteina ────────────────────
resumen <- df_long %>%
  group_by(glucoproteina, genotipo) %>%
  summarise(
    n       = n(),
    mediana = round(median(AUCa)),
    q1      = round(quantile(AUCa, 0.25)),
    q3      = round(quantile(AUCa, 0.75)),
    .groups = "drop"
  ) %>%
  mutate(resumen = sprintf("%d (RIC %d\u2013%d)", mediana, q1, q3))

cat("\n── Medianas (RIC) por grupo ──\n"); print(resumen)

stat_test <- df_long %>%
  group_by(glucoproteina) %>%
  rstatix::wilcox_test(AUCa ~ genotipo) %>%
  rstatix::add_significance("p")

cat("\n── Test de Wilcoxon ──\n"); print(stat_test)

# ── 5. Paletas (gB azul | gH rosa | gB/gH verde-rojo) ────────────────────────
cols_fill <- c(
  "gB_\u00danico"    = "#AED6F1", "gB_Mixto"    = "#4682B4",
  "gH_\u00danico"    = "pink",    "gH_Mixto"    = "hotpink3",
  "gB/gH_\u00danico" = "#4CAF7D", "gB/gH_Mixto" = "#E07070"
)
cols_jit <- c(
  "gB_\u00danico"    = "#2980B9",        "gB_Mixto"    = "#0D2B45",
  "gH_\u00danico"    = "palevioletred3", "gH_Mixto"    = "palevioletred4",
  "gB/gH_\u00danico" = "#2E7D52",        "gB/gH_Mixto" = "#B04040"
)

# ── 6. Límite del eje y posición de los corchetes ────────────────────────────
y_max <- 760
stat_test <- stat_test %>% mutate(y.position = 710)

# ── 7. Figura ────────────────────────────────────────────────────────────────
fig_auc <- ggplot(df_long, aes(x = genotipo, y = AUCa)) +
  geom_boxplot(aes(fill = cat_color),
               width = 0.45, color = "black", linewidth = 0.4,
               outlier.shape = NA, alpha = 0.85) +
  geom_jitter(aes(color = cat_color),
              width = 0.15, size = 1.8, alpha = 0.5,
              shape = 16, show.legend = FALSE) +
  facet_wrap(~ glucoproteina) +
  stat_pvalue_manual(
    stat_test,
    label        = "p.signif",
    tip.length   = 0.01,
    bracket.size = 0.6,
    size         = 5,
    color        = "black"
  ) +
  scale_fill_manual(values = cols_fill, guide = "none") +
  scale_color_manual(values = cols_jit, guide = "none") +
  scale_y_continuous(
    name   = expression("AUCa (log"[10]*"(UI/mL)\u00b7d\u00eda)"),
    breaks = seq(0, 700, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, y_max)) +
  labs(x = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position    = "none",
    panel.spacing.x    = unit(14, "pt"),
    axis.title.y       = element_text(size = 13, color = "black",
                                      margin = margin(r = 10)),
    axis.text.x        = element_text(size = 13, color = "black"),
    axis.text.y        = element_text(size = 12, color = "black"),
    axis.line          = element_line(color = "black", linewidth = 0.4),
    axis.ticks         = element_line(color = "black", linewidth = 0.3),
    strip.background   = element_blank(),
    strip.text         = element_text(size = 14, color = "black",
                                      face = "bold",
                                      margin = margin(t = 2, b = 6)),
    panel.grid         = element_blank(),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    plot.margin        = margin(15, 16, 10, 12)
  )

print(fig_auc)

cat("Pacientes con AUCa >", y_max, "(no representados):",
    sum(auc_pac$AUCa > y_max), "\n")

# ── 8. Guardar en SVG ─────────────────────────────────────────────────────────
# [guardado antiguo] ggsave(
# [guardado antiguo]   filename      = "~/Desktop/Figuras_tesis/Figura_36_AUCa.svg",
# [guardado antiguo]   plot          = fig_auc,
# [guardado antiguo]   device        = svglite::svglite,
# [guardado antiguo]   fix_text_size = FALSE,
# [guardado antiguo]   width         = 16,
# [guardado antiguo]   height        = 12,
# [guardado antiguo]   units         = "cm",
# [guardado antiguo]   bg            = "white"
# [guardado antiguo] )

## Requiere auc_pac del bloque [c] y la variable Ffinseg de pacientes

auc_anual <- auc_pac %>%
  left_join(pacientes %>% select(NHC, Ftx, Ffinseg), by = "NHC") %>%
  mutate(
    seg_anos  = as.numeric(as.Date(Ffinseg) - as.Date(Ftx)) / 365.25,
    AUCa_anio = AUCa / seg_anos
  ) %>%
  pivot_longer(c(gbm, ghm, gbghm), names_to = "gen", values_to = "mixto") %>%
  filter(!is.na(mixto)) %>%
  mutate(gen = factor(gen, levels = c("gbm", "ghm", "gbghm"),
                      labels = c("gB", "gH", "gB/gH")))

fmt_ric <- function(x, d = 0) sprintf(paste0("%.", d, "f (%.", d, "f\u2013%.", d, "f)"),
                                      median(x), quantile(x, 0.25), quantile(x, 0.75))

auc_anual %>%
  group_by(gen) %>%
  summarise(
    n_unico    = sum(mixto == 0),
    n_mixto    = sum(mixto == 1),
    AUCa_unico = fmt_ric(AUCa[mixto == 0]),
    AUCa_mixto = fmt_ric(AUCa[mixto == 1]),
    p_AUCa     = wilcox.test(AUCa ~ mixto, exact = FALSE)$p.value,
    anio_unico = fmt_ric(AUCa_anio[mixto == 0]),
    anio_mixto = fmt_ric(AUCa_anio[mixto == 1]),
    p_anio     = wilcox.test(AUCa_anio ~ mixto, exact = FALSE)$p.value,
    seg_unico  = fmt_ric(seg_anos[mixto == 0], 1),
    seg_mixto  = fmt_ric(seg_anos[mixto == 1], 1),
    p_seg      = wilcox.test(seg_anos ~ mixto, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(across(starts_with("p_"), ~ round(.x, 3))) %>%
  as.data.frame() %>%
  print()

cat("Pacientes con AUCa > 760 (no representados en la figura):",
    sum(auc_pac$AUCa > 760), "\n")

guardar_fig(fig_auc, "Figura_36", ancho = 16, alto = 12)

