# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.03. Anexo: Versión anterior de la Figura 16 (con n)
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Versión anterior de la Figura 16 (con n) ----
df <- muestras %>%
  mutate(
    anio   = factor(format(as.Date(Freg), "%Y")),
    gb     = rowSums(across(c(GB1, GB2, GB3, GB4)), na.rm = TRUE),
    gh     = rowSums(across(c(GH1, GH2)), na.rm = TRUE),
    gt_cat = case_when(
      gb >= 1 & gh >= 1 ~ "Genotipado completo",
      xor(gb >= 1, gh >= 1) ~ "Genotipado parcial",
      TRUE ~ "Genotipado negativo"
    ),
    gt_cat = factor(
      gt_cat,
      levels = c("Genotipado completo",
                 "Genotipado parcial",
                 "Genotipado negativo")))

totales <- df %>%
  count(anio, name = "total") %>%
  mutate(anio = factor(anio, levels = levels(df$anio)))

# ── Gráfico principal ──────────────────────────────────────────────────────────
graph_temporal <- ggplot(df, aes(x = anio, fill = gt_cat)) +
  geom_bar(width = 0.85) +
  geom_text(
    data = totales,
    aes(x = anio, y = total, label = total),
    inherit.aes = FALSE,
    vjust = -0.35, size = 4.5, fontface = "bold"
  ) +
  scale_fill_manual(values = c(
    "Genotipado completo" = "#4F94CD",
    "Genotipado parcial"  = "tan2",
    "Genotipado negativo" = "#CDC9C9"
  ), name = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(x = "", y = "Número de muestras") +
  theme_classic(base_size = 15) +
  theme(
    legend.position  = "top",
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 13),
    axis.text.y      = element_text(size = 13),
    axis.title.y     = element_text(size = 14),
    axis.line        = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
    legend.text      = element_text(size = 13),
    legend.key.size  = unit(11, "pt"),
    plot.margin      = margin(t = 10, r = 10, b = 0, l = 10))

# ── Franja de plataformas ──────────────────────────────────────────────────────
plataformas <- tibble::tribble(
  ~label,            ~xmin, ~xmax, ~color,
  "ARGENE\nSANGRE",   0.55,   4.45,  "#FF6B6B",  # 2011–2014 → pos 1–4
  "ALTONA\nSANGRE",   4.55,   6.45,  "#FF6B6B",  # 2015–2016 → pos 5–6
  "VERIS\nPLASMA",    6.55,   9.45,  "#FFD700",  # 2017–2019 → pos 7–9
  "ALTONA\nPLASMA",   9.55,  12.45,  "#FFD700"   # 2020–2022 → pos 10–12
) %>% mutate(xmid = (xmin + xmax) / 2)

franja <- ggplot(plataformas) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = color),
    alpha = 0.35, color = plataformas$color, linewidth = 0.5
  ) +
  geom_text(
    aes(x = xmid, y = 0.5, label = label),
    size = 3.8, fontface = "bold", lineheight = 0.9, color = "black"
  ) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0.5, 12.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1),      expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(t = 0, r = 10, b = 5, l = 10))

# ── Combinar ───────────────────────────────────────────────────────────────────
graph_temporal / franja + plot_layout(heights = c(10, 1.2))

