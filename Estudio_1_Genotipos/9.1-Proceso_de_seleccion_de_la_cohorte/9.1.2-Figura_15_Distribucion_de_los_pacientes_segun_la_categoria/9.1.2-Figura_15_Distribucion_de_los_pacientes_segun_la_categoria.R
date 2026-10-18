# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.1.2. Figura 15. Distribución de los pacientes según la categoría de DNAemia y el seroestatus D/R
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 15. Distribución de los pacientes según la categoría de DNAemia y el seroestatus D/R ----
datos_ser <- data.frame(
  seroestatus = rep(c("D+/R+", "D-/R+", "D+/R-", "D-/R-"), each = 3),
  grupo = rep(c("Siempre CV negativa",
                "CV positiva <1000 UI/mL",
                "CV positiva ≥1000 UI/mL"), 4),
  n = c(
    # D+/R+: neg, <1000, ≥1000
    84, 46, 135,
    # D-/R+: neg, <1000, ≥1000
    27,  5,  18,
    # D+/R-: neg, <1000, ≥1000
    4,  3,  22,
    # D-/R-: neg, <1000, ≥1000
    9,  1,   1
  )
) %>%
  group_by(seroestatus) %>%
  mutate(
    total = sum(n),
    pct   = n / total * 100
  ) %>%
  ungroup() %>%
  mutate(
    seroestatus = factor(seroestatus, levels = c("D+/R+", "D-/R+", "D+/R-", "D-/R-")),
    grupo = factor(grupo, levels = c(
      "Siempre CV negativa",
      "CV positiva <1000 UI/mL",
      "CV positiva ≥1000 UI/mL"))
  )

sero_high <- datos_ser %>%
  filter(grupo == "CV positiva ≥1000 UI/mL") %>%
  mutate(no = total - n) %>%
  select(seroestatus, si = n, no, total)

calc_p <- function(g1, g2) {
  d1 <- sero_high[sero_high$seroestatus == g1, ]
  d2 <- sero_high[sero_high$seroestatus == g2, ]
  m <- matrix(c(d1$si, d1$no, d2$si, d2$no), nrow = 2)
  fisher.test(m)$p.value
}

pares <- list(
  c("D+/R+", "D+/R-"),
  c("D-/R+", "D+/R-"),
  c("D+/R-", "D-/R-")
)

pvals <- sapply(pares, function(p) calc_p(p[1], p[2]))
print(data.frame(
  comparacion = sapply(pares, function(p) paste(p[1], "vs", p[2])),
  p_valor = round(pvals, 4)
))

stars_fn <- function(p) {
  if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else NA_character_
}

# Filtrar solo las parejas con diferencias significativas
sig <- !is.na(sapply(pvals, stars_fn))
sig_pares <- pares[sig]
sig_stars <- sapply(pvals[sig], stars_fn)

# Mapeo seroestatus -> posición exacta de la barra ≥1000 UI/mL en el eje X
# Con position_dodge(width = 0.8) y 3 categorías, la 3ª barra (≥1000) queda
# desplazada +0.8/3 ≈ 0.267 respecto al centro del grupo
dodge_offset <- 0.8 / 3
x_map <- setNames((1:4) + dodge_offset, c("D+/R+", "D-/R+", "D+/R-", "D-/R-"))
xmin_v <- sapply(sig_pares, function(p) x_map[p[1]])
xmax_v <- sapply(sig_pares, function(p) x_map[p[2]])

# Apilar los corchetes a alturas crecientes (los más cortos abajo)
orden <- order(xmax_v - xmin_v, xmin_v)
xmin_v <- xmin_v[orden]
xmax_v <- xmax_v[orden]
sig_stars <- sig_stars[orden]

y_base <- 90
y_step <- 7
y_pos <- y_base + (seq_along(sig_stars) - 1) * y_step

risk_labels <- c(
  "D+/R+" = "D+/R+\n(riesgo medio D+)\nn=265",
  "D-/R+" = "D-/R+\n(riesgo medio D-)\nn=50",
  "D+/R-" = "D+/R-\n(riesgo alto)\nn=29",
  "D-/R-" = "D-/R-\n(riesgo bajo)\nn=11"
)

brackets <- data.frame(
  x    = xmin_v,
  xend = xmax_v,
  y    = y_pos,
  star = sig_stars
)

fig_15 <- ggplot(datos_ser, aes(x = seroestatus, y = pct, fill = grupo)) +
  geom_col(width = 0.7, position = position_dodge(width = 0.8)) +
  geom_text(
    aes(label = paste0(round(pct), "%")),
    position = position_dodge(width = 0.8),
    vjust = -0.4, size = 4, fontface = "bold", color = "grey20"
  ) +
  # Línea horizontal del corchete
  geom_segment(
    data = brackets,
    aes(x = x, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.4, colour = "black"
  ) +
  # Patas verticales del corchete
  geom_segment(
    data = brackets,
    aes(x = x, xend = x, y = y - 1.2, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.4, colour = "black"
  ) +
  geom_segment(
    data = brackets,
    aes(x = xend, xend = xend, y = y - 1.2, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.4, colour = "black"
  ) +
  # Asteriscos
  geom_text(
    data = brackets,
    aes(x = (x + xend) / 2, y = y + 1.5, label = star),
    inherit.aes = FALSE,
    size = 5.5, fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Siempre CV negativa"     = "#95A5A6",
      "CV positiva <1000 UI/mL" = "#5DADE2",
      "CV positiva ≥1000 UI/mL" = "#C0392B"
    ),
    name = ""
  ) +
  scale_x_discrete(labels = risk_labels) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, 115),
    breaks = seq(0, 100, 25)
  ) +
  labs(x = "", y = "% pacientes") +
  theme_classic(base_size = 15) +
  theme(
    legend.position  = "top",
    legend.text      = element_text(size = 12),
    legend.key.size  = unit(12, "pt"),
    legend.spacing.x = unit(8, "pt"),
    axis.text.x      = element_text(size = 12, lineheight = 1),
    axis.text.y      = element_text(size = 13),
    axis.title       = element_text(size = 14),
    axis.line        = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
    plot.margin      = margin(10, 20, 10, 10)
  )
fig_15

guardar_fig(fig_15, "Figura_15", ancho = 16, alto = 12)


