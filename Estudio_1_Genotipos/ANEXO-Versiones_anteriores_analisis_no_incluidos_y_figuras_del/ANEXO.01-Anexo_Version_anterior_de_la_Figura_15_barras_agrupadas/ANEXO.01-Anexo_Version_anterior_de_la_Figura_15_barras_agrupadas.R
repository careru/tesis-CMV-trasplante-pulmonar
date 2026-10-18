# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.01. Anexo: Versión anterior de la Figura 15 (barras agrupadas)
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Versión anterior de la Figura 15 (barras agrupadas) ----
datos_ser <- data.frame(
  grupo = rep(c("Siempre CV negativa", "CV positiva <1000 UI/mL", "CV positiva ≥1000 UI/mL"), each = 4),
  seroestatus = rep(c("D+/R+", "D-/R+", "D+/R-", "D-/R-"), 3),
  n = c(
    84, 27, 4, 9,
    46,  5, 3, 1,
    135, 18, 22, 1
  )
) %>%
  group_by(grupo) %>%
  mutate(
    total = sum(n),
    pct   = n / total * 100
  ) %>%
  ungroup() %>%
  mutate(
    grupo = factor(grupo, levels = c(
      "CV positiva ≥1000 UI/mL",
      "CV positiva <1000 UI/mL",
      "Siempre CV negativa")),
    seroestatus = factor(seroestatus, levels = c("D+/R+", "D-/R+", "D+/R-", "D-/R-"))
  )

ggplot(datos_ser, aes(x = seroestatus, y = pct, fill = grupo)) +
  geom_col(width = 0.7, position = position_dodge(width = 0.8)) +
  geom_text(
    aes(label = paste0(round(pct), "%")),
    position = position_dodge(width = 0.8),
    vjust = -0.4, size = 4, fontface = "bold", color = "grey20"
  ) +
  scale_fill_manual(
    values = c(
      "CV positiva ≥1000 UI/mL" = "#C0392B",
      "CV positiva <1000 UI/mL" = "#E8A598",
      "Siempre CV negativa"     = "#5B9BD5"
    ),
    name = ""
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, ""),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(x = "", y = "% pacientes") +
  theme_classic(base_size = 15) +
  theme(
    legend.position  = "top",
    legend.text      = element_text(size = 12),
    legend.key.size  = unit(12, "pt"),
    legend.spacing.x = unit(8, "pt"),
    axis.text.x      = element_text(size = 13),
    axis.text.y      = element_text(size = 13),
    axis.title       = element_text(size = 14),
    axis.line        = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
    plot.margin      = margin(10, 20, 10, 10)
  )

