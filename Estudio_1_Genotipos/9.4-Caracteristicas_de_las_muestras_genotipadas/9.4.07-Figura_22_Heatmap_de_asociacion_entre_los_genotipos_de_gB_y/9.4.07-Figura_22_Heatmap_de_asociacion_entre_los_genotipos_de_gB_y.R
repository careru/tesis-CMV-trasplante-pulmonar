# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.4.07. Figura 22. Heatmap de asociación entre los genotipos de gB y gH
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 22. Heatmap de asociación entre los genotipos de gB y gH ----
# Preparación de datos
df_asoc <- df_gbgh %>%
  filter(!is.na(gB_combo), !is.na(gH_combo))

# Tabla de contingencia (frecuencia de combinaciones gB ↔ gH)
tabla_asoc <- df_asoc %>%
  count(gB_combo, gH_combo, name = "n") %>%
  pivot_wider(names_from = gH_combo, values_from = n, values_fill = 0) %>%
  tibble::column_to_rownames("gB_combo")

tabla_asoc

# Test estadístico de independencia (χ² y Cramér’s V)
m_asoc <- as.matrix(tabla_asoc)
test_chi <- suppressWarnings(chisq.test(m_asoc))

n <- sum(m_asoc)
phi2 <- test_chi$statistic / n
k <- min(nrow(m_asoc), ncol(m_asoc))
cramers_v <- sqrt(phi2 / (k - 1))

cat("=== Asociación gB ↔ gH ===\n")
cat("Chi-cuadrado:", round(test_chi$statistic, 3), "\n")
cat("Grados de libertad:", test_chi$parameter, "\n")
cat("p-valor:", round(test_chi$p.value, 4), "\n")
cat("Cramér’s V:", round(cramers_v, 3), "\n\n")

# Mapa de calor ordenado
tabla_larga <- df_asoc %>%
  count(gB_combo, gH_combo, name = "Frecuencia")

# Definir orden manual de los ejes
orden_gB <- c("gB1", "gB2", "gB3", "gB4",
              "gB1 + gB2", "gB1 + gB3", "gB1 + gB4",
              "gB2 + gB3", "gB2 + gB4", "gB3 + gB4",
              "gB1 + gB2 + gB3", "gB1 + gB2 + gB4", "gB1 + gB3 + gB4",
              "gB2 + gB3 + gB4")
orden_gH <- c("gH1", "gH2", "gH1 + gH2")

tabla_larga <- tabla_larga %>%
  mutate(
    gB_combo = factor(gB_combo, levels = orden_gB),
    gH_combo = factor(gH_combo, levels = orden_gH)
  )

fig_22 <- ggplot(tabla_larga, aes(x = gH_combo, y = gB_combo, fill = Frecuencia)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Frecuencia), color = "black", size = 3.5) +
  scale_fill_gradient(low = "#E0F3DB", high = "#0868AC") +
  scale_y_discrete(limits = rev(levels(tabla_larga$gB_combo))) +  # ← invierte el eje Y
  labs(x = "", y = "", title = "") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 11)
  )
fig_22

# Interpretación automática
if (test_chi$p.value > 0.05) {
  cat("No se detectan asociaciones estadísticamente significativas entre los genotipos gB y gH (p > 0.05).\n")
} else {
  cat("Se detecta asociación significativa entre gB y gH (p ≤ 0.05).\n")
}

cat("Fuerza de asociación (Cramér’s V): ",
    ifelse(cramers_v < 0.1, "muy débil",
           ifelse(cramers_v < 0.3, "moderada", "fuerte")), "\n")

guardar_fig(fig_22, "Figura_22", ancho = 12, alto = 14)

