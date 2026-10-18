# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.12. Anexo: Artículo: heatmap gB-gH en blanco y negro
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Artículo: heatmap gB-gH en blanco y negro ----
# 1) Datos
df_asoc <- df_gbgh %>%
  filter(!is.na(gB_combo), !is.na(gH_combo))

# 2) Orden de ejes
orden_gB <- c("gB1","gB2","gB3","gB4",
              "gB1 + gB2","gB1 + gB3","gB1 + gB4",
              "gB2 + gB3","gB2 + gB4","gB3 + gB4",
              "gB1 + gB2 + gB3","gB1 + gB2 + gB4","gB1 + gB3 + gB4",
              "gB2 + gB3 + gB4")
orden_gH <- c("gH1","gH2","gH1 + gH2")

# 3) Contingencia + Chi-cuadrado / Cramér's V
tabla_asoc <- df_asoc %>%
  count(gB_combo, gH_combo, name = "n") %>%
  pivot_wider(names_from = gH_combo, values_from = n, values_fill = 0) %>%
  tibble::column_to_rownames("gB_combo")

m_asoc    <- as.matrix(tabla_asoc)
test_chi  <- suppressWarnings(chisq.test(m_asoc))
n         <- sum(m_asoc)
k         <- min(nrow(m_asoc), ncol(m_asoc))
cramers_v <- sqrt((test_chi$statistic / n) / (k - 1))

cat(sprintf("Chi2 = %.2f | df = %d | p = %.4f | Cramer's V = %.3f\n",
            test_chi$statistic, test_chi$parameter, test_chi$p.value, cramers_v))

# 4) Formato largo
tabla_larga <- df_asoc %>%
  count(gB_combo, gH_combo, name = "Frecuencia") %>%
  mutate(
    gB_combo = factor(gB_combo, levels = orden_gB),
    gH_combo = factor(gH_combo, levels = orden_gH),
    txt_col  = ifelse(Frecuencia > 0.55 * max(Frecuencia), "white", "black")
  )

# 5) Etiqueta de asociación (p en cursiva y minúscula)
p_lab <- ifelse(test_chi$p.value < 0.001, "*p* < 0.001",
                paste0("*p* = ", formatC(test_chi$p.value, format = "f", digits = 3)))
v_lab <- sprintf("Cram\u00e9r's V = %.3f; %s", cramers_v, p_lab)

# 6) Heatmap (vertical: gH en x, gB en y)
p <- ggplot(tabla_larga, aes(x = gH_combo, y = gB_combo, fill = Frecuencia)) +
  geom_tile(color = "grey70", linewidth = 0.3) +
  geom_text(aes(label = Frecuencia, color = txt_col), size = 3.8) +
  scale_color_identity() +
  scale_fill_gradient(low = "grey95", high = "grey10", name = "Genotyped \nsamples") +
  scale_x_discrete(limits = orden_gH) +          # gH en horizontal
  scale_y_discrete(limits = rev(orden_gB)) +     # gB en vertical (gB1 arriba)
  labs(x = "gH genotype", y = "gB genotype", subtitle = v_lab) +
  coord_fixed() +
  theme_minimal(base_size = 13, base_family = "sans") +
  theme(
    plot.subtitle = element_markdown(size = 10, hjust = 0, margin = margin(b = 8)),
    axis.text.x   = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y   = element_text(color = "black"),
    axis.title.x  = element_text(size = 11, face = "bold", margin = margin(t = 12)),
    axis.title.y  = element_text(size = 11, face = "bold", margin = margin(r = 12)),
    panel.grid    = element_blank(),
    axis.ticks    = element_blank(),
    legend.key.width = unit(0.5, "cm"),
    legend.title  = element_text(size = 11),
    legend.text   = element_text(size = 11),
    plot.margin   = margin(t = 6, r = 10, b = 10, l = 10)
  )

# 7) Export en formato vertical (más alto que ancho)
ggsave("Figure_gB_gH_heatmap.tiff", p, width = 5.2, height = 6.6,
       dpi = 300, compression = "lzw")

p

