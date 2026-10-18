# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.6.3. Figura 50. Evolución de la carga viral y detección de MRA y VND
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 50. Evolución de la carga viral y detección de MRA y VND ----
p_5 <- figura_rtp2(variantes_rtp2_fig, tam_5,
                   y_mut_base = -0.1, y_mut_step = -0.5,
                   base_size = 16, tam_etq = 5.2, tam_num = 4.5) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
        axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 12)),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text  = element_text(size = 14),
        legend.key.size = unit(0.9, "cm"))
grid::grid.newpage(); print(p_5)
guardar_rtp2(p_5, "Figura_50", 42, 24)

