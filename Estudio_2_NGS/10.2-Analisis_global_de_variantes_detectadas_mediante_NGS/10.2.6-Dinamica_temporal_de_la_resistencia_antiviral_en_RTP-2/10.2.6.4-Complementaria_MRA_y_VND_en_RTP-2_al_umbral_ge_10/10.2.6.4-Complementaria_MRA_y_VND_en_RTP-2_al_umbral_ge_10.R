# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.6.4. Complementaria. MRA y VND en RTP-2 al umbral ≥ 10%
# Extracto del script completo (ver README.md del estudio)

## ---- Complementaria. MRA y VND en RTP-2 al umbral ≥ 10% ----
variantes_rtp2_10 <- variantes_rtp2_fig %>% filter(Prevalence >= 10)
p_10 <- figura_rtp2(variantes_rtp2_10, tam_10, y_mut_base = -0.8)
grid::grid.newpage(); print(p_10)
guardar_rtp2(p_10, "Complementaria_RTP2_MRA_VND_10pct", 35, 21)


