# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.6.2. Figura 49. Evolución de la carga viral y detección de MRA
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 49. Evolución de la carga viral y detección de MRA ----
mra_rtp2 <- variantes_rtp2 %>% filter(clase_var == "MRA")
cat("MRA en RTP-2:\n"); mra_rtp2 %>% select(dia, Protein, Mutation, Prevalence) %>% arrange(Protein, Mutation, dia) %>% print(n = Inf)

p_mra <- figura_rtp2(mra_rtp2, tam_5, y_mut_base = -0.4,
                     y_inf_trat = 6.28, y_sup_trat = 6.58, y_top_trat = 6.84, bar_h = 0.24, y_lim_sup = 7.25)
grid::grid.newpage(); print(p_mra)
guardar_rtp2(p_mra, "Figura_49", 32, 15)

