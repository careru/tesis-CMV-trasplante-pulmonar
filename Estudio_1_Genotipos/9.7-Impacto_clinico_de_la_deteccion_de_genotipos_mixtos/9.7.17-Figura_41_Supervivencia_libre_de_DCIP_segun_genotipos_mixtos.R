# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.7.17. Figura 41. Supervivencia libre de DCIP según genotipos mixtos
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 41. Supervivencia libre de DCIP según genotipos mixtos ----
fig_km_DCIP <- figura_km_2x2("DCIP", "F_DCIP", "Libre de DCIP", ylim = c(0.4, 1))
print(fig_km_DCIP)

## Exportar a SVG para una página A4 en Word (16 x 22 cm)
# [guardado antiguo] ggsave("~/Desktop/Figuras_tesis/Figura_42_KMext_DCIP_2x2.svg", fig_km_DCIP,
# [guardado antiguo]        device = svglite::svglite, fix_text_size = FALSE,
# [guardado antiguo]        width = 16, height = 22, units = "cm", bg = "white")

guardar_fig(fig_km_DCIP, "Figura_41", ancho = 16, alto = 22)

