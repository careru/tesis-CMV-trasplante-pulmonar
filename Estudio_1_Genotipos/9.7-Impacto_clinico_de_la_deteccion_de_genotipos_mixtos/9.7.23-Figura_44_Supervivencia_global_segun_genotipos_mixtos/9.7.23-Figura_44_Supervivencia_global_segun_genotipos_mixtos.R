# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.7.23. Figura 44. Supervivencia global según genotipos mixtos
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 44. Supervivencia global según genotipos mixtos ----
## Mismo estilo que la Figura 42 (funciones km_extendida / figura_km_2x2 de [7.12] [b2])
fig_km_exitus <- figura_km_2x2("exitus", "Fexitus", "Supervivencia global", ylim = c(0, 1))
print(fig_km_exitus)

## Exportar a SVG para una página A4 en Word (16 x 22 cm)
# [guardado antiguo] ggsave("~/Desktop/Figuras_tesis/Figura_KMext_mortalidad_2x2.svg", fig_km_exitus,
# [guardado antiguo]        device = svglite::svglite, fix_text_size = FALSE,
# [guardado antiguo]        width = 16, height = 22, units = "cm", bg = "white")

guardar_fig(fig_km_exitus, "Figura_44", ancho = 16, alto = 22)

