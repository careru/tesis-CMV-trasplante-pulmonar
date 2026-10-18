# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.7.21. Figura 43. Supervivencia libre de tumores sólidos no cutáneos según genotipos mixtos
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 43. Supervivencia libre de tumores sólidos no cutáneos según genotipos mixtos ----
## Mismo estilo que la Figura 42 (funciones km_extendida / figura_km_2x2 de [7.12] [b2])
fig_km_solido <- figura_km_2x2("solido", "F_onco", "Supervivencia libre de tumor sólido", ylim = c(0.4, 1))
print(fig_km_solido)

## Exportar a SVG para una página A4 en Word (16 x 22 cm)
# [guardado antiguo] ggsave("~/Desktop/Figuras_tesis/Figura_KMext_tumor_solido_2x2.svg", fig_km_solido,
# [guardado antiguo]        device = svglite::svglite, fix_text_size = FALSE,
# [guardado antiguo]        width = 16, height = 22, units = "cm", bg = "white")

guardar_fig(fig_km_solido, "Figura_43", ancho = 16, alto = 22)

# Tabla 22 (carcinoma pulmonar y relación temporal con los episodios de CMV):
# se elaboró revisando las historias clínicas; no se genera en este script.

