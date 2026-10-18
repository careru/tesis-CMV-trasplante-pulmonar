# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.26. Anexo: KM extendida de infecciones oportunistas (no incluida)
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: KM extendida de infecciones oportunistas (no incluida) ----
## Mismo estilo que la Figura 42 (funciones km_extendida / figura_km_2x2 de [7.12] [b2])
fig_km_inf_opo <- figura_km_2x2("inf_opo", "F_opo", "Supervivencia libre de infección oportunista", ylim = c(0.4, 1))
print(fig_km_inf_opo)

## Exportar a SVG para una página A4 en Word (16 x 22 cm)
ggsave("~/Desktop/Figuras_tesis/Figura_KMext_oportunistas_2x2.svg", fig_km_inf_opo,
       device = svglite::svglite, fix_text_size = FALSE,
       width = 16, height = 22, units = "cm", bg = "white")

