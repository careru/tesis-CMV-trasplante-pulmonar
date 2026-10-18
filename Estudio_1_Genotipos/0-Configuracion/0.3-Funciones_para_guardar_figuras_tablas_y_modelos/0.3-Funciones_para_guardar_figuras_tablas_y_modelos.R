# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 0.3. Funciones para guardar figuras, tablas y modelos
# Extracto del script completo (ver README.md del estudio)

## ---- Funciones para guardar figuras, tablas y modelos ----
guardar_fig <- function(p, nombre, ancho = 16, alto = 12) {
  ggsave(file.path(dir_salida, paste0(nombre, ".png")), p,
         width = ancho, height = alto, units = "cm", dpi = 300, bg = "white")
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(file.path(dir_salida, paste0(nombre, ".svg")), p,
           device = svglite::svglite, fix_text_size = FALSE,
           width = ancho, height = alto, units = "cm", bg = "white")
  }
  invisible(p)
}

guardar_tabla <- function(ft, nombre) {
  flextable::save_as_docx(ft, path = file.path(dir_salida, paste0(nombre, ".docx")))
  invisible(ft)
}

guardar_modelo <- function(modelo, nombre) {
  writeLines(c(capture.output(print(summary(modelo))), "",
               "IC 95%:", capture.output(print(confint(modelo)))),
             file.path(dir_salida, paste0(nombre, ".txt")))
  invisible(modelo)
}


