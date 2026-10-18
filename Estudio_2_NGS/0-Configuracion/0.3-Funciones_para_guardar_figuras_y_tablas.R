# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 0.3. Funciones para guardar figuras y tablas
# Extracto del script completo (ver README.md del estudio)

## ---- Funciones para guardar figuras y tablas ----
guardar_fig <- function(p, nombre, ancho = 16, alto = 12, dpi = 400) {
  ggsave(file.path(dir_salida, paste0(nombre, ".png")), p,
         width = ancho, height = alto, units = "cm", dpi = dpi, bg = "white")
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(file.path(dir_salida, paste0(nombre, ".svg")), p,
           device = svglite::svglite, fix_text_size = FALSE,
           width = ancho, height = alto, units = "cm", bg = "white")
  }
  invisible(p)
}

guardar_tabla <- function(ft, nombre, apaisada = FALSE) {
  if (apaisada) {
    doc <- read_docx() %>% body_add_flextable(ft) %>% body_end_section_landscape()
    print(doc, target = file.path(dir_salida, paste0(nombre, ".docx")))
  } else {
    flextable::save_as_docx(ft, path = file.path(dir_salida, paste0(nombre, ".docx")))
  }
  invisible(ft)
}


