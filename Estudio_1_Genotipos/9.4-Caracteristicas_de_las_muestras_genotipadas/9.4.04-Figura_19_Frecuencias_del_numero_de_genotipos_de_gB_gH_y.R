# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.4.04. Figura 19. Frecuencias del número de genotipos de gB, gH y gB+gH
# Extracto del script completo (ver README.md del estudio)

## ---- Figura 19. Frecuencias del número de genotipos de gB, gH y gB+gH ----
# Gráfico de base R: el bloque abre png() en la carpeta de salidas y lo cierra con dev.off().
## 1) Filtrar muestras con genotipado completo gB + gH
df_completo <- muestras_sel %>%
  mutate(
    gb_sum = rowSums(across(c(GB1, GB2, GB3, GB4), ~ replace_na(., 0))),
    gh_sum = rowSums(across(c(GH1, GH2),           ~ replace_na(., 0)))
  ) %>%
  filter(gb_sum > 0, gh_sum > 0)

## 1b) Porcentajes de genotipos mixtos, calculados desde los datos
pct_mixtos <- df_completo %>%
  summarise(
    gb   = 100 * mean(gb_sum > 1),
    gh   = 100 * mean(gh_sum > 1),
    gbgh = 100 * mean(gb_sum > 1 | gh_sum > 1)
  )

fmt_mix <- function(x) formatC(x, format = "f", digits = 1, decimal.mark = ",")

cat(sprintf("gB mixtos: %s%% | gH mixtos: %s%% | gB/gH mixtos: %s%%\n",
            fmt_mix(pct_mixtos$gb), fmt_mix(pct_mixtos$gh), fmt_mix(pct_mixtos$gbgh)))

## 2) Resumen gB
tabla_gb <- df_completo %>%
  mutate(n_gen = rowSums(across(GB1:GB4, ~ replace_na(., 0)))) %>%
  filter(n_gen %in% 1:3) %>%
  count(n_gen, name = "n") %>%
  mutate(
    tipo = case_when(
      n_gen == 1 ~ "1 genotipo",
      n_gen == 2 ~ "2 genotipos",
      n_gen == 3 ~ "3 genotipos"
    ),
    porcentaje = round(100 * n / sum(n), 1),
    etiqueta = paste0(tipo, "\n", n, " (", gsub("\\.", ",", porcentaje), "%)")
  )
Prop_gb <- tabla_gb$n
cols_gb <- c("#63B8FF", "dodgerblue3", "#27408B")[seq_along(Prop_gb)]
n_total_gb <- sum(Prop_gb)

## 3) Resumen gH
tabla_gh <- df_completo %>%
  mutate(n_gen = rowSums(across(GH1:GH2, ~ replace_na(., 0)))) %>%
  filter(n_gen %in% 1:2) %>%
  count(n_gen, name = "n") %>%
  mutate(
    tipo = ifelse(n_gen == 1, "1 genotipo", "2 genotipos"),
    porcentaje = round(100 * n / sum(n), 1),
    etiqueta = paste0(tipo, "\n", n, " (", gsub("\\.", ",", porcentaje), "%)")
  )
Prop_gh <- tabla_gh$n
cols_gh <- c("#DDA0DD", "#CD69C9")[seq_along(Prop_gh)]
n_total_gh <- sum(Prop_gh)

## 4) Resumen combinado gB + gH
tabla_gbgh <- tabla_grupos %>%
  mutate(
    grupo_label = case_when(
      grupo == "gBunico_gHunico" ~ "Genotipos únicos gB / gH",
      grupo == "gBmixto_gHunico" ~ "Genotipos mixtos gB / únicos gH",
      grupo == "gBunico_gHmixto" ~ "Genotipos únicos gB / mixtos gH",
      grupo == "gBmixto_gHmixto" ~ "Genotipos mixtos gB / gH"
    ),
    color = case_when(
      grupo == "gBunico_gHunico" ~ "#3CB371",
      grupo == "gBmixto_gHunico" ~ "#F08080",
      grupo == "gBunico_gHmixto" ~ "#CD5C5C",
      grupo == "gBmixto_gHmixto" ~ "#8B0000"
    ),
    etiqueta = paste0(grupo_label, "\n", n, " (", gsub("\\.", ",", porcentaje), "%)")
  )

## 5) Funciones de dibujo con más aire
plot_donut <- function(props, labels, cols, center_text) {
  pie(
    props,
    labels = labels,
    col = cols,
    border = "white",
    cex = 1.4,
    radius = 0.9,
    init.angle = 45,
    main = ""
  )
  symbols(0, 0,
          circles = 0.5,
          inches = FALSE,
          add = TRUE,
          bg = "white",
          fg = "white")
  text(0, 0, center_text, cex = 1.5, font = 2)
}

plot_donut_gbgh <- function(props, labels, cols, center_text) {
  pie(
    props,
    labels = labels,
    col = cols,
    border = "white",
    cex = 1.4,
    radius = 0.9,
    init.angle = 45,
    main = ""
  )
  symbols(0, 0,
          circles = 0.5,
          inches = FALSE,
          add = TRUE,
          bg = "white",
          fg = "white")
  text(0, 0, center_text, cex = 1.5, font = 2)
}

## 5b) Etiqueta con recuadro de color bajo el donut
box_label <- function(txt, col) {
  x0 <- 0; y0 <- -1.28
  w  <- strwidth(txt,  cex = 1.6, font = 2) + 0.15
  h  <- strheight(txt, cex = 1.6, font = 2) + 0.16
  rect(x0 - w/2, y0 - h/2, x0 + w/2, y0 + h/2,
       col = col, border = NA, xpd = NA)
  text(x0, y0, txt, col = "white", font = 2, cex = 1.5, xpd = NA)
}

## 6) Layout compuesto
layout(matrix(c(1, 2,
                3, 3), nrow = 2, byrow = TRUE),
       heights = c(1, 1))
par(mar = c(5, 2, 3, 2))

# (A) Donut gB
plot_donut(Prop_gb, tabla_gb$etiqueta, cols_gb,
           paste0("gB\nn = ", n_total_gb))
title(main = "A", adj = 0, cex.main = 3)
box_label(paste0("gB mixtos: ", fmt_mix(pct_mixtos$gb), "%"), "#27408B")

# (B) Donut gH
plot_donut(Prop_gh, tabla_gh$etiqueta, cols_gh,
           paste0("gH\nn = ", n_total_gh))
title(main = "B", adj = 0, cex.main = 3)
box_label(paste0("gH mixtos: ", fmt_mix(pct_mixtos$gh), "%"), "#CD69C9")

# (C) Donut combinado gB + gH
plot_donut_gbgh(tabla_gbgh$n, tabla_gbgh$etiqueta, tabla_gbgh$color,
                paste0("gB + gH\nn = ", sum(tabla_gbgh$n)))
title(main = "C", adj = 0.02, cex.main = 3)
box_label(paste0("gB/gH mixtos: ", fmt_mix(pct_mixtos$gbgh), "%"), "#8B0000")

## 8) Guardar para media pagina A4 vertical
png(file.path(dir_salida, "Figura_19.png"),
    width = 16.5, height = 12, units = "cm", res = 300, bg = "white",
    pointsize = 10)

layout(matrix(c(1, 2,
                3, 3), nrow = 2, byrow = TRUE),
       heights = c(1, 1.05))
par(mar = c(3.2, 1, 1.8, 1))

plot_donut(Prop_gb, tabla_gb$etiqueta, cols_gb,
           paste0("gB\nn = ", n_total_gb))
title(main = "A", adj = 0, cex.main = 1.6)
box_label(paste0("gB mixtos: ", fmt_mix(pct_mixtos$gb), "%"), "#27408B")

plot_donut(Prop_gh, tabla_gh$etiqueta, cols_gh,
           paste0("gH\nn = ", n_total_gh))
title(main = "B", adj = 0, cex.main = 1.6)
box_label(paste0("gH mixtos: ", fmt_mix(pct_mixtos$gh), "%"), "#CD69C9")

plot_donut_gbgh(tabla_gbgh$n, tabla_gbgh$etiqueta, tabla_gbgh$color,
                paste0("gB + gH\nn = ", sum(tabla_gbgh$n)))
title(main = "C", adj = 0.02, cex.main = 1.6)
box_label(paste0("gB/gH mixtos: ", fmt_mix(pct_mixtos$gbgh), "%"), "#8B0000")

dev.off()

