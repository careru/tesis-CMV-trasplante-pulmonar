# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.14. Anexo: Donuts por episodio (no incluida)
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Donuts por episodio (no incluida) ----
##
## dnaemias contiene una fila por episodio (n total = 234).
## gbsum y ghsum cuentan cuántos genotipos se detectaron en ese episodio.


## 1) Filtrar episodios con genotipado completo gB + gH ──────────────────────
df_completo <- dnaemias %>%
  mutate(across(c(gb1:gb4, gh1:gh2), ~ replace_na(., 0))) %>%
  filter(gbsum > 0, ghsum > 0)

## Helper para etiquetas con coma decimal ────────────────────────────────────
fmt_pct <- function(x) format(x, decimal.mark = ",", nsmall = 1, trim = TRUE)

## 2) Resumen gB ─────────────────────────────────────────────────────────────
tabla_gb <- df_completo %>%
  filter(gbsum %in% 1:3) %>%
  count(gbsum, name = "n") %>%
  mutate(
    tipo = case_when(
      gbsum == 1 ~ "1 genotipo",
      gbsum == 2 ~ "2 genotipos",
      gbsum == 3 ~ "3 genotipos"
    ),
    porcentaje = round(100 * n / sum(n), 1),
    etiqueta   = paste0(tipo, "\n", n, " (", fmt_pct(porcentaje), "%)")
  )

Prop_gb    <- tabla_gb$n
cols_gb    <- c("#63B8FF", "dodgerblue3", "#27408B")[seq_along(Prop_gb)]
n_total_gb <- sum(Prop_gb)

## 3) Resumen gH ─────────────────────────────────────────────────────────────
tabla_gh <- df_completo %>%
  filter(ghsum %in% 1:2) %>%
  count(ghsum, name = "n") %>%
  mutate(
    tipo       = ifelse(ghsum == 1, "1 genotipo", "2 genotipos"),
    porcentaje = round(100 * n / sum(n), 1),
    etiqueta   = paste0(tipo, "\n", n, " (", fmt_pct(porcentaje), "%)")
  )

Prop_gh    <- tabla_gh$n
cols_gh    <- c("#DDA0DD", "#CD69C9")[seq_along(Prop_gh)]
n_total_gh <- sum(Prop_gh)

## 4) Resumen combinado gB + gH ──────────────────────────────────────────────
tabla_gbgh <- df_completo %>%
  mutate(
    grupo = case_when(
      gbsum == 1 & ghsum == 1 ~ "gBunico_gHunico",
      gbsum >= 2 & ghsum == 1 ~ "gBmixto_gHunico",
      gbsum == 1 & ghsum >= 2 ~ "gBunico_gHmixto",
      gbsum >= 2 & ghsum >= 2 ~ "gBmixto_gHmixto"
    )
  ) %>%
  count(grupo, name = "n") %>%
  mutate(
    porcentaje  = round(100 * n / sum(n), 1),
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
    etiqueta = paste0(grupo_label, "\n", n, " (", fmt_pct(porcentaje), "%)")
  )

## 5) Función de dibujo del donut con etiquetas más espaciadas ──────────────
plot_donut <- function(props, labels, cols, center_text,
                       r_outer  = 0.55,   # radio exterior del donut
                       r_inner  = 0.32,   # radio interior (hueco central)
                       r_line1  = 0.60,   # inicio de la línea guía
                       r_line2  = 0.82,   # fin de la línea guía (más larga)
                       r_label  = 0.86,   # posición de la etiqueta
                       init_deg = 45) {
  
  total <- sum(props)
  cum   <- c(0, cumsum(props) / total)
  init  <- init_deg * pi / 180
  
  op <- par(xpd = NA)
  on.exit(par(op), add = TRUE)
  
  plot.new()
  plot.window(xlim = c(-1.4, 1.4), ylim = c(-1.05, 1.05), asp = 1)
  
  ## Sectores del donut
  for (i in seq_along(props)) {
    a1   <- init + 2 * pi * cum[i]
    a2   <- init + 2 * pi * cum[i + 1]
    aseq <- seq(a1, a2, length.out = 200)
    
    polygon(
      x      = c(r_outer * cos(aseq), r_inner * cos(rev(aseq))),
      y      = c(r_outer * sin(aseq), r_inner * sin(rev(aseq))),
      col    = cols[i],
      border = "white",
      lwd    = 1.5
    )
  }
  
  ## Texto central
  text(0, 0, center_text, cex = 1.3, font = 2)
  
  ## Leader lines + etiquetas externas
  for (i in seq_along(props)) {
    mid <- init + 2 * pi * (cum[i] + cum[i + 1]) / 2
    cs  <- cos(mid)
    sn  <- sin(mid)
    
    segments(
      x0  = r_line1 * cs, y0 = r_line1 * sn,
      x1  = r_line2 * cs, y1 = r_line2 * sn,
      col = "grey45",
      lwd = 0.8
    )
    
    h_adj <- if (cs >= 0) 0 else 1
    pad   <- 0.025 * sign(cs)   # separación extra entre línea y texto
    
    text(
      x      = r_label * cs + pad,
      y      = r_label * sn,
      labels = labels[i],
      cex    = 1.2,
      adj    = c(h_adj, 0.5)
    )
  }
}

## Alias para mantener la llamada del panel C
plot_donut_gbgh <- plot_donut

## 6) Layout compuesto ───────────────────────────────────────────────────────
layout(matrix(c(1, 2,
                3, 3), nrow = 2, byrow = TRUE),
       heights = c(1, 1))
par(mar = c(2, 2, 3, 2))

# (A) Donut gB
plot_donut(Prop_gb, tabla_gb$etiqueta, cols_gb, "gB")
title(main = "A", adj = 0)

# (B) Donut gH
plot_donut(Prop_gh, tabla_gh$etiqueta, cols_gh, "gH")
title(main = "B", adj = 0)

# (C) Donut combinado gB+gH
plot_donut_gbgh(tabla_gbgh$n, tabla_gbgh$etiqueta, tabla_gbgh$color, "gB+gH")
title(main = "C", adj = 0.02)

## 7) Guardar la figura ──────────────────────────────────────────────────────
donut_figura_compuesta_episodios <- recordPlot()
donut_figura_compuesta_episodios

