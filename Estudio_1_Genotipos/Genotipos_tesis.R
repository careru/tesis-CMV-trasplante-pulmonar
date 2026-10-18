################################################################################
#
#   TESIS DOCTORAL · CMV EN RECEPTORES DE TRASPLANTE PULMONAR
#   Estudio 1: genotipos de las glucoproteínas gB y gH
#
#   Carmen Alhena Reyes Ruiz · Hospital Universitario 12 de Octubre · UAM
#
#   Organización: una sección por apartado de Resultados (9.1 a 9.7) y, dentro
#   de cada una, un bloque por Figura o Tabla con la numeración de la tesis.
#   Todas las secciones se pliegan y despliegan (Alt+O / Alt+Shift+O en RStudio).
#   El script corre de principio a fin; las figuras y tablas se guardan en
#   'dir_salida' con su número de la tesis (Figura_27.png, Tabla_12.docx...).
#
#   Reorganizado a partir de "DEFINITIVO genotipos.R". El código de cada análisis
#   es el del script original; los cambios están marcados en los comentarios.
#
################################################################################


# ---- 0. Configuración ----

## ---- Rutas y archivos de datos ----
dir_datos  <- "~/Documents/CMV/ANALISIS"     # carpeta con los tres archivos de datos
dir_salida <- file.path(dir_datos, "salidas_tesis")

f_bd        <- "_BD CMV TXP CLINICA 22-08-24.xlsx"     # base clínica (pacientes, dnaemias...)
f_genotipos <- "Resultados_genotiposMIL_COMPLETO.xls"  # resultados de genotipado (618 muestras)
f_cargas    <- "DTA_CV_TXP.dta"                        # cargas virales seriadas

ejecutar_anexo <- FALSE   # TRUE para ejecutar también los bloques del anexo

setwd(dir_datos)
dir.create(dir_salida, showWarnings = FALSE)

## ---- Librerías ----
# MASS va primero para que dplyr::select no quede enmascarado.
library(MASS)
library(readxl)
library(haven)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(stringr)
library(lubridate)
library(scales)
library(ggplot2)
library(ggrepel)
library(ggtext)
library(ggpubr)
library(rstatix)
library(patchwork)
library(cowplot)
library(grid)
library(reshape2)
library(flextable)
library(officer)
library(survival)
library(survminer)
library(lme4)
library(lmerTest)
library(car)
library(dunn.test)
library(DiagrammeR)

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
recode <- dplyr::recode   # car::recode lo enmascara

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


# ---- 1. Datos comunes: AUC recalculada desde las cargas seriadas ----
## Trapecio sobre log10(CV) entre Fini y Ffin de cada episodio, con las fechas
## tal como están en dnaemias_sel. Negativas y "< LDC" = 0. Si no hay
## determinación en Ffin, se añade un 0 en esa fecha. Sustituye a la AUCCV de
## Stata, que sumaba el tramo entre la determinación previa y la primera positiva.
## Requiere DTA_CV_TXP.dta en el directorio de trabajo.

AUC_importada <- local({
  cargas <- read_dta(f_cargas) %>%
    zap_labels() %>%
    transmute(
      NHC   = as.numeric(NHC),
      fecha = as.Date(Freg),
      cv    = if_else(!is.na(CV) & CV > 0, as.numeric(CV), 0)
    ) %>%
    filter(!is.na(fecha)) %>%
    group_by(NHC, fecha) %>%
    summarise(cv = max(cv), .groups = "drop")
  
  episodios <- read_excel(f_bd,
                          sheet = "dnaemias_sel") %>%
    transmute(NHC  = as.numeric(NHC),
              Fini = as.Date(Fini),
              Ffin = as.Date(Ffin)) %>%
    filter(!is.na(Fini), !is.na(Ffin))
  
  auc_trapecio <- function(fecha, cv) {
    if (length(cv) < 2) return(NA_real_)
    y <- ifelse(cv > 0, pmax(log10(cv), 0), 0)
    t <- as.numeric(fecha)
    sum(diff(t) * (head(y, -1) + tail(y, -1)) / 2)
  }
  
  episodios %>%
    inner_join(cargas, by = "NHC", relationship = "many-to-many") %>%
    filter(fecha >= Fini, fecha <= Ffin) %>%
    group_by(NHC, Fini) %>%
    arrange(fecha, .by_group = TRUE) %>%
    summarise(
      AUC_importada = {
        f <- fecha; v <- cv
        if (!any(f == first(Ffin))) { f <- c(f, first(Ffin)); v <- c(v, 0) }
        auc_trapecio(f, v)
      },
      .groups = "drop"
    )
})

## Función para incorporar la AUC cada vez que se lee dnaemias_sel:
## AUCCV pasa a contener AUC_importada; la de Stata queda en AUCCV_stata.
incorporar_AUC <- function(df) {
  n_filas <- nrow(df)
  df <- df %>%
    select(-any_of(c("AUC_importada", "AUCCV_stata"))) %>%
    mutate(.NHC = as.numeric(NHC), .Fini = as.Date(Fini)) %>%
    left_join(AUC_importada %>% rename(.NHC = NHC, .Fini = Fini),
              by = c(".NHC", ".Fini"), relationship = "many-to-one") %>%
    select(-.NHC, -.Fini) %>%
    mutate(AUCCV_stata = AUCCV,
           AUCCV       = AUC_importada)
  stopifnot(nrow(df) == n_filas,
            sum(!is.na(df$AUCCV)) == sum(!is.na(df$AUCCV_stata)))
  df
}

## Comprobación (esperado: 220 episodios, mediana 73,0)
AUC_importada %>%
  summarise(episodios = n(),
            mediana   = median(AUC_importada),
            p25       = quantile(AUC_importada, 0.25),
            p75       = quantile(AUC_importada, 0.75))


# ---- 9.1. Proceso de selección de la cohorte ----

## ---- Figura 14. Proceso de selección de muestras y pacientes (datos del diagrama) ----
# La Figura 14 es un diagrama de flujo montado fuera de R; este bloque genera los recuentos.
datos_cv <- data.frame(
  categoria = c(
    "Alguna CV ≥ 1.000 UI/mL",
    "Alguna CV < 1.000 UI/mL",
    "Todas CV neg / <LDC",
    "Sin registros de CV"
  ),
  pct = c(47.0, 14.7, 33.2, 5.1)
) %>%
  mutate(categoria = factor(categoria, levels = categoria))

ggplot(datos_cv, aes(x = reorder(categoria, pct), y = pct, fill = categoria)) +
  geom_col(width = 0.6) +
  geom_text(
    aes(label = paste0(gsub("\\.", ",", format(pct, nsmall = 1)), "%")),
    hjust = -0.2, fontface = "bold", size = 4.5, color = "grey20"
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 55),
    expand = expansion(mult = c(0, 0)),
    labels = function(x) gsub("\\.", ",", x)
  ) +
  scale_fill_manual(
    values = c(
      "Alguna CV ≥ 1.000 UI/mL" = "#C0392B",
      "Alguna CV < 1.000 UI/mL" = "#E8A598",
      "Todas CV neg / <LDC"     = "#5B9BD5",
      "Sin registros de CV"     = "grey60"
    )
  ) +
  theme_classic(base_size = 15) +
  theme(
    legend.position  = "none",
    axis.text.x      = element_text(size = 13),
    axis.text.y      = element_text(size = 13),
    axis.title.x     = element_text(size = 14),
    axis.line        = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
    plot.margin      = margin(10, 20, 10, 10)
  ) +
  labs(x = "", y = "% pacientes")

## ---- Figura 15. Distribución de los pacientes según la categoría de DNAemia y el seroestatus D/R ----
datos_ser <- data.frame(
  seroestatus = rep(c("D+/R+", "D-/R+", "D+/R-", "D-/R-"), each = 3),
  grupo = rep(c("Siempre CV negativa",
                "CV positiva <1000 UI/mL",
                "CV positiva ≥1000 UI/mL"), 4),
  n = c(
    # D+/R+: neg, <1000, ≥1000
    84, 46, 135,
    # D-/R+: neg, <1000, ≥1000
    27,  5,  18,
    # D+/R-: neg, <1000, ≥1000
    4,  3,  22,
    # D-/R-: neg, <1000, ≥1000
    9,  1,   1
  )
) %>%
  group_by(seroestatus) %>%
  mutate(
    total = sum(n),
    pct   = n / total * 100
  ) %>%
  ungroup() %>%
  mutate(
    seroestatus = factor(seroestatus, levels = c("D+/R+", "D-/R+", "D+/R-", "D-/R-")),
    grupo = factor(grupo, levels = c(
      "Siempre CV negativa",
      "CV positiva <1000 UI/mL",
      "CV positiva ≥1000 UI/mL"))
  )

sero_high <- datos_ser %>%
  filter(grupo == "CV positiva ≥1000 UI/mL") %>%
  mutate(no = total - n) %>%
  select(seroestatus, si = n, no, total)

calc_p <- function(g1, g2) {
  d1 <- sero_high[sero_high$seroestatus == g1, ]
  d2 <- sero_high[sero_high$seroestatus == g2, ]
  m <- matrix(c(d1$si, d1$no, d2$si, d2$no), nrow = 2)
  fisher.test(m)$p.value
}

pares <- list(
  c("D+/R+", "D+/R-"),
  c("D-/R+", "D+/R-"),
  c("D+/R-", "D-/R-")
)

pvals <- sapply(pares, function(p) calc_p(p[1], p[2]))
print(data.frame(
  comparacion = sapply(pares, function(p) paste(p[1], "vs", p[2])),
  p_valor = round(pvals, 4)
))

stars_fn <- function(p) {
  if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else NA_character_
}

# Filtrar solo las parejas con diferencias significativas
sig <- !is.na(sapply(pvals, stars_fn))
sig_pares <- pares[sig]
sig_stars <- sapply(pvals[sig], stars_fn)

# Mapeo seroestatus -> posición exacta de la barra ≥1000 UI/mL en el eje X
# Con position_dodge(width = 0.8) y 3 categorías, la 3ª barra (≥1000) queda
# desplazada +0.8/3 ≈ 0.267 respecto al centro del grupo
dodge_offset <- 0.8 / 3
x_map <- setNames((1:4) + dodge_offset, c("D+/R+", "D-/R+", "D+/R-", "D-/R-"))
xmin_v <- sapply(sig_pares, function(p) x_map[p[1]])
xmax_v <- sapply(sig_pares, function(p) x_map[p[2]])

# Apilar los corchetes a alturas crecientes (los más cortos abajo)
orden <- order(xmax_v - xmin_v, xmin_v)
xmin_v <- xmin_v[orden]
xmax_v <- xmax_v[orden]
sig_stars <- sig_stars[orden]

y_base <- 90
y_step <- 7
y_pos <- y_base + (seq_along(sig_stars) - 1) * y_step

risk_labels <- c(
  "D+/R+" = "D+/R+\n(riesgo medio D+)\nn=265",
  "D-/R+" = "D-/R+\n(riesgo medio D-)\nn=50",
  "D+/R-" = "D+/R-\n(riesgo alto)\nn=29",
  "D-/R-" = "D-/R-\n(riesgo bajo)\nn=11"
)

brackets <- data.frame(
  x    = xmin_v,
  xend = xmax_v,
  y    = y_pos,
  star = sig_stars
)

fig_15 <- ggplot(datos_ser, aes(x = seroestatus, y = pct, fill = grupo)) +
  geom_col(width = 0.7, position = position_dodge(width = 0.8)) +
  geom_text(
    aes(label = paste0(round(pct), "%")),
    position = position_dodge(width = 0.8),
    vjust = -0.4, size = 4, fontface = "bold", color = "grey20"
  ) +
  # Línea horizontal del corchete
  geom_segment(
    data = brackets,
    aes(x = x, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.4, colour = "black"
  ) +
  # Patas verticales del corchete
  geom_segment(
    data = brackets,
    aes(x = x, xend = x, y = y - 1.2, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.4, colour = "black"
  ) +
  geom_segment(
    data = brackets,
    aes(x = xend, xend = xend, y = y - 1.2, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.4, colour = "black"
  ) +
  # Asteriscos
  geom_text(
    data = brackets,
    aes(x = (x + xend) / 2, y = y + 1.5, label = star),
    inherit.aes = FALSE,
    size = 5.5, fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Siempre CV negativa"     = "#95A5A6",
      "CV positiva <1000 UI/mL" = "#5DADE2",
      "CV positiva ≥1000 UI/mL" = "#C0392B"
    ),
    name = ""
  ) +
  scale_x_discrete(labels = risk_labels) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, 115),
    breaks = seq(0, 100, 25)
  ) +
  labs(x = "", y = "% pacientes") +
  theme_classic(base_size = 15) +
  theme(
    legend.position  = "top",
    legend.text      = element_text(size = 12),
    legend.key.size  = unit(12, "pt"),
    legend.spacing.x = unit(8, "pt"),
    axis.text.x      = element_text(size = 12, lineheight = 1),
    axis.text.y      = element_text(size = 13),
    axis.title       = element_text(size = 14),
    axis.line        = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
    plot.margin      = margin(10, 20, 10, 10)
  )
fig_15

guardar_fig(fig_15, "Figura_15", ancho = 16, alto = 12)


# ---- 9.2. Rendimiento del genotipado ----

## ---- Importación de las muestras genotipadas (n = 618) ----
muestras <- read_excel(f_genotipos,
                       col_types = c(
                         rep("numeric", 2), "date", "text", "numeric", "text",
                         rep("numeric", 14), "text", rep("skip", 3)))

length(unique(muestras$Nreg)) #692 muestras de 177 pacientes
length(unique(muestras$NHC))

muestras <- muestras %>% #quito 74 que no se encontraron / volumen insuf. (NE=.)
  filter(NE != ".")

head(muestras)
length(unique(muestras$Nreg)) #618 muestras de 169 pacientes
length(unique(muestras$NHC))

#Borrar entradas vacias  
muestras <- muestras %>% 
  filter(!is.na(Freg))

muestras %>%
  mutate(
    gb_sum = rowSums(across(c(GB1, GB2, GB3, GB4)), na.rm = TRUE),
    gh_sum = rowSums(across(c(GH1, GH2)), na.rm = TRUE),
    categoria = case_when(
      gb_sum > 0 & gh_sum > 0 ~ "Genotipado completo",
      gb_sum > 0 & gh_sum == 0 ~ "Parcial - Solo gB",
      gh_sum > 0 & gb_sum == 0 ~ "Parcial - solo gH",
      gb_sum == 0 & gh_sum == 0 ~ "Negativo gB y gH" )) %>%
  count(categoria) %>%
  mutate(porcentaje = round(100 * n / sum(n), 1))

## ---- Figura 16. Distribución anual de muestras genotipadas por técnica de carga viral ----
df <- muestras %>%
  mutate(
    anio   = factor(format(as.Date(Freg), "%Y")),
    gb     = rowSums(across(c(GB1, GB2, GB3, GB4)), na.rm = TRUE),
    gh     = rowSums(across(c(GH1, GH2)), na.rm = TRUE),
    gt_cat = case_when(
      gb >= 1 & gh >= 1 ~ "Genotipado completo",
      xor(gb >= 1, gh >= 1) ~ "Genotipado parcial",
      TRUE ~ "Genotipado negativo"
    ),
    gt_cat = factor(
      gt_cat,
      levels = c("Genotipado completo",
                 "Genotipado parcial",
                 "Genotipado negativo")))

totales <- df %>%
  count(anio, gt_cat, name = "n") %>%
  group_by(anio) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  filter(gt_cat == "Genotipado completo") %>%
  mutate(
    pct     = n / total * 100,
    pct_fmt = paste0(formatC(pct, format = "f", digits = 0, decimal.mark = ","), "%")
  ) %>%
  mutate(anio = factor(anio, levels = levels(df$anio)))

# ── Gráfico principal ──────────────────────────────────────────────────────────
graph_temporal <- ggplot(df, aes(x = anio, fill = gt_cat)) +
  geom_bar(width = 0.85) +
  geom_text(
    data = totales,
    aes(x = anio, y = total, label = pct_fmt),
    inherit.aes = FALSE,
    vjust = -0.35, size = 4.5, fontface = "bold"
  ) +
  scale_fill_manual(values = c(
    "Genotipado completo" = "#4F94CD",
    "Genotipado parcial"  = "tan2",
    "Genotipado negativo" = "#CDC9C9"
  ), name = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(x = "", y = "Número de muestras") +
  theme_classic(base_size = 15) +
  theme(
    legend.position  = "top",
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 13),
    axis.text.y      = element_text(size = 13),
    axis.title.y     = element_text(size = 14),
    axis.line        = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
    legend.text      = element_text(size = 13),
    legend.key.size  = unit(11, "pt"),
    plot.margin      = margin(t = 10, r = 10, b = 0, l = 10))

# ── Franja de plataformas ──────────────────────────────────────────────────────
plataformas <- tibble::tribble(
  ~label,            ~xmin, ~xmax, ~color,
  "ARGENE\nSANGRE",   0.55,   4.45,  "#FF6B6B",
  "ALTONA\nSANGRE",   4.55,   6.45,  "#FF6B6B",
  "VERIS\nPLASMA",    6.55,   9.45,  "#FFD700",
  "ALTONA\nPLASMA",   9.55,  12.45,  "#FFD700"
) %>% mutate(xmid = (xmin + xmax) / 2)

franja <- ggplot(plataformas) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = color),
    alpha = 0.35, color = plataformas$color, linewidth = 0.5
  ) +
  geom_text(
    aes(x = xmid, y = 0.5, label = label),
    size = 3.8, fontface = "bold", lineheight = 0.9, color = "black"
  ) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0.5, 12.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1),      expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(t = 0, r = 10, b = 5, l = 10))

# ── Combinar ───────────────────────────────────────────────────────────────────
fig_16 <- graph_temporal / franja + plot_layout(heights = c(10, 1.2))
fig_16

tabla_genotipos <- muestras %>%
  mutate(
    anio = factor(format(as.Date(Freg), "%Y")),
    gb = rowSums(across(c(GB1, GB2, GB3, GB4)), na.rm = TRUE),
    gh = rowSums(across(c(GH1, GH2)), na.rm = TRUE),
    gt_cat = case_when(
      gb >= 1 & gh >= 1 ~ "Genotipo completo",
      (gb >= 1 | gh >= 1) & !(gb >= 1 & gh >= 1) ~ "Genotipo parcial",
      TRUE ~ "Genotipo negativo"
    )
  ) %>%
  group_by(anio) %>%
  summarise(
    total      = n(),
    negativos  = sum(gt_cat == "Genotipo negativo"),
    parciales  = sum(gt_cat == "Genotipo parcial"),
    completos  = sum(gt_cat == "Genotipo completo"),
    perc_neg   = round(100 * negativos / total, 1),
    perc_par   = round(100 * parciales / total, 1),
    perc_comp  = round(100 * completos / total, 1),
    .groups    = "drop"
  )


tabla_genotipos

guardar_fig(fig_16, "Figura_16", ancho = 18, alto = 12)

## ---- Selección de las muestras con genotipado completo (n = 483) ----
muestras_sel <- muestras %>%
  mutate(
    gb = rowSums(across(c(GB1, GB2, GB3, GB4)), na.rm = TRUE),
    gh = rowSums(across(c(GH1, GH2)), na.rm = TRUE),
    gt_cat = case_when(
      gb >= 1 & gh >= 1 ~ "Genotipado completo",
      xor(gb >= 1, gh >= 1) ~ "Genotipado parcial",
      TRUE ~ "Genotipado negativo"
    )
  ) %>%
  filter(gt_cat == "Genotipado completo")

#Número de pacientes
n_distinct(muestras_sel$NHC)

# Número de muestras por paciente mínimo y maximo
muestras_por_paciente <- muestras_sel %>%
  count(NHC, name = "n_muestras")

min(muestras_por_paciente$n_muestras)
max(muestras_por_paciente$n_muestras)

# Mediana + IQR
quantile(muestras_por_paciente$n_muestras,
         probs = c(0.25, 0.5, 0.75),
         na.rm = TRUE)


# ---- 9.3. Características de la cohorte de estudio ----

## ---- Tabla 7. Características basales de la cohorte globales y por perfiles de riesgo ----
pacientes <- read_excel(f_bd,
                        sheet = "pacientes_sel")

pacientes <- pacientes %>%
  mutate(
    edad_receptor           = as.numeric(difftime(Ftx, Fnac, units = "days")) / 365.25,
    # Seguimiento hasta Ffinseg (fallecimiento, traslado o 31/12/2022), como en Métodos
    tiempo_seguimiento_anos = as.numeric(difftime(as.Date(Ffinseg), as.Date(Ftx),
                                                  units = "days")) / 365.25,
    serDR = factor(serDR, levels = c("D-/R+", "D+/R+", "D+/R-")),
    tac   = grepl("TACRO", IS_base),
    cicl  = grepl("CICLO", IS_base),
    mmf   = grepl("MICO",  IS_base),
    aza   = grepl("AZA",   IS_base)
  )

gl <- pacientes
dm <- pacientes %>% filter(serDR == "D-/R+")
pp <- pacientes %>% filter(serDR == "D+/R+")
pm <- pacientes %>% filter(serDR == "D+/R-")

# ── Helpers ────────────────────────────────────────────────────────────────────
med <- function(x) sprintf("%.1f ± %.1f", mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))
mdn <- function(x) sprintf("%.1f (%.1f–%.1f)",
                           median(x, na.rm=TRUE),
                           quantile(x, 0.25, na.rm=TRUE),
                           quantile(x, 0.75, na.rm=TRUE))
pct <- function(x, tot=length(x)) sprintf("%d (%.1f%%)",
                                          sum(x, na.rm=TRUE),
                                          100*sum(x, na.rm=TRUE)/tot)
fp  <- function(p) if(is.na(p)) "-" else if(p<0.001) "<0,001" else gsub("\\.", ",", sprintf("%.3f", p))

eb  <- function(g, cat) {
  n <- sum(g$ebase == cat, na.rm=TRUE)
  sprintf("%d (%.1f%%)", n, 100*n/nrow(g))
}

# ── p-valores ──────────────────────────────────────────────────────────────────
p_edad_rx  <- fp(summary(aov(edad_receptor           ~ serDR, data=pacientes))[[1]][["Pr(>F)"]][1])
p_edad_don <- fp(summary(aov(edad_don                ~ serDR, data=pacientes))[[1]][["Pr(>F)"]][1])
p_sex_rx   <- fp(chisq.test(table(pacientes$SexF,     pacientes$serDR))$p.value)
p_sex_don  <- fp(chisq.test(table(pacientes$SexF_don, pacientes$serDR))$p.value)
p_tipotx   <- fp(chisq.test(table(pacientes$tipotx,   pacientes$serDR))$p.value)
p_ebase    <- fp(fisher.test(table(pacientes$ebase,   pacientes$serDR),
                             simulate.p.value=TRUE, B=10000)$p.value)
p_seg      <- fp(kruskal.test(tiempo_seguimiento_anos ~ serDR, data=pacientes)$p.value)

p_tac  <- tryCatch(fp(chisq.test(table(pacientes$tac,  pacientes$serDR), simulate.p.value=TRUE)$p.value), error=function(e) "-")
p_cicl <- tryCatch(fp(fisher.test(table(pacientes$cicl, pacientes$serDR))$p.value), error=function(e) "-")
p_mmf  <- tryCatch(fp(chisq.test(table(pacientes$mmf,  pacientes$serDR), simulate.p.value=TRUE)$p.value), error=function(e) "-")
p_aza  <- tryCatch(fp(fisher.test(table(pacientes$aza,  pacientes$serDR))$p.value), error=function(e) "-")

# ── Tabla ──────────────────────────────────────────────────────────────────────
df <- data.frame(
  Caracteristica = c(
    "Edad del receptor (años), media ± DE",
    "Edad del donante (años), media ± DE",
    "Género masculino del receptor, n (%)",
    "Género masculino del donante, n (%)",
    "Tipo de trasplante, n (%)",
    "   Unilateral",
    "   Bilateral",
    "Enfermedad basal, n (%)",
    "   EPOC / enfisema",
    "   EPID",
    "   Fibrosis quística",
    "   Hipertensión pulmonar",
    "   Otros",
    "Inmunosupresión inicial, n (%)",
    "   Esteroides",
    "   Tacrolimus",
    "   Ciclosporina",
    "   Micofenolato mofetilo / ácido micofenólico",
    "   Azatioprina",
    "Tiempo de seguimiento (años), mediana (RIC)"
  ),
  Global = c(
    med(gl$edad_receptor), med(gl$edad_don),
    pct(gl$SexF==0), pct(gl$SexF_don==0),
    "", pct(gl$tipotx=="UNI"), pct(gl$tipotx=="BI"),
    "", eb(gl,"EPOC"), eb(gl,"EPID"), eb(gl,"FQ"), eb(gl,"HAP"), eb(gl,"OTROS"),
    "", pct(rep(TRUE,nrow(gl))), pct(gl$tac), pct(gl$cicl), pct(gl$mmf), pct(gl$aza),
    mdn(gl$tiempo_seguimiento_anos)
  ),
  DmRp = c(
    med(dm$edad_receptor), med(dm$edad_don),
    pct(dm$SexF==0, nrow(dm)), pct(dm$SexF_don==0, nrow(dm)),
    "", pct(dm$tipotx=="UNI", nrow(dm)), pct(dm$tipotx=="BI", nrow(dm)),
    "", eb(dm,"EPOC"), eb(dm,"EPID"), eb(dm,"FQ"), eb(dm,"HAP"), eb(dm,"OTROS"),
    "", pct(rep(TRUE,nrow(dm)),nrow(dm)), pct(dm$tac,nrow(dm)),
    pct(dm$cicl,nrow(dm)), pct(dm$mmf,nrow(dm)), pct(dm$aza,nrow(dm)),
    mdn(dm$tiempo_seguimiento_anos)
  ),
  DpRp = c(
    med(pp$edad_receptor), med(pp$edad_don),
    pct(pp$SexF==0, nrow(pp)), pct(pp$SexF_don==0, nrow(pp)),
    "", pct(pp$tipotx=="UNI", nrow(pp)), pct(pp$tipotx=="BI", nrow(pp)),
    "", eb(pp,"EPOC"), eb(pp,"EPID"), eb(pp,"FQ"), eb(pp,"HAP"), eb(pp,"OTROS"),
    "", pct(rep(TRUE,nrow(pp)),nrow(pp)), pct(pp$tac,nrow(pp)),
    pct(pp$cicl,nrow(pp)), pct(pp$mmf,nrow(pp)), pct(pp$aza,nrow(pp)),
    mdn(pp$tiempo_seguimiento_anos)
  ),
  DpRm = c(
    med(pm$edad_receptor), med(pm$edad_don),
    pct(pm$SexF==0, nrow(pm)), pct(pm$SexF_don==0, nrow(pm)),
    "", pct(pm$tipotx=="UNI", nrow(pm)), pct(pm$tipotx=="BI", nrow(pm)),
    "", eb(pm,"EPOC"), eb(pm,"EPID"), eb(pm,"FQ"), eb(pm,"HAP"), eb(pm,"OTROS"),
    "", pct(rep(TRUE,nrow(pm)),nrow(pm)), pct(pm$tac,nrow(pm)),
    pct(pm$cicl,nrow(pm)), pct(pm$mmf,nrow(pm)), pct(pm$aza,nrow(pm)),
    mdn(pm$tiempo_seguimiento_anos)
  ),
  p = c(
    p_edad_rx, p_edad_don, p_sex_rx, p_sex_don,
    p_tipotx, "", "",
    p_ebase, "", "", "", "", "",
    "", "-", p_tac, p_cicl, p_mmf, p_aza,
    p_seg
  ),
  stringsAsFactors = FALSE
)

colnames(df) <- c("Característica", "Global\n(n=146)", "D-/R+\n(n=15)",
                  "D+/R+\n(n=110)", "D+/R-\n(n=21)", "p")

# ── Flextable ──────────────────────────────────────────────────────────────────
apartados <- c("Tipo de trasplante, n (%)", "Enfermedad basal, n (%)",
               "Inmunosupresión inicial, n (%)")

ft_tabla_7 <- flextable(df) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  bold(i = ~ Característica %in% apartados, part = "body") %>%
  italic(i = ~ grepl("^   ", Característica), j = "Característica", part = "body") %>%
  align(j = "Característica", align = "left",   part = "all") %>%
  align(j = -1,               align = "center", part = "all") %>%
  fontsize(part = "all", size = 10) %>%
  padding(part = "all", padding = 3) %>%
  add_footer_lines("*p < 0,05. EPID: enfermedad pulmonar intersticial difusa; EPOC: enfermedad pulmonar obstructiva crónica. Variables categóricas comparadas mediante test exacto de Fisher con simulación de Monte Carlo (B=10.000). Cada fármaco inmunosupresor se analizó como variable dicotómica independiente.") %>%
  fontsize(part = "footer", size = 9) %>%
  autofit()

ft_tabla_7

guardar_tabla(ft_tabla_7, "Tabla_07")


# ---- 9.4. Características de las muestras genotipadas ----

## ---- Figura 17. Muestras genotipadas por año y número acumulado de receptores ----
pacientes <- read_excel(f_bd,
                        sheet = "pacientes_sel")

## 1) Muestras completamente genotipadas por año
df_comp <- muestras %>%
  mutate(
    anio = year(as.Date(Freg)),
    gb   = rowSums(across(c(GB1, GB2, GB3, GB4)), na.rm = TRUE),
    gh   = rowSums(across(c(GH1, GH2)), na.rm = TRUE),
    gt_cat = case_when(
      gb >= 1 & gh >= 1 ~ "Genotipado completo",
      xor(gb >= 1, gh >= 1) ~ "Genotipado parcial",
      TRUE ~ "Genotipado negativo"
    )
  ) %>%
  filter(gt_cat == "Genotipado completo") %>%
  count(anio, name = "n_muestras")

## 2) Trasplantes por año y acumulados a partir de `pacientes$Ftx`
tx_por_anio <- pacientes %>%
  mutate(anio = year(as_date(Ftx))) %>%
  count(anio, name = "n_tx") %>%
  arrange(anio) %>%
  mutate(tx_acum = cumsum(n_tx))

## 3) Unir tablas, asegurar años 2008–2022 y asignar metodología
anios_levels <- 2008:2022

metodologia_levels <- c(
  "ARGENE SANGRE (2008-2013)",
  "ALTONA SANGRE (2014-2016)",
  "VERIS PLASMA (2017-2019)",
  "ALTONA PLASMA (2020-2022)"
)

plot_df <- full_join(df_comp, tx_por_anio, by = "anio") %>%
  complete(anio = anios_levels) %>%
  arrange(anio) %>%
  mutate(
    anio_factor = factor(anio, levels = anios_levels),
    metodologia = case_when(
      anio %in% 2008:2013 ~ metodologia_levels[1],
      anio %in% 2014:2016 ~ metodologia_levels[2],
      anio %in% 2017:2019 ~ metodologia_levels[3],
      anio %in% 2020:2022 ~ metodologia_levels[4]
    ),
    metodologia = factor(metodologia, levels = metodologia_levels)
  )

## 4) Paleta: sangre en rojos, plasma en naranjas y amarillos
colores_metodologia <- c(
  "ARGENE SANGRE (2008-2013)" = "#8B3A3A",  # rojo granate
  "ALTONA SANGRE (2014-2016)" = "lightsalmon",  # rojo medio
  
  "VERIS PLASMA (2017-2019)"  = "darkorange3",  # naranja intenso
  "ALTONA PLASMA (2020-2022)" = "#FDD049"   # amarillo dorado
)

## 5) Ratio para el eje secundario
ratio <- max(plot_df$n_muestras, na.rm = TRUE) /
  max(plot_df$tx_acum, na.rm = TRUE)

## 6) Gráfico final
graph_acumulados <- ggplot(plot_df, aes(x = anio_factor)) +
  geom_col(aes(y = n_muestras, fill = metodologia),
           width = 0.8, na.rm = TRUE) +
  geom_text(aes(y = 0, label = n_muestras),
            vjust = 1.4, size = 4, fontface = "bold", na.rm = TRUE) +
  geom_line(aes(y = tx_acum * ratio, group = 1),
            linewidth = 0.8, na.rm = TRUE) +
  geom_point(aes(y = tx_acum * ratio),
             size = 2, na.rm = TRUE) +
  geom_text(aes(y = tx_acum * ratio, label = tx_acum),
            vjust = -1, size = 4, na.rm = TRUE) +
  scale_fill_manual(name = "Metodología:", values = colores_metodologia) +
  scale_y_continuous(
    name = "Número de muestras",
    expand = expansion(mult = c(0.05, 0.12)),
    sec.axis = sec_axis(~ . / ratio,
                        name = "Número de RTP (acumulado)")
  ) +
  coord_cartesian(clip = "off") +
  labs(x = "", y = "") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.line    = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks   = element_line(linewidth = 0.3, colour = "grey50")
  )
graph_acumulados

guardar_fig(graph_acumulados, "Figura_17", ancho = 18, alto = 12)

## ---- Figura 18. Muestras con DNAemia cuantificable por año y categoría de carga viral ----
DTA_CV_TXP <- read_dta(f_cargas)

# Años válidos
anios_todos <- DTA_CV_TXP %>%
  filter(!is.na(Freg)) %>%
  mutate(anio = as.numeric(format(as.Date(Freg), "%Y"))) %>%
  distinct(anio) %>%
  filter(!anio %in% c(2008, 2023)) %>%
  arrange(anio) %>%
  pull(anio)

# Limpieza + categorías CV
DTA_CV_TXP_cv2 <- DTA_CV_TXP %>%
  mutate(
    Resultado_interpretado = na_if(Resultado_interpretado, ""),
    anio = as.numeric(format(as.Date(Freg), "%Y"))
  ) %>%
  filter(
    !is.na(Freg),
    !is.na(Resultado_interpretado),
    !anio %in% c(2008, 2023)
  ) %>%
  mutate(
    anio = factor(anio, levels = anios_todos),
    cat_cv = case_when(
      Resultado_interpretado == "< LDC"                 ~ "CV < LDC",
      Resultado_interpretado == "POSITIVO" & CV < 1000  ~ "CV positivas < 1.000 UI/mL",
      Resultado_interpretado == "POSITIVO" & CV >= 1000 ~ "CV positivas ≥ 1.000 UI/mL",
      TRUE                                              ~ NA_character_
    )
  ) %>%
  filter(!is.na(cat_cv))

# Conteo
df_cv_cat <- DTA_CV_TXP_cv2 %>%
  group_by(anio, cat_cv) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(cat_cv = factor(cat_cv,
                         levels = c("CV < LDC",
                                    "CV positivas < 1.000 UI/mL",
                                    "CV positivas ≥ 1.000 UI/mL")
  ))


graph_cv_cat <- ggplot(df_cv_cat,
                       aes(x = anio, y = n, fill = cat_cv)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  
  scale_fill_manual(
    values = c(
      "CV < LDC"                     = "#FFC0CB",
      "CV positivas < 1.000 UI/mL"   = "#CD5555",
      "CV positivas ≥ 1.000 UI/mL"   = "#8B0000"
    ),
    name = "",
    labels = c(
      "CV < LDC",
      "CV positivas < 1.000 UI/mL",
      "CV positivas ≥ 1.000 UI/mL"
    )
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(
    x = "",
    y = "Número de muestras"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x     = element_text(angle = 45, hjust = 1),
    axis.line       = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks      = element_line(linewidth = 0.3, colour = "grey50"),
    legend.key.size = unit(9, "pt")
  )

graph_cv_cat

# Tabla con N y %
tabla_cv <- df_cv_cat %>%
  group_by(anio) %>%
  mutate(
    total_anio = sum(n),
    pct = round(100 * n / total_anio, 1)
  ) %>%
  ungroup() %>%
  arrange(anio, cat_cv)

# Mostrar en consola
tabla_cv


tabla_cv_wide <- tabla_cv %>%
  mutate(N_pct = paste0(n, " (", pct, "%)")) %>%
  select(anio, cat_cv, N_pct) %>%
  tidyr::pivot_wider(names_from = cat_cv, values_from = N_pct)

# Mostrar en consola
tabla_cv_wide

# Preparar datos para el gráfico (en %)
df_cv_cat <- DTA_CV_TXP_cv2 %>%
  group_by(anio, cat_cv) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(anio) %>%
  mutate(
    total_anio = sum(n),
    pct = round(100 * n / total_anio, 0)
  ) %>%
  ungroup()

pos_dodge <- position_dodge(width = 0.85)

graph_cv_cat <- ggplot(df_cv_cat,
                       aes(x = anio, y = pct, fill = cat_cv)) +
  geom_col(position = pos_dodge, width = 0.8) +
  
  # Etiquetas SOLO para CV ≥ 1.000 UI/mL (las rojas oscuras, a la derecha)
  geom_text(
    aes(
      label = ifelse(cat_cv == "CV ≥ 1.000 UI/mL",
                     paste0(round(pct), "%"),
                     "")
    ),
    position = pos_dodge,
    vjust = -0.5,
    size = 4
  ) +
  
  scale_fill_manual(
    values = c(
      "CV < 1.000 UI/mL" = "#FFAEB9",
      "CV ≥ 1.000 UI/mL" = "#8B0000"
    ),
    name = "",
    labels = c(
      "CV < 1.000 UI/mL",
      "CV ≥ 1.000 UI/mL"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.10)),
    limits = c(0, 100),
    labels = scales::label_number(accuracy = 1)
  ) +
  labs(
    x = "",
    y = "% muestras CV detectable"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    legend.text     = element_text(size = 12),
    axis.text       = element_text(size = 13),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 13),
    axis.title      = element_text(size = 14),
    axis.line       = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks      = element_line(linewidth = 0.3, colour = "grey50"),
    legend.key.size = unit(9, "pt")
  )

graph_cv_cat

guardar_fig(graph_cv_cat, "Figura_18", ancho = 18, alto = 12)

## ---- Resúmenes de frecuencias genotípicas por muestra (datos para las Figuras 19-21) ----
df_gb <- muestras_sel %>%
  mutate(
    gb_sum = rowSums(across(GB1:GB4, ~replace_na(., 0))),  # suma de GB1-GB4
    combo  = apply(select(., GB1:GB4) == 1, 1, function(r)
      paste(colnames(select(., GB1:GB4))[which(r)], collapse = "+"))
  ) %>%
  filter(gb_sum %in% 1:3)

tabla_resumen <- df_gb %>%
  count(gb_sum, name = "n") %>%
  mutate(
    tipo = factor(paste(gb_sum, "genotipo(s)"),
                  levels = c("1 genotipo(s)","2 genotipo(s)","3 genotipo(s)")),
    porcentaje = round(100 * n / sum(n), 1)
  ) %>%
  arrange(tipo) %>%
  select(tipo, n, porcentaje)

tabla_resumen

desglose <- df_gb %>%
  count(gb_sum, combo, name = "n") %>%
  group_by(gb_sum) %>%
  mutate(porcentaje = round(100 * n / sum(n), 1)) %>%
  arrange(gb_sum, desc(n)) %>%
  ungroup()

desglose_1 <- desglose %>% filter(gb_sum == 1) %>% select(combo, n, porcentaje)
desglose_2 <- desglose %>% filter(gb_sum == 2) %>% select(combo, n, porcentaje)
desglose_3 <- desglose %>% filter(gb_sum == 3) %>% select(combo, n, porcentaje)

desglose_1  # combinaciones cuando hay 1 genotipo
desglose_2  # combinaciones cuando hay 2 genotipos
desglose_3  # combinaciones cuando hay 3 genotipos

df_gh <- muestras_sel %>%
  mutate(
    gh_sum = rowSums(across(GH1:GH2, ~replace_na(., 0))),  # suma GH1–GH2
    combo  = apply(select(., GH1:GH2) == 1, 1, function(r)
      paste(colnames(select(., GH1:GH2))[which(r)], collapse = "+"))
  ) %>%
  filter(gh_sum %in% 1:2)

tabla_resumen_gh <- df_gh %>%
  count(gh_sum, name = "n") %>%
  mutate(
    tipo = factor(paste(gh_sum, "genotipo(s)"),
                  levels = c("1 genotipo(s)", "2 genotipo(s)")),
    porcentaje = round(100 * n / sum(n), 1)
  ) %>%
  arrange(tipo) %>%
  select(tipo, n, porcentaje)

tabla_resumen_gh

desglose_gh <- df_gh %>%
  count(gh_sum, combo, name = "n") %>%
  group_by(gh_sum) %>%
  mutate(porcentaje = round(100 * n / sum(n), 1)) %>%
  arrange(gh_sum, desc(n)) %>%
  ungroup()

desglose_1_gh <- desglose_gh %>%
  filter(gh_sum == 1) %>%
  select(combo, n, porcentaje)

desglose_2_gh <- desglose_gh %>%
  filter(gh_sum == 2) %>%
  select(combo, n, porcentaje)

desglose_1_gh   # combinaciones cuando hay 1 genotipo gH
desglose_2_gh   # combinaciones cuando hay 2 genotipos gH

# Partimos de `muestras_sel` (solo muestras con genotipado completo)

df_gbgh <- muestras_sel %>%
  mutate(
    # nº de genotipos presentes
    gb_sum = rowSums(across(GB1:GB4, ~ tidyr::replace_na(., 0))),
    gh_sum = rowSums(across(GH1:GH2, ~ tidyr::replace_na(., 0))),
    
    # combinaciones gB (formato gB1, gB1 + gB3, etc.)
    gB_combo = apply(select(., GB1:GB4) == 1, 1, function(r) {
      if (!any(r)) return(NA_character_)
      paste0("gB", which(r), collapse = " + ")
    }),
    
    # combinaciones gH (formato gH1, gH1 + gH2)
    gH_combo = apply(select(., GH1:GH2) == 1, 1, function(r) {
      if (!any(r)) return(NA_character_)
      paste0("gH", which(r), collapse = " + ")
    }),
    
    # tipo único/mixto por proteína
    gB_tipo = if_else(gb_sum == 1, "gBunico", "gBmixto"),
    gH_tipo = if_else(gh_sum == 1, "gHunico", "gHmixto"),
    
    # grupo combinado
    grupo = paste(gB_tipo, gH_tipo, sep = "_"),
    
    # combinación conjunta gB/gH
    combo_total = paste(gB_combo, gH_combo, sep = " / ")
  ) %>%
  # por seguridad, nos quedamos solo con filas con gB y gH definidos
  filter(!is.na(gB_combo), !is.na(gH_combo))

## 1) Tabla resumen de grupos (gBunico_gHunico, etc.)

niveles_grupo <- c(
  "gBunico_gHunico",
  "gBmixto_gHunico",
  "gBunico_gHmixto",
  "gBmixto_gHmixto"
)

tabla_grupos <- df_gbgh %>%
  count(grupo, name = "n") %>%
  mutate(
    grupo = factor(grupo, levels = niveles_grupo),
    porcentaje = round(100 * n / sum(n), 1)
  ) %>%
  arrange(grupo)

tabla_grupos

## Subdesglose gBunico_gHunico
desglose_gBunico_gHunico <- df_gbgh %>%
  filter(grupo == "gBunico_gHunico") %>%
  count(gB_combo, gH_combo, combo_total, name = "n") %>%
  mutate(
    porcentaje = round(100 * n / sum(n), 1)
  ) %>%
  arrange(desc(n))

desglose_gBunico_gHunico

## Subdesglose gBmixto_gHunico
desglose_gBmixto_gHunico <- df_gbgh %>%
  filter(grupo == "gBmixto_gHunico") %>%
  count(gB_combo, gH_combo, combo_total, name = "n") %>%
  mutate(
    porcentaje = round(100 * n / sum(n), 1)
  ) %>%
  arrange(desc(n))

desglose_gBmixto_gHunico

## Subdesglose gBunico_gHmixto
desglose_gBunico_gHmixto <- df_gbgh %>%
  filter(grupo == "gBunico_gHmixto") %>%
  count(gB_combo, gH_combo, combo_total, name = "n") %>%
  mutate(
    porcentaje = round(100 * n / sum(n), 1)
  ) %>%
  arrange(desc(n))

desglose_gBunico_gHmixto

## 5) Subdesglose gBmixto_gHmixto
desglose_gBmixto_gHmixto <- df_gbgh %>%
  filter(grupo == "gBmixto_gHmixto") %>%
  count(gB_combo, gH_combo, combo_total, name = "n") %>%
  mutate(
    porcentaje = round(100 * n / sum(n), 1)
  ) %>%
  arrange(desc(n))

desglose_gBmixto_gHmixto

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

## ---- Figura 20. Frecuencia de los genotipos de gB y gH a nivel de muestra ----
## Helper: porcentaje con coma decimal
fmt_pct <- function(x) paste0(formatC(x, format = "f", digits = 1, decimal.mark = ","), "%")

## 0) Bases comunes
df_gb <- muestras_sel %>%
  mutate(
    gb_sum = rowSums(across(GB1:GB4, ~ replace_na(., 0))),
    combo  = apply(select(., GB1:GB4) == 1, 1, function(r)
      paste0("gB", which(r), collapse = " + "))
  ) %>%
  filter(gb_sum %in% 1:3)

df_gh <- muestras_sel %>%
  mutate(
    gh_sum = rowSums(across(GH1:GH2, ~ replace_na(., 0))),
    combo  = apply(select(., GH1:GH2) == 1, 1, function(r)
      paste0("gH", which(r), collapse = " + "))
  ) %>%
  filter(gh_sum %in% 1:2)

## 1) Únicos gB (gB1-gB4)
tabla_uni_gb <- df_gb %>%
  filter(gb_sum == 1) %>%
  count(combo, name = "n") %>%
  arrange(n) %>%
  mutate(
    combo = factor(combo, levels = combo),
    pct   = round(100 * n / sum(n), 1),
    etiqueta = paste0(n, " (", fmt_pct(pct), ")")
  )
n_uni_gb <- sum(tabla_uni_gb$n)

cols_uni_gb <- c("gB1" = "#C6DBEF", "gB2" = "#6BAED6",
                 "gB3" = "#3182BD", "gB4" = "#08519C")

## 2) Únicos gH (gH1-gH2)
tabla_uni_gh <- df_gh %>%
  filter(gh_sum == 1) %>%
  count(combo, name = "n") %>%
  arrange(n) %>%
  mutate(
    combo = factor(combo, levels = combo),
    pct   = round(100 * n / sum(n), 1),
    etiqueta = paste0(n, " (", fmt_pct(pct), ")")
  )
n_uni_gh <- sum(tabla_uni_gh$n)

cols_uni_gh <- c("gH1" = "#DDA0DD", "gH2" = "#CD69C9")

## 3) Mixtos gB (todas las combinaciones, incl. n=1)
tabla_mix_gb <- df_gb %>%
  filter(gb_sum %in% 2:3) %>%
  count(combo, name = "n") %>%
  arrange(n) %>%
  mutate(
    combo = factor(combo, levels = combo),
    pct   = round(100 * n / sum(n), 1),
    etiqueta = paste0(n, " (", fmt_pct(pct), ")")
  )
n_mix_gb <- sum(tabla_mix_gb$n)

cols_mix_gb <- grDevices::hcl.colors(nrow(tabla_mix_gb), "Dark 3")
names(cols_mix_gb) <- as.character(tabla_mix_gb$combo)

## 4) Panel únicos gB
p_uni_gb <- ggplot(tabla_uni_gb, aes(x = n, y = combo, fill = combo)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = etiqueta), hjust = -0.15, size = 4.2) +
  scale_fill_manual(values = cols_uni_gb, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 1.30))) +
  labs(x = "", y = "",
       title = paste0("(A) Genotipos únicos gB (n = ", n_uni_gb, ")")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.margin  = margin(t = 6, r = 30, b = 6, l = 6),
    legend.position = "none",
    axis.text    = element_text(size = 11),
    axis.line    = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks   = element_line(linewidth = 0.3, colour = "grey50")
  )

## 5) Panel únicos gH
p_uni_gh <- ggplot(tabla_uni_gh, aes(x = n, y = combo, fill = combo)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = etiqueta), hjust = -0.15, size = 4.2) +
  scale_fill_manual(values = cols_uni_gh, guide = "none") +
  scale_x_continuous(
    breaks = c(0, 150, 300),
    limits = c(0, 500),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = "", y = "",
       title = paste0("(B) Genotipos únicos gH (n = ", n_uni_gh, ")")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.margin  = margin(t = 6, r = 30, b = 6, l = 6),
    legend.position = "none",
    axis.text    = element_text(size = 11),
    axis.line    = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks   = element_line(linewidth = 0.3, colour = "grey50")
  )

## 6) Panel mixtos gB
p_mix_gb <- ggplot(tabla_mix_gb, aes(x = n, y = combo, fill = combo)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = etiqueta), hjust = -0.15, size = 4.2) +
  scale_fill_manual(values = cols_mix_gb, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.30))) +
  labs(x = "Número de muestras", y = "",
       title = paste0("(C) Genotipos mixtos gB (n = ", n_mix_gb, ")")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0, size = 13,
                                margin = margin(b = 9, l = -35)),
    plot.margin  = margin(t = 6, r = 30, b = 6, l = 6),
    legend.position = "none",
    axis.text    = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    axis.line    = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks   = element_line(linewidth = 0.3, colour = "grey50")
  )

## 7) Composición: fila superior (gB | gH) y debajo mixtos gB a todo el ancho
panel_gb <- (p_uni_gb | p_uni_gh) / p_mix_gb +
  plot_layout(heights = c(max(nrow(tabla_uni_gb), nrow(tabla_uni_gh)),
                          nrow(tabla_mix_gb)))
panel_gb

guardar_fig(panel_gb, "Figura_20", ancho = 16, alto = 18)

## ---- Figura 21. Frecuencia de las combinaciones genotípicas gB/gH a nivel de muestra ----
## Helper: porcentaje con coma decimal
fmt_pct <- function(x) paste0(formatC(x, format = "f", digits = 1, decimal.mark = ","), "%")

## 0) Base gB+gH
df_gbgh <- muestras_sel %>%
  mutate(
    gb_sum = rowSums(across(GB1:GB4, ~ replace_na(., 0))),
    gh_sum = rowSums(across(GH1:GH2, ~ replace_na(., 0))),
    gB_combo = apply(select(., GB1:GB4) == 1, 1, function(r) {
      if (!any(r)) return(NA_character_)
      paste0("gB", which(r), collapse = " + ")
    }),
    gH_combo = apply(select(., GH1:GH2) == 1, 1, function(r) {
      if (!any(r)) return(NA_character_)
      paste0("gH", which(r), collapse = " + ")
    }),
    combo_total = paste(gB_combo, gH_combo, sep = " / "),
    # única solo si gB y gH son ambos únicos
    es_unica = (gb_sum == 1 & gh_sum == 1)
  ) %>%
  filter(!is.na(gB_combo), !is.na(gH_combo))

## 1) Combinaciones ÚNICAS (gB único / gH único)
tabla_uni <- df_gbgh %>%
  filter(es_unica) %>%
  count(combo_total, name = "n") %>%
  arrange(n) %>%
  mutate(
    combo_total = factor(combo_total, levels = combo_total),
    pct   = round(100 * n / sum(n), 1),
    etiqueta = paste0(n, " (", fmt_pct(pct), ")")
  )
n_uni <- sum(tabla_uni$n)

cols_uni <- colorRampPalette(c("#C7E9C0", "#00441B"))(nrow(tabla_uni))

## 2) Combinaciones MIXTAS (mezcla en gB, en gH o en ambas)
tabla_mix <- df_gbgh %>%
  filter(!es_unica) %>%
  count(combo_total, name = "n") %>%
  arrange(n) %>%
  mutate(
    combo_total = factor(combo_total, levels = combo_total),
    pct   = round(100 * n / sum(n), 1),
    etiqueta = paste0(n, " (", fmt_pct(pct), ")")
  )
n_mix <- sum(tabla_mix$n)

cols_mix <- colorRampPalette(c("#FCBBA1", "#67000D"))(nrow(tabla_mix))

## 3) Panel A: combinaciones únicas
p_uni <- ggplot(tabla_uni, aes(x = n, y = combo_total, fill = combo_total)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = etiqueta), hjust = -0.15, size = 4) +
  scale_fill_manual(values = cols_uni, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.5))) +
  labs(x = "", y = "",
       title = paste0("(A) Combinaciones únicas gB/gH (n = ", n_uni, ")")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0, size = 12,
                                margin = margin(b = 7, l = -30)),  
    plot.margin  = margin(t = 6, r = 30, b = 6, l = 6),
    legend.position = "none",
    axis.text    = element_text(size = 11),
    axis.line    = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks   = element_line(linewidth = 0.3, colour = "grey50")
  )

## 4) Panel B: combinaciones mixtas
p_mix <- ggplot(tabla_mix, aes(x = n, y = combo_total, fill = combo_total)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = etiqueta), hjust = -0.15, size = 4) +
  scale_fill_manual(values = cols_mix, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.50))) +
  labs(x = "Número de muestras", y = "",
       title = paste0("(B) Combinaciones mixtas gB/gH (n = ", n_mix, ")")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0, size = 12,
                                margin = margin(b = 7, l = -30)),
    plot.margin  = margin(t = 6, r = 30, b = 6, l = 6),
    legend.position = "none",
    axis.text    = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    axis.line    = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks   = element_line(linewidth = 0.3, colour = "grey50")
  )

## 5) Composición: únicas arriba, mixtas abajo (alturas según nº de barras)
panel_gbgh <- p_uni / p_mix +
  plot_layout(heights = c(nrow(tabla_uni), nrow(tabla_mix)))
panel_gbgh

guardar_fig(panel_gbgh, "Figura_21", ancho = 16, alto = 20)

## ---- Figura 22. Heatmap de asociación entre los genotipos de gB y gH ----
# Preparación de datos
df_asoc <- df_gbgh %>%
  filter(!is.na(gB_combo), !is.na(gH_combo))

# Tabla de contingencia (frecuencia de combinaciones gB ↔ gH)
tabla_asoc <- df_asoc %>%
  count(gB_combo, gH_combo, name = "n") %>%
  pivot_wider(names_from = gH_combo, values_from = n, values_fill = 0) %>%
  tibble::column_to_rownames("gB_combo")

tabla_asoc

# Test estadístico de independencia (χ² y Cramér’s V)
m_asoc <- as.matrix(tabla_asoc)
test_chi <- suppressWarnings(chisq.test(m_asoc))

n <- sum(m_asoc)
phi2 <- test_chi$statistic / n
k <- min(nrow(m_asoc), ncol(m_asoc))
cramers_v <- sqrt(phi2 / (k - 1))

cat("=== Asociación gB ↔ gH ===\n")
cat("Chi-cuadrado:", round(test_chi$statistic, 3), "\n")
cat("Grados de libertad:", test_chi$parameter, "\n")
cat("p-valor:", round(test_chi$p.value, 4), "\n")
cat("Cramér’s V:", round(cramers_v, 3), "\n\n")

# Mapa de calor ordenado
tabla_larga <- df_asoc %>%
  count(gB_combo, gH_combo, name = "Frecuencia")

# Definir orden manual de los ejes
orden_gB <- c("gB1", "gB2", "gB3", "gB4",
              "gB1 + gB2", "gB1 + gB3", "gB1 + gB4",
              "gB2 + gB3", "gB2 + gB4", "gB3 + gB4",
              "gB1 + gB2 + gB3", "gB1 + gB2 + gB4", "gB1 + gB3 + gB4",
              "gB2 + gB3 + gB4")
orden_gH <- c("gH1", "gH2", "gH1 + gH2")

tabla_larga <- tabla_larga %>%
  mutate(
    gB_combo = factor(gB_combo, levels = orden_gB),
    gH_combo = factor(gH_combo, levels = orden_gH)
  )

fig_22 <- ggplot(tabla_larga, aes(x = gH_combo, y = gB_combo, fill = Frecuencia)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Frecuencia), color = "black", size = 3.5) +
  scale_fill_gradient(low = "#E0F3DB", high = "#0868AC") +
  scale_y_discrete(limits = rev(levels(tabla_larga$gB_combo))) +  # ← invierte el eje Y
  labs(x = "", y = "", title = "") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 11)
  )
fig_22

# Interpretación automática
if (test_chi$p.value > 0.05) {
  cat("No se detectan asociaciones estadísticamente significativas entre los genotipos gB y gH (p > 0.05).\n")
} else {
  cat("Se detecta asociación significativa entre gB y gH (p ≤ 0.05).\n")
}

cat("Fuerza de asociación (Cramér’s V): ",
    ifelse(cramers_v < 0.1, "muy débil",
           ifelse(cramers_v < 0.3, "moderada", "fuerte")), "\n")

guardar_fig(fig_22, "Figura_22", ancho = 12, alto = 14)

## ---- Figura 23. Carga viral según genotipado completo, parcial o negativo ----
df <- muestras %>%
  mutate(
    gb = rowSums(across(c(GB1, GB2, GB3, GB4)), na.rm = TRUE),
    gh = rowSums(across(c(GH1, GH2)), na.rm = TRUE),
    gt_cat = case_when(
      gb >= 1 & gh >= 1 ~ "Genotipado completo",
      xor(gb >= 1, gh >= 1) ~ "Genotipado parcial",
      TRUE ~ "Genotipado negativo"
    ),
    gt_cat = factor(gt_cat, levels = c(
      "Genotipado completo",
      "Genotipado parcial",
      "Genotipado negativo"
    )),
    logCV = log10(CV)
  ) %>%
  drop_na(gt_cat, logCV)

comparaciones <- list(
  c("Genotipado completo", "Genotipado parcial"),
  c("Genotipado completo", "Genotipado negativo"),
  c("Genotipado parcial", "Genotipado negativo")
)

stats <- df %>%
  pairwise_wilcox_test(logCV ~ gt_cat, comparisons = comparaciones, p.adjust.method = "none") %>%
  add_significance("p") %>%
  add_y_position(fun = "max", step.increase = 0.15)

y_top <- max(stats$y.position, na.rm = TRUE) * 1.05

p <- ggplot(df, aes(x = gt_cat, y = logCV, fill = gt_cat)) +
  geom_boxplot(
    width = 0.50, color = "grey25", linewidth = 0.45,
    outlier.shape = 21, outlier.size = 1.8,
    outlier.fill = "white", outlier.color = "grey50",
    outlier.stroke = 0.4, outlier.alpha = 0.7
  ) +
  geom_jitter(
    aes(color = gt_cat),
    width = 0.12, size = 1.0, alpha = 0.35, show.legend = FALSE
  ) +
  scale_fill_manual(values = c(
    "Genotipado completo" = "#63B8FF",
    "Genotipado parcial"  = "#FFA54F",
    "Genotipado negativo" = "#BDBDBD"
  )) +
  scale_color_manual(values = c(
    "Genotipado completo" = "#1A6EBF",
    "Genotipado parcial"  = "#CC7A1A",
    "Genotipado negativo" = "#757575"
  )) +
  labs(x = NULL, y = "Carga viral CMV Log10 (UI/mL)") +
  scale_y_continuous(
    expand = expansion(mult = c(0.04, 0.12)),
    limits = c(NA, y_top)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position  = "none",
    panel.grid.major.y = element_line(color = "white", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_blank(),
    axis.line.x        = element_line(color = "grey30", linewidth = 0.5),
    axis.line.y        = element_line(color = "grey30", linewidth = 0.5),
    axis.ticks         = element_line(color = "grey30", linewidth = 0.4),
    axis.ticks.length  = unit(3, "pt"),
    axis.text.x        = element_text(size = 11, color = "grey20", margin = margin(t = 4)),
    axis.text.y        = element_text(size = 10, color = "grey30"),
    axis.title.y       = element_text(size = 11, color = "grey20", margin = margin(r = 8)),
    plot.margin        = margin(10, 15, 8, 8)
  ) +
  coord_cartesian(clip = "off") +
  stat_pvalue_manual(
    data = stats,
    label = "p.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.012,
    bracket.size = 0.4,
    size = 4.0,
    color = "grey30",
    inherit.aes = FALSE
  )

p

guardar_fig(p, "Figura_23", ancho = 14, alto = 12)

## ---- Figura 24. Carga viral según genotipos únicos y mixtos de gB, gH y gB+gH ----
## Base (solo genotipado completo gB y gH presentes)
dfc <- muestras_sel %>%
  mutate(
    gb_sum  = rowSums(across(GB1:GB4, ~ tidyr::replace_na(., 0))),
    gh_sum  = rowSums(across(GH1:GH2, ~ tidyr::replace_na(., 0))),
    CVlog10 = log10(CV)
  ) %>%
  filter(gb_sum >= 1, gh_sum >= 1) %>%
  drop_na(CV, CVlog10)

## Secciones
gb_ind <- dfc %>%
  filter(gb_sum == 1) %>%
  pivot_longer(GB1:GB4, names_to = "categoria_raw", values_to = "v") %>%
  filter(v == 1) %>%
  transmute(
    seccion   = "gB (individuales únicos)",
    categoria = paste0("g", substr(categoria_raw, 2, nchar(categoria_raw))),
    CV, CVlog10
  )

gh_ind <- dfc %>%
  filter(gh_sum == 1) %>%
  pivot_longer(GH1:GH2, names_to = "categoria_raw", values_to = "v") %>%
  filter(v == 1) %>%
  transmute(
    seccion   = "gH (individuales únicos)",
    categoria = paste0("g", substr(categoria_raw, 2, nchar(categoria_raw))),
    CV, CVlog10
  )

gb_cat <- dfc %>%
  transmute(
    seccion   = "gB (único vs mixto)",
    categoria = ifelse(gb_sum == 1, "gB único", "gB mixto"),
    CV, CVlog10
  )

gh_cat <- dfc %>%
  transmute(
    seccion   = "gH (único vs mixto)",
    categoria = ifelse(gh_sum == 1, "gH único", "gH mixto"),
    CV, CVlog10
  )

gbgh_cat <- dfc %>%
  transmute(
    seccion   = "gB+gH (único vs mixto)",
    categoria = ifelse(gb_sum == 1 & gh_sum == 1, "gB+gH único", "gB+gH mixto"),
    CV, CVlog10
  )

## Niveles
sec_levels <- c(
  "gB (individuales únicos)",
  "gH (individuales únicos)",
  "gB (único vs mixto)",
  "gH (único vs mixto)",
  "gB+gH (único vs mixto)"
)

cat_levels <- c(
  "gB1", "gB2", "gB3", "gB4",
  "gH1", "gH2",
  "gB único", "gB mixto",
  "gH único", "gH mixto",
  "gB+gH único", "gB+gH mixto"
)

df_plot <- bind_rows(gb_ind, gh_ind, gb_cat, gh_cat, gbgh_cat) %>%
  mutate(
    seccion   = factor(seccion, levels = sec_levels),
    categoria = factor(categoria, levels = cat_levels)
  )

## Paleta refinada
cols <- c(
  "gB1"         = "#5B9BD5",
  "gB2"         = "#5B9BD5",
  "gB3"         = "#5B9BD5",
  "gB4"         = "#5B9BD5",
  "gH1"         = "#D98FD6",
  "gH2"         = "#D98FD6",
  "gB único"    = "#4CAF7D",
  "gB mixto"    = "#E07070",
  "gH único"    = "#4CAF7D",
  "gH mixto"    = "#E07070",
  "gB+gH único" = "#4CAF7D",
  "gB+gH mixto" = "#E07070"
)

## Estadísticas
max_log_by_seccion <- function(data, sec_label) {
  data %>% filter(seccion == sec_label) %>% pull(CVlog10) %>% max(na.rm = TRUE)
}

st_gb_ind <- gb_ind %>%
  pairwise_wilcox_test(CV ~ categoria, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(seccion = sec_levels[1]) %>%
  arrange(group1, group2)

max_gb_log <- max_log_by_seccion(df_plot, sec_levels[1])
st_gb_ind$y.position <- max_gb_log + seq(0.12, by = 0.12, length.out = nrow(st_gb_ind))

st_gh_ind <- gh_ind %>%
  wilcox_test(CV ~ categoria) %>%
  add_significance("p") %>%
  mutate(seccion = sec_levels[2], y.position = max_log_by_seccion(df_plot, sec_levels[2]) + 0.12)

st_gb <- gb_cat %>%
  wilcox_test(CV ~ categoria) %>%
  add_significance("p") %>%
  mutate(seccion = sec_levels[3], y.position = max_log_by_seccion(df_plot, sec_levels[3]) + 0.12)

st_gh <- gh_cat %>%
  wilcox_test(CV ~ categoria) %>%
  add_significance("p") %>%
  mutate(seccion = sec_levels[4], y.position = max_log_by_seccion(df_plot, sec_levels[4]) + 0.12)

st_gbgh <- gbgh_cat %>%
  wilcox_test(CV ~ categoria) %>%
  add_significance("p") %>%
  mutate(seccion = sec_levels[5], y.position = max_log_by_seccion(df_plot, sec_levels[5]) + 0.12)

stats_all <- bind_rows(st_gb_ind, st_gh_ind, st_gb, st_gh, st_gbgh) %>%
  mutate(
    p.signif = dplyr::coalesce(p.signif, p.adj.signif),
    seccion  = factor(seccion, levels = sec_levels)
  )

stats_sig <- stats_all %>% filter(p.signif != "ns")

## Plot
p <- ggplot(df_plot, aes(categoria, CVlog10, fill = categoria)) +
  geom_boxplot(
    width          = 0.55,
    color          = "grey20",
    linewidth      = 0.42,
    outlier.shape  = 21,
    outlier.size   = 1.7,
    outlier.fill   = "white",
    outlier.color  = "grey50",
    outlier.stroke = 0.35,
    outlier.alpha  = 0.75
  ) +
  geom_jitter(
    aes(color = categoria),
    width = 0.10, size = 0.9, alpha = 0.30, show.legend = FALSE
  ) +
  facet_grid(
    ~ seccion,
    scales = "free_x",
    space  = "free_x",
    switch = "x",
    drop   = FALSE
  ) +
  scale_fill_manual(values = cols, breaks = names(cols)) +
  scale_color_manual(
    values = c(
      "gB1" = "#1A5FA8", "gB2" = "#1A5FA8", "gB3" = "#1A5FA8", "gB4" = "#1A5FA8",
      "gH1" = "#A04EA0", "gH2" = "#A04EA0",
      "gB único" = "#2E7D52", "gB mixto" = "#B04040",
      "gH único" = "#2E7D52", "gH mixto" = "#B04040",
      "gB+gH único" = "#2E7D52", "gB+gH mixto" = "#B04040"
    )
  ) +
  labs(x = NULL, y = "Carga viral CMV Log10 (UI/mL)") +
  scale_y_continuous(
    expand = expansion(mult = c(0.03, 0.12)),
    breaks = scales::pretty_breaks(n = 5)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "none",
    # Grid: solo horizontales suaves, sin verticales
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    # Cuadrante con borde fino y fondo muy levemente tintado
    panel.border       = element_rect(color = "grey55", fill = NA, linewidth = 0.45),
    panel.background   = element_rect(fill = "#FAFAFA", color = NA),
    # Strip labels debajo
    strip.placement    = "outside",
    strip.background   = element_rect(fill = "grey90", color = "grey70", linewidth = 0.35),
    strip.text         = element_text(size = 10.5, color = "grey20", face = "bold",
                                      margin = margin(t = 4, b = 4)),
    # Ejes
    axis.line          = element_blank(),   # el borde del panel ya hace de eje
    axis.ticks         = element_line(color = "grey50", linewidth = 0.35),
    axis.ticks.length  = unit(3, "pt"),
    axis.text.x        = element_text(size = 10, color = "grey20", margin = margin(t = 3)),
    axis.text.y        = element_text(size = 10, color = "grey30"),
    axis.title.y       = element_text(size = 11, color = "grey20", margin = margin(r = 8)),
    # Separación entre paneles
    panel.spacing.x    = grid::unit(10, "pt"),
    plot.margin        = margin(10, 15, 8, 8)
  ) +
  coord_cartesian(clip = "off") +
  stat_pvalue_manual(
    stats_sig,
    label        = "p.signif",
    xmin         = "group1",
    xmax         = "group2",
    y.position   = "y.position",
    tip.length   = 0.012,
    bracket.size = 0.42,
    size         = 3.8,
    color        = "grey25",
    inherit.aes  = FALSE
  )

p

## Resumen numérico
resumen_CV_gt <- df_plot %>%
  group_by(seccion, categoria) %>%
  summarise(
    n   = dplyr::n(),
    med = median(CV, na.rm = TRUE),
    p25 = quantile(CV, .25, na.rm = TRUE),
    p75 = quantile(CV, .75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(seccion, categoria) %>%
  mutate(
    med_fmt = scales::number(med, accuracy = 1, big.mark = ".", decimal.mark = ","),
    p25_fmt = scales::number(p25, accuracy = 1, big.mark = ".", decimal.mark = ","),
    p75_fmt = scales::number(p75, accuracy = 1, big.mark = ".", decimal.mark = ","),
    resumen_texto = paste0(
      "mediana de ", med_fmt,
      " UI/mL (IQR: ", p25_fmt, "–", p75_fmt,
      "; n = ", n, ")"
    )
  )

resumen_CV_gt

## Tabla de p-valores
tabla_pvalores <- stats_all %>%
  select(seccion, group1, group2, p, p.adj, p.signif) %>%
  mutate(
    p_lab       = ifelse(is.na(p), NA_character_,
                         formatC(p, format = "f", digits = 3, decimal.mark = ",")),
    p_adj_lab   = ifelse(is.na(p.adj), NA_character_,
                         formatC(p.adj, format = "f", digits = 3, decimal.mark = ",")),
    comparacion = paste0(group1, " vs ", group2)
  ) %>%
  select(seccion, comparacion, p_lab, p_adj_lab, p.signif)

tabla_pvalores

guardar_fig(p, "Figura_24", ancho = 22, alto = 12)

## ---- Figura 25. Valores de Ct según el genotipo de gB y gH ----
# Formato largo
df_ct_long <- muestras_sel %>%
  select(GB1CT, GB2CT, GB3CT, GB4CT, GH1CT, GH2CT) %>%
  pivot_longer(
    cols      = everything(),
    names_to  = "Genotipo",
    values_to = "CT"
  ) %>%
  mutate(
    Genotipo = factor(
      Genotipo,
      levels = c("GB1CT", "GB2CT", "GB3CT", "GB4CT", "GH1CT", "GH2CT")
    ),
    Etiqueta = recode(
      Genotipo,
      "GB1CT" = "gB1", "GB2CT" = "gB2", "GB3CT" = "gB3", "GB4CT" = "gB4",
      "GH1CT" = "gH1", "GH2CT" = "gH2"
    ),
    Glucoproteina = if_else(grepl("^gB", Etiqueta), "", "")
  ) %>%
  filter(!is.na(CT))

# Paletas
cols_fill <- c(
  "GB1CT" = "#5B9BD5",
  "GB2CT" = "#5B9BD5",
  "GB3CT" = "#5B9BD5",
  "GB4CT" = "#5B9BD5",
  "GH1CT" = "#D98FD6",
  "GH2CT" = "#D98FD6"
)

cols_jitter <- c(
  "GB1CT" = "#1A5FA8",
  "GB2CT" = "#1A5FA8",
  "GB3CT" = "#1A5FA8",
  "GB4CT" = "#1A5FA8",
  "GH1CT" = "#A04EA0",
  "GH2CT" = "#A04EA0"
)

y_max_plot <- max(df_ct_long$CT, na.rm = TRUE) + 1.5

# Plot
p_ct <- ggplot(df_ct_long, aes(x = Etiqueta, y = CT, fill = Genotipo)) +
  geom_boxplot(
    width          = 0.50,
    color          = "grey25",
    linewidth      = 0.45,
    outlier.shape  = 21,
    outlier.size   = 1.8,
    outlier.fill   = "white",
    outlier.color  = "grey50",
    outlier.stroke = 0.4,
    outlier.alpha  = 0.7
  ) +
  geom_jitter(
    aes(color = Genotipo),
    width = 0.12, size = 1.0, alpha = 0.35, show.legend = FALSE
  ) +
  facet_grid(
    ~ Glucoproteina,
    scales = "free_x",
    space  = "free_x"
  ) +
  scale_fill_manual(values  = cols_fill)  +
  scale_color_manual(values = cols_jitter) +
  scale_y_continuous(
    name   = "Valor Ct",
    limits = c(24, y_max_plot),
    breaks = seq(25, 40, 5),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(x = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "none",
    panel.grid.major.y = element_line(color = "white", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_blank(),
    # Línea de eje solo en el panel izquierdo para Y, y en ambos para X
    axis.line.x        = element_line(color = "grey30", linewidth = 0.5),
    axis.line.y        = element_line(color = "grey30", linewidth = 0.5),
    axis.ticks         = element_line(color = "grey30", linewidth = 0.4),
    axis.ticks.length  = unit(3, "pt"),
    axis.text.x        = element_text(size = 11, color = "grey20", margin = margin(t = 4)),
    axis.text.y        = element_text(size = 10, color = "grey30"),
    axis.title.y       = element_text(size = 11, color = "grey20", margin = margin(r = 8)),
    # Strip
    strip.background   = element_blank(),
    strip.text         = element_text(size = 12, face = "bold", color = "grey20",
                                      margin = margin(b = 6)),
    panel.spacing.x    = unit(20, "pt"),
    plot.margin        = margin(10, 15, 8, 8)
  ) +
  coord_cartesian(clip = "off")

p_ct

guardar_fig(p_ct, "Figura_25", ancho = 16, alto = 10)

## ---- Datos para el texto (9.4.4): Ct en genotipos únicos frente a mixtos ----
# ── 1. Construir dataset largo ────────────────────────────────────────────────
df_ct_comp <- bind_rows(
  muestras_sel %>%
    filter(GB1 == 1) %>%
    transmute(genotipo = "gB1", CT = GB1CT,
              tipo = if_else(gb == 1, "Único", "Mixto")),
  muestras_sel %>%
    filter(GB2 == 1) %>%
    transmute(genotipo = "gB2", CT = GB2CT,
              tipo = if_else(gb == 1, "Único", "Mixto")),
  muestras_sel %>%
    filter(GB3 == 1) %>%
    transmute(genotipo = "gB3", CT = GB3CT,
              tipo = if_else(gb == 1, "Único", "Mixto")),
  muestras_sel %>%
    filter(GB4 == 1) %>%
    transmute(genotipo = "gB4", CT = GB4CT,
              tipo = if_else(gb == 1, "Único", "Mixto")),
  muestras_sel %>%
    filter(GH1 == 1) %>%
    transmute(genotipo = "gH1", CT = GH1CT,
              tipo = if_else(gh == 1, "Único", "Mixto")),
  muestras_sel %>%
    filter(GH2 == 1) %>%
    transmute(genotipo = "gH2", CT = GH2CT,
              tipo = if_else(gh == 1, "Único", "Mixto"))
) %>%
  filter(!is.na(CT)) %>%
  mutate(
    genotipo = factor(genotipo,
                      levels = c("gB1", "gB2", "gB3", "gB4", "gH1", "gH2")),
    tipo = factor(tipo, levels = c("Único", "Mixto"))
  )

# ── 2. Resumen numérico ───────────────────────────────────────────────────────
cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  Ct por genotipo: único vs mixto (mediana [RIC], n)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

df_ct_comp %>%
  group_by(genotipo, tipo) %>%
  summarise(
    n       = n(),
    mediana = median(CT, na.rm = TRUE),
    p25     = quantile(CT, 0.25, na.rm = TRUE),
    p75     = quantile(CT, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(resumen = sprintf("%.1f (%.1f–%.1f), n=%d", mediana, p25, p75, n)) %>%
  select(genotipo, tipo, resumen) %>%
  pivot_wider(names_from = tipo, values_from = resumen) %>%
  print()

# ── 3. Tests de Wilcoxon por genotipo ─────────────────────────────────────────
cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  Wilcoxon único vs mixto por genotipo\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

stats_ct <- df_ct_comp %>%
  group_by(genotipo) %>%
  wilcox_test(CT ~ tipo) %>%
  add_significance("p")

print(stats_ct %>% select(genotipo, group1, group2, n1, n2, p, p.signif))

# ── 4. Posiciones y para las anotaciones ──────────────────────────────────────
y_pos <- df_ct_comp %>%
  group_by(genotipo) %>%
  summarise(y.position = max(CT, na.rm = TRUE) + 1.5, .groups = "drop")

stats_plot <- left_join(stats_ct, y_pos, by = "genotipo") %>%
  mutate(
    xmin = 1,
    xmax = 2
  ) %>%
  filter(p < 0.05)   # ← solo significativos

# ── 5. Figura ─────────────────────────────────────────────────────────────────
cols_tipo   <- c("Único" = "seagreen3", "Mixto" = "#E07070")
cols_jitter <- c("Único" = "seagreen4", "Mixto" = "#A03030")

figura_ct_comp <- ggplot(df_ct_comp,
                         aes(x = tipo, y = CT, fill = tipo)) +
  geom_boxplot(
    width         = 0.55,
    color         = "grey25",
    linewidth     = 0.45,
    outlier.shape = 21,
    outlier.size  = 1.8,
    outlier.fill  = "white",
    outlier.color = "grey50",
    outlier.alpha = 0.7
  ) +
  geom_jitter(
    aes(color = tipo),
    width = 0.12, size = 0.9, alpha = 0.35,
    show.legend = FALSE
  ) +
  stat_pvalue_manual(
    stats_plot,
    label        = "p.signif",
    xmin         = "xmin",
    xmax         = "xmax",
    y.position   = "y.position",
    tip.length   = 0.01,
    bracket.size = 0.45,
    size         = 4,
    inherit.aes  = FALSE
  ) +
  facet_wrap(~ genotipo, nrow = 1, scales = "free_y") +
  scale_fill_manual(values  = cols_tipo) +
  scale_color_manual(values = cols_jitter) +
  scale_y_continuous(
    name   = "Valor Ct",
    breaks = seq(25, 42, 5),
    expand = expansion(mult = c(0.02, 0.15))
  ) +
  labs(x = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "bottom",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_blank(),
    axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
    axis.line.y        = element_line(color = "grey60", linewidth = 0.4),
    axis.ticks         = element_line(color = "grey60", linewidth = 0.3),
    axis.text.x        = element_text(size = 10, color = "grey30"),
    axis.text.y        = element_text(size = 10, color = "grey30"),
    strip.background   = element_blank(),
    strip.text         = element_text(size = 12, face = "bold", color = "grey20"),
    panel.spacing.x    = unit(16, "pt"),
    plot.margin        = margin(10, 15, 8, 8)
  )

figura_ct_comp


# ---- 9.5. Episodios de DNAemia: distribución genotípica y cinética viral ----

## ---- Importación de los episodios (n = 234) con la AUC recalculada ----
# Única importación de dnaemias del script: incluye incorporar_AUC().
# (En el script antiguo, el bloque [4.1][a] repetido volvía a leerla sin la AUC recalculada.)
dnaemias <- read_excel(f_bd, 
                       sheet = "dnaemias_sel") %>%
  incorporar_AUC()

episodios_por_serostatus <- dnaemias %>%
  filter(!is.na(serDR)) %>%
  count(serDR) %>%
  rename(seroestatus = serDR, n_episodios = n)
episodios_por_serostatus

## ---- Tabla 8. Características clínicas y virológicas de los episodios por seroestatus D/R ----
## ── Preparación de datos ──────────────────────────────────────────────────

dnaemias_t <- dnaemias %>%
  mutate(
    Fini      = as.Date(Fini),
    Ffin      = as.Date(Ffin),
    FultCVpos = as.Date(FultCVpos),
    aclaramiento_dias = as.numeric(Ffin - Fini),
    no_tratado  = case_when(tto == 0 ~ 1L, tto == 1 ~ 0L, TRUE ~ NA_integer_),
    sintomatico = case_when(clinica == "VA" ~ 0L, !is.na(clinica) ~ 1L, TRUE ~ NA_integer_),
    alt_pp_bin  = if_else(!is.na(alt_pp) & alt_pp == 1, 1L, 0L)
  )

grupos <- c("Global", "D-/R+", "D+/R+", "D+/R-")

## ── Funciones auxiliares ──────────────────────────────────────────────────

fmt_iqr <- function(data, var, grp) {
  d <- if (grp == "Global") data else filter(data, serDR == grp)
  d <- filter(d, !is.na(.data[[var]]))
  if (nrow(d) == 0) return("\u2013")
  med <- median(d[[var]])
  q1  <- quantile(d[[var]], 0.25)
  q3  <- quantile(d[[var]], 0.75)
  sprintf("%.1f (%.1f\u2013%.1f)", med, q1, q3)
}

fmt_npct <- function(data, var, grp) {
  d <- if (grp == "Global") data else filter(data, serDR == grp)
  d <- filter(d, !is.na(.data[[var]]))
  if (nrow(d) == 0) return("\u2013")
  n   <- sum(d[[var]] == 1)
  tot <- nrow(d)
  if (n == 0) return("\u2013")
  sprintf("%d (%.1f)", n, 100 * n / tot)
}

fmt_p <- function(p) {
  if (is.na(p) || is.null(p)) return("")
  if (p < 0.001) return("< 0,001")
  sub("\\.", ",", sprintf("%.3f", p))
}

## Episodios por paciente (nivel paciente)
base_pac <- dnaemias_t %>% select(NHC, serDR, nepi_gtc) %>% distinct()

fmt_iqr_pac <- function(var, grp) {
  d <- if (grp == "Global") base_pac else filter(base_pac, serDR == grp)
  d <- filter(d, !is.na(.data[[var]]))
  med <- median(d[[var]])
  q1  <- quantile(d[[var]], 0.25)
  q3  <- quantile(d[[var]], 0.75)
  sprintf("%.1f (%.1f\u2013%.1f)", med, q1, q3)
}

## Episodios con >1 muestra genotipada
fmt_mas1 <- function(grp) {
  d   <- if (grp == "Global") dnaemias_t else filter(dnaemias_t, serDR == grp)
  n   <- sum(d$PCR_gt > 1, na.rm = TRUE)
  tot <- nrow(d)
  sprintf("%d (%.1f)", n, 100 * n / tot)
}

## alt_pp: NA → 0, excluir D-/R+
fmt_altpp <- function(grp) {
  if (grp == "D-/R+") return("\u2013")
  d   <- if (grp == "Global") dnaemias_t else filter(dnaemias_t, serDR == grp)
  n   <- sum(d$alt_pp_bin, na.rm = TRUE)
  tot <- nrow(d)
  if (n == 0) return("\u2013")
  sprintf("%d (%.1f)", n, 100 * n / tot)
}

## Resistencia: mostrar – si n = 0
fmt_resist <- function(grp) {
  d <- if (grp == "Global") dnaemias_t else filter(dnaemias_t, serDR == grp)
  n   <- sum(d$CMV_R == 1, na.rm = TRUE)
  tot <- nrow(d)
  if (n == 0) return("\u2013")
  sprintf("%d (%.1f)", n, 100 * n / tot)
}

## No tratado: mostrar – si n = 0
fmt_notrat <- function(grp) {
  d <- if (grp == "Global") dnaemias_t else filter(dnaemias_t, serDR == grp)
  d <- filter(d, !is.na(no_tratado))
  n   <- sum(d$no_tratado == 1)
  tot <- nrow(d)
  if (n == 0) return("\u2013")
  sprintf("%d (%.1f)", n, 100 * n / tot)
}

## ── p-values ──────────────────────────────────────────────────────────────

calc_p_cat <- function(data, var) {
  tab <- data %>%
    filter(!is.na(serDR), !is.na(.data[[var]]), serDR != "Global") %>%
    group_by(serDR) %>%
    summarise(ev = sum(.data[[var]] == 1), no = sum(.data[[var]] == 0), .groups = "drop") %>%
    column_to_rownames("serDR") %>% as.matrix()
  tryCatch(chisq.test(tab)$p.value, error = function(e) fisher.test(tab)$p.value)
}

calc_p_cont <- function(data, var) {
  d <- filter(data, !is.na(.data[[var]]), !is.na(serDR))
  tryCatch(
    kruskal.test(as.formula(paste(var, "~ serDR")), data = d)$p.value,
    error = function(e) NA_real_
  )
}

p_epi_pac  <- kruskal.test(nepi_gtc ~ serDR, data = base_pac)$p.value
p_pcr_gt   <- calc_p_cont(dnaemias_t, "PCR_gt")
p_mas1_pcr <- calc_p_cat(
  dnaemias_t %>% mutate(mas1 = as.integer(PCR_gt > 1)), "mas1"
)
p_primer   <- calc_p_cat(dnaemias_t, "primer_mil")

p_altpp <- {
  d <- filter(dnaemias_t, serDR %in% c("D+/R+", "D+/R-"))
  tab <- d %>%
    group_by(serDR) %>%
    summarise(ev = sum(alt_pp_bin), no = n() - sum(alt_pp_bin), .groups = "drop") %>%
    column_to_rownames("serDR") %>% as.matrix()
  tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_)
}

p_sint   <- calc_p_cat(dnaemias_t, "sintomatico")
p_ingr   <- calc_p_cat(dnaemias_t, "ingr")
p_ingcmv <- calc_p_cat(dnaemias_t, "ing_cmv")

p_notrat <- {
  d <- filter(dnaemias_t, !is.na(no_tratado), serDR %in% c("D-/R+", "D+/R+"))
  tab <- d %>%
    group_by(serDR) %>%
    summarise(ev = sum(no_tratado), no = n() - sum(no_tratado), .groups = "drop") %>%
    column_to_rownames("serDR") %>% as.matrix()
  tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_)
}

p_cvini  <- calc_p_cont(dnaemias_t, "cvini")
p_cvpico <- calc_p_cont(dnaemias_t, "cvpico")
p_aclar  <- calc_p_cont(filter(dnaemias_t, !is.na(aclaramiento_dias)), "aclaramiento_dias")
p_auc    <- calc_p_cont(filter(dnaemias_t, !is.na(AUCCV)), "AUCCV")

## ── Construcción de filas ─────────────────────────────────────────────────

make_row <- function(variable, vals, p_val, es_header = FALSE) {
  tibble(
    Variable = variable,
    Global   = vals[1],
    `D-/R+`  = vals[2],
    `D+/R+`  = vals[3],
    `D+/R-`  = vals[4],
    p        = if (es_header) "" else fmt_p(p_val),
    header   = es_header,
    p_raw    = if (es_header) NA_real_ else p_val
  )
}

tabla_r2 <- bind_rows(
  
  ## ── Sección 1: Distribución ──
  make_row("Distribución de episodios y muestras", rep("", 4), NA, es_header = TRUE),
  make_row(
    "Episodios por paciente, mediana (RIC)",
    sapply(grupos, fmt_iqr_pac, var = "nepi_gtc"),
    p_epi_pac
  ),
  make_row(
    "Muestras genotipadas por episodio, mediana (RIC)",
    sapply(grupos, function(g) fmt_iqr(dnaemias_t, "PCR_gt", g)),
    p_pcr_gt
  ),
  make_row(
    "Episodios con >1 muestra genotipada, n (%)",
    sapply(grupos, fmt_mas1),
    p_mas1_pcr
  ),
  
  ## ── Sección 2: Características clínicas ──
  make_row("Características clínicas", rep("", 4), NA, es_header = TRUE),
  make_row(
    "Primer episodio de CV \u2265 1.000 UI/mL, n (%)",
    sapply(grupos, function(g) fmt_npct(dnaemias_t, "primer_mil", g)),
    p_primer
  ),
  make_row(
    "Suspensión precoz de la profilaxis primaria, n (%)",
    sapply(grupos, fmt_altpp),
    p_altpp
  ),
  make_row(
    "Episodios sintomáticos, n (%)",
    sapply(grupos, function(g) fmt_npct(dnaemias_t, "sintomatico", g)),
    p_sint
  ),
  make_row(
    "Hospitalización durante el episodio, n (%)",
    sapply(grupos, function(g) fmt_npct(dnaemias_t, "ingr", g)),
    p_ingr
  ),
  make_row(
    "Hospitalización atribuible a CMV, n (%)",
    sapply(grupos, function(g) fmt_npct(dnaemias_t, "ing_cmv", g)),
    p_ingcmv
  ),
  make_row(
    "Resistencia antiviral demostrada, n (%)",
    sapply(grupos, fmt_resist),
    NA_real_
  ),
  make_row(
    "Episodio no tratado con antiviral, n (%)",
    sapply(grupos, fmt_notrat),
    p_notrat
  ),
  
  ## ── Sección 3: Cinética viral ──
  make_row("Características de cinética viral", rep("", 4), NA, es_header = TRUE),
  make_row(
    "CV inicial (UI/mL), mediana (RIC)",
    sapply(grupos, function(g) fmt_iqr(dnaemias_t, "cvini", g)),
    p_cvini
  ),
  make_row(
    "CV pico (UI/mL), mediana (RIC)",
    sapply(grupos, function(g) fmt_iqr(dnaemias_t, "cvpico", g)),
    p_cvpico
  ),
  make_row(
    "Tiempo a aclaramiento (días), mediana (RIC)",
    sapply(grupos, function(g) fmt_iqr(filter(dnaemias_t, !is.na(aclaramiento_dias)),
                                       "aclaramiento_dias", g)),
    p_aclar
  ),
  make_row(
    "AUC de exposición viral, mediana (RIC)",
    sapply(grupos, function(g) fmt_iqr(filter(dnaemias_t, !is.na(AUCCV)), "AUCCV", g)),
    p_auc
  )
)

## ── Flextable ─────────────────────────────────────────────────────────────

idx_headers <- which(tabla_r2$header)
idx_sig     <- which(!is.na(tabla_r2$p_raw) & tabla_r2$p_raw < 0.05)

ft_r2 <- tabla_r2 %>%
  select(Variable, Global, `D-/R+`, `D+/R+`, `D+/R-`, p) %>%
  flextable() %>%
  theme_booktabs() %>%
  set_header_labels(
    Variable = "Característica",
    Global   = "Global\nn = 234",
    `D-/R+`  = "D\u2212/R+\nn = 21",
    `D+/R+`  = "D+/R+\nn = 176",
    `D+/R-`  = "D+/R\u2212\nn = 37",
    p        = "p"
  ) %>%
  align(j = "Variable", align = "left",   part = "all") %>%
  align(j = c("Global","D-/R+","D+/R+","D+/R-","p"), align = "center", part = "all") %>%
  bold(i = idx_headers, j = "Variable", part = "body") %>%
  bg(i = idx_headers, bg = "#F2F2F2", part = "body") %>%
  bold(i = idx_sig, j = "p", part = "body") %>%
  padding(
    i = setdiff(seq_len(nrow(tabla_r2)), idx_headers),
    j = "Variable", padding.left = 14, part = "body"
  ) %>%
  bold(part = "header") %>%
  bg(bg = "#F2F2F2", part = "header") %>%
  add_footer_lines(values = c(
    "Los datos se presentan como mediana (RIC) o n (%), según corresponda.",
    "CV: carga viral; AUC: área bajo la curva de exposición viral acumulada.",
    "Comparaciones entre grupos: prueba de Kruskal-Wallis para variables continuas; ji-cuadrado o test exacto de Fisher para variables categóricas.",
    "– : ningún caso registrado en ese grupo.",
    "* p < 0,05 (valor destacado en negrita)."
  )) %>%
  fontsize(part = "footer", size = 8) %>%
  align(part = "footer", align = "left") %>%
  fontsize(part = "body",   size = 10) %>%
  fontsize(part = "header", size = 10) %>%
  autofit()

ft_r2

guardar_tabla(ft_r2, "Tabla_08")

## ---- Tablas 9, 10 y 11. Frecuencias genotípicas de gB, gH y gB+gH en muestras y episodios ----
## ── gB muestras ──────────────────────────────────────────────────────────

df_gb_m <- muestras_sel %>%
  mutate(
    gb_sum  = rowSums(across(GB1:GB4, ~ replace_na(., 0))),
    gB_combo = apply(select(., GB1:GB4) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gB", which(r), collapse = " + ")
    })
  ) %>%
  filter(gb_sum >= 1, !is.na(gB_combo))

n_m_total   <- nrow(df_gb_m)                          # 483
n_m_gb_uni  <- sum(df_gb_m$gb_sum == 1)               # únicos
n_m_gb_mix  <- sum(df_gb_m$gb_sum > 1)                # mixtos

freq_gb_uni_m <- df_gb_m %>%
  filter(gb_sum == 1) %>%
  count(gB_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_m_gb_uni, 1))

freq_gb_mix_m <- df_gb_m %>%
  filter(gb_sum > 1) %>%
  count(gB_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_m_gb_mix, 1)) %>%
  arrange(desc(n))

## ── gH muestras ──────────────────────────────────────────────────────────

df_gh_m <- muestras_sel %>%
  mutate(
    gh_sum   = rowSums(across(GH1:GH2, ~ replace_na(., 0))),
    gH_combo = apply(select(., GH1:GH2) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gH", which(r), collapse = " + ")
    })
  ) %>%
  filter(gh_sum >= 1, !is.na(gH_combo))

n_m_gh_uni <- sum(df_gh_m$gh_sum == 1)
n_m_gh_mix <- sum(df_gh_m$gh_sum > 1)

freq_gh_uni_m <- df_gh_m %>%
  filter(gh_sum == 1) %>%
  count(gH_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_m_gh_uni, 1))

freq_gh_mix_m <- df_gh_m %>%
  filter(gh_sum > 1) %>%
  count(gH_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_m_gh_mix, 1))

## ── gB+gH muestras ───────────────────────────────────────────────────────

df_gbgh_m <- muestras_sel %>%
  mutate(
    gb_sum   = rowSums(across(GB1:GB4, ~ replace_na(., 0))),
    gh_sum   = rowSums(across(GH1:GH2, ~ replace_na(., 0))),
    gB_combo = apply(select(., GB1:GB4) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gB", which(r), collapse = " + ")
    }),
    gH_combo = apply(select(., GH1:GH2) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gH", which(r), collapse = " + ")
    }),
    combo_total = paste(gB_combo, gH_combo, sep = " / "),
    mixto = gb_sum > 1 | gh_sum > 1
  ) %>%
  filter(!is.na(gB_combo), !is.na(gH_combo))

n_m_gbgh_uni <- sum(!df_gbgh_m$mixto)
n_m_gbgh_mix <- sum(df_gbgh_m$mixto)

freq_gbgh_uni_m <- df_gbgh_m %>%
  filter(!mixto) %>%
  count(combo_total, name = "n") %>%
  mutate(pct = round(100 * n / n_m_gbgh_uni, 1)) %>%
  arrange(desc(n))

freq_gbgh_mix_m <- df_gbgh_m %>%
  filter(mixto) %>%
  count(combo_total, name = "n") %>%
  mutate(pct = round(100 * n / n_m_gbgh_mix, 1)) %>%
  arrange(desc(n))


## ── gB episodios ─────────────────────────────────────────────────────────

df_gb_e <- dnaemias %>%
  mutate(
    gB_combo = apply(select(., gb1:gb4) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gB", which(r), collapse = " + ")
    })
  ) %>%
  filter(!is.na(gB_combo))

n_e_total  <- nrow(df_gb_e)                           # 234
n_e_gb_uni <- sum(df_gb_e$gbsum == 1)
n_e_gb_mix <- sum(df_gb_e$gbsum > 1)

freq_gb_uni_e <- df_gb_e %>%
  filter(gbsum == 1) %>%
  count(gB_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_e_gb_uni, 1))

freq_gb_mix_e <- df_gb_e %>%
  filter(gbsum > 1) %>%
  count(gB_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_e_gb_mix, 1)) %>%
  arrange(desc(n))

## ── gH episodios ─────────────────────────────────────────────────────────

df_gh_e <- dnaemias %>%
  mutate(
    gH_combo = apply(select(., gh1:gh2) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gH", which(r), collapse = " + ")
    })
  ) %>%
  filter(!is.na(gH_combo))

n_e_gh_uni <- sum(df_gh_e$ghsum == 1)
n_e_gh_mix <- sum(df_gh_e$ghsum > 1)

freq_gh_uni_e <- df_gh_e %>%
  filter(ghsum == 1) %>%
  count(gH_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_e_gh_uni, 1))

freq_gh_mix_e <- df_gh_e %>%
  filter(ghsum > 1) %>%
  count(gH_combo, name = "n") %>%
  mutate(pct = round(100 * n / n_e_gh_mix, 1))

## ── gB+gH episodios ──────────────────────────────────────────────────────

df_gbgh_e <- dnaemias %>%
  mutate(
    gB_combo = apply(select(., gb1:gb4) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gB", which(r), collapse = " + ")
    }),
    gH_combo = apply(select(., gh1:gh2) == 1, 1, function(r) {
      if (!any(r, na.rm = TRUE)) return(NA_character_)
      paste0("gH", which(r), collapse = " + ")
    }),
    combo_total = paste(gB_combo, gH_combo, sep = " / "),
    mixto = gbsum > 1 | ghsum > 1
  ) %>%
  filter(!is.na(gB_combo), !is.na(gH_combo))

n_e_gbgh_uni <- sum(!df_gbgh_e$mixto)
n_e_gbgh_mix <- sum(df_gbgh_e$mixto)

freq_gbgh_uni_e <- df_gbgh_e %>%
  filter(!mixto) %>%
  count(combo_total, name = "n") %>%
  mutate(pct = round(100 * n / n_e_gbgh_uni, 1)) %>%
  arrange(desc(n))

freq_gbgh_mix_e <- df_gbgh_e %>%
  filter(mixto) %>%
  count(combo_total, name = "n") %>%
  mutate(pct = round(100 * n / n_e_gbgh_mix, 1)) %>%
  arrange(desc(n))


fmt_np <- function(n, pct) sprintf("%d (%s)", n, format(pct, nsmall = 1))

hacer_ft <- function(df, titulo, nota) {
  
  idx_headers <- which(df$header)
  idx_indent  <- which(!df$header)
  
  ft <- df %>%
    select(Variable, Muestras, Episodios) %>%
    flextable() %>%
    theme_booktabs() %>%
    set_header_labels(
      Variable  = titulo,
      Muestras  = paste0("Muestras\n(n = ", n_m_total, ")"),
      Episodios = paste0("Episodios\n(n = ", n_e_total, ")")
    ) %>%
    align(j = "Variable",              align = "left",   part = "all") %>%
    align(j = c("Muestras","Episodios"), align = "center", part = "all") %>%
    bold(i = idx_headers, j = "Variable", part = "body") %>%
    bg(i = idx_headers, bg = "#F2F2F2", part = "body") %>%
    padding(i = idx_indent, j = "Variable", padding.left = 14, part = "body") %>%
    bold(part = "header") %>%
    bg(bg = "#F2F2F2", part = "header") %>%
    add_footer_lines(values = nota) %>%
    fontsize(part = "footer", size = 8) %>%
    align(part = "footer", align = "left") %>%
    fontsize(part = "body",   size = 10) %>%
    fontsize(part = "header", size = 10) %>%
    autofit()
  
  ft
}


combos_gb <- sort(union(freq_gb_uni_m$gB_combo, freq_gb_uni_e$gB_combo))
combos_gb_mix <- union(freq_gb_mix_m$gB_combo, freq_gb_mix_e$gB_combo)

get_m <- function(df, key) {
  r <- df[df[[1]] == key, ]
  if (nrow(r) == 0) return("\u2013")
  fmt_np(r$n, r$pct)
}

filas_r3 <- bind_rows(
  tibble(Variable = "Genotipos únicos gB, n (%)",
         Muestras  = fmt_np(n_m_gb_uni, round(100*n_m_gb_uni/n_m_total,1)),
         Episodios = fmt_np(n_e_gb_uni, round(100*n_e_gb_uni/n_e_total,1)),
         header = TRUE),
  
  purrr::map_dfr(combos_gb, function(k) {
    tibble(Variable  = k,
           Muestras  = get_m(freq_gb_uni_m, k),
           Episodios = get_m(freq_gb_uni_e, k),
           header    = FALSE)
  }),
  
  tibble(Variable  = "Genotipos mixtos gB, n (%)",
         Muestras  = fmt_np(n_m_gb_mix, round(100*n_m_gb_mix/n_m_total,1)),
         Episodios = fmt_np(n_e_gb_mix, round(100*n_e_gb_mix/n_e_total,1)),
         header    = TRUE),
  
  purrr::map_dfr(combos_gb_mix, function(k) {
    tibble(Variable  = k,
           Muestras  = get_m(rename(freq_gb_mix_m, gB_combo=gB_combo), k),
           Episodios = get_m(rename(freq_gb_mix_e, gB_combo=gB_combo), k),
           header    = FALSE)
  })
)

ft_r3 <- hacer_ft(
  filas_r3,
  "Genotipo gB",
  c("Los porcentajes de gB1\u2013gB4 se calculan respecto al total de genotipos únicos.",
    "Los porcentajes de genotipos mixtos se calculan respecto al total de genotipos mixtos.")
)

ft_r3


filas_r4 <- bind_rows(
  tibble(Variable  = "Genotipos únicos gH, n (%)",
         Muestras  = fmt_np(n_m_gh_uni, round(100*n_m_gh_uni/n_m_total,1)),
         Episodios = fmt_np(n_e_gh_uni, round(100*n_e_gh_uni/n_e_total,1)),
         header    = TRUE),
  
  purrr::map_dfr(sort(union(freq_gh_uni_m$gH_combo, freq_gh_uni_e$gH_combo)), function(k) {
    tibble(Variable  = k,
           Muestras  = get_m(freq_gh_uni_m, k),
           Episodios = get_m(freq_gh_uni_e, k),
           header    = FALSE)
  }),
  
  tibble(Variable  = "Genotipos mixtos gH, n (%)",
         Muestras  = fmt_np(n_m_gh_mix, round(100*n_m_gh_mix/n_m_total,1)),
         Episodios = fmt_np(n_e_gh_mix, round(100*n_e_gh_mix/n_e_total,1)),
         header    = TRUE),
  
  purrr::map_dfr(union(freq_gh_mix_m$gH_combo, freq_gh_mix_e$gH_combo), function(k) {
    tibble(Variable  = k,
           Muestras  = get_m(freq_gh_mix_m, k),
           Episodios = get_m(freq_gh_mix_e, k),
           header    = FALSE)
  })
)

ft_r4 <- hacer_ft(
  filas_r4,
  "Genotipo gH",
  c("Los porcentajes de gH1 y gH2 se calculan respecto al total de genotipos únicos.",
    "Los porcentajes de genotipos mixtos se calculan respecto al total de genotipos mixtos.")
)

ft_r4


combos_gbgh_uni  <- union(freq_gbgh_uni_m$combo_total,  freq_gbgh_uni_e$combo_total)
combos_gbgh_mix  <- union(freq_gbgh_mix_m$combo_total,  freq_gbgh_mix_e$combo_total)

## Ordenar por frecuencia en muestras (desc)
combos_gbgh_uni <- freq_gbgh_uni_m %>%
  full_join(freq_gbgh_uni_e, by = "combo_total", suffix = c("_m","_e")) %>%
  arrange(desc(replace_na(n_m, 0))) %>%
  pull(combo_total)

combos_gbgh_mix <- freq_gbgh_mix_m %>%
  full_join(freq_gbgh_mix_e, by = "combo_total", suffix = c("_m","_e")) %>%
  arrange(desc(replace_na(n_m, 0))) %>%
  pull(combo_total)

filas_r5 <- bind_rows(
  tibble(Variable  = "Genotipos únicos gB / gH, n (%)",
         Muestras  = fmt_np(n_m_gbgh_uni, round(100*n_m_gbgh_uni/n_m_total,1)),
         Episodios = fmt_np(n_e_gbgh_uni, round(100*n_e_gbgh_uni/n_e_total,1)),
         header    = TRUE),
  
  purrr::map_dfr(combos_gbgh_uni, function(k) {
    tibble(Variable  = k,
           Muestras  = get_m(freq_gbgh_uni_m, k),
           Episodios = get_m(freq_gbgh_uni_e, k),
           header    = FALSE)
  }),
  
  tibble(Variable  = "Genotipos mixtos gB / gH, n (%)",
         Muestras  = fmt_np(n_m_gbgh_mix, round(100*n_m_gbgh_mix/n_m_total,1)),
         Episodios = fmt_np(n_e_gbgh_mix, round(100*n_e_gbgh_mix/n_e_total,1)),
         header    = TRUE),
  
  purrr::map_dfr(combos_gbgh_mix, function(k) {
    tibble(Variable  = k,
           Muestras  = get_m(freq_gbgh_mix_m, k),
           Episodios = get_m(freq_gbgh_mix_e, k),
           header    = FALSE)
  })
)

ft_r5 <- hacer_ft(
  filas_r5,
  "Combinación gB / gH",
  c("Los porcentajes de genotipos únicos se calculan respecto al total de genotipos únicos.",
    "Los porcentajes de genotipos mixtos se calculan respecto al total de genotipos mixtos.")
)

ft_r3

ft_r4

ft_r5


  

guardar_tabla(ft_r3, "Tabla_09")

guardar_tabla(ft_r4, "Tabla_10")

guardar_tabla(ft_r5, "Tabla_11")

## ---- Datos para el texto (9.5.2): frecuencias de mixtos por seroestatus ----
# 1) Preparar datos (mismos 3 serDR que estabas usando)
df <- dnaemias %>%
  filter(serDR %in% c("D-/R+", "D+/R+", "D+/R-")) %>%
  mutate(
    gb_mixto = (gbmix == 1),
    
    # gH mixto: si existe 'ghmix' lo usa; si no, lo deja en NA
    gh_mixto = if ("ghmix" %in% names(.)) (ghmix == 1) else NA,
    
    # gB+gH mixto: si existe 'gbghmix' lo usa; si no, si existe ghmix calcula ambos mixtos
    gbgh_mixto = if ("gbghmix" %in% names(.)) {
      (gbghmix == 1)
    } else if ("ghmix" %in% names(.)) {
      (gbmix == 1 & ghmix == 1)
    } else {
      NA
    }
  )

# 2) Función resumen por grupo
resumen_por_grupo <- function(data) {
  data %>%
    summarise(
      episodios = n(),
      gb_mix_n  = sum(gb_mixto, na.rm = TRUE),
      gb_mix_p  = round(100 * gb_mix_n / episodios, 1),
      
      gh_mix_n  = if (all(is.na(gh_mixto))) NA_integer_ else sum(gh_mixto, na.rm = TRUE),
      gh_mix_p  = if (all(is.na(gh_mixto))) NA_real_    else round(100 * gh_mix_n / episodios, 1),
      
      gbgh_mix_n = if (all(is.na(gbgh_mixto))) NA_integer_ else sum(gbgh_mixto, na.rm = TRUE),
      gbgh_mix_p = if (all(is.na(gbgh_mixto))) NA_real_    else round(100 * gbgh_mix_n / episodios, 1)
    )
}

# 3) Tabla por serDR + fila Global
tab <- df %>%
  group_by(serDR) %>%
  group_modify(~ resumen_por_grupo(.x)) %>%
  ungroup()

tab_global <- resumen_por_grupo(df) %>%
  mutate(serDR = "Global")

tab_final <- bind_rows(tab, tab_global) %>%
  relocate(serDR)

# 4) Formato "n (%)" y flextable final
tab_final_fmt <- tab_final %>%
  transmute(
    serDR,
    `Episodios (n)` = episodios,
    `gB mixto`      = ifelse(is.na(gb_mix_n),  NA, sprintf("%d (%.1f%%)", gb_mix_n, gb_mix_p)),
    `gH mixto`      = ifelse(is.na(gh_mix_n),  NA, sprintf("%d (%.1f%%)", gh_mix_n, gh_mix_p)),
    `gB+gH mixto`   = ifelse(is.na(gbgh_mix_n), NA, sprintf("%d (%.1f%%)", gbgh_mix_n, gbgh_mix_p))
  )

ft <- flextable(tab_final_fmt) %>%
  autofit()

ft


# 1) Preparar datos
df <- dnaemias %>%
  filter(serDR %in% c("D-/R+", "D+/R+", "D+/R-")) %>%
  mutate(
    gb_mixto = (gbmix == 1),
    gh_mixto = if ("ghmix" %in% names(.)) (ghmix == 1) else NA,
    gbgh_mixto = if ("gbghmix" %in% names(.)) {
      (gbghmix == 1)
    } else if ("ghmix" %in% names(.)) {
      (gbmix == 1 & ghmix == 1)
    } else {
      NA
    }
  )

# 2) Calcular frecuencias por grupo
tab_fig <- df %>%
  group_by(serDR) %>%
  summarise(
    episodios = n(),
    `gB mixto` = round(100 * sum(gb_mixto, na.rm = TRUE) / episodios, 1),
    `gH mixto` = if (all(is.na(gh_mixto))) NA_real_ 
    else round(100 * sum(gh_mixto, na.rm = TRUE) / episodios, 1),
    `gB/gH mixto` = if (all(is.na(gbgh_mixto))) NA_real_ 
    else round(100 * sum(gbgh_mixto, na.rm = TRUE) / episodios, 1)
  ) %>%
  ungroup()

# 3) Convertir a formato largo para ggplot
tab_long <- tab_fig %>%
  pivot_longer(
    cols = c(`gB mixto`, `gH mixto`, `gB/gH mixto`),
    names_to = "tipo_genotipo",
    values_to = "porcentaje"
  ) %>%
  # Ordenar grupos serológicos
  mutate(
    serDR = factor(serDR, levels = c("D+/R-", "D+/R+", "D-/R+")),
    tipo_genotipo = factor(tipo_genotipo, levels = c("gB mixto", "gH mixto", "gB/gH mixto"))
  )

# 4) Crear figura
figura <- ggplot(tab_long, aes(x = serDR, y = porcentaje, fill = tipo_genotipo)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    aes(label = sprintf("%.1f%%", porcentaje)),
    position = position_dodge(width = 0.8),
    vjust = -0.5,
    size = 3.5
  ) +
  scale_fill_manual(
    values = c("gB mixto" = "#E69F00", "gH mixto" = "#56B4E9", "gB/gH mixto" = "#009E73"),
    name = "Tipo de genotipo"
  ) +
  labs(
    x = "Perfil serológico D/R",
    y = "Frecuencia (%)",
    title = "Frecuencia de genotipos mixtos por episodio según perfil serológico"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.position = "bottom",
    panel.grid.major.x = element_blank()
  ) +
  ylim(0, 45)

# Mostrar figura
print(figura)

# Opcional: guardar figura
# ggsave("figura_genotipos_mixtos.png", figura, width = 8, height = 6, dpi = 300)

## ---- Figura 26. Frecuencia de genotipos mixtos por episodio, global y según perfil D/R ----
df <- dnaemias %>%
  filter(
    serDR %in% c("D-/R+", "D+/R+", "D+/R-"),
    primer_mil == 1
  ) %>%
  mutate(
    gb_mixto = gbmix == 1,
    gh_mixto = if ("ghmix" %in% names(.)) ghmix == 1 else NA,
    gbgh_mixto = if ("gbghmix" %in% names(.)) {
      gbghmix == 1
    } else if ("ghmix" %in% names(.)) {
      gbmix == 1 & ghmix == 1
    } else NA
  )

tab <- df %>%
  group_by(serDR) %>%
  summarise(
    episodios = n(),
    gb_mix = sprintf("%d (%.1f%%)",
                     sum(gb_mixto, na.rm = TRUE),
                     100 * sum(gb_mixto, na.rm = TRUE) / episodios),
    
    gh_mix = if (all(is.na(gh_mixto))) NA else
      sprintf("%d (%.1f%%)",
              sum(gh_mixto, na.rm = TRUE),
              100 * sum(gh_mixto, na.rm = TRUE) / episodios),
    
    gbgh_mix = if (all(is.na(gbgh_mixto))) NA else
      sprintf("%d (%.1f%%)",
              sum(gbgh_mixto, na.rm = TRUE),
              100 * sum(gbgh_mixto, na.rm = TRUE) / episodios)
  ) %>%
  ungroup()

# fila global
tab_global <- df %>%
  summarise(
    serDR = "Global",
    episodios = n(),
    gb_mix = sprintf("%d (%.1f%%)",
                     sum(gb_mixto, na.rm = TRUE),
                     100 * sum(gb_mixto, na.rm = TRUE) / episodios),
    gh_mix = if (all(is.na(gh_mixto))) NA else
      sprintf("%d (%.1f%%)",
              sum(gh_mixto, na.rm = TRUE),
              100 * sum(gh_mixto, na.rm = TRUE) / episodios),
    gbgh_mix = if (all(is.na(gbgh_mixto))) NA else
      sprintf("%d (%.1f%%)",
              sum(gbgh_mixto, na.rm = TRUE),
              100 * sum(gbgh_mixto, na.rm = TRUE) / episodios)
  )

tab_final <- bind_rows(tab, tab_global)

flextable(tab_final) %>%
  autofit()

# ── 1. Etiquetas del eje x con n= ─────────────────────────────────────────────
n_total <- nrow(df)
n_DR    <- sum(df$serDR == "D+/R-",  na.rm = TRUE)
n_DpR   <- sum(df$serDR == "D+/R+",  na.rm = TRUE)
n_DmR   <- sum(df$serDR == "D-/R+",  na.rm = TRUE)

etiq_x <- c(
  "Global" = paste0("Global\n(n = ", n_total, ")"),
  "D+/R-"  = paste0("D+/R-\n(n = ", n_DR,    ")"),
  "D+/R+"  = paste0("D+/R+\n(n = ", n_DpR,   ")"),
  "D-/R+"  = paste0("D-/R+\n(n = ", n_DmR,   ")")
)

# ── 2a. Tabla por seroestatus ─────────────────────────────────────────────────
tab_ser <- df %>%
  group_by(serDR) %>%
  summarise(
    episodios     = n(),
    `gB mixto`    = 100 * sum(gb_mixto,   na.rm = TRUE) / episodios,
    `gH mixto`    = if (all(is.na(gh_mixto)))    NA_real_ else
      100 * sum(gh_mixto,   na.rm = TRUE) / episodios,
    `gB/gH mixto` = if (all(is.na(gbgh_mixto))) NA_real_ else
      100 * sum(gbgh_mixto, na.rm = TRUE) / episodios,
    .groups = "drop"
  )

# ── 2b. Fila global ───────────────────────────────────────────────────────────
tab_global <- df %>%
  summarise(
    serDR         = "Global",
    episodios     = n(),
    `gB mixto`    = 100 * sum(gb_mixto,   na.rm = TRUE) / episodios,
    `gH mixto`    = if (all(is.na(gh_mixto)))    NA_real_ else
      100 * sum(gh_mixto,   na.rm = TRUE) / episodios,
    `gB/gH mixto` = if (all(is.na(gbgh_mixto))) NA_real_ else
      100 * sum(gbgh_mixto, na.rm = TRUE) / episodios
  )

# ── 2c. Combinar y pivotar ────────────────────────────────────────────────────
tab_fig <- bind_rows(tab_global, tab_ser) %>%
  pivot_longer(
    cols      = c(`gB mixto`, `gH mixto`, `gB/gH mixto`),
    names_to  = "tipo_genotipo",
    values_to = "porcentaje"
  ) %>%
  filter(!is.na(porcentaje)) %>%
  mutate(
    serDR = factor(serDR, levels = c("Global", "D+/R-", "D+/R+", "D-/R+")),
    tipo_genotipo = factor(tipo_genotipo,
                           levels = c("gB mixto", "gH mixto", "gB/gH mixto"))
  )

# ── 3. Figura ─────────────────────────────────────────────────────────────────
figura_R15 <- ggplot(tab_fig,
                     aes(x = serDR, y = porcentaje, fill = tipo_genotipo)) +
  geom_col(
    position  = position_dodge(width = 0.8),
    width     = 0.7,
    color     = "white",
    linewidth = 0.2
  ) +
  geom_vline(
    xintercept = 1.5,           # separador visual Global | por seroestatus
    color      = "grey75",
    linewidth  = 0.4,
    linetype   = "dashed"
  ) +
  geom_text(
    aes(label = paste0(formatC(porcentaje, digits = 0, format = "f",
                               decimal.mark = ","), "%")),
    position  = position_dodge(width = 0.75),
    vjust     = -0.5,
    size      = 4,
    color     = "grey20",
    na.rm     = TRUE
  ) +
  scale_fill_manual(
    values = c(
      "gB mixto"    = "#7EC0EE",
      "gH mixto"    = "#EEA9B8",
      "gB/gH mixto" = "tomato3"
    ),
    name = NULL
  ) +
  scale_x_discrete(labels = etiq_x) +
  scale_y_continuous(
    limits = c(0, 50),
    expand = expansion(mult = c(0, 0.08)),
    breaks = seq(0, 50, 10),
    labels = scales::label_number(suffix = "%", decimal.mark = ",")
  ) +
  labs(x = "", y = "Frecuencia (%)") +
  theme_minimal(base_size = 13) +
  theme(
    panel.border       = element_blank(),
    axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
    axis.line.y        = element_line(color = "grey60", linewidth = 0.4),
    axis.ticks.x       = element_blank(),
    axis.ticks.y       = element_line(color = "grey60", linewidth = 0.3),
    axis.text          = element_text(color = "grey30", size = 11),
    axis.title         = element_text(color = "grey20", size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "white", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    legend.position    = "top",
    legend.text        = element_text(size = 11),
    legend.key.size    = unit(0.45, "cm"),
    legend.spacing.x   = unit(0.4, "cm"),
    plot.margin        = margin(t = 8, r = 10, b = 4, l = 4)
  )

figura_R15

guardar_fig(figura_R15, "Figura_26", ancho = 16, alto = 11)

## ---- Tabla 12. Cinética viral según genotipo individual (genotipos únicos de gB y gH) ----
# ── 1. Función genérica ───────────────────────────────────────────────────────
resumen_bloque <- function(data, geno_var, geno_levels, geno_labels,
                           vars, var_labels, tipo_comparacion) {
  map_dfr(seq_along(vars), function(i) {
    v   <- vars[i]
    lab <- var_labels[i]
    
    celdas <- map_chr(geno_levels, function(g) {
      x <- data %>% filter(.data[[geno_var]] == g) %>%
        pull(.data[[v]]) %>% na.omit()
      if (length(x) == 0) return("–")
      sprintf("%s (%s–%s)\nn = %d",
              formatC(median(x),          format = "f", digits = 0),
              formatC(quantile(x, 0.25),  format = "f", digits = 0),
              formatC(quantile(x, 0.75),  format = "f", digits = 0),
              length(x))
    })
    
    grupos <- map(geno_levels, function(g)
      data %>% filter(.data[[geno_var]] == g) %>%
        pull(.data[[v]]) %>% na.omit())
    
    p_val <- tryCatch({
      if (tipo_comparacion == "kruskal") kruskal.test(grupos)$p.value
      else wilcox.test(grupos[[1]], grupos[[2]], exact = FALSE)$p.value
    }, error = function(e) NA_real_)
    
    p_fmt <- if (is.na(p_val)) "" else
      if (p_val < 0.001) "<0,001" else
        sub("\\.", ",", sprintf("%.3f", p_val))
    
    row        <- as_tibble_row(setNames(celdas, geno_labels))
    row$Parámetro <- lab
    row$p         <- p_fmt
    row$p_raw     <- p_val
    row
  })
}

# ── 2. Preparar datos ─────────────────────────────────────────────────────────
dnaemias <- dnaemias %>%
  mutate(
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini)
  )

vars     <- c("cvini", "cvpico", "aclaramiento_dias", "AUCCV")
var_labs <- c("CV inicial (UI/mL)", "CV pico (UI/mL)",
              "Aclaramiento (días)", "AUC viral (log10(UI/mL)·día)")

# ── 3. Bloque gB ──────────────────────────────────────────────────────────────
df_gb <- dnaemias %>%
  filter(gbsum == 1) %>%
  mutate(gb_tipo = case_when(
    gb1 == 1 ~ "gB1", gb2 == 1 ~ "gB2",
    gb3 == 1 ~ "gB3", gb4 == 1 ~ "gB4"
  ))

bloque_gb <- resumen_bloque(
  data             = df_gb,
  geno_var         = "gb_tipo",
  geno_levels      = c("gB1","gB2","gB3","gB4"),
  geno_labels      = c("gB1","gB2","gB3","gB4"),
  vars             = vars,
  var_labels       = var_labs,
  tipo_comparacion = "kruskal"
)

# ── 4. Bloque gH ──────────────────────────────────────────────────────────────
df_gh <- dnaemias %>%
  filter(ghsum == 1) %>%
  mutate(gh_tipo = case_when(
    gh1 == 1 ~ "gH1", gh2 == 1 ~ "gH2"
  ))

bloque_gh_raw <- resumen_bloque(
  data             = df_gh,
  geno_var         = "gh_tipo",
  geno_levels      = c("gH1","gH2"),
  geno_labels      = c("gH1","gH2"),
  vars             = vars,
  var_labels       = var_labs,
  tipo_comparacion = "wilcox"
)

# Añadir columnas gB3/gB4 vacías y renombrar gH1/gH2 → gB1/gB2
# (el flextable usará los nombres de columna, las cabeceras se personalizan)
bloque_gh <- bloque_gh_raw %>%
  rename(gB1 = gH1, gB2 = gH2) %>%
  mutate(gB3 = NA_character_, gB4 = NA_character_) %>%
  select(Parámetro, gB1, gB2, gB3, gB4, p, p_raw)

# ── 5. Filas cabecera ─────────────────────────────────────────────────────────
cab <- function(txt) tibble(
  Parámetro = txt,
  gB1 = NA_character_, gB2 = NA_character_,
  gB3 = NA_character_, gB4 = NA_character_,
  p   = NA_character_,  p_raw = NA_real_
)

# ── 6. Tabla completa ─────────────────────────────────────────────────────────
tabla_final <- bind_rows(
  cab("Glucoproteína gB"),
  bloque_gb,
  cab("Glucoproteína gH"),
  bloque_gh
)

# ── 7. Índices para formato ───────────────────────────────────────────────────
idx_cab <- which(tabla_final$Parámetro %in%
                   c("Glucoproteína gB", "Glucoproteína gH"))
idx_sig <- which(!is.na(tabla_final$p_raw) & tabla_final$p_raw < 0.05)
idx_gh  <- which(tabla_final$Parámetro == "Glucoproteína gH") +
  seq_len(length(vars))   # filas de datos de gH

# ── 8. Flextable ──────────────────────────────────────────────────────────────
ft <- tabla_final %>%
  select(-p_raw) %>%
  flextable() %>%
  # cabecera con nombres correctos por sección
  set_header_labels(
    Parámetro = "Parámetro",
    gB1 = "gB1 / gH1",
    gB2 = "gB2 / gH2",
    gB3 = "gB3",
    gB4 = "gB4",
    p   = "p"
  ) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  # cabeceras de sección
  bold(i = idx_cab, j = 1, part = "body") %>%
  bg(i = idx_cab,   bg = "#F2F2F2", part = "body") %>%
  # p significativas en negrita
  bold(i = idx_sig, j = "p", part = "body") %>%
  # alineación
  align(j = "Parámetro", align = "left",   part = "all") %>%
  align(j = c("gB1","gB2","gB3","gB4","p"),
        align = "center", part = "all") %>%
  # sangría en filas de datos
  padding(i = setdiff(seq_len(nrow(tabla_final)), idx_cab),
          j = 1, padding.left = 14, part = "body") %>%
  # gris para celdas vacías (gB3/gB4 en filas de gH)
  bg(i = idx_gh, j = c("gB3","gB4"), bg = "#F8F8F8", part = "body") %>%
  autofit() %>%
  add_footer_lines(values = c(
    "Datos presentados como mediana (RIC). CV: carga viral; AUC: área bajo la curva viral.",
    "Solo se incluyen episodios con genotipo único de cada glucoproteína.",
    "Glucoproteína gB: comparación entre gB1, gB2, gB3 y gB4 (Kruskal-Wallis).",
    "Glucoproteína gH: comparación entre gH1 y gH2 (U de Mann-Whitney).",
    "Las columnas gB3 y gB4 no aplican a gH (–). Los valores de p < 0,05 en negrita."
  )) %>%
  fontsize(part = "footer", size = 8) %>%
  align(part = "footer", align = "left")

ft

guardar_tabla(ft, "Tabla_12")

## ---- Figura 27. Parámetros de cinética viral entre los genotipos únicos de gB ----
# ── 1. Datos ──────────────────────────────────────────────────────────────────
df_gb <- dnaemias %>%
  mutate(
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini),
    cvini_log         = log10(cvini),
    cvpico_log        = log10(cvpico)
  ) %>%
  filter(gbsum == 1) %>%
  mutate(
    gb_tipo = factor(case_when(
      gb1 == 1 ~ "gB1", gb2 == 1 ~ "gB2",
      gb3 == 1 ~ "gB3", gb4 == 1 ~ "gB4"
    ), levels = c("gB1","gB2","gB3","gB4"))
  )

# ── 2. Comparaciones post-hoc significativas ──────────────────────────────────
comp_aclar <- tibble(
  group1       = c("gB1", "gB2"),
  group2       = c("gB3", "gB3"),
  p.adj.signif = c("**", "**"),
  y.position   = c(160, 175)
)

comp_auc <- tibble(
  group1       = c("gB1", "gB2"),
  group2       = c("gB3", "gB3"),
  p.adj.signif = c("**", "**"),
  y.position   = c(320, 355)
)

comp_vacio <- tibble(
  group1       = character(0),
  group2       = character(0),
  p.adj.signif = character(0),
  y.position   = numeric(0)
)

# ── 3. Colores: paleta de azules equilibrada (gB3 no tan oscuro) ───────────────
cols_gb <- c(
  "gB1" = "#D6EAF8",
  "gB2" = "#85C1E9",
  "gB3" = "#3498DB",
  "gB4" = "#A9CCE3"
)

cols_jitter <- c(
  "gB1" = "#2E86C1",
  "gB2" = "#2471A3",
  "gB3" = "#1B4F72",
  "gB4" = "#5499C7"
)

# ── 4. Función para construir cada panel ──────────────────────────────────────
make_panel <- function(data, var, ylab, comparaciones,
                       cap_vis = NULL, ylim_max = NULL) {
  
  d <- data %>%
    select(gb_tipo, valor = !!sym(var)) %>%
    filter(!is.na(valor))
  
  if (!is.null(cap_vis)) d <- d %>% filter(valor <= cap_vis)
  
  p <- ggplot(d, aes(x = gb_tipo, y = valor, fill = gb_tipo)) +
    geom_boxplot(
      width         = 0.55,
      color         = "grey25",
      linewidth     = 0.45,
      outlier.shape = 21,
      outlier.size  = 1.8,
      outlier.fill  = "white",
      outlier.color = "grey50",
      outlier.alpha = 0.7
    ) +
    # Mediana resaltada en blanco por encima del relleno
    stat_summary(
      fun = median, geom = "crossbar",
      width = 0.55, fatten = 0,
      color = "white", linewidth = 0.7
    ) +
    geom_jitter(
      aes(color = gb_tipo),
      width = 0.12, size = 0.9, alpha = 0.35,
      show.legend = FALSE
    ) +
    scale_fill_manual(values  = cols_gb) +
    scale_color_manual(values = cols_jitter) +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = 15) +
    theme(
      legend.position    = "none",
      axis.line          = element_line(color = "black", linewidth = 0.5),
      axis.ticks         = element_line(color = "black", linewidth = 0.4),
      axis.text.x        = element_text(size = 14, color = "black"),
      axis.text.y        = element_text(size = 12, color = "black"),
      axis.title.y       = element_text(size = 14, color = "black", face = "bold"),
      plot.margin        = margin(8, 12, 6, 6)
    )
  
  if (!is.null(ylim_max)) {
    p <- p + scale_y_continuous(
      name   = ylab,
      limits = c(NA, ylim_max),
      expand = expansion(mult = c(0.02, 0.05))
    )
  } else {
    p <- p + scale_y_continuous(
      name   = ylab,
      expand = expansion(mult = c(0.02, 0.15))
    )
  }
  
  if (nrow(comparaciones) > 0) {
    p <- p + stat_pvalue_manual(
      comparaciones,
      label        = "p.adj.signif",
      xmin         = "group1",
      xmax         = "group2",
      y.position   = "y.position",
      tip.length   = 0.01,
      bracket.size = 0.45,
      size         = 6,
      inherit.aes  = FALSE
    )
  }
  
  p
}

# ── 5. Construir los 4 paneles ────────────────────────────────────────────────
p_cvini <- make_panel(
  data          = df_gb,
  var           = "cvini_log",
  ylab          = expression(bold(log[10]*" CV inicial (UI/mL)")),
  comparaciones = comp_vacio
)

p_cvpico <- make_panel(
  data          = df_gb,
  var           = "cvpico_log",
  ylab          = expression(bold(log[10]*" CV pico (UI/mL)")),
  comparaciones = comp_vacio
)

p_aclar <- make_panel(
  data          = df_gb,
  var           = "aclaramiento_dias",
  ylab          = "Aclaramiento (d\u00edas)",
  comparaciones = comp_aclar,
  cap_vis       = 200,
  ylim_max      = 200
)

p_auc <- make_panel(
  data          = df_gb,
  var           = "AUCCV",
  ylab          = expression(bold("AUC (log"[10]*"(UI/mL)\u00b7d\u00eda)")),
  comparaciones = comp_auc,
  cap_vis       = 400,
  ylim_max      = 400
)

# ── 6. Combinar 2x2 ──────────────────────────────────────────────────────────
figura_gb_cinetica <- (p_cvini | p_cvpico) / (p_aclar | p_auc)

figura_gb_cinetica


# ── 7. Guardar en SVG ─────────────────────────────────────────────────────────
# [guardado antiguo] ggsave(
# [guardado antiguo]   filename      = "~/Desktop/Figuras_tesis/Figura_27_cinetica_gB.svg",
# [guardado antiguo]   plot          = figura_gb_cinetica,
# [guardado antiguo]   device        = svglite::svglite,
# [guardado antiguo]   fix_text_size = FALSE,
# [guardado antiguo]   width         = 16,
# [guardado antiguo]   height        = 15,
# [guardado antiguo]   units         = "cm",
# [guardado antiguo]   bg            = "white"
# [guardado antiguo] )

# ── Preparar datos ────────────────────────────────────────────────────────────
df_gb <- dnaemias %>%
  mutate(
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini)
  ) %>%
  filter(gbsum == 1) %>%
  mutate(gb_tipo = case_when(
    gb1 == 1 ~ "gB1", gb2 == 1 ~ "gB2",
    gb3 == 1 ~ "gB3", gb4 == 1 ~ "gB4"
  ))

# ── Post-hoc Dunn — Aclaramiento ──────────────────────────────────────────────
dunn.test(
  x      = df_gb$aclaramiento_dias,
  g      = df_gb$gb_tipo,
  method = "bonferroni",
  kw     = TRUE,
  label  = TRUE
)

# ── Post-hoc Dunn — AUC ───────────────────────────────────────────────────────
dunn.test(
  x      = df_gb$AUCCV,
  g      = df_gb$gb_tipo,
  method = "bonferroni",
  kw     = TRUE,
  label  = TRUE
)

guardar_fig(figura_gb_cinetica, "Figura_27", ancho = 16, alto = 15)

## ---- Tabla 13. Modelo de regresión lineal del aclaramiento en episodios con genotipo único de gB ----
# Preparación, análisis univariante (Figura R19, no incluida) y selección backward;
# la Tabla 13 corresponde al modelo clínico con variables elegidas a priori (modelo_clinico_gb3).
#  APARTADO 4.5.2 — IMPACTO DE LOS GENOTIPOS gB Y gH EN LA CINÉTICA VIRAL
#
#  Tesis CMV en RTP — Carmen, HU12O
#
#  Contenido del script:
#    1. Preparación de datos (df_gb3_uni)
#    2. Análisis descriptivo por genotipo  → Tabla R6, Figura R18
#    3. Análisis univariante (confusores)  → Figura R19
#    4. Análisis multivariante (backward)  → Tabla R7
#    5. Análisis de sensibilidad (modelo mixto)
#    6. Análisis exploratorio gB-gH (retaguardia)
#    7. Análisis de sensibilidad sin gB4 (retaguardia)
#
#  Outputs:
#    Tabla R6   → tabla_R6_cinetica_genotipos.docx (o pegada en Word)
#    Figura R18 → figura_R18_boxplots_cinetica_gB.png
#    Figura R19 → figura_R19_univariante_gB3.png
#    Tabla R7   → tabla_R7_multivariante_gB3.docx (o pegada en Word)
#
# ════════════════════════════════════════════════════════════════════════════


# ════════════════════════════════════════════════════════════════════════════
#  1. PREPARACIÓN DE DATOS
# ════════════════════════════════════════════════════════════════════════════

df_gb3_uni <- dnaemias %>%
  left_join(pacientes %>%
              select(NHC, edad_tx, tipotx, ebase,
                     SexF, edad_don),
            by = "NHC") %>%
  filter(gbsum == 1) %>%
  mutate(
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini),
    Tdesdetx_meses    = as.numeric(Fini - as.Date(Ftx)) / 30.44,
    sintomatico       = case_when(
      clinica == "VA"  ~ 0L,
      !is.na(clinica)  ~ 1L,
      TRUE             ~ NA_integer_
    ),
    risk_f       = factor(
      ifelse(serDR == "D+/R-", "Alto riesgo", "Riesgo intermedio"),
      levels = c("Riesgo intermedio", "Alto riesgo")
    ),
    gb3_f        = factor(gb3,         levels = c(0,1), labels = c("No gB3","gB3")),
    ingr_f       = factor(ingr,        levels = c(0,1), labels = c("No","Sí")),
    ing_cmv_f    = factor(ing_cmv,     levels = c(0,1), labels = c("No","Sí")),
    sintom_f     = factor(sintomatico, levels = c(0,1), labels = c("No","Sí")),
    profilaxis_f = factor(antiviral,   levels = c(0,1), labels = c("No","Sí")),
    tto_f        = factor(tto,         levels = c(0,1), labels = c("No","Sí")),
    tipotx_f     = factor(tipotx,      levels = c("UNI","BI")),
    ebase_f      = factor(ebase),
    sexo_f       = factor(SexF, levels = c(0,1), labels = c("Hombre","Mujer"))
  ) %>%
  filter(!is.na(aclaramiento_dias), !is.na(gb3_f))

cat(sprintf("\n── N episodios mono-gB: %d │ N pacientes: %d ──\n",
            nrow(df_gb3_uni), n_distinct(df_gb3_uni$NHC)))


# ════════════════════════════════════════════════════════════════════════════
#  2. ANÁLISIS DESCRIPTIVO POR GENOTIPO
#     → Tabla R6 (medianas e IQR por genotipo)
#     → Figura R18 (boxplots, código no incluido aquí)
# ════════════════════════════════════════════════════════════════════════════

# 2.1 Construir df ampliado con genotipo gB de 4 niveles ────────────────────
df_gb_4niv <- dnaemias %>%
  filter(gbsum == 1) %>%
  mutate(
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini),
    gb_tipo = case_when(
      gb1 == 1 ~ "gB1",
      gb2 == 1 ~ "gB2",
      gb3 == 1 ~ "gB3",
      gb4 == 1 ~ "gB4"
    )
  ) %>%
  filter(!is.na(gb_tipo))

# 2.2 Resumen Tabla R6: cinética por genotipo gB ────────────────────────────
cat("\n── Tabla R6: Parámetros cinéticos por genotipo gB ──\n")
df_gb_4niv %>%
  group_by(gb_tipo) %>%
  summarise(
    n_cv         = sum(!is.na(cvini)),
    cv_inicial   = sprintf("%.0f (%.0f–%.0f)",
                           median(cvini, na.rm=TRUE),
                           quantile(cvini, 0.25, na.rm=TRUE),
                           quantile(cvini, 0.75, na.rm=TRUE)),
    cv_pico      = sprintf("%.0f (%.0f–%.0f)",
                           median(cvpico, na.rm=TRUE),
                           quantile(cvpico, 0.25, na.rm=TRUE),
                           quantile(cvpico, 0.75, na.rm=TRUE)),
    n_ac         = sum(!is.na(aclaramiento_dias)),
    aclaramiento = sprintf("%.0f (%.0f–%.0f)",
                           median(aclaramiento_dias, na.rm=TRUE),
                           quantile(aclaramiento_dias, 0.25, na.rm=TRUE),
                           quantile(aclaramiento_dias, 0.75, na.rm=TRUE)),
    n_auc        = sum(!is.na(AUCCV)),
    auc_viral    = sprintf("%.0f (%.0f–%.0f)",
                           median(AUCCV, na.rm=TRUE),
                           quantile(AUCCV, 0.25, na.rm=TRUE),
                           quantile(AUCCV, 0.75, na.rm=TRUE)),
    .groups      = "drop"
  ) %>%
  print()

# 2.3 Tests globales Kruskal-Wallis (4 genotipos gB) ────────────────────────
cat("\n── Tests Kruskal-Wallis (4 genotipos gB) ──\n")
cat(sprintf("CV inicial:    p = %.3f\n",
            kruskal.test(cvini ~ gb_tipo, data = df_gb_4niv)$p.value))
cat(sprintf("CV pico:       p = %.3f\n",
            kruskal.test(cvpico ~ gb_tipo, data = df_gb_4niv)$p.value))
cat(sprintf("Aclaramiento:  p = %.3f\n",
            kruskal.test(aclaramiento_dias ~ gb_tipo, data = df_gb_4niv)$p.value))
cat(sprintf("AUC viral:     p = %.3f\n",
            kruskal.test(AUCCV ~ gb_tipo, data = df_gb_4niv)$p.value))

# 2.4 Post-hoc Dunn con corrección Bonferroni ───────────────────────────────
if (requireNamespace("FSA", quietly = TRUE)) {
  cat("\n── Post-hoc Dunn (Bonferroni) - aclaramiento ──\n")
  print(FSA::dunnTest(aclaramiento_dias ~ gb_tipo,
                      data = df_gb_4niv, method = "bonferroni"))
  
  cat("\n── Post-hoc Dunn (Bonferroni) - AUC ──\n")
  print(FSA::dunnTest(AUCCV ~ gb_tipo,
                      data = df_gb_4niv, method = "bonferroni"))
} else {
  message("Instalar FSA para post-hoc Dunn: install.packages('FSA')")
}


# ════════════════════════════════════════════════════════════════════════════
#  3. ANÁLISIS UNIVARIANTE DE CONFUSORES
#     → Figura R19 (forest plot tipo tabla con cajas por dominio)
# ════════════════════════════════════════════════════════════════════════════

# 3.1 Correlaciones de Spearman (variables continuas) ───────────────────────
cor_vars <- list(
  "Edad receptor (años)"     = "edad_tx",
  "Edad donante (años)"      = "edad_don",
  "CV inicial (UI/mL)"       = "cvini",
  "CV pico (UI/mL)"          = "cvpico",
  "Tiempo desde TXP (meses)" = "Tdesdetx_meses"
)

res_cor <- map_dfr(names(cor_vars), function(lab) {
  v <- cor_vars[[lab]]
  x <- as.numeric(df_gb3_uni[[v]])
  if (v %in% c("cvini","cvpico")) x <- log10(x + 1)
  idx <- !is.na(x) & !is.na(df_gb3_uni$aclaramiento_dias)
  test <- cor.test(x[idx], df_gb3_uni$aclaramiento_dias[idx],
                   method = "spearman", exact = FALSE)
  tibble(Variable = lab, p_val = test$p.value,
         label_est = sprintf("\u03C1 = %+.3f", test$estimate))
})

# 3.2 Wilcoxon (categóricas 2 niveles) ──────────────────────────────────────
cat_vars <- list(
  "Genotipo gB3"                      = "gb3_f",
  "Sexo (Mujer vs Hombre)"            = "sexo_f",
  "D+/R- (vs R+)"                     = "risk_f",
  "Hospitalización cualquier causa"   = "ingr_f",
  "Hospitalización por CMV"           = "ing_cmv_f",
  "Enfermedad sintomática"            = "sintom_f",
  "Profilaxis al inicio del episodio" = "profilaxis_f",
  "Tratamiento antiviral"             = "tto_f",
  "Tipo de trasplante (BI vs UNI)"    = "tipotx_f"
)

res_cat <- map_dfr(names(cat_vars), function(lab) {
  v  <- cat_vars[[lab]]
  g  <- df_gb3_uni[[v]]
  ac <- df_gb3_uni$aclaramiento_dias
  d0 <- ac[g == levels(g)[1] & !is.na(g)]
  d1 <- ac[g == levels(g)[2] & !is.na(g)]
  test <- wilcox.test(d1, d0, exact = FALSE)
  dif  <- median(d1, na.rm = TRUE) - median(d0, na.rm = TRUE)
  tibble(Variable = lab, p_val = test$p.value,
         label_est = sprintf("\u0394 = %+.0f d\u00edas", dif))
})

# 3.3 Kruskal-Wallis (categórica multinivel: enfermedad de base) ────────────
kw_eb <- kruskal.test(aclaramiento_dias ~ ebase_f, data = df_gb3_uni)
res_eb <- tibble(Variable = "Enfermedad de base",
                 p_val = kw_eb$p.value,
                 label_est = "(no calculable)")

# 3.4 Combinar y categorizar ────────────────────────────────────────────────
orden_vars <- c(
  # Demográficas
  "Edad receptor (años)",
  "Edad donante (años)",
  "Sexo (Mujer vs Hombre)",
  "Tipo de trasplante (BI vs UNI)",
  "Enfermedad de base",
  "D+/R- (vs R+)",
  # Cinéticas
  "Tiempo desde TXP (meses)",
  "CV inicial (UI/mL)",
  "CV pico (UI/mL)",
  "Genotipo gB3",
  # Clínicas
  "Enfermedad sintomática",
  "Profilaxis al inicio del episodio",
  "Tratamiento antiviral",
  "Hospitalización cualquier causa",
  "Hospitalización por CMV"
)

grupos_map <- c(
  "Edad receptor (años)"              = "Demográficas",
  "Edad donante (años)"               = "Demográficas",
  "Sexo (Mujer vs Hombre)"            = "Demográficas",
  "Tipo de trasplante (BI vs UNI)"    = "Demográficas",
  "Enfermedad de base"                = "Demográficas",
  "D+/R- (vs R+)"                     = "Demográficas",
  "Tiempo desde TXP (meses)"          = "Cinéticas",
  "CV inicial (UI/mL)"                = "Cinéticas",
  "CV pico (UI/mL)"                   = "Cinéticas",
  "Genotipo gB3"                      = "Cinéticas",
  "Enfermedad sintomática"            = "Clínicas",
  "Profilaxis al inicio del episodio" = "Clínicas",
  "Tratamiento antiviral"             = "Clínicas",
  "Hospitalización cualquier causa"   = "Clínicas",
  "Hospitalización por CMV"           = "Clínicas"
)

res_all <- bind_rows(res_cor, res_cat, res_eb) %>%
  mutate(
    sig   = p_val < 0.05,
    cand  = p_val < 0.20 & !sig,
    p_fmt = case_when(
      p_val < 0.001 ~ "<0,001",
      TRUE          ~ sub("\\.", ",", sprintf("%.3f", p_val))
    ),
    grupo    = factor(grupos_map[as.character(Variable)],
                      levels = c("Demográficas", "Cinéticas", "Clínicas")),
    Variable = factor(Variable, levels = rev(orden_vars)),
    color_cat = factor(case_when(
      sig  ~ "Significativa (p < 0,05)",
      cand ~ "Tendencia (p < 0,20)",
      TRUE ~ "No asociada"
    ), levels = c("Significativa (p < 0,05)",
                  "Tendencia (p < 0,20)",
                  "No asociada"))
  )

cat("\n── Univariante gB3: resultados ordenados por p ──\n")
print(res_all %>% arrange(p_val) %>%
        select(Variable, label_est, p_fmt), n = Inf)


# 3.5 Figura R19 ─────────────────────────────────────────────────────────────
X_EST_COL  <- -0.05
X_PLOT_MIN <-  0.00
X_PLOT_MAX <-  1.10
X_P_COL    <-  1.12
X_LEFT     <- -0.40
X_RIGHT    <-  1.32

cabeceras <- tibble(
  grupo = factor("Demográficas",
                 levels = c("Demográficas", "Cinéticas", "Clínicas")),
  Variable = factor("Edad receptor (años)",
                    levels = levels(res_all$Variable))
)

figura_R19_univariante_gb3 <- ggplot(res_all, aes(x = p_val, y = Variable)) +
  annotate("rect",
           xmin = X_PLOT_MIN, xmax = X_PLOT_MAX,
           ymin = -Inf, ymax = Inf,
           fill = "grey97", color = NA) +
  geom_vline(xintercept = 0.05, color = "firebrick3",
             linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  geom_vline(xintercept = 0.20, color = "grey55",
             linetype = "dashed", linewidth = 0.4, alpha = 0.6) +
  geom_text(aes(x = X_EST_COL, label = label_est),
            hjust = 1, size = 3.6, color = "grey25",
            show.legend = FALSE) +
  geom_point(aes(color = color_cat), size = 3.0) +
  geom_text(aes(x = X_P_COL, label = p_fmt, color = color_cat),
            hjust = 0, size = 3.6, fontface = "bold",
            show.legend = FALSE) +
  geom_text(data = cabeceras,
            aes(x = X_EST_COL, y = Inf, label = "Estimador"),
            hjust = 1, vjust = -0.5, size = 3.8,
            fontface = "bold", color = "grey15",
            inherit.aes = FALSE) +
  geom_text(data = cabeceras,
            aes(x = (X_PLOT_MIN + X_PLOT_MAX) / 2, y = Inf,
                label = "valor p"),
            hjust = 0.5, vjust = -0.5, size = 3.8,
            fontface = "bold", color = "grey15",
            inherit.aes = FALSE) +
  facet_grid(grupo ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  scale_x_continuous(
    name   = NULL,
    limits = c(X_LEFT, X_RIGHT),
    breaks = c(0.05, 0.20, 0.50, 1.00),
    labels = c("0,05", "0,20", "0,50", "1,00"),
    expand = c(0, 0)
  ) +
  scale_color_manual(
    name   = NULL,
    values = c(
      "Significativa (p < 0,05)" = "firebrick3",
      "Tendencia (p < 0,20)"     = "#E08E45",
      "No asociada"              = "grey55"
    )
  ) +
  coord_cartesian(clip = "off") +
  labs(y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position      = "bottom",
    legend.text          = element_text(size = 10),
    legend.key.size      = unit(10, "pt"),
    legend.margin        = margin(t = 4, b = 0),
    legend.box.spacing   = unit(2, "pt"),
    panel.grid.major.y   = element_line(color = "grey95", linewidth = 0.3),
    panel.grid.major.x   = element_blank(),
    panel.grid.minor     = element_blank(),
    axis.line.x          = element_line(color = "grey60", linewidth = 0.4),
    axis.ticks.x         = element_line(color = "grey60", linewidth = 0.3),
    axis.ticks.length.x  = unit(3, "pt"),
    axis.text.y          = element_text(size = 11, color = "grey15"),
    axis.text.x          = element_text(size = 10, color = "grey40"),
    strip.placement      = "outside",
    strip.text.y.left    = element_text(angle = 90, size = 11.5,
                                        face = "bold", color = "white",
                                        margin = margin(r = 5, l = 5)),
    strip.background.y   = element_rect(fill = "grey35", color = NA),
    panel.spacing.y      = unit(8, "pt"),
    panel.border         = element_rect(color = "grey80",
                                        fill = NA, linewidth = 0.4),
    plot.margin          = margin(28, 14, 8, 8)
  )

figura_R19_univariante_gb3

# [guardado antiguo] ggsave("~/Desktop/Figuras_tesis/figura_R19_univariante_gB3.png",
# [guardado antiguo]        figura_R19_univariante_gb3,
# [guardado antiguo]        width = 9, height = 6, dpi = 320, bg = "white")


# ════════════════════════════════════════════════════════════════════════════
#  4. ANÁLISIS MULTIVARIANTE (regresión lineal con backward)
#     → Tabla R7
# ════════════════════════════════════════════════════════════════════════════

# 4.1 Construir df con variables candidatas (p < 0,20) ──────────────────────
df_multi_gb3 <- df_gb3_uni %>%
  select(NHC, aclaramiento_dias,
         gb3_f, risk_f, cvini, cvpico, Tdesdetx_meses,
         sintom_f, profilaxis_f, ing_cmv_f, ingr_f, tto_f) %>%
  filter(complete.cases(.))

cat(sprintf("\n── N para multivariante (datos completos): %d episodios / %d pacientes ──\n",
            nrow(df_multi_gb3), n_distinct(df_multi_gb3$NHC)))

# 4.2 Modelo completo ────────────────────────────────────────────────────────
modelo_completo_gb3 <- lm(
  aclaramiento_dias ~ gb3_f + risk_f + log10(cvini + 1) + log10(cvpico + 1) +
    Tdesdetx_meses + sintom_f + profilaxis_f + ing_cmv_f + ingr_f,
  data = df_multi_gb3
)

cat("\n── Modelo completo (todas las candidatas) ──\n")
print(summary(modelo_completo_gb3))

# 4.3 Selección backward (Akaike) ───────────────────────────────────────────
modelo_backward_gb3 <- step(modelo_completo_gb3, direction = "backward",
                            trace = FALSE)

cat("\n══════════════════════════════════════════════\n")
cat("  MODELO FINAL (post-backward)\n")
cat("══════════════════════════════════════════════\n")
print(summary(modelo_backward_gb3))

cat("\n── IC 95% del modelo final ──\n")
print(confint(modelo_backward_gb3))

cat("\n── VIF (multicolinealidad) ──\n")
print(car::vif(modelo_backward_gb3))


# ════════════════════════════════════════════════════════════════════════════
#  5. ANÁLISIS DE SENSIBILIDAD: modelo lineal mixto (paciente aleatorio)
# ════════════════════════════════════════════════════════════════════════════

modelo_mixto_gb3 <- lmer(
  aclaramiento_dias ~ gb3_f + log10(cvini + 1) + log10(cvpico + 1) +
    profilaxis_f + ing_cmv_f + (1 | NHC),
  data = df_multi_gb3
)

cat("\n══════════════════════════════════════════════\n")
cat("  MODELO MIXTO: paciente como efecto aleatorio\n")
cat("══════════════════════════════════════════════\n")
print(summary(modelo_mixto_gb3))

cat("\n── IC 95% (Wald) ──\n")
print(confint(modelo_mixto_gb3, method = "Wald"))

vc  <- as.data.frame(VarCorr(modelo_mixto_gb3))
icc <- vc$vcov[1] / sum(vc$vcov)
cat(sprintf("\n── ICC (varianza atribuible al paciente): %.3f ──\n", icc))

# 5.1 Comparación lineal vs mixto ───────────────────────────────────────────
coef_lin <- coef(modelo_backward_gb3)["gb3_fgB3"]
ci_lin   <- confint(modelo_backward_gb3)["gb3_fgB3", ]
p_lin    <- summary(modelo_backward_gb3)$coefficients["gb3_fgB3", "Pr(>|t|)"]

coef_mix <- fixef(modelo_mixto_gb3)["gb3_fgB3"]
ci_mix   <- confint(modelo_mixto_gb3, method = "Wald")["gb3_fgB3", ]
p_mix    <- summary(modelo_mixto_gb3)$coefficients["gb3_fgB3", "Pr(>|t|)"]

cat("\n── Comparación LINEAL vs MIXTO (efecto gB3) ──\n")
cat(sprintf("  Lineal:  %+.1f días (IC 95%%: %+.1f a %+.1f), p = %.4f\n",
            coef_lin, ci_lin[1], ci_lin[2], p_lin))
cat(sprintf("  Mixto:   %+.1f días (IC 95%%: %+.1f a %+.1f), p = %.4f\n",
            coef_mix, ci_mix[1], ci_mix[2], p_mix))


# ════════════════════════════════════════════════════════════════════════════
#  6. ANÁLISIS EXPLORATORIO: combinaciones gB-gH (RETAGUARDIA)
# ════════════════════════════════════════════════════════════════════════════

df_gbgh_combo <- dnaemias %>%
  filter(gbsum == 1, ghsum == 1) %>%
  mutate(
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini),
    gb_tipo = case_when(
      gb1 == 1 ~ "gB1", gb2 == 1 ~ "gB2",
      gb3 == 1 ~ "gB3", gb4 == 1 ~ "gB4"
    ),
    gh_tipo = case_when(
      gh1 == 1 ~ "gH1", gh2 == 1 ~ "gH2"
    ),
    combo = paste0(gb_tipo, "-", gh_tipo)
  ) %>%
  filter(!is.na(aclaramiento_dias), !is.na(combo))

cat("\n── Combinaciones gB-gH (8 grupos) ──\n")
print(
  df_gbgh_combo %>%
    group_by(combo) %>%
    summarise(
      n       = n(),
      mediana = median(aclaramiento_dias),
      q1      = quantile(aclaramiento_dias, 0.25),
      q3      = quantile(aclaramiento_dias, 0.75),
      .groups = "drop"
    ) %>% arrange(desc(mediana))
)

cat("\n── Kruskal-Wallis global (8 combinaciones) ──\n")
print(kruskal.test(aclaramiento_dias ~ combo, data = df_gbgh_combo))

# 6.1 Foco en gB3: ¿efecto del gH? ──────────────────────────────────────────
df_gb3_combo <- df_gbgh_combo %>% filter(gb_tipo == "gB3")

cat("\n── Dentro de gB3: gH1 vs gH2 ──\n")
print(df_gb3_combo %>%
        group_by(combo) %>%
        summarise(n = n(),
                  mediana = median(aclaramiento_dias),
                  .groups = "drop"))

cat("\n── Wilcoxon gB3-gH1 vs gB3-gH2 ──\n")
print(wilcox.test(aclaramiento_dias ~ combo,
                  data = df_gb3_combo, exact = FALSE))


# ════════════════════════════════════════════════════════════════════════════
#  7. ANÁLISIS DE SENSIBILIDAD: excluyendo gB4 (RETAGUARDIA)
# ════════════════════════════════════════════════════════════════════════════

df_gb_no_gb4 <- df_gbgh_combo %>% filter(gb_tipo != "gB4")

cat("\n── Kruskal-Wallis sin gB4 (3 grupos: gB1, gB2, gB3) ──\n")
print(kruskal.test(aclaramiento_dias ~ gb_tipo, data = df_gb_no_gb4))


# ════════════════════════════════════════════════════════════════════════════
#  FIN DEL SCRIPT
# ════════════════════════════════════════════════════════════════════════════

# ── Modelo clínico gB3 (variables elegidas a priori) ──────────────────────────
df_clin_gb3 <- df_gb3_uni %>%
  select(NHC, aclaramiento_dias,
         gb3_f, risk_f, cvini, cvpico, profilaxis_f) %>%
  filter(complete.cases(.))

cat(sprintf("\n── N modelo clínico: %d episodios / %d pacientes ──\n",
            nrow(df_clin_gb3), n_distinct(df_clin_gb3$NHC)))

modelo_clinico_gb3 <- lm(
  aclaramiento_dias ~ gb3_f + risk_f + log10(cvini + 1) + log10(cvpico + 1) +
    profilaxis_f,
  data = df_clin_gb3
)

cat("\n══════════════════════════════════════════════\n")
cat("  MODELO CLÍNICO gB3 (variables a priori)\n")
cat("══════════════════════════════════════════════\n")
print(summary(modelo_clinico_gb3))

cat("\n── IC 95% ──\n")
print(confint(modelo_clinico_gb3))

cat("\n── VIF (multicolinealidad) ──\n")
print(car::vif(modelo_clinico_gb3))

# Efecto ajustado de gB3 en una línea
coef_gb3 <- coef(modelo_clinico_gb3)["gb3_fgB3"]
ci_gb3   <- confint(modelo_clinico_gb3)["gb3_fgB3", ]
p_gb3    <- summary(modelo_clinico_gb3)$coefficients["gb3_fgB3", "Pr(>|t|)"]

cat(sprintf("\n── Efecto ajustado gB3: %+.1f días (IC 95%%: %+.1f a %+.1f), p = %.4f ──\n",
            coef_gb3, ci_gb3[1], ci_gb3[2], p_gb3))


modelo_mixto_clinico <- lmer(
  aclaramiento_dias ~ gb3_f + risk_f + log10(cvini + 1) + log10(cvpico + 1) +
    profilaxis_f + (1 | NHC),
  data = df_clin_gb3
)

print(summary(modelo_mixto_clinico))
print(confint(modelo_mixto_clinico, method = "Wald")["gb3_fgB3", ])

vc  <- as.data.frame(VarCorr(modelo_mixto_clinico))
icc <- vc$vcov[1] / sum(vc$vcov)
cat(sprintf("\nICC: %.3f\n", icc))

guardar_modelo(modelo_clinico_gb3, "Tabla_13_modelo")

## ---- Tabla 14. Cinética viral según presencia de genotipos mixtos de gB, gH y gB/gH ----
# ── 1. Calcular aclaramiento si no existe ─────────────────────────────────────
dnaemias <- dnaemias %>%
  mutate(
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini)
  )

# ── 2. Funciones auxiliares ───────────────────────────────────────────────────
fmt_med <- function(x) {
  n   <- sum(!is.na(x))
  med <- median(x, na.rm = TRUE)
  q1  <- quantile(x, 0.25, na.rm = TRUE)
  q3  <- quantile(x, 0.75, na.rm = TRUE)
  sprintf("%s (%s\u2013%s)\n(n=%d)",
          formatC(med, format="f", digits=0),
          formatC(q1,  format="f", digits=0),
          formatC(q3,  format="f", digits=0), n)
}

fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0,001*")
  txt <- sprintf("%.3f", p)
  txt <- sub("\\.", ",", txt)
  if (p < 0.05) paste0(txt, "*") else txt
}

wilcox_p <- function(data, var, grupo) {
  d <- data %>% filter(!is.na(.data[[var]]), !is.na(.data[[grupo]]))
  if (nrow(d) < 4) return(NA_real_)
  tryCatch(
    wilcox.test(as.formula(paste(var, "~", grupo)), data = d)$p.value,
    error = function(e) NA_real_
  )
}

# ── 3. Función que construye 4 filas para un grupo genotípico ─────────────────
bloque <- function(data, var_mix) {
  d0 <- data %>% filter(.data[[var_mix]] == 0)
  d1 <- data %>% filter(.data[[var_mix]] == 1)
  
  vars  <- c("cvini", "cvpico", "aclaramiento_dias", "AUCCV")
  labs  <- c("CV inicial (UI/mL)", "CV pico (UI/mL)",
             "Aclaramiento (días)", "AUC viral")
  
  tibble(
    Parámetro = labs,
    `Genotipo único` = sapply(vars, function(v) fmt_med(d0[[v]])),
    `Genotipo mixto`  = sapply(vars, function(v) fmt_med(d1[[v]])),
    p = sapply(vars, function(v) fmt_p(wilcox_p(data, v, var_mix)))
  )
}

# ── 4. Construir tabla completa ───────────────────────────────────────────────
cabecera <- function(txt) tibble(
  Parámetro        = txt,
  `Genotipo único` = NA_character_,
  `Genotipo mixto` = NA_character_,
  p                = NA_character_
)

tabla_cin <- bind_rows(
  cabecera("Glucoproteína gB"),
  bloque(dnaemias, "gbmix"),
  cabecera("Glucoproteína gH"),
  bloque(dnaemias, "ghmix"),
  cabecera("Glucoproteínas gB y gH"),
  bloque(dnaemias, "gbghmix")
)

# ── 5. Flextable ──────────────────────────────────────────────────────────────
# Índices de filas cabecera
idx_cab <- which(tabla_cin$Parámetro %in%
                   c("Glucoproteína gB", "Glucoproteína gH",
                     "Glucoproteínas gB y gH"))

ft_cin <- flextable(tabla_cin) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  bold(i = idx_cab, j = 1, part = "body") %>%
  bg(i = idx_cab, bg = "#F2F2F2", part = "body") %>%
  align(j = "Parámetro", align = "left",   part = "all") %>%
  align(j = c("Genotipo único", "Genotipo mixto", "p"),
        align = "center", part = "all") %>%
  padding(i = setdiff(seq_len(nrow(tabla_cin)), idx_cab),
          j = 1, padding.left = 14, part = "body") %>%
  autofit() %>%
  add_footer_lines(values = c(
    "Los datos se presentan como mediana (RIC).",
    "CV: carga viral; AUC: área bajo la curva.",
    "Comparaciones mediante la prueba U de Mann-Whitney.",
    "*p < 0,05."
  )) %>%
  fontsize(part = "footer", size = 8) %>%
  align(part = "footer", align = "left")

ft_cin

guardar_tabla(ft_cin, "Tabla_14")

## ---- Datos para el texto (9.5.3): cinética por seroestatus, todos los episodios y primer episodio ----
# ── 1. Datos ──────────────────────────────────────────────────────────────────
df <- dnaemias

vars <- c("cvini", "cvpico", "diasdur_vir", "AUCCV")
var_labels <- c(
  cvini       = "CV inicial",
  cvpico      = "CV pico",
  diasdur_vir = "Aclaramiento",
  AUCCV       = "AUC viral"
)

# ── 2. Función de cálculo ─────────────────────────────────────────────────────
calcular_p <- function(mixvar, mix_label) {
  map_dfr(vars, function(v) {
    d0 <- df %>% filter(.data[[mixvar]] == 0) %>% pull(.data[[v]])
    d1 <- df %>% filter(.data[[mixvar]] == 1) %>% pull(.data[[v]])
    p_global <- tryCatch(
      wilcox.test(d0, d1, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
    grupos_res <- map_dfr(c("D-/R+","D+/R+","D+/R-"), function(grupo) {
      dg <- df %>% filter(serDR == grupo)
      d0g <- dg %>% filter(.data[[mixvar]] == 0) %>% pull(.data[[v]])
      d1g <- dg %>% filter(.data[[mixvar]] == 1) %>% pull(.data[[v]])
      p <- if (length(na.omit(d0g)) > 0 & length(na.omit(d1g)) > 0)
        tryCatch(wilcox.test(d0g, d1g, exact = FALSE)$p.value,
                 error = function(e) NA_real_)
      else NA_real_
      tibble(serDR = grupo, p = p)
    }) %>%
      pivot_wider(names_from = serDR, values_from = p,
                  names_prefix = "p_")
    bind_cols(
      tibble(Glucoproteína = mix_label, Variable = var_labels[v], p_global),
      grupos_res
    )
  })
}

# ── 3. Calcular y combinar ────────────────────────────────────────────────────
tabla_raw <- bind_rows(
  calcular_p("gbmix",   "gB"),
  calcular_p("ghmix",   "gH"),
  calcular_p("gbghmix", "gB+gH")
)

# ── 4. Formatear p-valores ────────────────────────────────────────────────────
fmt_p <- function(p) {
  if (is.na(p))    return("")
  if (p < 0.001)   return("<0,001")
  txt <- sub("\\.", ",", sprintf("%.3f", p))
  if (p < 0.05) paste0(txt, "*") else txt
}

tabla_fmt <- tabla_raw %>%
  mutate(across(where(is.numeric), ~ sapply(., fmt_p)))

# ── 5. Añadir filas cabecera por glucoproteína ────────────────────────────────
cab <- function(txt) tibble(
  Glucoproteína = txt,
  Variable      = NA_character_,
  p_global      = NA_character_,
  `p_D-/R+`     = NA_character_,
  `p_D+/R+`     = NA_character_,
  `p_D+/R-`     = NA_character_
)

tabla_con_cab <- bind_rows(
  cab("Glucoproteína gB"),
  tabla_fmt %>% filter(Glucoproteína == "gB"),
  cab("Glucoproteína gH"),
  tabla_fmt %>% filter(Glucoproteína == "gH"),
  cab("Glucoproteínas gB+gH"),
  tabla_fmt %>% filter(Glucoproteína == "gB+gH")
)

# Recalcular índices sobre la tabla con cabeceras
idx_cab <- which(is.na(tabla_con_cab$Variable))

# Recalcular índices sig sobre tabla_raw reordenada
tabla_raw_ord <- bind_rows(
  tabla_raw %>% filter(Glucoproteína == "gB"),
  tabla_raw %>% filter(Glucoproteína == "gH"),
  tabla_raw %>% filter(Glucoproteína == "gB+gH")
)

# Posiciones en tabla_con_cab (desplazadas por las filas cabecera)
offset <- function(idx_raw, idx_cab_vec) {
  sapply(idx_raw, function(i) {
    i + sum(idx_cab_vec <= (i + sum(idx_cab_vec < i)))
  })
}

idx_sig_global <- which(!is.na(tabla_raw_ord$p_global) & tabla_raw_ord$p_global < 0.05)
idx_sig_DmR    <- which(!is.na(tabla_raw_ord$`p_D-/R+`) & tabla_raw_ord$`p_D-/R+` < 0.05)
idx_sig_DpR    <- which(!is.na(tabla_raw_ord$`p_D+/R+`) & tabla_raw_ord$`p_D+/R+` < 0.05)
idx_sig_DpRm   <- which(!is.na(tabla_raw_ord$`p_D+/R-`) & tabla_raw_ord$`p_D+/R-` < 0.05)

# Convertir índices raw a índices en tabla con cabeceras
to_table_idx <- function(raw_idx) {
  # Cada bloque de 4 filas tiene 1 cabecera antes
  raw_idx + ceiling(raw_idx / 4)
}

ft <- tabla_con_cab %>%
  select(-Glucoproteína) %>%
  flextable() %>%
  set_header_labels(
    Variable  = "Variable",
    p_global  = "Global",
    `p_D-/R+` = "D-/R+",
    `p_D+/R+` = "D+/R+",
    `p_D+/R-` = "D+/R-"
  ) %>%
  add_header_row(
    values     = c("", "p (único vs mixto)"),
    colwidths  = c(1, 4)
  ) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  bold(i = idx_cab, j = 1, part = "body") %>%
  bg(i = idx_cab, bg = "#F2F2F2", part = "body") %>%
  padding(i = setdiff(seq_len(nrow(tabla_con_cab)), idx_cab),
          j = 1, padding.left = 14, part = "body") %>%
  bold(i = to_table_idx(idx_sig_global), j = "p_global",  part = "body") %>%
  bold(i = to_table_idx(idx_sig_DmR),    j = "p_D-/R+",  part = "body") %>%
  bold(i = to_table_idx(idx_sig_DpR),    j = "p_D+/R+",  part = "body") %>%
  bold(i = to_table_idx(idx_sig_DpRm),   j = "p_D+/R-",  part = "body") %>%
  align(j = "Variable", align = "left",   part = "all") %>%
  align(j = c("p_global","p_D-/R+","p_D+/R+","p_D+/R-"),
        align = "center", part = "all") %>%
  autofit() %>%
  add_footer_lines(values = c(
    "Comparaciones mediante la prueba U de Mann-Whitney (único vs mixto).",
    "Los p-valores se muestran para el análisis global y estratificado por perfil serológico D/R.",
    "*p < 0,05."
  )) %>%
  fontsize(part = "footer", size = 8) %>%
  align(part = "footer", align = "left")

ft

# ── 1. Datos — solo primer episodio por paciente ──────────────────────────────
df <- dnaemias %>% filter(primer_mil == 1)

vars <- c("cvini", "cvpico", "diasdur_vir", "AUCCV")
var_labels <- c(
  cvini       = "CV inicial",
  cvpico      = "CV pico",
  diasdur_vir = "Aclaramiento",
  AUCCV       = "AUC viral"
)

# ── 2. Función de cálculo ─────────────────────────────────────────────────────
calcular_p <- function(mixvar, mix_label) {
  map_dfr(vars, function(v) {
    d0 <- df %>% filter(.data[[mixvar]] == 0) %>% pull(.data[[v]])
    d1 <- df %>% filter(.data[[mixvar]] == 1) %>% pull(.data[[v]])
    p_global <- tryCatch(
      wilcox.test(d0, d1, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
    grupos_res <- map_dfr(c("D-/R+","D+/R+","D+/R-"), function(grupo) {
      dg <- df %>% filter(serDR == grupo)
      d0g <- dg %>% filter(.data[[mixvar]] == 0) %>% pull(.data[[v]])
      d1g <- dg %>% filter(.data[[mixvar]] == 1) %>% pull(.data[[v]])
      p <- if (length(na.omit(d0g)) > 0 & length(na.omit(d1g)) > 0)
        tryCatch(wilcox.test(d0g, d1g, exact = FALSE)$p.value,
                 error = function(e) NA_real_)
      else NA_real_
      tibble(serDR = grupo, p = p)
    }) %>%
      pivot_wider(names_from = serDR, values_from = p,
                  names_prefix = "p_")
    bind_cols(
      tibble(Glucoproteína = mix_label, Variable = var_labels[v], p_global),
      grupos_res
    )
  })
}

# ── 3. Calcular y combinar ────────────────────────────────────────────────────
tabla_raw <- bind_rows(
  calcular_p("gbmix",   "gB"),
  calcular_p("ghmix",   "gH"),
  calcular_p("gbghmix", "gB+gH")
)

# ── 4. Formatear p-valores ────────────────────────────────────────────────────
fmt_p <- function(p) {
  if (is.na(p))    return("")
  if (p < 0.001)   return("<0,001")
  txt <- sub("\\.", ",", sprintf("%.3f", p))
  if (p < 0.05) paste0(txt, "*") else txt
}

tabla_fmt <- tabla_raw %>%
  mutate(across(where(is.numeric), ~ sapply(., fmt_p)))

# ── 5. Añadir filas cabecera por glucoproteína ────────────────────────────────
cab <- function(txt) tibble(
  Glucoproteína = txt,
  Variable      = NA_character_,
  p_global      = NA_character_,
  `p_D-/R+`     = NA_character_,
  `p_D+/R+`     = NA_character_,
  `p_D+/R-`     = NA_character_
)

tabla_con_cab <- bind_rows(
  cab("Glucoproteína gB"),
  tabla_fmt %>% filter(Glucoproteína == "gB"),
  cab("Glucoproteína gH"),
  tabla_fmt %>% filter(Glucoproteína == "gH"),
  cab("Glucoproteínas gB+gH"),
  tabla_fmt %>% filter(Glucoproteína == "gB+gH")
)

# Recalcular índices sobre la tabla con cabeceras
idx_cab <- which(is.na(tabla_con_cab$Variable))

# Recalcular índices sig sobre tabla_raw reordenada
tabla_raw_ord <- bind_rows(
  tabla_raw %>% filter(Glucoproteína == "gB"),
  tabla_raw %>% filter(Glucoproteína == "gH"),
  tabla_raw %>% filter(Glucoproteína == "gB+gH")
)

# Posiciones en tabla_con_cab (desplazadas por las filas cabecera)
offset <- function(idx_raw, idx_cab_vec) {
  sapply(idx_raw, function(i) {
    i + sum(idx_cab_vec <= (i + sum(idx_cab_vec < i)))
  })
}

idx_sig_global <- which(!is.na(tabla_raw_ord$p_global) & tabla_raw_ord$p_global < 0.05)
idx_sig_DmR    <- which(!is.na(tabla_raw_ord$`p_D-/R+`) & tabla_raw_ord$`p_D-/R+` < 0.05)
idx_sig_DpR    <- which(!is.na(tabla_raw_ord$`p_D+/R+`) & tabla_raw_ord$`p_D+/R+` < 0.05)
idx_sig_DpRm   <- which(!is.na(tabla_raw_ord$`p_D+/R-`) & tabla_raw_ord$`p_D+/R-` < 0.05)

# Convertir índices raw a índices en tabla con cabeceras
to_table_idx <- function(raw_idx) {
  # Cada bloque de 4 filas tiene 1 cabecera antes
  raw_idx + ceiling(raw_idx / 4)
}

ft <- tabla_con_cab %>%
  select(-Glucoproteína) %>%
  flextable() %>%
  set_header_labels(
    Variable  = "Variable",
    p_global  = "Global",
    `p_D-/R+` = "D-/R+",
    `p_D+/R+` = "D+/R+",
    `p_D+/R-` = "D+/R-"
  ) %>%
  add_header_row(
    values     = c("", "PRIMER EPISODIO p (único vs mixto)"),
    colwidths  = c(1, 4)
  ) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  bold(i = idx_cab, j = 1, part = "body") %>%
  bg(i = idx_cab, bg = "#F2F2F2", part = "body") %>%
  padding(i = setdiff(seq_len(nrow(tabla_con_cab)), idx_cab),
          j = 1, padding.left = 14, part = "body") %>%
  bold(i = to_table_idx(idx_sig_global), j = "p_global",  part = "body") %>%
  bold(i = to_table_idx(idx_sig_DmR),    j = "p_D-/R+",  part = "body") %>%
  bold(i = to_table_idx(idx_sig_DpR),    j = "p_D+/R+",  part = "body") %>%
  bold(i = to_table_idx(idx_sig_DpRm),   j = "p_D+/R-",  part = "body") %>%
  align(j = "Variable", align = "left",   part = "all") %>%
  align(j = c("p_global","p_D-/R+","p_D+/R+","p_D+/R-"),
        align = "center", part = "all") %>%
  autofit() %>%
  add_footer_lines(values = c(
    "Comparaciones mediante la prueba U de Mann-Whitney (único vs mixto).",
    "Los p-valores se muestran para el análisis global y estratificado por perfil serológico D/R.",
    "*p < 0,05."
  )) %>%
  fontsize(part = "footer", size = 8) %>%
  align(part = "footer", align = "left")

ft

## ---- Figura 28. Cinética viral según genotipos mixtos, estratificada por perfil serológico ----
# ── Variables de control ──────────────────────────────────────────────────────
id_patient <- "NHC"
cv_y_max   <- 7

## Etiquetas de las filas (orden de arriba abajo, como en la tesis)
## Se escriben como expresiones plotmath (log[10] = subíndice) y se
## interpretan con label_parsed en facet_grid()
## atop() reparte cada etiqueta en dos líneas para que quepa en la cabecera
lab_acl  <- 'atop(bold("Aclaramiento"), bold("(días)"))'
lab_auc  <- 'atop(bold("AUC"), bold("(log"[10]*"(UI/mL)·día)"))'
lab_cvi  <- 'atop(bold("log"[10]*" CV inicial"), bold("(UI/mL)"))'
lab_cvp  <- 'atop(bold("log"[10]*" CV pico"), bold("(UI/mL)"))'
orden_var <- c(lab_acl, lab_auc, lab_cvi, lab_cvp)

df_sub <- dnaemias

vars_to_pivot <- c("diasdur_vir", "AUCCV", "cvpico", "cvini")

# ── Formato largo ─────────────────────────────────────────────────────────────
## case_when() en lugar de recode(): car::recode() enmascara a dplyr::recode()
df_long <- df_sub %>%
  mutate(
    gB   = ifelse(gbmix   == 1, "gB mixto",    "gB único"),
    gH   = ifelse(ghmix   == 1, "gH mixto",    "gH único"),
    gBGH = ifelse(gbghmix == 1, "gB+gH mixto", "gB+gH único")
  ) %>%
  pivot_longer(cols = c(gB, gH, gBGH),
               names_to  = "tipo_gen",
               values_to = "genotipo") %>%
  select(serDR, tipo_gen, genotipo,
         !!sym(id_patient), all_of(vars_to_pivot)) %>%
  pivot_longer(cols      = all_of(vars_to_pivot),
               names_to  = "variable",
               values_to = "valor_raw") %>%
  mutate(
    valor = case_when(
      variable %in% c("cvpico", "cvini") & !is.na(valor_raw) ~
        log10(as.numeric(valor_raw) + 1),
      TRUE ~ as.numeric(valor_raw)
    ),
    variable = case_when(
      variable == "diasdur_vir" ~ lab_acl,
      variable == "AUCCV"       ~ lab_auc,
      variable == "cvini"       ~ lab_cvi,
      variable == "cvpico"      ~ lab_cvp
    ),
    serDR    = factor(serDR, levels = c("D-/R+", "D+/R+", "D+/R-")),
    genotipo = factor(genotipo, levels = c(
      "gB único",    "gB mixto",
      "gH único",    "gH mixto",
      "gB+gH único", "gB+gH mixto"
    ))
  )

# ── Posiciones X ──────────────────────────────────────────────────────────────
x_map <- c(
  "gB único"    = 1,   "gB mixto"    = 2,
  "gH único"    = 3.3, "gH mixto"    = 4.3,
  "gB+gH único" = 5.6, "gB+gH mixto" = 6.6
)
df_plot <- df_long %>%
  mutate(x_pos = x_map[as.character(genotipo)])

# ── Límites visuales por panel ────────────────────────────────────────────────
caps_dyn <- tibble(
  variable = orden_var,
  cap_vis  = c(100, 400, cv_y_max, cv_y_max)
)

df_plot_vis <- df_plot %>%
  left_join(caps_dyn, by = "variable") %>%
  filter(is.na(valor) | valor <= cap_vis) %>%
  mutate(variable = factor(variable, levels = orden_var))

# ── Comparaciones Wilcoxon sobre datos SIN filtrar (solo significativas) ──────
comparaciones <- df_plot %>%
  group_by(variable, serDR, tipo_gen) %>%
  wilcox_test(valor ~ genotipo) %>%
  add_significance("p") %>%
  ungroup() %>%
  filter(p.signif != "ns") %>%
  mutate(
    group1 = case_when(
      tipo_gen == "gB"   ~ "gB único",
      tipo_gen == "gH"   ~ "gH único",
      tipo_gen == "gBGH" ~ "gB+gH único"
    ),
    group2 = case_when(
      tipo_gen == "gB"   ~ "gB mixto",
      tipo_gen == "gH"   ~ "gH mixto",
      tipo_gen == "gBGH" ~ "gB+gH mixto"
    ),
    xmin = x_map[group1],
    xmax = x_map[group2]
  )

base_positions <- df_plot_vis %>%
  mutate(variable = as.character(variable)) %>%
  group_by(variable, serDR) %>%
  summarise(maxval  = max(valor,   na.rm = TRUE),
            cap_vis = first(cap_vis), .groups = "drop") %>%
  mutate(base = pmin(maxval * 1.12, cap_vis * 0.90))

comp_pos <- comparaciones %>%
  left_join(base_positions, by = c("variable", "serDR")) %>%
  group_by(variable, serDR) %>%
  arrange(p) %>%
  mutate(idx = row_number(),
         y   = base + (idx - 1) * 0.05 * cap_vis) %>%
  ungroup() %>%
  mutate(variable = factor(variable, levels = orden_var))

print(comparaciones %>% select(variable, serDR, tipo_gen, n1, n2, p, p.signif))

# ── Etiquetas serológicas con n ───────────────────────────────────────────────
counts_serDR <- df_sub %>%
  group_by(serDR) %>%
  summarise(
    episodios = n(),
    pacientes = n_distinct(.data[[id_patient]]),
    .groups   = "drop"
  ) %>%
  mutate(label = paste0(
    "<b>", serDR, "</b><br>",
    "<span style='font-size:10pt'>",
    pacientes, " pac. / ", episodios, " ep.</span>"
  ))

serdr_labels        <- counts_serDR$label
names(serdr_labels) <- counts_serDR$serDR

# ── Colores ───────────────────────────────────────────────────────────────────
colores <- c(
  "gB único"    = "#63B8FF",
  "gB mixto"    = "dodgerblue3",
  "gH único"    = "#DDA0DD",
  "gH mixto"    = "#CD69C9",
  "gB+gH único" = "#3CB371",
  "gB+gH mixto" = "#F08080"
)

# ── Figura ────────────────────────────────────────────────────────────────────
p_global <- ggplot(df_plot_vis,
                   aes(x = x_pos, y = valor, fill = genotipo)) +
  geom_boxplot(
    width         = 0.55,
    color         = "grey25",
    linewidth     = 0.4,
    outlier.shape = 21,
    outlier.size  = 1.4,
    outlier.fill  = "white",
    outlier.color = "grey50",
    outlier.alpha = 0.6
  ) +
  geom_vline(xintercept = c(2.65, 4.95),
             color = "grey80", linewidth = 0.35, linetype = "dashed") +
  facet_grid(variable ~ serDR, scales = "free_y",
             labeller = labeller(serDR    = as_labeller(serdr_labels),
                                 variable = label_parsed)) +
  geom_text(data = comp_pos,
            aes(x = (xmin + xmax) / 2, y = y, label = p.signif),
            inherit.aes = FALSE, size = 3.8, vjust = 0, color = "grey20") +
  geom_segment(data = comp_pos,
               aes(x    = xmin, xend = xmax,
                   y    = y - 0.012 * cap_vis,
                   yend = y - 0.012 * cap_vis),
               inherit.aes = FALSE,
               linewidth = 0.45, color = "grey30") +
  scale_fill_manual(values = colores) +
  scale_x_continuous(
    breaks = unname(x_map),
    labels = names(x_map),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position    = "none",
    panel.border       = element_rect(color = "grey70", fill = NA,
                                      linewidth = 0.4),
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.line          = element_blank(),
    axis.ticks         = element_line(color = "grey60", linewidth = 0.3),
    axis.text.x        = element_text(size = 10, color = "grey25",
                                      angle = 45, hjust = 1),
    axis.text.y        = element_text(size = 12, color = "grey30"),
    strip.background   = element_rect(fill = "grey95", color = "grey70",
                                      linewidth = 0.4),
    strip.text.x       = element_markdown(size = 13, face = "bold",
                                          color = "grey20"),
    strip.text.y       = element_text(size = 11, face = "bold",
                                      color = "grey20"),
    panel.spacing.x    = unit(6,  "pt"),
    panel.spacing.y    = unit(10, "pt"),
    plot.margin        = margin(8, 10, 6, 6)
  )

p_global

# ── Guardar en SVG ────────────────────────────────────────────────────────────
# [guardado antiguo] ggsave(
# [guardado antiguo]   filename      = "~/Desktop/Figuras_tesis/Figura_28_cinetica_serologia.svg",
# [guardado antiguo]   plot          = p_global,
# [guardado antiguo]   device        = svglite::svglite,
# [guardado antiguo]   fix_text_size = FALSE,
# [guardado antiguo]   width         = 16,
# [guardado antiguo]   height        = 22,
# [guardado antiguo]   units         = "cm",
# [guardado antiguo]   bg            = "white"
# [guardado antiguo] )

guardar_fig(p_global, "Figura_28", ancho = 16, alto = 22)

## ---- Tabla 15. Modelo multivariante del aclaramiento en el análisis de genotipos mixtos gB/gH ----
# Análisis univariante de confusores y modelos previos; la Tabla 15 corresponde a modelo_clinico_mix.
# ── 1. Preparar datos ─────────────────────────────────────────────────────────
df_uni <- dnaemias %>%
  left_join(pacientes %>% select(NHC, edad_tx, tipotx, ebase), by = "NHC") %>%
  mutate(
    risk_f = factor(
      ifelse(serDR == "D+/R-", "Alto riesgo", "Riesgo intermedio"),
      levels = c("Riesgo intermedio", "Alto riesgo")
    ),
    Fini              = as.Date(Fini),
    Ffin              = as.Date(Ffin),
    aclaramiento_dias = as.numeric(Ffin - Fini),
    Tdesdetx_meses    = as.numeric(Fini - as.Date(Ftx)) / 30.44,
    sintomatico       = case_when(
      clinica == "VA"  ~ 0L,
      !is.na(clinica)  ~ 1L,
      TRUE             ~ NA_integer_
    ),
    gbghmix_f    = factor(gbghmix,  levels = c(0,1), labels = c("Único","Mixto")),
    ingr_f       = factor(ingr,     levels = c(0,1), labels = c("No","Sí")),
    ing_cmv_f    = factor(ing_cmv,  levels = c(0,1), labels = c("No","Sí")),
    sintom_f     = factor(sintomatico, levels = c(0,1), labels = c("No","Sí")),
    profilaxis_f = factor(antiviral, levels = c(0,1), labels = c("No","Sí")),
    tto_f        = factor(tto,      levels = c(0,1), labels = c("No","Sí")),
    cmvr_f       = factor(CMV_R,    levels = c(0,1), labels = c("No","Sí")),
    tipotx_f     = factor(tipotx,   levels = c("UNI","BI")),
    ebase_f      = factor(ebase)
  ) %>%
  filter(!is.na(aclaramiento_dias))

# ── 2. Correlaciones de Spearman (variables continuas) ────────────────────────
cor_vars <- list(
  "Edad receptor (años)"        = "edad_tx",
  "CV inicial (UI/mL)"          = "cvini",
  "CV pico (UI/mL)"             = "cvpico",
  "Tiempo desde TXP (meses)"     = "Tdesdetx_meses"
)

res_cor <- map_dfr(names(cor_vars), function(lab) {
  v <- cor_vars[[lab]]
  x <- as.numeric(df_uni[[v]])
  if (v %in% c("cvini","cvpico")) x <- log10(x + 1)
  idx  <- !is.na(x) & !is.na(df_uni$aclaramiento_dias)
  x_ok <- x[idx]
  y_ok <- df_uni$aclaramiento_dias[idx]
  if (length(x_ok) < 5) return(tibble(Variable = lab, tipo = "continua",
                                      estimador = NA_real_, p_val = NA_real_,
                                      label_est = "n insuf."))
  test <- tryCatch(
    cor.test(x_ok, y_ok, method = "spearman", exact = FALSE),
    error = function(e) NULL
  )
  if (is.null(test)) return(tibble(Variable = lab, tipo = "continua",
                                   estimador = NA_real_, p_val = NA_real_,
                                   label_est = "error"))
  tibble(
    Variable  = lab,
    tipo      = "continua",
    estimador = round(test$estimate, 3),
    p_val     = test$p.value,
    label_est = sprintf("rho = %+.3f", test$estimate)
  )
})

# ── 3. Comparaciones de grupos (variables categóricas 2 grupos) ───────────────
cat_vars <- list(
  "Genotipos mixtos gB/gH"            = "gbghmix_f",
  "Perfil de riesgo serológico D/R"   = "risk_f",      
  "Hospitalización cualquier causa"   = "ingr_f",
  "Hospitalización por CMV"           = "ing_cmv_f",
  "Enfermedad sintomática"            = "sintom_f",
  "Profilaxis al inicio del episodio" = "profilaxis_f",
  "Tratamiento antiviral"             = "tto_f",
  "Tipo de trasplante (BI vs UNI)"    = "tipotx_f"
)

res_cat <- map_dfr(names(cat_vars), function(lab) {
  v  <- cat_vars[[lab]]
  g  <- df_uni[[v]]
  ac <- df_uni$aclaramiento_dias
  d0 <- ac[g == levels(g)[1] & !is.na(g)]
  d1 <- ac[g == levels(g)[2] & !is.na(g)]
  test <- wilcox.test(d1, d0, exact = FALSE)
  dif  <- median(d1, na.rm = TRUE) - median(d0, na.rm = TRUE)
  tibble(
    Variable  = lab,
    tipo      = "categorica",
    estimador = round(dif, 1),
    p_val     = test$p.value,
    label_est = sprintf("\u0394 = %+.1f d\u00edas", dif)
  )
})

# ── 4. Kruskal-Wallis (variables categóricas >2 grupos) ───────────────────────
kw_eb  <- kruskal.test(aclaramiento_dias ~ ebase_f, data = df_uni)

res_eb <- tibble(
  Variable  = "Enfermedad de base",
  tipo      = "categorica",
  estimador = NA_real_,
  p_val     = kw_eb$p.value,
  label_est = "Kruskal-Wallis"
)

res_eb <- tibble(
  Variable  = "Enfermedad de base",
  tipo      = "categorica",
  estimador = NA_real_,
  p_val     = kw_eb$p.value,
  label_est = "Kruskal-Wallis"
)

# ── 5. Combinar y preparar para figura ────────────────────────────────────────
orden_vars <- c(
  "Edad receptor (años)",
  "Tipo de trasplante (BI vs UNI)",
  "Enfermedad de base",
  "Perfil de riesgo serológico D/R",    # <-- corregido
  "Tiempo desde TXP (meses)",
  "CV inicial (UI/mL)",
  "CV pico (UI/mL)",
  "Genotipos mixtos gB/gH",
  "Enfermedad sintomática",
  "Profilaxis al inicio del episodio",
  "Tratamiento antiviral",
  "Hospitalización cualquier causa",
  "Hospitalización por CMV"
)

grupos_map <- c(
  "Edad receptor (años)"              = "Demográficas",
  "Tipo de trasplante (BI vs UNI)"    = "Demográficas",
  "Enfermedad de base"                = "Demográficas",
  "Perfil de riesgo serológico D/R"   = "Demográficas",   
  "Tiempo desde TXP (meses)"          = "Cinéticas",
  "CV inicial (UI/mL)"                = "Cinéticas",
  "CV pico (UI/mL)"                   = "Cinéticas",
  "Genotipos mixtos gB/gH"            = "Cinéticas",
  "Enfermedad sintomática"            = "Clínicas",
  "Profilaxis al inicio del episodio" = "Clínicas",
  "Tratamiento antiviral"             = "Clínicas",
  "Hospitalización cualquier causa"   = "Clínicas",
  "Hospitalización por CMV"           = "Clínicas"
)

res_all <- bind_rows(res_cor, res_cat, res_eb) %>%
  mutate(
    sig     = p_val < 0.05,
    p_fmt   = case_when(
      p_val < 0.001 ~ "<0,001",
      TRUE          ~ sub("\\.", ",", sprintf("%.3f", p_val))
    ),
    p_label = paste0("p = ", p_fmt),
    grupo   = factor(grupos_map[as.character(Variable)],
                     levels = c("Demográficas", "Cinéticas", "Clínicas")),
    Variable = factor(Variable, levels = rev(orden_vars))
  )

# ── 6. Figura ─────────────────────────────────────────────────────────────────
figura_univariante <- ggplot(res_all, aes(x = p_val, y = Variable)) +
  geom_vline(
    xintercept = 0.05,
    color = "firebrick3",
    linetype = "dashed",
    linewidth = 0.5,
    alpha = 0.5
  ) +
  geom_vline(xintercept = 0.20, color = "grey65",
             linetype = "dashed", linewidth = 0.4) +
  geom_point(aes(color = sig, size = sig)) +
  geom_text(aes(label = p_label, color = sig),
            hjust = -0.15, size = 4) +
  geom_text(aes(x = 1.1, label = label_est),
            hjust = 1, size = 4, color = "grey40") +
  facet_grid(grupo ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(
    name   = "valor p",
    limits = c(0, 1.15),
    breaks = c(0.05, 0.20, 0.50, 1.00),
    labels = c("0,05", "0,20", "0,50", "1,00")
  ) +
  scale_color_manual(values = c("TRUE" = "firebrick3", "FALSE" = "grey50")) +
  scale_size_manual(values  = c("TRUE" = 3.5,          "FALSE" = 2.5)) +
  labs(y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "none",
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.3),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
    axis.text.y        = element_text(size = 12, color = "grey20"),
    axis.text.x        = element_text(size = 12, color = "grey40"),
    axis.title.x       = element_text(size = 11, color = "grey20"),
    strip.text.y       = element_text(angle = -90, size = 12,
                                      face = "bold", color = "grey20",
                                      margin = margin(l = 6, r = 6)),
    strip.background   = element_rect(fill = "grey93", color = NA),
    panel.spacing.y    = unit(0.6, "lines"),
    plot.margin        = margin(8, 60, 8, 8)
  )

figura_univariante

# Variables candidatas (p < 0,20 en univariante)
df_multi <- df_uni %>%
  filter(
    !is.na(aclaramiento_dias),
    !is.na(gbghmix_f),
    !is.na(risk_f),
    !is.na(cvini),
    !is.na(cvpico),
    !is.na(Tdesdetx_meses),
    !is.na(sintom_f),
    !is.na(profilaxis_f),
    !is.na(ing_cmv_f)
  )

modelo_completo <- lm(
  aclaramiento_dias ~ gbghmix_f + risk_f + log10(cvini + 1) + log10(cvpico + 1) +
    Tdesdetx_meses + sintom_f + profilaxis_f + ing_cmv_f,
  data = df_multi
)

summary(modelo_completo)

# ── Selección backward ────────────────────────────────────────────────────────
modelo_backward <- step(modelo_completo, direction = "backward", trace = TRUE)
summary(modelo_backward)

# ── Intervalos de confianza y VIF ─────────────────────────────────────────────

confint(modelo_backward, level = 0.95)
vif(modelo_backward)


confint(modelo_backward, level = 0.95)

# ── 1. Datos completos para las variables del modelo clínico ──────────────────
df_clin_mix <- df_uni %>%
  select(NHC, aclaramiento_dias,
         gbghmix_f, risk_f, cvini, cvpico, profilaxis_f) %>%
  filter(complete.cases(.))

cat(sprintf("\n── N modelo clínico mixtos: %d episodios / %d pacientes ──\n",
            nrow(df_clin_mix), n_distinct(df_clin_mix$NHC)))

# ── 2. Modelo lineal múltiple ─────────────────────────────────────────────────
modelo_clinico_mix <- lm(
  aclaramiento_dias ~ gbghmix_f + risk_f + log10(cvini + 1) + log10(cvpico + 1) +
    profilaxis_f,
  data = df_clin_mix
)

cat("\n══════════════════════════════════════════════\n")
cat("  MODELO CLÍNICO MIXTOS gB+gH (variables a priori)\n")
cat("══════════════════════════════════════════════\n")
print(summary(modelo_clinico_mix))

cat("\n── IC 95% ──\n")
print(confint(modelo_clinico_mix))

cat("\n── VIF (multicolinealidad) ──\n")
print(car::vif(modelo_clinico_mix))

# ── 3. Efecto ajustado de los mixtos en una línea ─────────────────────────────
# El nombre del coeficiente depende del nivel: comprobamos cuál es
coef_name <- grep("gbghmix_f", names(coef(modelo_clinico_mix)), value = TRUE)
cat("\n── Coeficiente de interés:", coef_name, "──\n")

coef_mix <- coef(modelo_clinico_mix)[coef_name]
ci_mix   <- confint(modelo_clinico_mix)[coef_name, ]
p_mix    <- summary(modelo_clinico_mix)$coefficients[coef_name, "Pr(>|t|)"]

cat(sprintf("\n── Efecto ajustado mixtos gB+gH: %+.1f días (IC 95%%: %+.1f a %+.1f), p = %.4f ──\n",
            coef_mix, ci_mix[1], ci_mix[2], p_mix))

# ── 4. Análisis de sensibilidad: modelo mixto (paciente aleatorio) ────────────
modelo_mixto_mix <- lmer(
  aclaramiento_dias ~ gbghmix_f + risk_f + log10(cvini + 1) + log10(cvpico + 1) +
    profilaxis_f + (1 | NHC),
  data = df_clin_mix
)

cat("\n══════════════════════════════════════════════\n")
cat("  MODELO MIXTO: paciente como efecto aleatorio\n")
cat("══════════════════════════════════════════════\n")
print(summary(modelo_mixto_mix))

cat("\n── IC 95% (Wald) del efecto mixtos ──\n")
print(confint(modelo_mixto_mix, method = "Wald")[coef_name, ])

vc  <- as.data.frame(VarCorr(modelo_mixto_mix))
icc <- vc$vcov[1] / sum(vc$vcov)
cat(sprintf("\n── ICC (varianza atribuible al paciente): %.3f ──\n", icc))

guardar_modelo(modelo_clinico_mix, "Tabla_15_modelo")

## ---- Figura 29. Aclaramiento viral según presencia de gB3, en genotipos únicos o mixtos de gB ----
# ── 1. Crear variables de clasificación ───────────────────────────────────────
df_gb3 <- df_uni %>%
  mutate(
    # gB3 presente en el episodio (independientemente de si es único o mixto)
    gb3_presente = as.integer(gb3 == 1),
    
    # Cuatro grupos: único/mixto x gB3 sí/no
    grupo_gb3 = case_when(
      gbmix == 0 & gb3 == 0 ~ "Único sin gB3",
      gbmix == 0 & gb3 == 1 ~ "Único con gB3",
      gbmix == 1 & gb3 == 0 ~ "Mixto sin gB3",
      gbmix == 1 & gb3 == 1 ~ "Mixto con gB3",
      TRUE ~ NA_character_
    ),
    grupo_gb3 = factor(grupo_gb3,
                       levels = c("Único sin gB3", "Único con gB3",
                                  "Mixto sin gB3", "Mixto con gB3"))
  ) %>%
  filter(!is.na(aclaramiento_dias), !is.na(grupo_gb3))

# ── 2. Análisis A: entre episodios mixtos, gB3 vs no gB3 ─────────────────────
df_mixtos <- df_gb3 %>% filter(gbmix == 1)

cat("\n── Análisis A: episodios mixtos, gB3 vs no gB3 ──\n")
df_mixtos %>%
  group_by(gb3_presente) %>%
  summarise(
    n        = n(),
    mediana  = median(aclaramiento_dias, na.rm = TRUE),
    q1       = quantile(aclaramiento_dias, 0.25, na.rm = TRUE),
    q3       = quantile(aclaramiento_dias, 0.75, na.rm = TRUE)
  ) %>% print()

wilcox.test(aclaramiento_dias ~ gb3_presente, data = df_mixtos, exact = FALSE)

# ── 3. Análisis B: cuatro grupos ──────────────────────────────────────────────
cat("\n── Análisis B: cuatro grupos ──\n")
df_gb3 %>%
  group_by(grupo_gb3) %>%
  summarise(
    n        = n(),
    mediana  = median(aclaramiento_dias, na.rm = TRUE),
    q1       = quantile(aclaramiento_dias, 0.25, na.rm = TRUE),
    q3       = quantile(aclaramiento_dias, 0.75, na.rm = TRUE)
  ) %>% print()

kruskal.test(aclaramiento_dias ~ grupo_gb3, data = df_gb3)

# Post-hoc si Kruskal significativo
dunn.test(df_gb3$aclaramiento_dias, df_gb3$grupo_gb3,
          method = "bonferroni", kw = TRUE)

etiquetas_n <- c(
  "Único sin gB3" = paste0("Único sin gB3\n(n=", sum(df_gb3$grupo_gb3 == "Único sin gB3", na.rm=TRUE), ")"),
  "Único con gB3" = paste0("Único con gB3\n(n=", sum(df_gb3$grupo_gb3 == "Único con gB3", na.rm=TRUE), ")"),
  "Mixto sin gB3" = paste0("Mixto sin gB3\n(n=", sum(df_gb3$grupo_gb3 == "Mixto sin gB3", na.rm=TRUE), ")"),
  "Mixto con gB3" = paste0("Mixto con gB3\n(n=", sum(df_gb3$grupo_gb3 == "Mixto con gB3", na.rm=TRUE), ")")
)

figura_gb3 <- ggplot(
  df_gb3,
  aes(
    x = grupo_gb3,
    y = aclaramiento_dias,
    fill = ifelse(grepl("con gB3", grupo_gb3), "con gB3", "sin gB3")
  )
) +
  geom_boxplot(
    alpha         = 0.75,
    outlier.shape = NA,      # oculta outliers nativos del boxplot
    width         = 0.55,
    color         = "grey25",
    linewidth     = 0.4
  ) +
  # Puntos individuales solo hasta 150
  geom_jitter(
    data        = df_gb3 %>% filter(aclaramiento_dias <= 150),
    aes(x = grupo_gb3, y = aclaramiento_dias,
        fill = ifelse(grepl("con gB3", grupo_gb3), "con gB3", "sin gB3")),
    shape       = 21,
    size        = 1.5,
    color       = "grey40",
    alpha       = 0.6,
    width       = 0.15,
    height      = 0,
    inherit.aes = FALSE
  ) + 
  geom_vline(
    xintercept = 2.5,
    linetype   = "dashed",
    color      = "grey70",
    linewidth  = 0.4
  ) +
  # Barras de significación
  annotate("segment", x = 1,   xend = 2, y = 170, yend = 170, linewidth = 0.4) +
  annotate("segment", x = 1,   xend = 1, y = 165, yend = 170, linewidth = 0.4) +
  annotate("segment", x = 2,   xend = 2, y = 165, yend = 170, linewidth = 0.4) +
  annotate("text",    x = 1.5, y = 174,  label = "***",  size = 5) +
  annotate("segment", x = 1,   xend = 4, y = 185, yend = 185, linewidth = 0.4) +
  annotate("segment", x = 1,   xend = 1, y = 180, yend = 185, linewidth = 0.4) +
  annotate("segment", x = 4,   xend = 4, y = 180, yend = 185, linewidth = 0.4) +
  annotate("text",    x = 2.5, y = 189,  label = "****", size = 5) +
  # Etiquetas de bloque
  annotate("text", x = 1.5, y = 208, label = "Genotipo único",
           size = 4, color = "grey30", fontface = "bold") +
  annotate("text", x = 3.5, y = 208, label = "Genotipo mixto",
           size = 4, color = "grey30", fontface = "bold") +
  scale_fill_manual(
    values = c(
      "sin gB3" = "#CDC0B0",
      "con gB3" = "#4F94CD"
    )
  ) +
  scale_x_discrete(labels = etiquetas_n) +
  scale_y_continuous(
    name   = "Tiempo hasta aclaramiento viral (días)",
    breaks = seq(0, 150, by = 50)
  ) +
  coord_cartesian(ylim = c(0, 215), clip = "off") +
  labs(x = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35),
    axis.line.y        = element_line(color = "grey60", linewidth = 0.4),
    axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
    axis.text.x        = element_text(size = 11, color = "grey20"),
    axis.text.y        = element_text(size = 11, color = "grey40"),
    axis.title.y       = element_text(size = 11, color = "grey20"),
    plot.margin        = margin(8, 8, 8, 8)
  )

figura_gb3

guardar_fig(figura_gb3, "Figura_29", ancho = 14, alto = 11)


# ---- 9.6. Diversidad y evolución de los genotipos ----

## ---- Preparación: muestras genotipadas por episodio ----
# ── 1. Cargar muestras ────────────────────────────────────────────────────────
muestras <- read_excel(f_bd,
                       sheet = "genotipado",
                       col_types = c("numeric", "numeric", "date", "numeric",
                                     "numeric", "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "text"))

# ── 1b. Filtrar solo muestras con genotipado completo ─────────────────────────
muestras <- muestras %>%
  mutate(
    gb = rowSums(across(c(GB1, GB2, GB3, GB4)), na.rm = TRUE),
    gh = rowSums(across(c(GH1, GH2)), na.rm = TRUE)
  ) %>%
  filter(gb >= 1 & gh >= 1)

cat(sprintf("Muestras con genotipado completo: %d\n", nrow(muestras)))
cat(sprintf("Pacientes: %d\n", n_distinct(muestras$NHC)))

# ── 2. Asignar cada muestra a su episodio ─────────────────────────────────────
muestras_ep <- muestras %>%
  inner_join(
    dnaemias %>%
      select(NHC, Fini, Ffin, gbghmix) %>%
      mutate(Fini = as.Date(Fini), Ffin = as.Date(Ffin)),
    by = "NHC",
    relationship = "many-to-many"
  ) %>%
  filter(
    Freg >= Fini &
      (is.na(Ffin) | Freg <= Ffin)
  )

cat(sprintf("Muestras asignadas: %d\n", nrow(muestras_ep)))
cat(sprintf("Episodios representados: %d\n",
            n_distinct(muestras_ep %>% select(NHC, Fini))))

# ── 3. Episodios con ≥2 muestras ──────────────────────────────────────────────
ep_n_muestras <- muestras_ep %>%
  group_by(NHC, Fini, Ffin, gbghmix) %>%
  summarise(n_muestras = n(), .groups = "drop")

n_total_ep <- 234
n_intra_ep <- ep_n_muestras %>% filter(n_muestras >= 2) %>% nrow()

cat(sprintf("Episodios con ≥2 muestras: %d (%.1f%%)\n",
            n_intra_ep, 100 * n_intra_ep / n_total_ep))

ep_intra <- ep_n_muestras %>% filter(n_muestras >= 2)

cat(sprintf("Mediana muestras por episodio: %.0f (RIQ %.0f-%.0f)\n",
            median(ep_intra$n_muestras),
            quantile(ep_intra$n_muestras, 0.25),
            quantile(ep_intra$n_muestras, 0.75)))

# ── 4. Clasificar episodios únicos vs mixtos ──────────────────────────────────
cat("\n── Clasificación episodios con ≥2 muestras ──\n")
ep_intra %>%
  mutate(clase = ifelse(gbghmix == 0, "Único", "Mixto")) %>%
  count(clase) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()

# ── 5. Variabilidad intraepisodio ─────────────────────────────────────────────
muestras_intra <- muestras_ep %>%
  semi_join(ep_intra, by = c("NHC", "Fini")) %>%
  mutate(patron = paste(GB1, GB2, GB3, GB4, GH1, GH2, sep = "-"))

variabilidad_ep <- muestras_intra %>%
  group_by(NHC, Fini, gbghmix) %>%
  summarise(
    n_muestras = n(),
    n_patrones = n_distinct(patron),
    .groups    = "drop"
  ) %>%
  mutate(
    clase    = ifelse(gbghmix == 0, "Único", "Mixto"),
    variable = ifelse(n_patrones > 1, "Variable", "Estable")
  )

cat("\n── Variabilidad por clase ──\n")
variabilidad_ep %>%
  group_by(clase, variable) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(clase) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()

# ── 6. Episodios mixtos variables ─────────────────────────────────────────────
mixtos_var <- variabilidad_ep %>%
  filter(clase == "Mixto", variable == "Variable")

cat(sprintf("\n── Episodios mixtos variables: %d episodios, %d pacientes ──\n",
            nrow(mixtos_var), n_distinct(mixtos_var$NHC)))

mixtos_var %>%
  left_join(pacientes %>% select(NHC, serDR), by = "NHC") %>%
  distinct(NHC, serDR) %>%
  count(serDR) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()

## ---- Datos para la Figura 30 (variabilidad intraepisodio; la figura de la tesis es una infografía) ----
df_fig <- variabilidad_ep %>%
  group_by(clase, variable) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(clase) %>%
  mutate(
    total = sum(n),
    pct   = 100 * n / total,
    label = paste0(n, " (", sub("\\.", ",", round(pct, 1)), "%)")
  ) %>%
  ungroup() %>%
  mutate(
    clase    = factor(clase, levels = c("Único", "Mixto"),
                      labels = c("Genotipo único gB/gH", "Genotipo mixto gB/gH")),
    variable = factor(variable, levels = c("Variable", "Estable"))
  )

figura_R21 <- ggplot(df_fig, aes(x = clase, y = pct, fill = variable)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4, fontface = "bold", color = "white"
  ) +
  scale_fill_manual(
    values = c("Estable"  = "#2E7D32",
               "Variable" = "#E65100"),
    labels = c("Estable"  = "Patrón dominante",
               "Variable" = "Genotipo fluctuante"),
    name   = ""
  ) +
  scale_y_continuous(
    name   = "% de episodios",
    limits = c(0, 105),
    breaks = seq(0, 100, 25),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "top",
    legend.text        = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
    axis.line.y        = element_line(color = "grey60", linewidth = 0.4),
    axis.text.x        = element_text(size = 12, color = "grey20"),
    axis.text.y        = element_text(size = 11, color = "grey40"),
    axis.title.y       = element_text(size = 11, color = "grey20"),
    plot.margin        = margin(8, 8, 8, 8)
  )

figura_R21

## ---- Tabla 16. Cambios intraepisodio en episodios variables con genotipos mixtos gB/gH ----
muestras_mv <- muestras_ep %>%
  semi_join(mixtos_var, by = c("NHC", "Fini")) %>%
  mutate(patron = paste(GB1, GB2, GB3, GB4, GH1, GH2, sep = "-")) %>%
  arrange(NHC, Fini, Freg)

# Número de cambios por episodio
cambios_ep <- muestras_mv %>%
  group_by(NHC, Fini) %>%
  arrange(Freg, .by_group = TRUE) %>%
  mutate(cambio = patron != lag(patron)) %>%
  summarise(
    n_muestras = n(),
    n_cambios  = sum(cambio, na.rm = TRUE),
    .groups    = "drop"
  )

cat("── Número de cambios por episodio ──\n")
cat(sprintf("Rango: %d – %d\n", min(cambios_ep$n_cambios), max(cambios_ep$n_cambios)))
cat(sprintf("Mediana: %.0f (RIQ %.0f–%.0f)\n",
            median(cambios_ep$n_cambios),
            quantile(cambios_ep$n_cambios, 0.25),
            quantile(cambios_ep$n_cambios, 0.75)))

# Estabilidad por genotipo
genotipos <- c("GB1","GB2","GB3","GB4","GH1","GH2")
nombres   <- c("gB1","gB2","gB3","gB4","gH1","gH2")

estabilidad_gt <- muestras_mv %>%
  group_by(NHC, Fini) %>%
  summarise(
    across(all_of(genotipos), ~ {
      presente <- sum(. == 1, na.rm = TRUE)
      n_tot    <- n()
      if (presente == 0) "ausente"
      else if (presente == n_tot) "estable"
      else "variable"
    }),
    .groups = "drop"
  )

cat("\n── Frecuencia genotipos estables ──\n")
for (i in seq_along(genotipos)) {
  g <- genotipos[i]
  n_est <- sum(estabilidad_gt[[g]] == "estable")
  cat(sprintf("%s estable: %d/35 (%.1f%%)\n", nombres[i], n_est, 100*n_est/35))
}

cat("\n── Frecuencia genotipos variables ──\n")
for (i in seq_along(genotipos)) {
  g <- genotipos[i]
  n_var <- sum(estabilidad_gt[[g]] == "variable")
  cat(sprintf("%s variable: %d/35 (%.1f%%)\n", nombres[i], n_var, 100*n_var/35))
}

estabilidad_gt <- estabilidad_gt %>%
  mutate(tiene_estable = rowSums(across(all_of(genotipos), ~ . == "estable")) > 0)

cat(sprintf("\n── Episodios con ≥1 genotipo estable: %d/35 (%.1f%%)\n",
            sum(estabilidad_gt$tiene_estable),
            100 * mean(estabilidad_gt$tiene_estable)))

cat(sprintf("── Episodios con cambio completo: %d\n",
            sum(!estabilidad_gt$tiene_estable)))

tabla_mv <- estabilidad_gt %>%
  select(-any_of("n_muestras")) %>%
  left_join(
    cambios_ep %>% select(NHC, Fini, n_cambios, n_muestras),
    by = c("NHC", "Fini")
  ) %>%
  left_join(pacientes %>% select(NHC, serDR), by = "NHC") %>%
  mutate(
    estables = pmap_chr(
      list(GB1, GB2, GB3, GB4, GH1, GH2),
      function(gb1, gb2, gb3, gb4, gh1, gh2) {
        gt <- c(
          if (gb1 == "estable") "gB1",
          if (gb2 == "estable") "gB2",
          if (gb3 == "estable") "gB3",
          if (gb4 == "estable") "gB4",
          if (gh1 == "estable") "gH1",
          if (gh2 == "estable") "gH2"
        )
        if (length(gt) == 0) "—" else paste(gt, collapse = "/")
      }
    ),
    variables = pmap_chr(
      list(GB1, GB2, GB3, GB4, GH1, GH2),
      function(gb1, gb2, gb3, gb4, gh1, gh2) {
        gt <- c(
          if (gb1 == "variable") "gB1",
          if (gb2 == "variable") "gB2",
          if (gb3 == "variable") "gB3",
          if (gb4 == "variable") "gB4",
          if (gh1 == "variable") "gH1",
          if (gh2 == "variable") "gH2"
        )
        if (length(gt) == 0) "—" else paste(gt, collapse = "/")
      }
    )
  ) %>%
  arrange(NHC, Fini) %>%
  mutate(
    pac_num   = as.integer(factor(NHC, levels = unique(NHC))),
    ep_num    = row_number(),
    pac_label = paste0(pac_num, " ", serDR)
  ) %>%
  select(ep_num, pac_label, n_muestras, estables, variables, n_cambios)

tabla_ft <- tabla_mv %>%
  rename(
    "Episodio"              = ep_num,
    "Paciente"              = pac_label,
    "N muestras"            = n_muestras,
    "Patrón dominante"      = estables,
    "Genotipos fluctuantes" = variables,
    "N cambios"             = n_cambios
  )

ft <- flextable(tabla_ft) %>%
  width(j = "Episodio",              width = 0.7) %>%
  width(j = "Paciente",              width = 1.6) %>%
  width(j = "N muestras",            width = 0.9) %>%
  width(j = "Patrón dominante",      width = 2.0) %>%
  width(j = "Genotipos fluctuantes", width = 2.0) %>%
  width(j = "N cambios",             width = 0.9) %>%
  align(align = "center", part = "all") %>%
  align(j = c("Patrón dominante", "Genotipos fluctuantes"),
        align = "left", part = "body") %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header") %>%
  bg(bg = "white", part = "all") %>%
  bold(i = 1, part = "body") %>%
  color(i = 1, color = "firebrick3", part = "body") %>%
  border_remove() %>%
  hline_top(part = "header",
            border = fp_border(color = "black", width = 1.2)) %>%
  hline(part = "header",
        border = fp_border(color = "black", width = 0.8)) %>%
  hline_bottom(part = "body",
               border = fp_border(color = "black", width = 1.2)) %>%
  hline(
    i = which(diff(as.integer(factor(tabla_ft$Paciente))) != 0),
    border = fp_border(color = "grey50", width = 0.5),
    part = "body"
  ) %>%
  padding(padding.top = 2, padding.bottom = 2,
          padding.left = 4, padding.right = 4, part = "all") %>%
  add_footer_lines(values = c(
    "Patrón dominante: genotipos presentes en todas las muestras del episodio.",
    "Genotipos fluctuantes: genotipos con aparición o desaparición entre muestras consecutivas.",
    "El episodio en rojo corresponde al único caso con cambio completo del patrón dominante.",
    "N cambios: número de transiciones de patrón genotípico entre muestras consecutivas."
  )) %>%
  fontsize(size = 8, part = "footer") %>%
  color(color = "grey30", part = "footer") %>%
  set_caption(caption = "Tabla R8. Episodios con variabilidad intraepisodio en genotipos mixtos gB/gH.") %>%
  set_table_properties(layout = "autofit")

ft

guardar_tabla(ft, "Tabla_16")

## ---- Figura 31. Distancia al percentil 75 del Ct en genotipos estables y fluctuantes ----
ct_cols  <- c("GB1CT","GB2CT","GB3CT","GB4CT","GH1CT","GH2CT")
gt_cols  <- c("GB1","GB2","GB3","GB4","GH1","GH2")
gt_names <- c("gB1","gB2","gB3","gB4","gH1","gH2")

gt_largo <- muestras_mv %>%
  select(NHC, Fini, Freg, all_of(gt_cols)) %>%
  pivot_longer(cols = all_of(gt_cols),
               names_to = "genotipo_var", values_to = "presente") %>%
  filter(presente == 1) %>%
  mutate(genotipo = gt_names[match(genotipo_var, gt_cols)])

ct_largo <- muestras_mv %>%
  select(NHC, Fini, Freg, all_of(ct_cols)) %>%
  pivot_longer(cols = all_of(ct_cols),
               names_to = "ct_var", values_to = "Ct") %>%
  mutate(genotipo_var = gt_cols[match(ct_var, ct_cols)])

det_largo <- gt_largo %>%
  left_join(ct_largo %>% select(NHC, Fini, Freg, genotipo_var, Ct),
            by = c("NHC", "Fini", "Freg", "genotipo_var")) %>%
  filter(!is.na(Ct))

estabilidad_largo <- estabilidad_gt %>%
  select(NHC, Fini, all_of(gt_cols)) %>%
  pivot_longer(cols = all_of(gt_cols),
               names_to = "genotipo_var", values_to = "estabilidad") %>%
  filter(estabilidad != "ausente")

det_largo <- det_largo %>%
  left_join(estabilidad_largo, by = c("NHC", "Fini", "genotipo_var")) %>%
  filter(!is.na(estabilidad))

# P75 por genotipo sobre todas las muestras
todas_ct <- muestras %>%
  select(NHC, Freg, all_of(gt_cols), all_of(ct_cols)) %>%
  pivot_longer(cols = all_of(gt_cols),
               names_to = "genotipo_var", values_to = "presente") %>%
  filter(presente == 1) %>%
  left_join(
    muestras %>%
      select(NHC, Freg, all_of(ct_cols)) %>%
      pivot_longer(cols = all_of(ct_cols),
                   names_to = "ct_var", values_to = "Ct") %>%
      mutate(genotipo_var = gt_cols[match(ct_var, ct_cols)]) %>%
      select(NHC, Freg, genotipo_var, Ct),
    by = c("NHC", "Freg", "genotipo_var")
  ) %>%
  filter(!is.na(Ct))

p75_genotipo <- todas_ct %>%
  group_by(genotipo_var) %>%
  summarise(P75 = quantile(Ct, 0.75, na.rm = TRUE), .groups = "drop")

det_largo <- det_largo %>%
  left_join(p75_genotipo, by = "genotipo_var") %>%
  mutate(dist_p75 = P75 - Ct)

cat("\n── Distancia al P75: estable vs fluctuante ──\n")
det_largo %>%
  group_by(estabilidad) %>%
  summarise(
    n             = n(),
    mediana       = median(dist_p75, na.rm = TRUE),
    q1            = quantile(dist_p75, 0.25, na.rm = TRUE),
    q3            = quantile(dist_p75, 0.75, na.rm = TRUE),
    pct_sobre_p75 = 100 * mean(dist_p75 < 0, na.rm = TRUE)
  ) %>% print()

wilcox.test(dist_p75 ~ estabilidad, data = det_largo, exact = FALSE)

etiq <- det_largo %>%
  group_by(estabilidad) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = paste0(ifelse(estabilidad == "estable",
                               "Genotipos estables", "Genotipos fluctuantes"),
                        "\n(n=", n, ")"))

det_fig <- det_largo %>%
  mutate(
    grupo = ifelse(estabilidad == "estable",
                   "Patrón estable", "Genotipos fluctuantes"),
    grupo = factor(grupo, levels = c("Patrón estable", "Genotipos fluctuantes"))
  )

etiq_x <- etiq %>%
  mutate(grupo = ifelse(estabilidad == "estable",
                        "Patrón estable", "Genotipos fluctuantes"),
         grupo = factor(grupo, levels = c("Patrón estable", "Genotipos fluctuantes")))

figura_R22 <- ggplot(det_fig, aes(x = grupo, y = dist_p75, fill = grupo)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "firebrick3", linewidth = 0.6) +
  annotate("text", x = 2.42, y = 0.8, label = "P75",
           size = 5, color = "firebrick3", fontface = "bold") +
  geom_jitter(aes(color = grupo), width = 0.15, size = 0.9,
              alpha = 0.35, show.legend = FALSE) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA, width = 0.5,
               color = "grey15", linewidth = 0.7) +
  # Mediana resaltada en blanco
  stat_summary(fun = median, geom = "crossbar",
               width = 0.5, fatten = 0,
               color = "white", linewidth = 0.7) +
  annotate("segment", x = 1, xend = 2, y = 9,   yend = 9,
           color = "black", linewidth = 0.5) +
  annotate("segment", x = 1, xend = 1, y = 8.5, yend = 9,
           color = "black", linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2, y = 8.5, yend = 9,
           color = "black", linewidth = 0.5) +
  annotate("text", x = 1.5, y = 9.7, label = "***",
           size = 7, color = "black", fontface = "bold") +
  scale_x_discrete(labels = setNames(etiq_x$label, etiq_x$grupo)) +
  scale_fill_manual(
    values = c("Patrón estable"        = "#2E7D32",
               "Genotipos fluctuantes" = "#E65100")
  ) +
  scale_color_manual(
    values = c("Patrón estable"        = "#2E7D32",
               "Genotipos fluctuantes" = "#E65100")
  ) +
  scale_y_continuous(
    name   = "Distancia al P75 del genotipo (ciclos Ct)",
    breaks = seq(-6, 10, by = 2)
  ) +
  coord_cartesian(ylim = c(-7, 11)) +
  labs(x = NULL) +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "none",
    axis.line       = element_line(color = "black", linewidth = 0.6),
    axis.ticks      = element_line(color = "black", linewidth = 0.5),
    axis.text.x     = element_text(size = 14, color = "black"),
    axis.text.y     = element_text(size = 13, color = "black"),
    axis.title.y    = element_text(size = 14, color = "black", face = "bold"),
    plot.margin     = margin(8, 8, 8, 8)
  )

figura_R22

guardar_fig(figura_R22, "Figura_31", ancho = 12, alto = 13)

## ---- Preparación: clasificación interepisodio ----
# ── 1. Pacientes con >1 episodio genotipado ───────────────────────────────────
ep_por_paciente <- dnaemias %>%
  filter(!is.na(gbghmix)) %>%
  mutate(Fini = as.Date(Fini)) %>%
  group_by(NHC) %>%
  summarise(n_episodios = n(), .groups = "drop") %>%
  filter(n_episodios > 1)

cat(sprintf("Pacientes con >1 episodio: %d\n", nrow(ep_por_paciente)))
cat(sprintf("Total episodios: %d\n", sum(ep_por_paciente$n_episodios)))

cat("\n── Distribución por número de episodios ──\n")
ep_por_paciente %>%
  count(n_episodios) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()

# ── 2. Episodios de pacientes con >1 episodio ─────────────────────────────────
ep_multi <- dnaemias %>%
  filter(!is.na(gbghmix)) %>%
  mutate(Fini = as.Date(Fini)) %>%
  semi_join(ep_por_paciente, by = "NHC")

cat("\n── Clasificación único/mixto ──\n")
ep_multi %>%
  mutate(clase = ifelse(gbghmix == 0, "Único", "Mixto")) %>%
  count(clase) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()

# ── 3. Clasificar pacientes según comportamiento interepisodio ────────────────
patron_paciente <- ep_multi %>%
  group_by(NHC) %>%
  summarise(
    n_episodios = n(),
    n_unicos    = sum(gbghmix == 0),
    n_mixtos    = sum(gbghmix == 1),
    .groups     = "drop"
  ) %>%
  mutate(
    patron = case_when(
      n_unicos > 0 & n_mixtos == 0 ~ "Siempre único",
      n_unicos == 0 & n_mixtos > 0 ~ "Siempre mixto",
      TRUE                          ~ "Mixtos y únicos"
    )
  )

cat("\n── Clasificación pacientes ──\n")
patron_paciente %>%
  count(patron) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()

cat("\n── Episodios por grupo ──\n")
patron_paciente %>%
  group_by(patron) %>%
  summarise(n_episodios = sum(n_episodios), .groups = "drop") %>%
  mutate(pct = 100 * n_episodios / sum(n_episodios)) %>%
  print()

# ── 4. Perfil serológico por grupo ────────────────────────────────────────────
cat("\n── Perfil serológico por grupo ──\n")
patron_paciente %>%
  left_join(pacientes %>% select(NHC, serDR), by = "NHC") %>%
  group_by(patron, serDR) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(patron) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()

# ── 5. Tiempo entre episodios por grupo ───────────────────────────────────────
cat("\n── Tiempo entre episodios (Ffin → Fini siguiente, días) ──\n")
ep_multi %>%
  mutate(Ffin = as.Date(Ffin)) %>%
  arrange(NHC, Fini) %>%
  group_by(NHC) %>%
  mutate(tiempo_entre = as.numeric(Fini - lag(Ffin))) %>%
  filter(!is.na(tiempo_entre), tiempo_entre >= 0) %>%
  left_join(patron_paciente %>% select(NHC, patron), by = "NHC") %>%
  group_by(patron) %>%
  summarise(
    mediana = median(tiempo_entre, na.rm = TRUE),
    q1      = quantile(tiempo_entre, 0.25, na.rm = TRUE),
    q3      = quantile(tiempo_entre, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

# ── 6. Análisis genotipos estables/variables interepisodio ────────────────────
pacientes_cambio <- patron_paciente %>%
  filter(patron == "Mixtos y únicos") %>%
  pull(NHC)

ep_cambio <- ep_multi %>%
  filter(NHC %in% pacientes_cambio) %>%
  arrange(NHC, Fini)

genotipos <- c("gb1","gb2","gb3","gb4","gh1","gh2")
nombres   <- c("gB1","gB2","gB3","gB4","gH1","gH2")

estab_inter <- ep_cambio %>%
  group_by(NHC) %>%
  summarise(
    n_ep = n(),
    across(all_of(genotipos), ~ {
      presente <- sum(. == 1, na.rm = TRUE)
      n_tot    <- n()
      if (presente == 0) "ausente"
      else if (presente == n_tot) "estable"
      else "variable"
    }),
    .groups = "drop"
  ) %>%
  mutate(
    tiene_estable = rowSums(across(all_of(genotipos), ~ . == "estable")) > 0
  )

cat(sprintf("\n── Pacientes con ≥1 genotipo estable: %d/21 (%.1f%%)\n",
            sum(estab_inter$tiene_estable),
            100 * mean(estab_inter$tiene_estable)))

cat(sprintf("── Pacientes sin ningún genotipo estable: %d\n",
            sum(!estab_inter$tiene_estable)))

cat("\n── Frecuencia genotipos estables ──\n")
for (i in seq_along(genotipos)) {
  n_est <- sum(estab_inter[[genotipos[i]]] == "estable")
  cat(sprintf("%s estable: %d/21 (%.1f%%)\n", nombres[i], n_est, 100*n_est/21))
}

cat("\n── Frecuencia genotipos variables ──\n")
for (i in seq_along(genotipos)) {
  n_var <- sum(estab_inter[[genotipos[i]]] == "variable")
  cat(sprintf("%s variable: %d/21 (%.1f%%)\n", nombres[i], n_var, 100*n_var/21))
}

# ── 7. Número de cambios entre episodios consecutivos ─────────────────────────
cambios_inter <- ep_cambio %>%
  arrange(NHC, Fini) %>%
  group_by(NHC) %>%
  mutate(
    patron_ep = paste(gb1, gb2, gb3, gb4, gh1, gh2, sep = "-"),
    cambio    = patron_ep != lag(patron_ep)
  ) %>%
  summarise(
    n_cambios = sum(cambio, na.rm = TRUE),
    .groups   = "drop"
  )

## ---- Tabla 17. Patrón de cambio de genotipos entre episodios ----
# ── 8. Construir tabla ────────────────────────────────────────────────────────
tabla_inter <- estab_inter %>%
  left_join(cambios_inter, by = "NHC") %>%
  left_join(pacientes %>% select(NHC, serDR), by = "NHC") %>%
  mutate(
    estables = pmap_chr(
      list(gb1, gb2, gb3, gb4, gh1, gh2),
      function(gb1, gb2, gb3, gb4, gh1, gh2) {
        gt <- c(
          if (gb1 == "estable") "gB1",
          if (gb2 == "estable") "gB2",
          if (gb3 == "estable") "gB3",
          if (gb4 == "estable") "gB4",
          if (gh1 == "estable") "gH1",
          if (gh2 == "estable") "gH2"
        )
        if (length(gt) == 0) "—" else paste(gt, collapse = "/")
      }
    ),
    variables = pmap_chr(
      list(gb1, gb2, gb3, gb4, gh1, gh2),
      function(gb1, gb2, gb3, gb4, gh1, gh2) {
        gt <- c(
          if (gb1 == "variable") "gB1",
          if (gb2 == "variable") "gB2",
          if (gb3 == "variable") "gB3",
          if (gb4 == "variable") "gB4",
          if (gh1 == "variable") "gH1",
          if (gh2 == "variable") "gH2"
        )
        if (length(gt) == 0) "—" else paste(gt, collapse = "/")
      }
    )
  ) %>%
  arrange(desc(n_ep), desc(n_cambios)) %>%
  mutate(pac_num = row_number()) %>%
  select(pac_num, serDR, n_ep, estables, variables, n_cambios)

# ── 9. Flextable ──────────────────────────────────────────────────────────────
tabla_inter_ft <- tabla_inter %>%
  rename(
    "Paciente"              = pac_num,
    "Serología D/R"         = serDR,
    "N episodios"           = n_ep,
    "Patrón dominante"      = estables,
    "Genotipos fluctuantes" = variables,
    "N cambios"             = n_cambios
  )

fila_sin_estable <- which(tabla_inter_ft$`Patrón dominante` == "—")

ft_inter <- flextable(tabla_inter_ft) %>%
  width(j = "Paciente",              width = 0.7) %>%
  width(j = "Serología D/R",         width = 1.2) %>%
  width(j = "N episodios",           width = 0.9) %>%
  width(j = "Patrón dominante",      width = 2.0) %>%
  width(j = "Genotipos fluctuantes", width = 2.2) %>%
  width(j = "N cambios",             width = 0.9) %>%
  align(align = "center", part = "all") %>%
  align(j = c("Patrón dominante", "Genotipos fluctuantes"),
        align = "left", part = "body") %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header") %>%
  bg(bg = "white", part = "all") %>%
  bold(i = fila_sin_estable, part = "body") %>%
  color(i = fila_sin_estable, color = "firebrick3", part = "body") %>%
  border_remove() %>%
  hline_top(part = "header",
            border = fp_border(color = "black", width = 1.2)) %>%
  hline(part = "header",
        border = fp_border(color = "black", width = 0.8)) %>%
  hline_bottom(part = "body",
               border = fp_border(color = "black", width = 1.2)) %>%
  padding(padding.top = 2, padding.bottom = 2,
          padding.left = 4, padding.right = 4, part = "all") %>%
  add_footer_lines(values = c(
    "Patrón dominante: genotipos presentes en todos los episodios del paciente.",
    "Genotipos fluctuantes: genotipos con aparición o desaparición entre episodios.",
    "El paciente en rojo corresponde al único caso sin ningún genotipo estable entre episodios.",
    "N cambios: número de transiciones de patrón genotípico entre episodios consecutivos."
  )) %>%
  fontsize(size = 8, part = "footer") %>%
  color(color = "grey30", part = "footer") %>%
  set_caption(
    caption = "Tabla R9. Patrón de cambio de genotipos gB/gH en pacientes con infecciones múltiples que presentaron variabilidad interepisodio."
  ) %>%
  set_table_properties(layout = "autofit")

ft_inter

guardar_tabla(ft_inter, "Tabla_17")

## ---- Datos para las Figuras 30 y 32 (diagrama de flujo intra e interepisodio) ----
# Las Figuras 30 y 32 de la tesis son infografías; este diagrama recoge sus recuentos.
grViz("
digraph flujo {

  graph [layout = dot, rankdir = TB,
         fontname = 'Arial', splines = ortho,
         nodesep = 0.6, ranksep = 0.8]

  node [shape = rectangle, style = filled,
        fontname = 'Arial', fontsize = 19,
        width = 2.0, height = 0.7, margin = '0.18,0.12']

  edge [color = '#999999', arrowsize = 0.5,
        arrowhead = vee]

subgraph cluster_A {
  label = < <B>A. Variabilidad intraepisodio</B> >
  fontsize = 25
  labeljust = 'c'
  margin = 40
  color = white

    A1 [label = '234 episodios\ngenotipado gB/gH',     fillcolor = '#607D8B', fontcolor = white]
    A2 [label = '131 (56,0%)\n1 sola muestra',         fillcolor = '#90A4AE', fontcolor = white]
    A3 [label = '103 (44,0%)\n≥2 muestras',            fillcolor = '#607D8B', fontcolor = white]
    A4 [label = '61 (59,2%)\ngenotipo único',           fillcolor = '#1565C0', fontcolor = white]
    A5 [label = '42 (40,8%)\ngenotipo mixto',           fillcolor = '#BF360C', fontcolor = white]
    A6 [label = '61 (100%)\npatrón estable',            fillcolor = '#2E7D32', fontcolor = white]
    A7 [label = '7 (16,7%)\nestable',                   fillcolor = '#2E7D32', fontcolor = white]
    A8 [label = '35 (83,3%)\nfluctuante',               fillcolor = '#E65100', fontcolor = white]
    A9 [label = '34 (97,1%)\n≥1 dominante',             fillcolor = '#2E7D32', fontcolor = white]
    A10[label = '1 (2,9%)\ncambio completo',            fillcolor = '#7B1FA2', fontcolor = white]

    A1 -> A2
    A1 -> A3
    A3 -> A4
    A3 -> A5
    A4 -> A6
    A5 -> A7
    A5 -> A8
    A8 -> A9
    A8 -> A10

    {rank = same; A2; A3}
    {rank = same; A4; A5}
    {rank = same; A6; A7; A8}
    {rank = same; A9; A10}
  }

  subgraph cluster_B {
    label = < <B>B. Variabilidad interepisodio</B> >
    fontsize = 25
    labeljust = 'c'
      margin = 40
    color = white

    B1 [label = '146 pacientes\ngenotipado gB/gH',     fillcolor = '#607D8B', fontcolor = white]
    B2 [label = '90 (61,6%)\n1 solo episodio',         fillcolor = '#90A4AE', fontcolor = white]
    B3 [label = '56 (38,4%)\n>1 episodio',             fillcolor = '#607D8B', fontcolor = white]
    B4 [label = '35 (62,5%)\nmantienen patrón',        fillcolor = '#1565C0', fontcolor = white]
    B5 [label = '21 (37,5%)\ncambian patrón',          fillcolor = '#BF360C', fontcolor = white]
    B6 [label = '26 (74,3%)\nsiempre único',           fillcolor = '#1565C0', fontcolor = white]
    B7 [label = '9 (25,7%)\nsiempre mixto',            fillcolor = '#BF360C', fontcolor = white]
    B8 [label = '20 (95,2%)\n≥1 dominante',            fillcolor = '#2E7D32', fontcolor = white]
    B9 [label = '1 (4,8%)\ncambio completo',           fillcolor = '#7B1FA2', fontcolor = white]

    B1 -> B2
    B1 -> B3
    B3 -> B4
    B3 -> B5
    B4 -> B6
    B4 -> B7
    B5 -> B8
    B5 -> B9

    {rank = same; B2; B3}
    {rank = same; B4; B5}
    {rank = same; B6; B7; B8; B9}
  }
}
")


# ---- 9.7. Impacto clínico de la detección de genotipos mixtos ----

## ---- Tabla 18. Características basales de los pacientes según genotipos mixtos gB/gH ----
# Este bloque vuelve a leer pacientes_sel: es la base de pacientes del resto de 9.7.
pacientes <- read_excel(f_bd,
                        sheet = "pacientes_sel")

# ── 1. Preparar datos ─────────────────────────────────────────────────────────
pacientes_tab <- pacientes %>%
  mutate(
    Ftx         = as.Date(Ftx),
    Ffinseg     = as.Date(Ffinseg),
    t_seg_anios = as.numeric(Ffinseg - Ftx) / 365.25,
    grupo       = factor(gbghm, levels = c(0,1),
                         labels = c("Genotipos únicos", "Genotipos mixtos"))
  )

g0  <- pacientes_tab %>% filter(grupo == "Genotipos únicos")
g1  <- pacientes_tab %>% filter(grupo == "Genotipos mixtos")
tot <- pacientes_tab
n0  <- nrow(g0); n1 <- nrow(g1); nt <- nrow(tot)

# ── 2. Funciones auxiliares ───────────────────────────────────────────────────
med_ric <- function(x) {
  x <- x[!is.na(x)]
  sprintf("%.1f (%.1f\u2013%.1f)",
          median(x), quantile(x, 0.25), quantile(x, 0.75))
}

n_pct <- function(x, n_total) {
  sprintf("%d (%.1f%%)", sum(x, na.rm = TRUE),
          100 * sum(x, na.rm = TRUE) / n_total)
}

fmt_p <- function(p) {
  ifelse(p < 0.001, "<0,001",
         sub("\\.", ",", sprintf("%.3f", p)))
}

p_cont <- function(x, g) {
  fmt_p(wilcox.test(x ~ g, exact = FALSE)$p.value)
}

p_cat <- function(x, g) {
  tab <- table(x, g)
  test <- tryCatch(
    fisher.test(tab, simulate.p.value = TRUE, B = 10000),
    error = function(e) chisq.test(tab)
  )
  fmt_p(test$p.value)
}

# ── 3. Construir filas ────────────────────────────────────────────────────────
filas <- bind_rows(
  
  tibble(variable = "n",
         total = as.character(nt),
         unicos = as.character(n0),
         mixtos = as.character(n1),
         p = ""),
  
  # ── Demográficas ──────────────────────────────────────────────────────────
  tibble(variable = "CARACTERÍSTICAS DEMOGRÁFICAS",
         total = "", unicos = "", mixtos = "", p = ""),
  
  tibble(variable = "Edad receptor (años), mediana (RIC)",
         total  = med_ric(tot$edad_tx),
         unicos = med_ric(g0$edad_tx),
         mixtos = med_ric(g1$edad_tx),
         p = p_cont(pacientes_tab$edad_tx, pacientes_tab$grupo)),
  
  tibble(variable = "Edad donante (años), mediana (RIC)",
         total  = med_ric(tot$edad_don),
         unicos = med_ric(g0$edad_don),
         mixtos = med_ric(g1$edad_don),
         p = p_cont(pacientes_tab$edad_don, pacientes_tab$grupo)),
  
  tibble(variable = "Sexo femenino receptor, n (%)",
         total  = n_pct(tot$SexF == 1, nt),
         unicos = n_pct(g0$SexF  == 1, n0),
         mixtos = n_pct(g1$SexF  == 1, n1),
         p = p_cat(pacientes_tab$SexF, pacientes_tab$grupo)),
  
  tibble(variable = "Sexo femenino donante, n (%)",
         total  = n_pct(tot$SexF_don == 1, nt),
         unicos = n_pct(g0$SexF_don  == 1, n0),
         mixtos = n_pct(g1$SexF_don  == 1, n1),
         p = p_cat(pacientes_tab$SexF_don, pacientes_tab$grupo)),
  
  # ── Seguimiento ──────────────────────────────────────────────────────────
  tibble(variable = "CARACTERÍSTICAS DE SEGUIMIENTO",
         total = "", unicos = "", mixtos = "", p = ""),
  
  tibble(variable = "Tiempo de seguimiento (años), mediana (RIC)",
         total  = med_ric(tot$t_seg_anios),
         unicos = med_ric(g0$t_seg_anios),
         mixtos = med_ric(g1$t_seg_anios),
         p = p_cont(pacientes_tab$t_seg_anios, pacientes_tab$grupo)),
  
  tibble(variable = "Episodios CMV totales, mediana (RIC)",
         total  = med_ric(tot$nepi_tot),
         unicos = med_ric(g0$nepi_tot),
         mixtos = med_ric(g1$nepi_tot),
         p = p_cont(pacientes_tab$nepi_tot, pacientes_tab$grupo)),
  
  tibble(variable = "Episodios CMV genotipados, mediana (RIC)",
         total  = med_ric(tot$nepi_gtc),
         unicos = med_ric(g0$nepi_gtc),
         mixtos = med_ric(g1$nepi_gtc),
         p = p_cont(pacientes_tab$nepi_gtc, pacientes_tab$grupo)),
  
  tibble(variable = "Determinaciones CV totales, mediana (RIC)",
         total  = med_ric(tot$Npcrs),
         unicos = med_ric(g0$Npcrs),
         mixtos = med_ric(g1$Npcrs),
         p = p_cont(pacientes_tab$Npcrs, pacientes_tab$grupo)),
  
  tibble(variable = "Antigenemias totales, mediana (RIC)",
         total  = med_ric(tot$Nag),
         unicos = med_ric(g0$Nag),
         mixtos = med_ric(g1$Nag),
         p = p_cont(pacientes_tab$Nag, pacientes_tab$grupo)),
  
  # ── Tipo de trasplante ────────────────────────────────────────────────────
  tibble(variable = "Tipo de trasplante, n (%)",
         total = "", unicos = "", mixtos = "",
         p = p_cat(pacientes_tab$tipotx, pacientes_tab$grupo)),
  
  tibble(variable = "   Unipulmonar",
         total  = n_pct(tot$tipotx == "UNI", nt),
         unicos = n_pct(g0$tipotx  == "UNI", n0),
         mixtos = n_pct(g1$tipotx  == "UNI", n1), p = ""),
  
  tibble(variable = "   Bipulmonar",
         total  = n_pct(tot$tipotx == "BI", nt),
         unicos = n_pct(g0$tipotx  == "BI", n0),
         mixtos = n_pct(g1$tipotx  == "BI", n1), p = ""),
  
  # ── Enfermedad de base ────────────────────────────────────────────────────
  tibble(variable = "Enfermedad de base, n (%)",
         total = "", unicos = "", mixtos = "",
         p = p_cat(pacientes_tab$ebase, pacientes_tab$grupo)),
  
  tibble(variable = "   EPID",
         total  = n_pct(tot$ebase == "EPID", nt),
         unicos = n_pct(g0$ebase  == "EPID", n0),
         mixtos = n_pct(g1$ebase  == "EPID", n1), p = ""),
  
  tibble(variable = "   EPOC",
         total  = n_pct(tot$ebase == "EPOC", nt),
         unicos = n_pct(g0$ebase  == "EPOC", n0),
         mixtos = n_pct(g1$ebase  == "EPOC", n1), p = ""),
  
  tibble(variable = "   FPI",
         total  = n_pct(tot$ebase == "FPI", nt),
         unicos = n_pct(g0$ebase  == "FPI", n0),
         mixtos = n_pct(g1$ebase  == "FPI", n1), p = ""),
  
  tibble(variable = "   FQ",
         total  = n_pct(tot$ebase == "FQ", nt),
         unicos = n_pct(g0$ebase  == "FQ", n0),
         mixtos = n_pct(g1$ebase  == "FQ", n1), p = ""),
  
  tibble(variable = "   HAP",
         total  = n_pct(tot$ebase == "HAP", nt),
         unicos = n_pct(g0$ebase  == "HAP", n0),
         mixtos = n_pct(g1$ebase  == "HAP", n1), p = ""),
  
  tibble(variable = "   Otros",
         total  = n_pct(!tot$ebase %in% c("EPID","EPOC","FPI","FQ","HAP"), nt),
         unicos = n_pct(!g0$ebase  %in% c("EPID","EPOC","FPI","FQ","HAP"), n0),
         mixtos = n_pct(!g1$ebase  %in% c("EPID","EPOC","FPI","FQ","HAP"), n1),
         p = ""),
  
  # ── Perfil serológico ─────────────────────────────────────────────────────
  tibble(variable = "Perfil serológico D/R, n (%)",
         total = "", unicos = "", mixtos = "",
         p = p_cat(pacientes_tab$serDR, pacientes_tab$grupo)),
  
  tibble(variable = "   D-/R+",
         total  = n_pct(tot$serDR == "D-/R+", nt),
         unicos = n_pct(g0$serDR  == "D-/R+", n0),
         mixtos = n_pct(g1$serDR  == "D-/R+", n1), p = ""),
  
  tibble(variable = "   D+/R+",
         total  = n_pct(tot$serDR == "D+/R+", nt),
         unicos = n_pct(g0$serDR  == "D+/R+", n0),
         mixtos = n_pct(g1$serDR  == "D+/R+", n1), p = ""),
  
  tibble(variable = "   D+/R-",
         total  = n_pct(tot$serDR == "D+/R-", nt),
         unicos = n_pct(g0$serDR  == "D+/R-", n0),
         mixtos = n_pct(g1$serDR  == "D+/R-", n1), p = ""),
  
  # ── Inmunosupresión ───────────────────────────────────────────────────────
  tibble(variable = "INMUNOSUPRESIÓN BASAL, n (%)",
         total = "", unicos = "", mixtos = "", p = ""),
  
  tibble(variable = "   Tacrolimus",
         total  = n_pct(grepl("TACRO", tot$IS_base), nt),
         unicos = n_pct(grepl("TACRO", g0$IS_base),  n0),
         mixtos = n_pct(grepl("TACRO", g1$IS_base),  n1),
         p = p_cat(grepl("TACRO", pacientes_tab$IS_base),
                   pacientes_tab$grupo)),
  
  tibble(variable = "   Ciclosporina",
         total  = n_pct(grepl("CICLO", tot$IS_base), nt),
         unicos = n_pct(grepl("CICLO", g0$IS_base),  n0),
         mixtos = n_pct(grepl("CICLO", g1$IS_base),  n1),
         p = p_cat(grepl("CICLO", pacientes_tab$IS_base),
                   pacientes_tab$grupo)),
  
  tibble(variable = "   Micofenolato",
         total  = n_pct(grepl("MICO", tot$IS_base), nt),
         unicos = n_pct(grepl("MICO", g0$IS_base),  n0),
         mixtos = n_pct(grepl("MICO", g1$IS_base),  n1),
         p = p_cat(grepl("MICO", pacientes_tab$IS_base),
                   pacientes_tab$grupo)),
  
  tibble(variable = "   Azatioprina",
         total  = n_pct(grepl("AZA", tot$IS_base), nt),
         unicos = n_pct(grepl("AZA", g0$IS_base),  n0),
         mixtos = n_pct(grepl("AZA", g1$IS_base),  n1),
         p = p_cat(grepl("AZA", pacientes_tab$IS_base),
                   pacientes_tab$grupo)),
  
  tibble(variable = "   Esteroides",
         total  = "146 (100%)",
         unicos = "95 (100%)",
         mixtos = "51 (100%)",
         p = "\u2014")
)

# ── 4. Flextable ──────────────────────────────────────────────────────────────
filas_seccion <- which(filas$variable %in% c(
  "CARACTERÍSTICAS DEMOGRÁFICAS",
  "CARACTERÍSTICAS DE SEGUIMIENTO",
  "INMUNOSUPRESIÓN BASAL, n (%)"
))

filas_padre <- which(filas$variable %in% c(
  "Tipo de trasplante, n (%)",
  "Enfermedad de base, n (%)",
  "Perfil serológico D/R, n (%)"
))

colnames(filas) <- c("Variable",
                     "Total\n(n = 146)",
                     "Genotipos únicos\n(n = 95)",
                     "Genotipos mixtos\n(n = 51)",
                     "p")

# p significativas para resaltar
p_sig <- c("<0,001")

ft <- flextable(filas) %>%
  width(j = "Variable",                   width = 4.5) %>%
  width(j = "Total\n(n = 146)",           width = 2.0) %>%
  width(j = "Genotipos únicos\n(n = 95)", width = 2.0) %>%
  width(j = "Genotipos mixtos\n(n = 51)", width = 2.0) %>%
  width(j = "p",                          width = 0.9) %>%
  align(align = "left",   part = "all") %>%
  align(j = c(2,3,4,5), align = "center", part = "all") %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header") %>%
  bg(bg = "white", part = "all") %>%
  bold(i = filas_seccion, part = "body") %>%
  bg(i = filas_seccion, bg = "#F0F0F0", part = "body") %>%
  italic(i = filas_padre, part = "body") %>%
  color(i = ~ p %in% p_sig, j = "p",
        color = "firebrick3", part = "body") %>%
  bold(i  = ~ p %in% p_sig, j = "p", part = "body") %>%
  border_remove() %>%
  hline_top(part = "header",
            border = fp_border(color = "black", width = 1.2)) %>%
  hline(part = "header",
        border = fp_border(color = "black", width = 0.8)) %>%
  hline_bottom(part = "body",
               border = fp_border(color = "black", width = 1.2)) %>%
  hline(i = filas_seccion - 1,
        border = fp_border(color = "grey70", width = 0.4),
        part = "body") %>%
  padding(padding.top = 2, padding.bottom = 2,
          padding.left = 4, padding.right = 4, part = "all") %>%
  add_footer_lines(values = c(
    "RIC: rango intercuartílico; EPID: enfermedad pulmonar intersticial difusa; EPOC: enfermedad pulmonar obstructiva crónica; FPI: fibrosis pulmonar idiopática; FQ: fibrosis quística; HAP: hipertensión arterial pulmonar.",
    "Variables continuas: test de Wilcoxon. Variables categóricas: test exacto de Fisher o chi-cuadrado.",
    "Valores de p en rojo indican diferencia estadísticamente significativa (p < 0,05)."
  )) %>%
  fontsize(size = 8, part = "footer") %>%
  color(color = "grey30", part = "footer") %>%
  set_caption(
    caption = "Tabla R11. Características basales de los pacientes según presencia de genotipos mixtos gB/gH."
  ) %>%
  set_table_properties(layout = "autofit")

ft

guardar_tabla(ft, "Tabla_18")

## ---- Tabla 19. Asociación entre genotipos mixtos y enfermedad sintomática por CMV ----
# ── 1. Preparar datos a nivel paciente ────────────────────────────────────────
pacientes_sintom <- pacientes %>%
  mutate(
    sintomatico = as.integer(vira == 0),
    gbghm_f     = factor(gbghm, levels = c(0,1), labels = c("Único","Mixto")),
    gbm_f       = factor(gbm,   levels = c(0,1), labels = c("Único","Mixto")),
    ghm_f       = factor(ghm,   levels = c(0,1), labels = c("Único","Mixto")),
    serDR_f     = factor(serDR)
  )

cat("── Distribución sintomático vs grupo ──\n")
pacientes_sintom %>%
  count(sintomatico, gbghm_f) %>%
  group_by(gbghm_f) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  print()

# ── 2. Funciones ──────────────────────────────────────────────────────────────
calc_or <- function(var_mix, data, label) {
  m_crudo <- glm(sintomatico ~ get(var_mix),
                 data = data, family = binomial)
  or_c <- exp(coef(m_crudo)[2])
  ci_c <- exp(confint.default(m_crudo)[2,])
  p_c  <- summary(m_crudo)$coefficients[2,4]
  
  m_ajus <- glm(sintomatico ~ get(var_mix) + serDR_f,
                data = data, family = binomial)
  or_a <- exp(coef(m_ajus)[2])
  ci_a <- exp(confint.default(m_ajus)[2,])
  p_a  <- summary(m_ajus)$coefficients[2,4]
  
  cat(sprintf("\n── %s ──\n", label))
  cat(sprintf("OR cruda:    %.2f (%.2f–%.2f)  p = %.3f\n",
              or_c, ci_c[1], ci_c[2], p_c))
  cat(sprintf("OR ajustada: %.2f (%.2f–%.2f)  p = %.3f\n",
              or_a, ci_a[1], ci_a[2], p_a))
  
  tibble(
    genotipo    = label,
    or_cruda    = sprintf("%.2f (%.2f\u2013%.2f)", or_c, ci_c[1], ci_c[2]),
    p_cruda     = ifelse(p_c < 0.001, "<0,001",
                         sub("\\.", ",", sprintf("%.3f", p_c))),
    or_ajustada = sprintf("%.2f (%.2f\u2013%.2f)", or_a, ci_a[1], ci_a[2]),
    p_ajustada  = ifelse(p_a < 0.001, "<0,001",
                         sub("\\.", ",", sprintf("%.3f", p_a)))
  )
}

calc_or_mixto <- function(var_mix, data, label) {
  m <- tryCatch(
    glmer(sintomatico ~ get(var_mix) + (1|NHC),
          data = data, family = binomial,
          control = glmerControl(optimizer = "bobyqa")),
    error = function(e) NULL
  )
  if (is.null(m)) {
    cat(sprintf("\n── %s: error en modelo mixto ──\n", label))
    return(tibble(genotipo = label, or_mixto = "error", p_mixto = "error"))
  }
  coefs <- summary(m)$coefficients
  or    <- exp(coefs[2, 1])
  se    <- coefs[2, 2]
  ci_lo <- exp(coefs[2, 1] - 1.96 * se)
  ci_hi <- exp(coefs[2, 1] + 1.96 * se)
  p     <- coefs[2, 4]
  
  cat(sprintf("\n── %s (modelo mixto) ──\n", label))
  cat(sprintf("OR: %.2f (%.2f–%.2f)  p = %.3f\n", or, ci_lo, ci_hi, p))
  
  tibble(
    genotipo = label,
    or_mixto = sprintf("%.2f (%.2f\u2013%.2f)", or, ci_lo, ci_hi),
    p_mixto  = ifelse(p < 0.001, "<0,001",
                      sub("\\.", ",", sprintf("%.3f", p)))
  )
}

freq_gt <- function(var, label, data) {
  tab    <- table(data[[var]], data$sintomatico)
  test   <- fisher.test(tab, simulate.p.value = TRUE, B = 10000)
  n0     <- sum(data$sintomatico == 0)
  n1     <- sum(data$sintomatico == 1)
  n_var0 <- sum(data[[var]] == 1 & data$sintomatico == 0, na.rm = TRUE)
  n_var1 <- sum(data[[var]] == 1 & data$sintomatico == 1, na.rm = TRUE)
  tibble(
    genotipo = label,
    asintom  = sprintf("%d (%.1f%%)", n_var0, 100*n_var0/n0),
    sintom   = sprintf("%d (%.1f%%)", n_var1, 100*n_var1/n1),
    p        = ifelse(test$p.value < 0.001, "<0,001",
                      sub("\\.", ",", sprintf("%.3f", test$p.value)))
  )
}

# ── 3. OR cruda y ajustada a nivel paciente ───────────────────────────────────
res_or <- bind_rows(
  calc_or("gbm_f",   pacientes_sintom, "gB mixtos"),
  calc_or("ghm_f",   pacientes_sintom, "gH mixtos"),
  calc_or("gbghm_f", pacientes_sintom, "gB+gH mixtos")
)

cat("\n── Tabla OR resumen ──\n")
print(res_or)

# ── 4. Preparar datos a nivel episodio ────────────────────────────────────────
dnaemias_sintom <- dnaemias %>%
  mutate(
    sintomatico = as.integer(!is.na(clinica) & clinica != "VA"),
    gbghm_f     = factor(gbghmix, levels = c(0,1), labels = c("Único","Mixto")),
    gbm_f       = factor(gbmix,   levels = c(0,1), labels = c("Único","Mixto")),
    ghm_f       = factor(ghmix,   levels = c(0,1), labels = c("Único","Mixto"))
  ) %>%
  left_join(
    pacientes %>% select(NHC, serDR_pac = serDR),
    by = "NHC"
  ) %>%
  mutate(serDR_f = factor(serDR_pac))

cat(sprintf("\nTotal episodios: %d\n", nrow(dnaemias_sintom)))
cat(sprintf("Sintomáticos: %d\n",     sum(dnaemias_sintom$sintomatico)))

# ── 5. Modelo mixto por episodio ──────────────────────────────────────────────
res_mixto <- bind_rows(
  calc_or_mixto("gbm_f",   dnaemias_sintom, "gB mixtos"),
  calc_or_mixto("ghm_f",   dnaemias_sintom, "gH mixtos"),
  calc_or_mixto("gbghm_f", dnaemias_sintom, "gB+gH mixtos")
)

cat("\n── Tabla OR mixto resumen ──\n")
print(res_mixto)

# ── 6. Frecuencias genotípicas por episodio ───────────────────────────────────
res_freq <- bind_rows(
  freq_gt("gbmix",   "gB mixtos",    dnaemias_sintom),
  freq_gt("ghmix",   "gH mixtos",    dnaemias_sintom),
  freq_gt("gbghmix", "gB+gH mixtos", dnaemias_sintom)
)

cat("\n── Frecuencias genotípicas ──\n")
print(res_freq)


# ── Combinar resultados ───────────────────────────────────────────────────────
tabla_sintom <- res_or %>%
  left_join(res_mixto, by = "genotipo") %>%
  rename(
    "Genotipos mixtos"          = genotipo,
    "OR cruda\n(IC 95%)"        = or_cruda,
    "p"                         = p_cruda,
    "OR ajustada*\n(IC 95%)"    = or_ajustada,
    "p "                        = p_ajustada,
    "OR por episodio**\n(IC 95%)"= or_mixto,
    "p  "                       = p_mixto
  )

ft_sintom <- flextable(tabla_sintom) %>%
  width(j = "Genotipos mixtos",           width = 1.6) %>%
  width(j = "OR cruda\n(IC 95%)",         width = 2.0) %>%
  width(j = "p",                          width = 0.8) %>%
  width(j = "OR ajustada*\n(IC 95%)",     width = 2.0) %>%
  width(j = "p ",                         width = 0.8) %>%
  width(j = "OR por episodio**\n(IC 95%)",width = 2.0) %>%
  width(j = "p  ",                        width = 0.8) %>%
  align(align = "center", part = "all") %>%
  align(j = "Genotipos mixtos", align = "left", part = "all") %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header") %>%
  bg(bg = "white", part = "all") %>%
  bg(i = seq(2, nrow(tabla_sintom), 2), bg = "#F7F7F7", part = "body") %>%
  border_remove() %>%
  hline_top(part = "header",
            border = fp_border(color = "black", width = 1.2)) %>%
  hline(part = "header",
        border = fp_border(color = "black", width = 0.8)) %>%
  hline_bottom(part = "body",
               border = fp_border(color = "black", width = 1.2)) %>%
  padding(padding.top = 3, padding.bottom = 3,
          padding.left = 4, padding.right = 4, part = "all") %>%
  add_footer_lines(values = c(
    "*Ajustada por perfil serológico D/R.",
    "**Modelo de regresión logística mixta con efecto aleatorio por paciente."
  )) %>%
  fontsize(size = 8, part = "footer") %>%
  color(color = "grey30", part = "footer") %>%
  set_caption(
    caption = "Tabla R12. Asociación entre genotipos mixtos gB/gH y enfermedad sintomática por CMV."
  ) %>%
  set_table_properties(layout = "autofit")

ft_sintom

guardar_tabla(ft_sintom, "Tabla_19")

## ---- Figura 33. Tiempo desde el trasplante hasta la primera DNAemia genotipada ----
## Panel A: Genotipos mixtos a lo largo del seguimiento (nivel paciente)
##          -> variables gbm/ghm/gbghm desde tabla `pacientes`
## Panel B: Genotipos mixtos en el primer episodio (nivel episodio)
##          -> variables gbmix/ghmix/gbghmix desde tabla `dnaemias`
## Filtro: primer_ev = 1
## Orden de paneles: D+/R- -> D+/R+ -> D-/R+


DODGE <- 0.7

#  Paletas (compartidas)
cols_fill <- c(
  "gB_\u00danico"     = "#AED6F1",  "gB_Mixto"     = "#4682B4",
  "gH_\u00danico"     = "pink",     "gH_Mixto"     = "hotpink3",
  "gB+gH_\u00danico"  = "#4CAF7D",  "gB+gH_Mixto"  = "#E07070"
)

cols_jit <- c(
  "gB_\u00danico"     = "#2980B9",  "gB_Mixto"     = "#0D2B45",
  "gH_\u00danico"     = "pink3",    "gH_Mixto"     = "palevioletred4",
  "gB+gH_\u00danico"  = "#2E7D52",  "gB+gH_Mixto"  = "#B04040"
)

df_prof <- tibble::tibble(
  serDR = factor(c("D+/R-", "D+/R+", "D-/R+"),
                 levels = c("D+/R-", "D+/R+", "D-/R+")),
  prof_meses = c(12, 6, 6)
)

#  FUNCION CONSTRUCTORA DEL PANEL
construir_panel <- function(df_long, etiquetas_riesgo, titulo_panel) {
  
  # Tests significativos del panel
  stat_test <- df_long %>%
    dplyr::group_by(serDR, glucoproteina) %>%
    rstatix::wilcox_test(meses_tx_inicio ~ genotipo) %>%
    rstatix::add_significance("p") %>%
    dplyr::ungroup() %>%
    dplyr::filter(p.signif != "ns") %>%
    dplyr::mutate(
      x_pos      = as.numeric(glucoproteina),
      xmin       = x_pos - DODGE / 4,
      xmax       = x_pos + DODGE / 4,
      y.position = 22
    )
  
  ggplot(df_long,
         aes(x = glucoproteina, y = meses_tx_inicio,
             fill = cat_color,
             group = interaction(glucoproteina, genotipo))) +
    geom_hline(data = df_prof,
               aes(yintercept = prof_meses),
               color = "grey55", linetype = "dashed", linewidth = 0.5) +
    geom_boxplot(
      width          = 0.55,
      position       = position_dodge(width = DODGE),
      color          = "grey25",
      linewidth      = 0.4,
      outlier.shape  = 21,
      outlier.size   = 1.6,
      outlier.fill   = "white",
      outlier.color  = "grey50",
      outlier.stroke = 0.3,
      outlier.alpha  = 0.7
    ) +
    geom_jitter(
      aes(color = cat_color),
      position    = position_jitterdodge(jitter.width = 0.10,
                                         dodge.width  = DODGE),
      size        = 1.0,
      alpha       = 0.4,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = cols_fill, guide = "none") +
    scale_color_manual(values = cols_jit, guide = "none") +
    facet_wrap(~ serDR, labeller = labeller(serDR = etiquetas_riesgo)) +
    {
      if (nrow(stat_test) > 0)
        ggpubr::stat_pvalue_manual(
          stat_test,
          label        = "p.signif",
          xmin         = "xmin",
          xmax         = "xmax",
          y.position   = "y.position",
          tip.length   = 0.01,
          bracket.size = 0.5,
          size         = 6,
          color        = "grey20",
          inherit.aes  = FALSE
        )
    } +
    labs(
      x = NULL,
      y = "Meses desde TX a 1er episodio",
      title = titulo_panel
    ) +
    scale_y_continuous(
      breaks = seq(0, 24, 6),
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    coord_cartesian(ylim = c(0, 24)) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position    = "none",
      panel.border       = element_blank(),
      panel.spacing.x    = unit(15, "pt"),
      axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
      axis.line.y        = element_line(color = "grey60", linewidth = 0.4),
      axis.ticks         = element_line(color = "grey60", linewidth = 0.3),
      axis.text.x        = element_text(size = 12, color = "grey20"),
      axis.text.y        = element_text(size = 11, color = "grey30"),
      axis.title.y       = element_text(size = 12, color = "grey15",
                                        margin = margin(r = 8)),
      plot.title         = element_text(size = 13, color = "grey15",
                                        face = "bold",
                                        margin = margin(b = 8)),
      strip.background   = element_rect(fill = "grey94",
                                        color = "grey70", linewidth = 0.4),
      strip.text         = element_text(size = 12, color = "grey15",
                                        face = "bold",
                                        margin = margin(t = 5, b = 5)),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.margin        = margin(10, 12, 6, 10)
    )
}

#  PANEL A - Nivel PACIENTE (gbm, ghm, gbghm desde tabla pacientes)
df_fig_A <- dnaemias %>%
  dplyr::filter(primer_ev == 1) %>%
  dplyr::mutate(
    Ftx              = as.Date(Ftx),
    Fini             = as.Date(Fini),
    meses_tx_inicio  = as.numeric(Fini - Ftx) / 30.44
  ) %>%
  dplyr::select(-any_of(c("gbm", "ghm", "gbghm"))) %>%
  dplyr::left_join(pacientes %>% dplyr::select(NHC, gbm, ghm, gbghm), by = "NHC") %>%
  dplyr::mutate(
    serDR = factor(serDR, levels = c("D+/R-", "D+/R+", "D-/R+"))
  ) %>%
  dplyr::filter(!is.na(serDR), !is.na(meses_tx_inicio))

df_long_A <- df_fig_A %>%
  dplyr::select(NHC, serDR, meses_tx_inicio, gbm, ghm, gbghm) %>%
  tidyr::pivot_longer(
    cols      = c(gbm, ghm, gbghm),
    names_to  = "glucoproteina",
    values_to = "mixto"
  ) %>%
  dplyr::filter(!is.na(mixto)) %>%
  dplyr::mutate(
    glucoproteina = factor(glucoproteina,
                           levels = c("gbm", "ghm", "gbghm"),
                           labels = c("gB", "gH", "gB+gH")),
    genotipo      = factor(ifelse(mixto == 1, "Mixto", "\u00danico"),
                           levels = c("\u00danico", "Mixto")),
    cat_color = paste(glucoproteina, genotipo, sep = "_"),
    cat_color = factor(cat_color,
                       levels = c("gB_\u00danico",   "gB_Mixto",
                                  "gH_\u00danico",   "gH_Mixto",
                                  "gB+gH_\u00danico","gB+gH_Mixto"))
  )

counts_A <- df_fig_A %>%
  dplyr::group_by(serDR) %>%
  dplyr::summarise(n = dplyr::n_distinct(NHC), .groups = "drop")

etiquetas_A <- setNames(
  sprintf("%s (n = %d)", counts_A$serDR, counts_A$n),
  as.character(counts_A$serDR)
)

panel_A <- construir_panel(
  df_long          = df_long_A,
  etiquetas_riesgo = etiquetas_A,
  titulo_panel     = "A. Genotipos únicos/mixtos a lo largo del seguimiento"
)

#  PANEL B - Nivel EPISODIO (gbmix, ghmix, gbghmix desde tabla dnaemias)
df_fig_B <- dnaemias %>%
  dplyr::filter(primer_ev == 1) %>%
  dplyr::mutate(
    Ftx              = as.Date(Ftx),
    Fini             = as.Date(Fini),
    meses_tx_inicio  = as.numeric(Fini - Ftx) / 30.44,
    serDR = factor(serDR, levels = c("D+/R-", "D+/R+", "D-/R+"))
  ) %>%
  dplyr::filter(!is.na(serDR), !is.na(meses_tx_inicio))

df_long_B <- df_fig_B %>%
  dplyr::select(NHC, serDR, meses_tx_inicio, gbmix, ghmix, gbghmix) %>%
  tidyr::pivot_longer(
    cols      = c(gbmix, ghmix, gbghmix),
    names_to  = "glucoproteina",
    values_to = "mixto"
  ) %>%
  dplyr::filter(!is.na(mixto)) %>%
  dplyr::mutate(
    glucoproteina = factor(glucoproteina,
                           levels = c("gbmix", "ghmix", "gbghmix"),
                           labels = c("gB", "gH", "gB+gH")),
    genotipo      = factor(ifelse(mixto == 1, "Mixto", "\u00danico"),
                           levels = c("\u00danico", "Mixto")),
    cat_color = paste(glucoproteina, genotipo, sep = "_"),
    cat_color = factor(cat_color,
                       levels = c("gB_\u00danico",   "gB_Mixto",
                                  "gH_\u00danico",   "gH_Mixto",
                                  "gB+gH_\u00danico","gB+gH_Mixto"))
  )

counts_B <- df_fig_B %>%
  dplyr::group_by(serDR) %>%
  dplyr::summarise(n = dplyr::n_distinct(NHC), .groups = "drop")

etiquetas_B <- setNames(
  sprintf("%s (n = %d)", counts_B$serDR, counts_B$n),
  as.character(counts_B$serDR)
)

panel_B <- construir_panel(
  df_long          = df_long_B,
  etiquetas_riesgo = etiquetas_B,
  titulo_panel     = "B. Genotipos únicos/mixtos en el primer episodio"
)

#  COMBINAR PANELES A y B
figura_R20_doble <- panel_A / panel_B +
  patchwork::plot_layout(heights = c(1, 1))

figura_R20_doble

# [guardado antiguo] ggsave("~/Desktop/Figuras_tesis/figura_R20_doble_paciente_vs_episodio.png",
# [guardado antiguo]        figura_R20_doble,
# [guardado antiguo]        width = 10, height = 9, dpi = 320, bg = "white")

guardar_fig(figura_R20_doble, "Figura_33", ancho = 25.4, alto = 22.9)

## ---- Datos para el texto (9.7.2): tiempo hasta la primera DNAemia ----
df_base <- dnaemias %>%
  dplyr::filter(primer_ev == 1) %>%
  dplyr::mutate(
    Ftx              = as.Date(Ftx),
    Fini             = as.Date(Fini),
    meses_tx_inicio  = as.numeric(Fini - Ftx) / 30.44,
    serDR = factor(serDR, levels = c("D+/R-", "D+/R+", "D-/R+"))
  ) %>%
  dplyr::filter(!is.na(serDR), !is.na(meses_tx_inicio))

df_A <- df_base %>%
  dplyr::select(-any_of(c("gbm", "ghm", "gbghm"))) %>%
  dplyr::left_join(pacientes %>% dplyr::select(NHC, gbm, ghm, gbghm),
                   by = "NHC")

df_B <- df_base

fmt_mediana_ric <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("NA")
  med <- median(x)
  ric <- quantile(x, probs = c(0.25, 0.75))
  sprintf("%.1f (RIC %.1f-%.1f)", med, ric[1], ric[2])
}

# >>> FUNCION CORREGIDA: extraer vectores directamente para wilcox.test <<<
comparar_grupos <- function(datos, var_mixto, serologia, etiqueta) {
  sub <- datos %>%
    dplyr::filter(serDR == serologia) %>%
    dplyr::filter(!is.na(.data[[var_mixto]])) %>%
    dplyr::filter(!is.na(meses_tx_inicio))
  
  if (nrow(sub) == 0) {
    cat(sprintf("  %s: sin datos\n", etiqueta))
    return(invisible(NULL))
  }
  
  # Extraer vectores directamente con [[ ]]
  vals     <- sub$meses_tx_inicio
  grupos   <- sub[[var_mixto]]
  
  mixtos <- vals[grupos == 1]
  unicos <- vals[grupos == 0]
  
  cat(sprintf("  %s:\n", etiqueta))
  cat(sprintf("    Mixtos  (n=%d): %s\n",
              length(mixtos), fmt_mediana_ric(mixtos)))
  cat(sprintf("    Unicos  (n=%d): %s\n",
              length(unicos), fmt_mediana_ric(unicos)))
  
  # >>> CORRECCION: pasar vectores directamente, no formula <<<
  if (length(mixtos) >= 1 && length(unicos) >= 1) {
    p <- tryCatch(
      stats::wilcox.test(x = mixtos, y = unicos, exact = FALSE)$p.value,
      error = function(e) {
        cat(sprintf("    ERROR en wilcox: %s\n", e$message))
        NA
      }
    )
    if (!is.na(p)) {
      if (p < 0.001) {
        cat("    p < 0.001\n")
      } else {
        cat(sprintf("    p = %.3f\n", p))
      }
    } else {
      cat("    p = NA\n")
    }
  }
}

for (grupo in c("D+/R+", "D-/R+", "D+/R-")) {
  cat("\n========================================\n")
  cat(sprintf("GRUPO %s\n", grupo))
  cat("========================================\n")
  
  n_grupo <- df_base %>%
    dplyr::filter(serDR == grupo) %>%
    dplyr::pull(NHC) %>%
    dplyr::n_distinct()
  cat(sprintf("N pacientes: %d\n\n", n_grupo))
  
  med_global <- df_base %>%
    dplyr::filter(serDR == grupo) %>%
    dplyr::pull(meses_tx_inicio)
  cat(sprintf("Mediana global: %s\n\n", fmt_mediana_ric(med_global)))
  
  cat("ABORDAJE A (nivel paciente):\n")
  comparar_grupos(df_A, "gbm",   grupo, "gB")
  comparar_grupos(df_A, "ghm",   grupo, "gH")
  comparar_grupos(df_A, "gbghm", grupo, "gB+gH")
  
  cat("\nABORDAJE B (nivel episodio):\n")
  comparar_grupos(df_B, "gbmix",   grupo, "gB")
  comparar_grupos(df_B, "ghmix",   grupo, "gH")
  comparar_grupos(df_B, "gbghmix", grupo, "gB+gH")
}

cat("\n========================================\n")
cat("DISTRIBUCION TEMPORAL D+/R+\n")
cat("========================================\n")

dist_dpr <- df_base %>%
  dplyr::filter(serDR == "D+/R+") %>%
  dplyr::mutate(tras_6m = meses_tx_inicio > 6) %>%
  dplyr::count(tras_6m) %>%
  dplyr::mutate(pct = round(100 * n / sum(n), 1))
print(dist_dpr)

## ---- Figura 34. Distribución del número de episodios por paciente según genotipos mixtos ----
# ── 1. Preparar datos ────────────────────────────────────────────────────────
df_dist <- pacientes %>%
  filter(!is.na(nepi_tot), !is.na(gbghm)) %>%
  mutate(
    genotipo = factor(ifelse(gbghm == 1, "Mixto", "Único"),
                      levels = c("Único", "Mixto")),
    cat_ep = factor(case_when(
      nepi_tot == 1 ~ "1",
      nepi_tot == 2 ~ "2",
      nepi_tot == 3 ~ "3",
      nepi_tot == 4 ~ "4",
      nepi_tot >= 5 ~ "\u22655"
    ), levels = c("1", "2", "3", "4", "\u22655"))
  )

# Conteos por grupo y categoría
df_pir <- df_dist %>%
  count(genotipo, cat_ep) %>%
  group_by(genotipo) %>%
  mutate(
    pct       = 100 * n / sum(n),
    pct_signo = ifelse(genotipo == "Único", -pct, pct),
    etiqueta  = sprintf("%d (%.1f%%)", n, pct) %>%
      sub("\\.", ",", .)
  ) %>%
  ungroup()

# Mostrar conteos por consola para confirmar
cat("--- Distribución por número de episodios ---\n")
print(df_pir)

# ── 2. Paleta coherente con el resto del documento ──────────────────────────
cols_pir <- c("Único" = "#4CAF7D", "Mixto" = "#E07070")

n_unico <- sum(df_dist$genotipo == "Único")
n_mixto <- sum(df_dist$genotipo == "Mixto")

# ── 3. Figura pirámide ──────────────────────────────────────────────────────
fig_pir <- ggplot(df_pir, aes(x = pct_signo, y = cat_ep, fill = genotipo)) +
  geom_col(width = 0.7, color = "grey25", linewidth = 0.3) +
  geom_text(
    aes(label = etiqueta,
        hjust = ifelse(genotipo == "Único", 1.05, -0.05)),
    size = 3.4, color = "grey15"
  ) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.4) +
  scale_x_continuous(
    name   = "Porcentaje de pacientes (%)",
    breaks = seq(-50, 50, 10),
    labels = function(x) paste0(abs(x), "%"),
    limits = c(-55, 55),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = cols_pir,
                    labels = c(sprintf("Genotipos \u00fanicos (n = %d)", n_unico),
                               sprintf("Genotipos mixtos (n = %d)", n_mixto))) +
  labs(
    y    = "N episodios por paciente",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "bottom",
    legend.text        = element_text(size = 11),
    legend.key.size    = unit(12, "pt"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
    axis.text.y        = element_text(size = 12, color = "grey15", face = "bold"),
    axis.text.x        = element_text(size = 10, color = "grey40"),
    axis.title.x       = element_text(size = 11, color = "grey20",
                                      margin = margin(t = 6)),
    axis.title.y       = element_text(size = 11, color = "grey20",
                                      margin = margin(r = 8)),
    plot.margin        = margin(10, 14, 8, 8)
  )

fig_pir

# [guardado antiguo] ggsave("~/Desktop/Figuras_tesis/figura_R35bis_distribucion_episodios.png",
# [guardado antiguo]        fig_pir, width = 6, height = 3, dpi = 320, bg = "white")


df_test <- pacientes %>%
  filter(!is.na(nepi_tot), !is.na(gbghm)) %>%
  mutate(genotipo = factor(ifelse(gbghm == 1, "Mixto", "Único"),
                           levels = c("Único", "Mixto")))

# Mann-Whitney
test_mw <- wilcox.test(nepi_tot ~ genotipo, data = df_test, exact = FALSE)
cat("Mann-Whitney p =", format.pval(test_mw$p.value, digits = 3), "\n")

# Chi-cuadrado de la distribución categórica (≥4 vs <4 episodios)
df_test <- df_test %>% mutate(cat_4 = ifelse(nepi_tot >= 4, "≥4", "<4"))
test_chi <- chisq.test(table(df_test$genotipo, df_test$cat_4))
cat("Chi-cuadrado p =", format.pval(test_chi$p.value, digits = 3), "\n")
print(table(df_test$genotipo, df_test$cat_4))


# ── 1. AUCCV acumulada por paciente ──────────────────────────────────────────
# Limpieza previa de columnas si existen de ejecuciones anteriores
pacientes <- pacientes %>%
  select(-any_of(c("AUCCV_AC", "n_episodios_AUC",
                   "AUCCV_AC.x", "AUCCV_AC.y",
                   "n_episodios_AUC.x", "n_episodios_AUC.y")))

auc_paciente <- dnaemias %>%
  group_by(NHC) %>%
  summarise(
    AUCCV_AC        = sum(AUCCV, na.rm = TRUE),
    n_episodios_AUC = sum(!is.na(AUCCV)),
    .groups         = "drop"
  )

pacientes <- pacientes %>%
  left_join(auc_paciente, by = "NHC")

# ── 2. Datos para comparaciones ──────────────────────────────────────────────
df_comp <- pacientes %>%
  filter(!is.na(n_episodios_AUC), n_episodios_AUC > 0) %>%
  mutate(
    gbm_lab   = factor(ifelse(gbm   == 1, "Mixto", "Único"),
                       levels = c("Único", "Mixto")),
    ghm_lab   = factor(ifelse(ghm   == 1, "Mixto", "Único"),
                       levels = c("Único", "Mixto")),
    gbghm_lab = factor(ifelse(gbghm == 1, "Mixto", "Único"),
                       levels = c("Único", "Mixto"))
  )

# ── 3. Paletas por glucoproteína ─────────────────────────────────────────────
# gB: azul claro / azul oscuro
cols_fill_gb <- c("Único" = "#AED6F1", "Mixto" = "#1A5276")
cols_jit_gb  <- c("Único" = "#2980B9", "Mixto" = "#0D2B45")

# gH: rosa / rosa oscuro
cols_fill_gh <- c("Único" = "pink",    "Mixto" = "hotpink3")
cols_jit_gh  <- c("Único" = "pink3",   "Mixto" = "palevioletred4")

# gB+gH: verde / rojo
cols_fill_gbgh <- c("Único" = "#4CAF7D", "Mixto" = "#E07070")
cols_jit_gbgh  <- c("Único" = "#2E7D52", "Mixto" = "#B04040")

# ── 4. Función para construir cada panel ─────────────────────────────────────
# Límite Y fijo en 2000; asteriscos en vez de p
make_panel_paciente <- function(data, var_grupo, titulo,
                                cols_fill, cols_jit,
                                ylim_max = 1000) {
  
  d <- data %>%
    select(grupo = all_of(var_grupo), AUCCV_AC) %>%
    filter(!is.na(AUCCV_AC), !is.na(grupo))
  
  # Test sobre datos completos (sin filtrar por cap visual)
  stat <- d %>%
    rstatix::wilcox_test(AUCCV_AC ~ grupo) %>%
    rstatix::add_significance("p") %>%
    mutate(y.position = ylim_max * 0.95)
  
  # Datos con cap visual aplicado
  d_vis <- d %>% filter(AUCCV_AC <= ylim_max)
  
  ggplot(d_vis, aes(x = grupo, y = AUCCV_AC, fill = grupo)) +
    geom_boxplot(
      width         = 0.55,
      color         = "grey25",
      linewidth     = 0.45,
      outlier.shape = 21,
      outlier.size  = 1.8,
      outlier.fill  = "white",
      outlier.color = "grey50",
      outlier.alpha = 0.7
    ) +
    geom_jitter(
      aes(color = grupo),
      width = 0.12, size = 1.1, alpha = 0.45,
      show.legend = FALSE
    ) +
    scale_fill_manual(values  = cols_fill) +
    scale_color_manual(values = cols_jit) +
    stat_pvalue_manual(
      stat,
      label        = "p.signif",
      tip.length   = 0.01,
      bracket.size = 0.45,
      size         = 5,
      inherit.aes  = FALSE
    ) +
    labs(x = NULL,
         y = "AUCCV acumulada (UI/mL\u00b7d\u00eda)",
         title = titulo) +
    scale_y_continuous(
      limits = c(NA, ylim_max),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position    = "none",
      panel.border       = element_blank(),
      axis.line.x        = element_line(color = "grey60", linewidth = 0.4),
      axis.line.y        = element_line(color = "grey60", linewidth = 0.4),
      axis.ticks         = element_line(color = "grey60", linewidth = 0.3),
      axis.text.x        = element_text(size = 11, color = "grey25"),
      axis.text.y        = element_text(size = 10, color = "grey30"),
      axis.title.y       = element_text(size = 10, color = "grey20"),
      plot.title         = element_text(size = 13, color = "grey15",
                                        face = "bold", hjust = 0.5),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.margin        = margin(8, 10, 6, 6)
    )
}

# ── 5. Construir los 3 paneles ───────────────────────────────────────────────
p_gb <- make_panel_paciente(
  df_comp, "gbm_lab", "gB",
  cols_fill = cols_fill_gb, cols_jit = cols_jit_gb
)

p_gh <- make_panel_paciente(
  df_comp, "ghm_lab", "gH",
  cols_fill = cols_fill_gh, cols_jit = cols_jit_gh
)

p_gbgh <- make_panel_paciente(
  df_comp, "gbghm_lab", "gB/gH",
  cols_fill = cols_fill_gbgh, cols_jit = cols_jit_gbgh
)

# ── 6. Combinar y mostrar ────────────────────────────────────────────────────
figura_AUCCV_AC <- p_gb | p_gh | p_gbgh

figura_AUCCV_AC

# ── 7. Guardar ───────────────────────────────────────────────────────────────
# ggsave("~/Desktop/Figuras_tesis/figura_R22_AUCCV_AC_mixto.png",
#        figura_AUCCV_AC,
#        width = 10, height = 4, dpi = 320, bg = "white")

guardar_fig(fig_pir, "Figura_34", ancho = 15.2, alto = 7.6)

## ---- Figura 35. Tasa de episodios de DNAemia en pacientes con genotipos únicos y mixtos ----
# Incluye los modelos de regresión binomial negativa (IRR) citados en el texto.
# Resolver conflictos
select <- dplyr::select
if (exists("offset")) rm(offset)

# ── 1. Preparar datos IRR ─────────────────────────────────────────────────────
df_irr <- pacientes %>%
  dplyr::mutate(
    Ftx_dt       = as.Date(Ftx),
    Ffinseg_dt   = as.Date(Ffinseg),
    persona_anos = as.numeric(Ffinseg_dt - Ftx_dt) / 365.25,
    risk_f       = factor(risk, levels = c("MEDIO", "ALTO")),
    gbghm_grupo  = factor(gbghm, levels = c(0,1),
                          labels = c("Genotipos únicos",
                                     "Genotipos mixtos"))
  ) %>%
  dplyr::filter(persona_anos > 0)

# ── 2. Estadísticos descriptivos ──────────────────────────────────────────────
cat("── Total episodios CMV ──\n")
cat(sprintf("n = %d\n", sum(pacientes$nepi_tot, na.rm = TRUE)))

cat("\n── Tasas por grupo ──\n")
df_irr %>%
  dplyr::group_by(gbghm_grupo) %>%
  dplyr::summarise(
    n               = n(),
    total_episodios = sum(nepi_tot, na.rm = TRUE),
    total_pa        = round(sum(persona_anos, na.rm = TRUE), 1),
    tasa_pa         = round(sum(nepi_tot) / sum(persona_anos), 3),
    media_ep        = round(mean(nepi_tot, na.rm = TRUE), 2),
    mediana_ep      = median(nepi_tot, na.rm = TRUE),
    .groups         = "drop"
  ) %>% print()

# ── 3. Sobredispersión ────────────────────────────────────────────────────────
cat("\n── Sobredispersión ──\n")
df_irr %>%
  dplyr::group_by(gbghm_grupo) %>%
  dplyr::summarise(
    media    = round(mean(nepi_tot), 2),
    varianza = round(var(nepi_tot), 2),
    ratio    = round(var(nepi_tot) / mean(nepi_tot), 2),
    .groups  = "drop"
  ) %>% print()

m_pois <- glm(nepi_tot ~ gbghm_grupo + offset(log(persona_anos)),
              family = poisson, data = df_irr)
cat(sprintf("Poisson deviance/df = %.2f\n",
            m_pois$deviance / m_pois$df.residual))

# ── 4. Modelos binomial negativa ──────────────────────────────────────────────
cat("\n── IRR crudo (binomial negativa) ──\n")
m_nb <- glm.nb(nepi_tot ~ gbghm_grupo + offset(log(persona_anos)),
               data = df_irr)
print(round(exp(cbind(IRR = coef(m_nb), confint(m_nb))), 3))
cat("\np-valor gbghm:\n")
print(round(summary(m_nb)$coefficients, 3))

cat("\n── IRR ajustado (binomial negativa) ──\n")
m_nb_adj <- glm.nb(nepi_tot ~ gbghm_grupo + edad_tx + risk_f +
                     offset(log(persona_anos)), data = df_irr)
print(round(exp(cbind(IRR = coef(m_nb_adj), confint(m_nb_adj))), 3))
cat("\nCoeficientes con p-valores:\n")
print(round(summary(m_nb_adj)$coefficients, 3))

select <- dplyr::select

# ── 1. Preparar datos base (a nivel paciente) ────────────────────────────────
df_irr3 <- pacientes %>%
  dplyr::mutate(
    Ftx_dt       = as.Date(Ftx),
    Ffinseg_dt   = as.Date(Ffinseg),
    persona_anos = as.numeric(Ffinseg_dt - Ftx_dt) / 365.25,
    risk_f       = factor(risk, levels = c("MEDIO", "ALTO"))
  ) %>%
  dplyr::filter(persona_anos > 0)

# ── 2. Funcion: IRR ajustado para una variable de grupo ──────────────────────
calc_irr <- function(data, var_grupo) {
  data$grupo <- data[[var_grupo]]
  m <- glm.nb(nepi_tot ~ grupo + edad_tx + risk_f + offset(log(persona_anos)),
              data = data)
  irr <- round(exp(coef(m)[2]), 2)
  ci  <- round(exp(confint.default(m)[2, ]), 2)   # Wald (rapido y estable)
  p   <- summary(m)$coefficients[2, 4]
  list(irr = irr, ci_lo = ci[1], ci_hi = ci[2], p = p)
}

# ── 3. Calcular IRR para las 3 glucoproteinas ────────────────────────────────
vars   <- c(gbm = "gbm", ghm = "ghm", gbghm = "gbghm")
labs_g <- c(gbm = "gB", ghm = "gH", gbghm = "gB/gH")

res_irr <- lapply(names(vars), function(v) {
  d <- df_irr3 %>%
    dplyr::mutate(grupo = factor(.data[[v]], levels = c(0, 1),
                                 labels = c("Únicos", "Mixtos"))) %>%
    dplyr::filter(!is.na(grupo))
  r <- calc_irr(d %>% dplyr::mutate(!!v := grupo), v)
  data.frame(
    glucoproteina = labs_g[[v]],
    irr = r$irr, ci_lo = r$ci_lo, ci_hi = r$ci_hi, p = r$p
  )
}) %>% bind_rows() %>%
  dplyr::mutate(
    glucoproteina = factor(glucoproteina, levels = c("gB", "gH", "gB/gH")),
    label_irr = sprintf("IRR = %.2f (%.2f\u2013%.2f)\np = %s",
                        irr, ci_lo, ci_hi,
                        ifelse(p < 0.001, "<0,001",
                               sub("\\.", ",", sprintf("%.3f", p))))
  )

cat("── IRR ajustados ──\n"); print(res_irr)

# ── 4. Datos en formato largo: paciente x glucoproteina ──────────────────────
df_long3 <- df_irr3 %>%
  dplyr::select(NHC, nepi_tot, persona_anos, gbm, ghm, gbghm) %>%
  tidyr::pivot_longer(c(gbm, ghm, gbghm),
                      names_to = "gen", values_to = "mixto") %>%
  dplyr::filter(!is.na(mixto)) %>%
  dplyr::mutate(
    glucoproteina = factor(gen, levels = c("gbm", "ghm", "gbghm"),
                           labels = c("gB", "gH", "gB/gH")),
    grupo     = factor(ifelse(mixto == 1, "Mixtos", "Únicos"),
                       levels = c("Únicos", "Mixtos")),
    cat_color = paste(glucoproteina,
                      ifelse(mixto == 1, "Mixto", "Unico"), sep = "_")
  )

# ── 5. Tasas (ep/paciente-año) por panel y grupo ─────────────────────────────
tasas3 <- df_long3 %>%
  dplyr::group_by(glucoproteina, grupo) %>%
  dplyr::summarise(
    tasa = round(sum(nepi_tot) / sum(persona_anos), 2),
    .groups = "drop"
  )

# ── 6. n por panel y grupo (para etiquetas del eje) ──────────────────────────
n_grupo <- df_long3 %>%
  dplyr::count(glucoproteina, grupo)

# ── 7. Paletas (gB azul | gH rosa | gB/gH verde-rojo) ────────────────────────
cols_fill <- c(
  "gB_Unico"    = "#AED6F1", "gB_Mixto"    = "#4682B4",
  "gH_Unico"    = "pink",    "gH_Mixto"    = "hotpink3",
  "gB/gH_Unico" = "darkseagreen3", "gB/gH_Mixto" = "#c0392b"
)
cols_jit <- c(
  "gB_Unico"    = "#2980B9",        "gB_Mixto"    = "#0D2B45",
  "gH_Unico"    = "palevioletred3", "gH_Mixto"    = "palevioletred4",
  "gB/gH_Unico" = "#698B69",        "gB/gH_Mixto" = "#922b21"
)

# ── 8. Figura ────────────────────────────────────────────────────────────────
figura_R35_triple <- ggplot(df_long3,
                            aes(x = grupo, y = nepi_tot)) +
  geom_boxplot(aes(fill = cat_color),
               width = 0.45, outlier.shape = NA, alpha = 0.8,
               color = "black", linewidth = 0.6) +
  geom_jitter(aes(color = cat_color), width = 0.15, size = 1.8,
              alpha = 0.6, shape = 16, show.legend = FALSE) +
  # Etiqueta de tasa sobre cada caja
  geom_label(data = tasas3,
             aes(x = grupo, y = 10,
                 label = paste0(gsub("\\.", ",", tasa))),
             fill = "white", color = "black", size = 4,
             label.size = 0.25, fontface = "bold") +
  # Barra de significacion + IRR por panel
  
  geom_text(data = res_irr,
            aes(x = 1.5, y = 13, label = label_irr),
            size = 4, color = "black", fontface = "italic",
            inherit.aes = FALSE) +
  facet_wrap(~ glucoproteina) +
  scale_fill_manual(values = cols_fill, guide = "none") +
  scale_color_manual(values = cols_jit, guide = "none") +
  scale_y_continuous(
    name   = "Número de episodios de DNAemia",
    breaks = seq(0, 12, 2),
    limits = c(0, 15),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    panel.spacing.x    = unit(14, "pt"),
    axis.title.y       = element_text(size = 13, color = "black",
                                      margin = margin(r = 10)),
    axis.text.x        = element_text(size = 12, color = "black"),
    axis.text.y        = element_text(size = 12, color = "black"),
    axis.line          = element_line(color = "black", linewidth = 0.4),
    axis.ticks         = element_line(color = "black", linewidth = 0.3),
    strip.background   = element_blank(),
    strip.text         = element_text(size = 14, color = "black",
                                      face = "bold",
                                      margin = margin(t = 2, b = 6)),
    panel.grid         = element_blank(),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    plot.margin        = margin(15, 16, 10, 12)
  )

print(figura_R35_triple)

# ════════════════════════════════════════════════════════════════════════════
#  DATOS NUMERICOS — Figura R35 triple (episodios CMV por gB, gH, gB/gH)
# ════════════════════════════════════════════════════════════════════════════

# ── 1. n, tasa, media, mediana (RIC) por panel y grupo ───────────────────────
cat("\n========== DESCRIPTIVO POR GRUPO ==========\n")
tabla_desc <- df_long3 %>%
  dplyr::group_by(glucoproteina, grupo) %>%
  dplyr::summarise(
    n          = dplyr::n(),
    episodios  = sum(nepi_tot),
    persona_a  = round(sum(persona_anos), 1),
    tasa_pa    = round(sum(nepi_tot) / sum(persona_anos), 2),
    media      = round(mean(nepi_tot), 2),
    mediana    = median(nepi_tot),
    q1         = quantile(nepi_tot, 0.25),
    q3         = quantile(nepi_tot, 0.75),
    .groups    = "drop"
  ) %>%
  dplyr::mutate(
    `mediana (RIC)` = sprintf("%g (%g\u2013%g)", mediana, q1, q3)
  )
print(tabla_desc, n = Inf)

# ── 2. IRR ajustados con IC y p (formato tesis) ──────────────────────────────
cat("\n========== IRR AJUSTADOS (NB: edad + risk) ==========\n")
res_irr %>%
  dplyr::mutate(
    `IRR (IC95%)` = sprintf("%.2f (%.2f\u2013%.2f)", irr, ci_lo, ci_hi),
    `p`           = ifelse(p < 0.001, "<0,001",
                           sub("\\.", ",", sprintf("%.3f", p)))
  ) %>%
  dplyr::select(glucoproteina, `IRR (IC95%)`, `p`) %>%
  print()

# ── 3. p crudas de comparacion de distribucion (Wilcoxon) por panel ──────────
#     (util si quieres contrastar tasa vs distribucion de episodios)
cat("\n========== WILCOXON nepi_tot ~ grupo (crudo) ==========\n")
df_long3 %>%
  dplyr::group_by(glucoproteina) %>%
  dplyr::summarise(
    p_raw = wilcox.test(nepi_tot ~ grupo, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_fmt = ifelse(p_raw < 0.001, "<0,001",
                   sub("\\.", ",", sprintf("%.3f", p_raw)))
  ) %>%
  print()

guardar_fig(figura_R35_triple, "Figura_35", ancho = 18, alto = 10)

## ---- Figura 36. AUC acumulada de exposición viral por paciente ----
# ── 1. Episodios: todos los de dnaemias, excluyendo los sin aclaramiento ──────
n_inicial <- nrow(dnaemias)

df_ep <- dnaemias %>%
  filter(!is.na(aclaramiento_dias), !is.na(AUCCV))

n_excluidos <- n_inicial - nrow(df_ep)

cat(sprintf("Episodios totales en dnaemias: %d\n", n_inicial))
cat(sprintf("Excluidos (sin aclaramiento/AUC): %d\n", n_excluidos))
cat(sprintf("Episodios analizados: %d\n", nrow(df_ep)))

# ── 2. AUC acumulada por paciente (suma de AUCCV de sus episodios) ────────────
auc_pac <- df_ep %>%
  group_by(NHC) %>%
  summarise(AUCa = sum(AUCCV, na.rm = TRUE), .groups = "drop") %>%
  left_join(pacientes %>% select(NHC, gbm, ghm, gbghm), by = "NHC")

cat(sprintf("Pacientes con AUC acumulada: %d\n", nrow(auc_pac)))

# ── 3. Formato largo: una fila por paciente x glucoproteina ──────────────────
df_long <- auc_pac %>%
  pivot_longer(c(gbm, ghm, gbghm),
               names_to = "gen", values_to = "mixto") %>%
  filter(!is.na(mixto)) %>%
  mutate(
    glucoproteina = factor(gen,
                           levels = c("gbm", "ghm", "gbghm"),
                           labels = c("gB", "gH", "gB/gH")),
    genotipo  = factor(ifelse(mixto == 1, "Mixtos", "\u00danicos"),
                       levels = c("\u00danicos", "Mixtos")),
    cat_color = paste(glucoproteina, ifelse(mixto == 1, "Mixto", "\u00danico"),
                      sep = "_")
  )

# ── 4. Medianas, RIC y test de Wilcoxon por glucoproteina ────────────────────
resumen <- df_long %>%
  group_by(glucoproteina, genotipo) %>%
  summarise(
    n       = n(),
    mediana = round(median(AUCa)),
    q1      = round(quantile(AUCa, 0.25)),
    q3      = round(quantile(AUCa, 0.75)),
    .groups = "drop"
  ) %>%
  mutate(resumen = sprintf("%d (RIC %d\u2013%d)", mediana, q1, q3))

cat("\n── Medianas (RIC) por grupo ──\n"); print(resumen)

stat_test <- df_long %>%
  group_by(glucoproteina) %>%
  rstatix::wilcox_test(AUCa ~ genotipo) %>%
  rstatix::add_significance("p")

cat("\n── Test de Wilcoxon ──\n"); print(stat_test)

# ── 5. Paletas (gB azul | gH rosa | gB/gH verde-rojo) ────────────────────────
cols_fill <- c(
  "gB_\u00danico"    = "#AED6F1", "gB_Mixto"    = "#4682B4",
  "gH_\u00danico"    = "pink",    "gH_Mixto"    = "hotpink3",
  "gB/gH_\u00danico" = "#4CAF7D", "gB/gH_Mixto" = "#E07070"
)
cols_jit <- c(
  "gB_\u00danico"    = "#2980B9",        "gB_Mixto"    = "#0D2B45",
  "gH_\u00danico"    = "palevioletred3", "gH_Mixto"    = "palevioletred4",
  "gB/gH_\u00danico" = "#2E7D52",        "gB/gH_Mixto" = "#B04040"
)

# ── 6. Límite del eje y posición de los corchetes ────────────────────────────
y_max <- 760
stat_test <- stat_test %>% mutate(y.position = 710)

# ── 7. Figura ────────────────────────────────────────────────────────────────
fig_auc <- ggplot(df_long, aes(x = genotipo, y = AUCa)) +
  geom_boxplot(aes(fill = cat_color),
               width = 0.45, color = "black", linewidth = 0.4,
               outlier.shape = NA, alpha = 0.85) +
  geom_jitter(aes(color = cat_color),
              width = 0.15, size = 1.8, alpha = 0.5,
              shape = 16, show.legend = FALSE) +
  facet_wrap(~ glucoproteina) +
  stat_pvalue_manual(
    stat_test,
    label        = "p.signif",
    tip.length   = 0.01,
    bracket.size = 0.6,
    size         = 5,
    color        = "black"
  ) +
  scale_fill_manual(values = cols_fill, guide = "none") +
  scale_color_manual(values = cols_jit, guide = "none") +
  scale_y_continuous(
    name   = expression("AUCa (log"[10]*"(UI/mL)\u00b7d\u00eda)"),
    breaks = seq(0, 700, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, y_max)) +
  labs(x = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position    = "none",
    panel.spacing.x    = unit(14, "pt"),
    axis.title.y       = element_text(size = 13, color = "black",
                                      margin = margin(r = 10)),
    axis.text.x        = element_text(size = 13, color = "black"),
    axis.text.y        = element_text(size = 12, color = "black"),
    axis.line          = element_line(color = "black", linewidth = 0.4),
    axis.ticks         = element_line(color = "black", linewidth = 0.3),
    strip.background   = element_blank(),
    strip.text         = element_text(size = 14, color = "black",
                                      face = "bold",
                                      margin = margin(t = 2, b = 6)),
    panel.grid         = element_blank(),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    plot.margin        = margin(15, 16, 10, 12)
  )

print(fig_auc)

cat("Pacientes con AUCa >", y_max, "(no representados):",
    sum(auc_pac$AUCa > y_max), "\n")

# ── 8. Guardar en SVG ─────────────────────────────────────────────────────────
# [guardado antiguo] ggsave(
# [guardado antiguo]   filename      = "~/Desktop/Figuras_tesis/Figura_36_AUCa.svg",
# [guardado antiguo]   plot          = fig_auc,
# [guardado antiguo]   device        = svglite::svglite,
# [guardado antiguo]   fix_text_size = FALSE,
# [guardado antiguo]   width         = 16,
# [guardado antiguo]   height        = 12,
# [guardado antiguo]   units         = "cm",
# [guardado antiguo]   bg            = "white"
# [guardado antiguo] )

## Requiere auc_pac del bloque [c] y la variable Ffinseg de pacientes

auc_anual <- auc_pac %>%
  left_join(pacientes %>% select(NHC, Ftx, Ffinseg), by = "NHC") %>%
  mutate(
    seg_anos  = as.numeric(as.Date(Ffinseg) - as.Date(Ftx)) / 365.25,
    AUCa_anio = AUCa / seg_anos
  ) %>%
  pivot_longer(c(gbm, ghm, gbghm), names_to = "gen", values_to = "mixto") %>%
  filter(!is.na(mixto)) %>%
  mutate(gen = factor(gen, levels = c("gbm", "ghm", "gbghm"),
                      labels = c("gB", "gH", "gB/gH")))

fmt_ric <- function(x, d = 0) sprintf(paste0("%.", d, "f (%.", d, "f\u2013%.", d, "f)"),
                                      median(x), quantile(x, 0.25), quantile(x, 0.75))

auc_anual %>%
  group_by(gen) %>%
  summarise(
    n_unico    = sum(mixto == 0),
    n_mixto    = sum(mixto == 1),
    AUCa_unico = fmt_ric(AUCa[mixto == 0]),
    AUCa_mixto = fmt_ric(AUCa[mixto == 1]),
    p_AUCa     = wilcox.test(AUCa ~ mixto, exact = FALSE)$p.value,
    anio_unico = fmt_ric(AUCa_anio[mixto == 0]),
    anio_mixto = fmt_ric(AUCa_anio[mixto == 1]),
    p_anio     = wilcox.test(AUCa_anio ~ mixto, exact = FALSE)$p.value,
    seg_unico  = fmt_ric(seg_anos[mixto == 0], 1),
    seg_mixto  = fmt_ric(seg_anos[mixto == 1], 1),
    p_seg      = wilcox.test(seg_anos ~ mixto, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(across(starts_with("p_"), ~ round(.x, 3))) %>%
  as.data.frame() %>%
  print()

cat("Pacientes con AUCa > 760 (no representados en la figura):",
    sum(auc_pac$AUCa > 760), "\n")

guardar_fig(fig_auc, "Figura_36", ancho = 16, alto = 12)

## ---- Tabla 20. Recurrencia de episodios según el patrón genotípico ----


# ── 1. Identificar primer episodio por paciente ─────────────────────────────
# Definición: el episodio con Fini más cercano a Ftx (mínimo Fini por NHC)

primer_episodio <- dnaemias %>%
  mutate(
    Fini = as.Date(Fini),
    Ffin = as.Date(Ffin)
  ) %>%
  group_by(NHC) %>%
  arrange(Fini, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    Ftx              = as.Date(Ftx),
    dias_tx_a_inicio = as.numeric(Fini - Ftx)
  )

# ── 2. Excluir pacientes sin Ffin documentado ────────────────────────────────
# (4 fallecidos antes de aclaramiento + 3 trasladados = 7 excluidos)
primer_ep_evaluables <- primer_episodio %>%
  filter(!is.na(Ffin))

cat(sprintf("── Pacientes evaluables: %d ──\n", nrow(primer_ep_evaluables)))

# ── 3. Identificar segundo episodio por paciente ────────────────────────────
segundo_episodio <- dnaemias %>%
  mutate(
    Fini = as.Date(Fini),
    Ffin = as.Date(Ffin)
  ) %>%
  inner_join(
    primer_ep_evaluables %>% select(NHC, Ffin_primer = Ffin),
    by = "NHC"
  ) %>%
  filter(Fini > Ffin_primer) %>%
  group_by(NHC) %>%
  arrange(Fini, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(NHC,
         Fini_segundo    = Fini,
         Ffin_segundo    = Ffin,
         gbghmix_segundo = gbghmix,
         gbmix_segundo   = gbmix,
         ghmix_segundo   = ghmix,
         cvini_segundo   = cvini,
         clinica_segundo = clinica)

# ── 4. Construir base df_recur ──────────────────────────────────────────────
df_recur <- primer_ep_evaluables %>%
  select(NHC,
         Ftx,
         Fini_primer    = Fini,
         Ffin_primer    = Ffin,
         gbghmix_primer = gbghmix,
         gbmix_primer   = gbmix,
         ghmix_primer   = ghmix,
         cvini_primer   = cvini,
         cvpico_primer  = cvpico,
         dias_tx_a_inicio) %>%
  left_join(segundo_episodio, by = "NHC") %>%
  left_join(
    pacientes %>%
      mutate(
        Ffinseg = as.Date(Ffinseg),
        Fexitus = as.Date(Fexitus)
      ) %>%
      select(NHC, Ffinseg, exitus, Fexitus, serDR),
    by = "NHC"
  )

# ── 5. Calcular tiempos y eventos para análisis ──────────────────────────────
df_recur <- df_recur %>%
  mutate(
    fecha_fin_seg = pmin(Ffinseg, Fexitus, na.rm = TRUE),
    dias_a_recur  = as.numeric(Fini_segundo - Ffin_primer),
    dias_seguim   = as.numeric(fecha_fin_seg - Ffin_primer),
    
    # Recurrencia a 6 meses
    recur_180 = case_when(
      !is.na(dias_a_recur) & dias_a_recur <= 180 ~ 1L,
      TRUE                                       ~ 0L
    ),
    tiempo_km = case_when(
      recur_180 == 1     ~ dias_a_recur,
      dias_seguim < 180  ~ dias_seguim,
      TRUE               ~ 180
    ),
    evento_km = recur_180,
    
    # Recurrencia a 12 meses
    recur_365 = case_when(
      !is.na(dias_a_recur) & dias_a_recur <= 365 ~ 1L,
      TRUE                                       ~ 0L
    ),
    tiempo_km_365 = case_when(
      recur_365 == 1     ~ dias_a_recur,
      dias_seguim < 365  ~ dias_seguim,
      TRUE               ~ 365
    ),
    evento_km_365 = recur_365,
    
    # Recurrencia en seguimiento completo
    recur_global = as.integer(!is.na(Fini_segundo)),
    tiempo_completo = case_when(
      recur_global == 1 ~ dias_a_recur,
      TRUE              ~ dias_seguim
    ),
    
    # Factores de exposición
    exp_grupo = factor(gbghmix_primer, levels = c(0, 1),
                       labels = c("Único", "Mixto")),
    gb_grupo  = factor(gbmix_primer,   levels = c(0, 1),
                       labels = c("gB único", "gB mixto")),
    gh_grupo  = factor(ghmix_primer,   levels = c(0, 1),
                       labels = c("gH único", "gH mixto")),
    serDR_f   = factor(serDR)
  )

# ── 5b. Excepción: aclaramiento posterior al cierre del seguimiento ──────────
# El paciente "NHC_anonimizado" se aclaró después del cierre administrativo
# del seguimiento (31/12/2022). Se incluye con tiempo 0 (censurado en el
# aclaramiento) en lugar de un tiempo negativo.
df_recur <- df_recur %>%
  mutate(across(c(tiempo_km, tiempo_km_365, tiempo_completo), ~ pmax(.x, 0)))

cat("Tiempos negativos tras la corrección:", sum(df_recur$tiempo_km_365 < 0), "\n")

# Objetos Surv
surv_180      <- Surv(time = df_recur$tiempo_km,     event = df_recur$evento_km)
surv_365      <- Surv(time = df_recur$tiempo_km_365, event = df_recur$evento_km_365)
surv_completo <- Surv(time = df_recur$tiempo_completo, event = df_recur$recur_global)

# ── 6 meses ─────────────────────────────────────────────────────────────────
cox_gbgh_6m     <- coxph(surv_180 ~ exp_grupo, data = df_recur)
cox_gb_6m       <- coxph(surv_180 ~ gb_grupo,  data = df_recur)
cox_gh_6m       <- coxph(surv_180 ~ gh_grupo,  data = df_recur)
cox_gbgh_6m_adj <- coxph(surv_180 ~ exp_grupo + serDR_f + cvpico_primer,
                         data = df_recur)
cox_gb_6m_adj   <- coxph(surv_180 ~ gb_grupo  + serDR_f + cvpico_primer,
                         data = df_recur)
cox_gh_6m_adj   <- coxph(surv_180 ~ gh_grupo  + serDR_f + cvpico_primer,
                         data = df_recur)

# ── 12 meses ────────────────────────────────────────────────────────────────
cox_gbgh_12m     <- coxph(surv_365 ~ exp_grupo, data = df_recur)
cox_gb_12m       <- coxph(surv_365 ~ gb_grupo,  data = df_recur)
cox_gh_12m       <- coxph(surv_365 ~ gh_grupo,  data = df_recur)
cox_gbgh_12m_adj <- coxph(surv_365 ~ exp_grupo + serDR_f + cvpico_primer,
                          data = df_recur)
cox_gb_12m_adj   <- coxph(surv_365 ~ gb_grupo  + serDR_f + cvpico_primer,
                          data = df_recur)
cox_gh_12m_adj   <- coxph(surv_365 ~ gh_grupo  + serDR_f + cvpico_primer,
                          data = df_recur)

# ── Seguimiento completo ────────────────────────────────────────────────────
cox_gbgh_sc     <- coxph(surv_completo ~ exp_grupo, data = df_recur)
cox_gb_sc       <- coxph(surv_completo ~ gb_grupo,  data = df_recur)
cox_gh_sc       <- coxph(surv_completo ~ gh_grupo,  data = df_recur)
cox_gbgh_sc_adj <- coxph(surv_completo ~ exp_grupo + serDR_f + cvpico_primer,
                         data = df_recur)
cox_gb_sc_adj   <- coxph(surv_completo ~ gb_grupo  + serDR_f + cvpico_primer,
                         data = df_recur)
cox_gh_sc_adj   <- coxph(surv_completo ~ gh_grupo  + serDR_f + cvpico_primer,
                         data = df_recur)

# ── Funciones auxiliares ────────────────────────────────────────────────────
extraer_hr_fmt <- function(cox) {
  hr <- exp(coef(cox)[1])
  ci <- exp(confint(cox)[1, ])
  p  <- summary(cox)$coefficients[1, 5]
  list(
    hr_ci = sprintf("%.2f (%.2f\u2013%.2f)", hr, ci[1], ci[2]),
    p     = ifelse(p < 0.001, "<0,001",
                   sub("\\.", ",", sprintf("%.3f", p)))
  )
}

fmt_n_pct <- function(n, N) {
  sub("\\.", ",", sprintf("%d (%.1f%%)", n, 100 * n / N))
}

# ── Construir filas ─────────────────────────────────────────────────────────
filas_R13 <- bind_rows(
  
  # Sección 1: 6 meses
  tibble(variable = "RECURRENCIA A 6 MESES (n = 33)",
         unicos = "", mixtos = "", hr_crudo = "", p_crudo = "",
         hr_adj = "", p_adj = ""),
  
  tibble(variable = "gB+gH",
         unicos   = fmt_n_pct(21, 99),
         mixtos   = fmt_n_pct(12, 40),
         hr_crudo = extraer_hr_fmt(cox_gbgh_6m)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gbgh_6m)$p,
         hr_adj   = extraer_hr_fmt(cox_gbgh_6m_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gbgh_6m_adj)$p),
  
  tibble(variable = "gB",
         unicos   = fmt_n_pct(25, 108),
         mixtos   = fmt_n_pct(8, 31),
         hr_crudo = extraer_hr_fmt(cox_gb_6m)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gb_6m)$p,
         hr_adj   = extraer_hr_fmt(cox_gb_6m_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gb_6m_adj)$p),
  
  tibble(variable = "gH",
         unicos   = fmt_n_pct(24, 112),
         mixtos   = fmt_n_pct(9, 27),
         hr_crudo = extraer_hr_fmt(cox_gh_6m)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gh_6m)$p,
         hr_adj   = extraer_hr_fmt(cox_gh_6m_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gh_6m_adj)$p),
  
  # Sección 2: 12 meses
  tibble(variable = "RECURRENCIA A 12 MESES (n = 42)",
         unicos = "", mixtos = "", hr_crudo = "", p_crudo = "",
         hr_adj = "", p_adj = ""),
  
  tibble(variable = "gB+gH",
         unicos   = fmt_n_pct(25, 99),
         mixtos   = fmt_n_pct(17, 40),
         hr_crudo = extraer_hr_fmt(cox_gbgh_12m)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gbgh_12m)$p,
         hr_adj   = extraer_hr_fmt(cox_gbgh_12m_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gbgh_12m_adj)$p),
  
  tibble(variable = "gB",
         unicos   = fmt_n_pct(30, 108),
         mixtos   = fmt_n_pct(12, 31),
         hr_crudo = extraer_hr_fmt(cox_gb_12m)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gb_12m)$p,
         hr_adj   = extraer_hr_fmt(cox_gb_12m_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gb_12m_adj)$p),
  
  tibble(variable = "gH",
         unicos   = fmt_n_pct(29, 112),
         mixtos   = fmt_n_pct(13, 27),
         hr_crudo = extraer_hr_fmt(cox_gh_12m)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gh_12m)$p,
         hr_adj   = extraer_hr_fmt(cox_gh_12m_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gh_12m_adj)$p),
  
  # Sección 3: Seguimiento completo
  tibble(variable = "RECURRENCIA EN SEGUIMIENTO COMPLETO (n = 56)",
         unicos = "", mixtos = "", hr_crudo = "", p_crudo = "",
         hr_adj = "", p_adj = ""),
  
  tibble(variable = "gB+gH",
         unicos   = fmt_n_pct(35, 99),
         mixtos   = fmt_n_pct(21, 40),
         hr_crudo = extraer_hr_fmt(cox_gbgh_sc)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gbgh_sc)$p,
         hr_adj   = extraer_hr_fmt(cox_gbgh_sc_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gbgh_sc_adj)$p),
  
  tibble(variable = "gB",
         unicos   = fmt_n_pct(42, 108),
         mixtos   = fmt_n_pct(14, 31),
         hr_crudo = extraer_hr_fmt(cox_gb_sc)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gb_sc)$p,
         hr_adj   = extraer_hr_fmt(cox_gb_sc_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gb_sc_adj)$p),
  
  tibble(variable = "gH",
         unicos   = fmt_n_pct(40, 112),
         mixtos   = fmt_n_pct(16, 27),
         hr_crudo = extraer_hr_fmt(cox_gh_sc)$hr_ci,
         p_crudo  = extraer_hr_fmt(cox_gh_sc)$p,
         hr_adj   = extraer_hr_fmt(cox_gh_sc_adj)$hr_ci,
         p_adj    = extraer_hr_fmt(cox_gh_sc_adj)$p)
)

colnames(filas_R13) <- c("Genotipos",
                         "Únicos\nn (%)",
                         "Mixtos\nn (%)",
                         "HR cruda\n(IC 95%)",
                         "p",
                         "HR ajustada*\n(IC 95%)",
                         "p ")

# ── Filas especiales por índice ─────────────────────────────────────────────
filas_seccion_R13 <- c(1, 5, 9)

es_p_sig <- function(x) {
  num <- suppressWarnings(as.numeric(sub(",", ".", x)))
  !is.na(num) & num < 0.05
}

filas_p_sig_crudo <- which(es_p_sig(filas_R13[[5]]))
filas_p_sig_adj   <- which(es_p_sig(filas_R13[[7]]))

# ── Construir flextable ─────────────────────────────────────────────────────
ft_R13 <- flextable(filas_R13) %>%
  width(j = 1, width = 2.0) %>%
  width(j = 2, width = 1.4) %>%
  width(j = 3, width = 1.4) %>%
  width(j = 4, width = 1.8) %>%
  width(j = 5, width = 0.9) %>%
  width(j = 6, width = 1.8) %>%
  width(j = 7, width = 0.9) %>%
  align(align = "left", part = "all") %>%
  align(j = c(2, 3, 4, 5, 6, 7), align = "center", part = "all") %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header") %>%
  bg(bg = "white", part = "all") %>%
  bold(i = filas_seccion_R13, part = "body") %>%
  bg(i = filas_seccion_R13, bg = "#F0F0F0", part = "body") %>%
  border_remove() %>%
  hline_top(part = "header",
            border = fp_border(color = "black", width = 1.2)) %>%
  hline(part = "header",
        border = fp_border(color = "black", width = 0.8)) %>%
  hline_bottom(part = "body",
               border = fp_border(color = "black", width = 1.2)) %>%
  hline(i = c(4, 8),
        border = fp_border(color = "grey70", width = 0.4),
        part = "body") %>%
  padding(padding.top = 2, padding.bottom = 2,
          padding.left = 4, padding.right = 4, part = "all") %>%
  add_footer_lines(values = c(
    "*Modelo de regresión de Cox ajustado por perfil serológico D/R y carga viral pico del primer episodio.",
    "Las HR comparan el grupo mixto frente al grupo único, considerado como referencia.",
    "Valores de p en rojo indican diferencia estadísticamente significativa (p < 0,05)."
  )) %>%
  fontsize(size = 8, part = "footer") %>%
  color(color = "grey30", part = "footer") %>%
  set_caption(
    caption = "Tabla R13. Recurrencia de DNAemia por CMV según el patrón genotípico del primer episodio en tres ventanas temporales."
  )

if (length(filas_p_sig_crudo) > 0) {
  ft_R13 <- ft_R13 %>%
    color(i = filas_p_sig_crudo, j = 5, color = "firebrick3", part = "body") %>%
    bold( i = filas_p_sig_crudo, j = 5, part = "body")
}
if (length(filas_p_sig_adj) > 0) {
  ft_R13 <- ft_R13 %>%
    color(i = filas_p_sig_adj, j = 7, color = "firebrick3", part = "body") %>%
    bold( i = filas_p_sig_adj, j = 7, part = "body")
}

ft_R13 <- ft_R13 %>%
  set_table_properties(layout = "autofit")

ft_R13

guardar_tabla(ft_R13, "Tabla_20")

## ---- Figura 37. Probabilidad libre de recurrencia a 12 meses según el patrón genotípico ----
## Mismo estilo que las Figuras 41-44. Requiere df_recur de [7.12] [a].
## Valor p: modelo de Cox ajustado por perfil serológico D/R y CV pico (Tabla 20).

pal_rec <- c("Genotipos únicos" = "#2E8B57", "Genotipos mixtos" = "firebrick3")

km_recur <- function(var_grupo, titulo, ylab = "Libre de recurrencia",
                     ylim = c(0.4, 1), tiempos = seq(0, 12, 3), base = 12) {
  
  d <- df_recur %>%
    transmute(
      meses  = tiempo_km_365 * 12 / 365,   # 365 días = 12 meses exactos
      evento = evento_km_365,
      grupo  = factor(ifelse(.data[[var_grupo]] == 1, "Genotipos mixtos", "Genotipos únicos"),
                      levels = names(pal_rec)),
      serDR_f, cvpico_primer
    ) %>%
    filter(!is.na(grupo), !is.na(meses))
  
  # Valor p del modelo de Cox ajustado (mismo modelo que la Tabla 20)
  m <- coxph(Surv(meses, evento) ~ grupo + serDR_f + cvpico_primer, data = d)
  p <- summary(m)$coefficients["grupoGenotipos mixtos", "Pr(>|z|)"]
  etiqueta_p <- paste0("p = ", ifelse(p < 0.001, "<0,001", sub("\\.", ",", sprintf("%.3f", p))))
  
  # Curvas de Kaplan-Meier
  sf  <- survfit(Surv(meses, evento) ~ grupo, data = d)
  cur <- data.frame(t = sf$time, s = sf$surv,
                    grupo = rep(names(sf$strata), sf$strata)) %>%
    mutate(grupo = factor(sub("grupo=", "", grupo), levels = names(pal_rec)))
  cur <- bind_rows(data.frame(t = 0, s = 1, grupo = factor(names(pal_rec), levels = names(pal_rec))),
                   cur) %>%
    arrange(grupo, t)
  
  xmax    <- max(tiempos)
  ybreaks <- seq(ylim[1], ylim[2], by = if (diff(ylim) > 0.5) 0.2 else 0.1)
  
  g <- ggplot(cur, aes(t, s, color = grupo)) +
    geom_step(linewidth = 1.1, lineend = "round") +
    scale_color_manual(values = pal_rec, name = NULL) +
    scale_x_continuous(breaks = tiempos, limits = c(0, xmax),
                       expand = expansion(mult = c(0.08, 0.03))) +
    scale_y_continuous(limits = ylim, breaks = ybreaks,
                       labels = function(x) sub("\\.", ",", x),
                       expand = expansion(mult = c(0, 0.02))) +
    annotate("text", x = xmax, y = ylim[2] - 0.02 * diff(ylim), label = etiqueta_p,
             hjust = 1, vjust = 1, size = base * 0.3) +
    labs(x = "Meses desde el aclaramiento", y = ylab, title = titulo) +
    theme_classic(base_size = base) +
    theme(plot.title         = element_text(size = base, face = "bold", hjust = 0.5),
          legend.position    = "none",
          panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
          axis.line          = element_line(color = "black", linewidth = 0.5),
          axis.ticks         = element_line(color = "black", linewidth = 0.5),
          axis.ticks.length  = unit(3, "pt"),
          axis.text          = element_text(color = "black", size = base * 0.95),
          axis.title         = element_text(size = base),
          plot.margin        = margin(4, 8, 2, 4))
  
  # Pacientes en riesgo (tiempo >= t)
  riesgo <- lapply(tiempos, function(tt) {
    d %>% filter(meses >= tt) %>% count(grupo, .drop = FALSE) %>% mutate(t = tt)
  }) %>% bind_rows() %>%
    complete(grupo, t = tiempos, fill = list(n = 0))
  
  tab <- ggplot(riesgo, aes(t, grupo, label = n, color = grupo)) +
    geom_text(size = base * 0.32) +
    scale_color_manual(values = pal_rec) +
    scale_x_continuous(limits = c(0, xmax), breaks = tiempos,
                       expand = expansion(mult = c(0.08, 0.03))) +
    scale_y_discrete(limits = rev(names(pal_rec)), labels = c("Mixtos", "Únicos")) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL, title = "Pacientes en riesgo") +
    theme_void(base_size = base) +
    theme(plot.title      = element_text(size = base * 0.8, hjust = 0),
          legend.position = "none",
          axis.text.y     = element_text(size = base * 0.8, hjust = 1, color = "black",
                                         margin = margin(r = 6)),
          plot.margin     = margin(0, 8, 4, 4))
  
  cat(sprintf("%s: n = %d, eventos = %d, %s\n", titulo, nrow(d), sum(d$evento), etiqueta_p))
  (g / tab) + plot_layout(heights = c(4, 1))
}

leyenda_rec <- function(base = 12) {
  g <- ggplot(data.frame(x = c(1, 2, 1, 2), y = c(1, 1, 2, 2),
                         grupo = factor(rep(names(pal_rec), each = 2), levels = names(pal_rec))),
              aes(x, y, color = grupo, group = grupo)) +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = pal_rec, name = NULL) +
    guides(color = guide_legend(ncol = 1)) +
    theme_void() +
    theme(legend.position  = "bottom",
          legend.text      = element_text(size = base * 0.95),
          legend.key.width = unit(24, "pt"))
  patchwork::wrap_elements(full = cowplot::get_legend(g))
}

b_gb   <- km_recur("gbmix_primer",   "gB")
b_gh   <- km_recur("ghmix_primer",   "gH")
b_gbgh <- km_recur("gbghmix_primer", "gB/gH")

fig_km_recur12 <- (b_gb | b_gh) / (b_gbgh | leyenda_rec()) +
  plot_layout(widths = c(1, 1), heights = c(1, 1)) &
  theme(plot.background = element_rect(fill = "white", color = NA))

print(fig_km_recur12)

## Exportar a SVG para una página A4 en Word (16 x 22 cm)
# [guardado antiguo] ggsave("~/Desktop/Figuras_tesis/Figura_37_KM_recurrencia_12m.svg", fig_km_recur12,
# [guardado antiguo]        device = svglite::svglite, fix_text_size = FALSE,
# [guardado antiguo]        width = 16, height = 20, units = "cm", bg = "white")

guardar_fig(fig_km_recur12, "Figura_37", ancho = 16, alto = 20)

## ---- Preparación común: modelos con exposición dependiente del tiempo y curvas KM extendidas ----
## EJECUTAR ESTA SECCIÓN ANTES que [7.4], [7.5], [7.6], [7.9] y [7.11].
##
## Requiere `pacientes` y `dnaemias` ya cargados (mismas hojas que el resto del script).
## Compara, para cada clasificación (gB, gH, gB/gH) y desenlace:
##   1. Fijo          : modelo actual de la tesis (exposición fija desde el trasplante)
##   2. TD            : exposición dependiente del tiempo (0 hasta la primera DNAemia mixta, 1 después)
##   3. TD + desfase  : la exposición empieza a contar 6 meses después de la primera mezcla
##   4. TD + n_epi    : TD ajustado por número acumulado de episodios genotipados (capacidad de detección)
##   5. Landmark 1 a  : estado de exposición fijado a 1 año; se excluyen eventos y seguimiento previos
##   6. Landmark 2 a  : ídem a 2 años
##   7. Fijo, >=2 epi : modelo fijo restringido a pacientes con >=2 episodios genotipados
##   8. Tres estados  : sin genotipar / unico demostrado / mixto, todos dependientes del tiempo
##   9. Entrada diferida: cada paciente entra en riesgo en su primer episodio genotipado


fmt_num <- function(x, d = 2) sub("\\.", ",", sprintf(paste0("%.", d, "f"), x))
fmt_p   <- function(p) ifelse(p < 0.001, "<0,001", fmt_num(p, 3))

# ── 0. Fechas de primera DNAemia mixta por clasificación ─────────────────────
mix_flags <- c(gbm = "gbmix", ghm = "ghmix", gbghm = "gbghmix")
stopifnot(all(mix_flags %in% names(dnaemias)))

epis <- dnaemias %>%
  mutate(Fini = as.Date(Fini)) %>%
  filter(!is.na(Fini)) %>%
  mutate(genotipado = (!is.na(gbsum) & gbsum >= 1) | (!is.na(ghsum) & ghsum >= 1),
         valido_gbm   = !is.na(gbsum) & gbsum >= 1,
         valido_ghm   = !is.na(ghsum) & ghsum >= 1,
         valido_gbghm = valido_gbm & valido_ghm)

primera_mix <- Reduce(function(a, b) full_join(a, b, by = "NHC"),
                      lapply(names(mix_flags), function(v) {
                        epis %>%
                          filter(.data[[mix_flags[[v]]]] == 1) %>%
                          group_by(NHC) %>%
                          summarise(!!paste0("Fmix_", v) := min(Fini), .groups = "drop")
                      }))

n_geno <- epis %>% filter(genotipado) %>% count(NHC, name = "n_geno")

# ── 1. Base de pacientes ──────────────────────────────────────────────────────
base <- pacientes %>%
  transmute(
    NHC,
    Ftx     = as.Date(Ftx),
    Ffinseg = as.Date(Ffinseg),
    edad_tx,
    risk_f  = factor(risk, levels = c("MEDIO", "ALTO")),
    gbm, ghm, gbghm,
    F_onco  = as.Date(F_onco),  onco,
    solido  = as.integer(onco == 1 & onco_tipo == "SOL"),
    F_RA1   = as.Date(F_RA1),   RA,
    F_DCIP  = as.Date(F_DCIP),  DCIP,
    Fexitus = as.Date(Fexitus), exitus,
    F_opo   = as.Date(F_opo),   inf_opo,
    serDR   = factor(serDR)
  ) %>%
  left_join(primera_mix, by = "NHC") %>%
  left_join(n_geno,      by = "NHC") %>%
  mutate(n_geno = ifelse(is.na(n_geno), 0L, n_geno))


# Excepción: el paciente "NHC_anonimizado" se aclaró después del cierre
# administrativo del seguimiento (31/12/2022). Se incluye con tiempo 0
# (censurado en el aclaramiento) en lugar de un tiempo negativo.
df_recur <- df_recur %>%
  mutate(across(c(tiempo_km, tiempo_km_365, tiempo_completo), ~ pmax(.x, 0)))
# Comprobación de coherencia: etiqueta fija vs existencia de fecha de mezcla
cat("── Coherencia etiqueta fija / fecha de primera mezcla ──\n")
for (v in names(mix_flags)) {
  cat(v, ":\n")
  print(table(fijo = base[[v]], con_fecha = !is.na(base[[paste0("Fmix_", v)]]),
              useNA = "ifany"))
}

# ── 2. Función genérica ───────────────────────────────────────────────────────
analisis_td <- function(clasif, flag_evento, fecha_evento, etiqueta,
                        covars = c("edad_tx", "risk_f"),
                        lag_anios = 0.5, landmarks = c(1, 2),
                        solo_principal = TRUE) {
  ## solo_principal = TRUE imprime solo el modelo elegido (entrada diferida);
  ## FALSE imprime ademas todos los modelos de sensibilidad.
  
  aj <- paste(covars, collapse = " + ")
  f  <- function(lhs, expo) as.formula(paste(lhs, "~", expo, "+", aj))
  
  d <- base %>%
    mutate(
      ev    = as.integer(!is.na(.data[[flag_evento]]) & .data[[flag_evento]] == 1),
      F_fin = if_else(ev == 1, .data[[fecha_evento]], Ffinseg),
      t_fin = as.numeric(F_fin - Ftx) / 365.25,
      t_mix = as.numeric(.data[[paste0("Fmix_", clasif)]] - Ftx) / 365.25,
      fijo  = .data[[clasif]]
    ) %>%
    filter(!is.na(fijo), !is.na(t_fin))
  n_t0 <- sum(d$t_fin <= 0)
  d <- d %>% filter(t_fin > 0) %>%
    select(NHC, all_of(covars), fijo, ev, t_fin, t_mix, n_geno)
  
  extrae <- function(m, term) {
    s  <- summary(m)$coefficients
    ci <- exp(confint(m))
    data.frame(HR = exp(s[term, "coef"]), lo = ci[term, 1], hi = ci[term, 2],
               p = s[term, ncol(s)])
  }
  
  fila <- function(nombre, m, term, n, n_ev, n_ev_exp) {
    r <- extrae(m, term)
    data.frame(Modelo = nombre, n = n, eventos = n_ev, eventos_expuesto = n_ev_exp,
               HR = fmt_num(r$HR), IC95 = paste0(fmt_num(r$lo), " a ", fmt_num(r$hi)),
               p = fmt_p(r$p), stringsAsFactors = FALSE)
  }
  
  # 1. Fijo (reproduce la tesis)
  m_fijo <- coxph(f("Surv(t_fin, ev)", "fijo"), data = d)
  res <- fila("Fijo desde trasplante (tesis)", m_fijo, "fijo",
              nrow(d), sum(d$ev), sum(d$ev[d$fijo == 1]))
  
  # 2. Dependiente del tiempo
  d_id <- d %>% select(NHC, all_of(covars), fijo, t_fin, ev)
  d_td <- tmerge(d_id, d_id, id = NHC, ev_td = event(t_fin, ev))
  d_td <- tmerge(d_td, d %>% filter(!is.na(t_mix)), id = NHC, mix_td = tdc(t_mix))
  d_td$mix_td[is.na(d_td$mix_td)] <- 0
  m_td <- coxph(f("Surv(tstart, tstop, ev_td)", "mix_td"),
                data = d_td, cluster = NHC)
  res <- rbind(res, fila("Dependiente del tiempo", m_td, "mix_td",
                         nrow(d), sum(d_td$ev_td), sum(d_td$ev_td[d_td$mix_td == 1])))
  
  # 3. TD con desfase
  d_lag <- tmerge(d_id, d_id, id = NHC, ev_td = event(t_fin, ev))
  d_lag <- tmerge(d_lag, d %>% filter(!is.na(t_mix)) %>% mutate(t_lag = t_mix + lag_anios),
                  id = NHC, mix_lag = tdc(t_lag))
  d_lag$mix_lag[is.na(d_lag$mix_lag)] <- 0
  m_lag <- coxph(f("Surv(tstart, tstop, ev_td)", "mix_lag"),
                 data = d_lag, cluster = NHC)
  res <- rbind(res, fila(sprintf("Dependiente del tiempo, desfase %d meses", round(lag_anios * 12)),
                         m_lag, "mix_lag",
                         nrow(d), sum(d_lag$ev_td), sum(d_lag$ev_td[d_lag$mix_lag == 1])))
  
  # 4. TD ajustado por episodios genotipados acumulados
  epi_t <- epis %>% filter(genotipado) %>% select(NHC, Fini) %>%
    inner_join(d %>% select(NHC, t_fin), by = "NHC") %>%
    inner_join(base %>% select(NHC, Ftx), by = "NHC") %>%
    mutate(t_epi = as.numeric(Fini - Ftx) / 365.25) %>%
    filter(t_epi > 0, t_epi <= t_fin) %>%
    select(NHC, t_epi)
  d_epi <- tmerge(d_td, epi_t, id = NHC, n_epi = cumtdc(t_epi))
  d_epi$n_epi[is.na(d_epi$n_epi)] <- 0
  m_epi <- coxph(f("Surv(tstart, tstop, ev_td)", "mix_td + n_epi"),
                 data = d_epi, cluster = NHC)
  res <- rbind(res, fila("Dependiente del tiempo + n episodios genotipados", m_epi, "mix_td",
                         nrow(d), sum(d_epi$ev_td), sum(d_epi$ev_td[d_epi$mix_td == 1])))
  
  # 5-6. Landmark
  for (L in landmarks) {
    dL <- d %>% filter(t_fin > L) %>%
      mutate(exp_L = as.integer(!is.na(t_mix) & t_mix <= L), t_L = t_fin - L)
    if (sum(dL$exp_L) >= 3 && sum(dL$ev) >= 3) {
      m_L <- coxph(f("Surv(t_L, ev)", "exp_L"), data = dL)
      res <- rbind(res, fila(sprintf("Landmark %d año(s)", L), m_L, "exp_L",
                             nrow(dL), sum(dL$ev), sum(dL$ev[dL$exp_L == 1])))
    }
  }
  
  # 7. Fijo restringido a >=2 episodios genotipados
  d2 <- d %>% filter(n_geno >= 2)
  if (sum(d2$ev) >= 3) {
    m_2 <- coxph(f("Surv(t_fin, ev)", "fijo"), data = d2)
    res <- rbind(res, fila("Fijo, solo pacientes con >=2 episodios genotipados", m_2, "fijo",
                           nrow(d2), sum(d2$ev), sum(d2$ev[d2$fijo == 1])))
  }
  
  
  # 8. Tres estados dependientes del tiempo:
  #    0 = sin episodio genotipado todavía; 1 = genotipado, solo únicos; 2 = mezcla detectada
  t_geno1 <- epis %>% filter(.data[[paste0("valido_", clasif)]]) %>% select(NHC, Fini) %>%
    inner_join(base %>% select(NHC, Ftx), by = "NHC") %>%
    mutate(t_g = as.numeric(Fini - Ftx) / 365.25) %>%
    group_by(NHC) %>% summarise(t_g = min(t_g), .groups = "drop") %>%
    filter(t_g > 0)
  d_3 <- tmerge(d_id, d_id, id = NHC, ev_td = event(t_fin, ev))
  d_3 <- tmerge(d_3, t_geno1, id = NHC, geno_td = tdc(t_g))
  d_3 <- tmerge(d_3, d %>% filter(!is.na(t_mix)), id = NHC, mix_td = tdc(t_mix))
  d_3$geno_td[is.na(d_3$geno_td)] <- 0; d_3$mix_td[is.na(d_3$mix_td)] <- 0
  d_3$estado <- factor(ifelse(d_3$mix_td == 1, "mixto",
                              ifelse(d_3$geno_td == 1, "unico", "sin_genotipar")),
                       levels = c("unico", "sin_genotipar", "mixto"))
  m_3 <- coxph(f("Surv(tstart, tstop, ev_td)", "estado"), data = d_3, cluster = NHC)
  res <- rbind(res, fila("Tres estados: mixto frente a unico demostrado", m_3, "estadomixto",
                         nrow(d), sum(d_3$ev_td), sum(d_3$ev_td[d_3$estado == "mixto"])))
  res <- rbind(res, fila("Tres estados: sin genotipar frente a unico demostrado", m_3, "estadosin_genotipar",
                         nrow(d), sum(d_3$ev_td), sum(d_3$ev_td[d_3$estado == "sin_genotipar"])))
  
  # 9. ANALISIS PRINCIPAL. Entrada diferida: el paciente entra en el conjunto de riesgo en la
  #    fecha de su primer episodio con genotipado valido para el gen estudiado, si esta libre
  #    del desenlace. El tiempo cero sigue siendo el trasplante (no se reinicia el reloj).
  #    Quien ya es mixto en esa primera muestra entra como mixto; el resto pasa a mixto en la
  #    primera deteccion de mezcla y permanece en esa categoria.
  d_9 <- d_td %>% inner_join(t_geno1, by = "NHC") %>%
    filter(tstop > t_g) %>%
    mutate(tstart = pmax(tstart, t_g))
  if (sum(d_9$ev_td) >= 3) {
    m_9 <- coxph(f("Surv(tstart, tstop, ev_td)", "mix_td"), data = d_9, cluster = NHC)
    res <- rbind(res, fila("Entrada diferida en el primer genotipado (PRINCIPAL)", m_9, "mix_td",
                           length(unique(d_9$NHC)), sum(d_9$ev_td), sum(d_9$ev_td[d_9$mix_td == 1])))
  }
  
  # Orden de presentacion: principal (entrada diferida) primero
  orden <- c("Entrada diferida en el primer genotipado (PRINCIPAL)",
             "Dependiente del tiempo",
             "Tres estados: mixto frente a unico demostrado",
             "Tres estados: sin genotipar frente a unico demostrado",
             "Dependiente del tiempo, desfase 6 meses",
             "Dependiente del tiempo + n episodios genotipados",
             "Landmark 1 año(s)", "Landmark 2 año(s)",
             "Fijo, solo pacientes con >=2 episodios genotipados",
             "Fijo desde trasplante (tesis)")
  res <- res[order(match(res$Modelo, orden)), ]
  
  cat(sprintf("\n══════ %s ══════\n", etiqueta))
  if (n_t0 > 0) cat(sprintf("Excluidos %d paciente(s) con evento en la fecha del trasplante (tiempo 0)\n", n_t0))
  cat(sprintf("Eventos con mezcla posterior al evento (mal clasificados en el modelo fijo): %d de %d\n",
              sum(d$ev == 1 & d$fijo == 1 & !is.na(d$t_mix) & d$t_mix >= d$t_fin),
              sum(d$ev == 1 & d$fijo == 1)))
  if (solo_principal) print(res[grepl("PRINCIPAL", res$Modelo), ], row.names = FALSE)
  else                print(res, row.names = FALSE)
  invisible(list(res = res, d_principal = d_9, modelo_principal = if (exists("m_9")) m_9 else NULL))
}

##
## Con exposición dependiente del tiempo no procede una KM por grupo fijo. Se usa la
## KM extendida (Simon y Makuch 1984; Snapinn et al. 2005): cada paciente aporta
## tiempo al grupo "único" hasta su primera mezcla y al grupo "mixto" después.
## Se dibuja sobre la población del análisis principal (entrada diferida en el
## primer genotipado válido del gen), con el tiempo desde el trasplante en el eje X.


pal_km <- c("Genotipos únicos" = "#2E8B57", "Genotipos mixtos" = "firebrick3")

km_extendida <- function(clasif, flag_evento, fecha_evento, titulo,
                         ylab = "Probabilidad libre de evento",
                         ylim = c(0, 1), mostrar_y = TRUE, tiempos_riesgo = seq(0, 12, 2),
                         covars = c("edad_tx", "risk_f")) {
  
  capture.output(a <- suppressWarnings(analisis_td(clasif, flag_evento, fecha_evento, titulo,
                                                   covars = covars, solo_principal = TRUE)))
  d9 <- a$d_principal %>%
    mutate(grupo = factor(ifelse(mix_td == 1, "Genotipos mixtos", "Genotipos únicos"),
                          levels = names(pal_km)))
  sf <- survfit(Surv(tstart, tstop, ev_td) ~ grupo, data = d9)
  
  # Curvas
  cur <- data.frame(t = sf$time, s = sf$surv, lo = sf$lower, hi = sf$upper,
                    grupo = rep(names(sf$strata), sf$strata)) %>%
    mutate(grupo = factor(sub("grupo=", "", grupo), levels = names(pal_km)))
  ini <- cur %>% group_by(grupo) %>%
    summarise(t = min(d9$tstart[d9$grupo == first(grupo)]), s = 1, lo = 1, hi = 1,
              .groups = "drop") %>% select(t, s, lo, hi, grupo)
  cur <- bind_rows(ini, cur) %>% arrange(grupo, t)
  
  # HR del modelo principal para la anotación
  r <- a$res[grepl("PRINCIPAL", a$res$Modelo), ]
  etiqueta_hr <- sprintf("p = %s", r$p)
  
  xmax <- max(d9$tstop, na.rm = TRUE)
  p <- ggplot(cur, aes(t, s, color = grupo, fill = grupo)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.08, color = NA) +
    geom_step(linewidth = 0.8) +
    scale_color_manual(values = pal_km, name = NULL) +
    scale_fill_manual(values = pal_km, name = NULL) +
    scale_x_continuous(breaks = seq(0, floor(xmax), 2), limits = c(0, xmax),
                       expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = ylim, labels = function(x) sub("\\.", ",", x)) +
    annotate("text", x = 0.2, y = ylim[1] + 0.06 * diff(ylim), label = etiqueta_hr,
             hjust = 0, size = 3.3) +
    labs(x = "Tiempo desde trasplante (años)", y = if (mostrar_y) ylab else NULL,
         title = titulo) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          legend.position = "none",
          axis.text = element_text(color = "black"))
  
  # Pacientes en riesgo por estado y tiempo
  riesgo <- lapply(tiempos_riesgo, function(tt) {
    d9 %>% filter(tstart <= tt, tstop > tt) %>% count(grupo, .drop = FALSE) %>%
      mutate(t = tt)
  }) %>% bind_rows() %>%
    tidyr::complete(grupo, t = tiempos_riesgo, fill = list(n = 0))
  tab <- ggplot(riesgo, aes(t, grupo, label = n, color = grupo)) +
    geom_text(size = 3.3) +
    scale_color_manual(values = pal_km) +
    scale_x_continuous(limits = c(0, xmax), breaks = seq(0, floor(xmax), 2),
                       expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_discrete(limits = rev(names(pal_km))) +
    labs(x = NULL, y = NULL, title = "Pacientes en riesgo") +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(size = 10), legend.position = "none",
          axis.text.y = element_text(size = 9, hjust = 1, color = "black"))
  
  list(plot = p, tabla = tab, bloque = (p / tab) + plot_layout(heights = c(4, 1)),
       res = a$res)
}

leyenda_km <- function() {
  g <- ggplot(data.frame(x = c(1, 2, 1, 2), y = c(1, 1, 2, 2),
                         grupo = factor(rep(names(pal_km), each = 2), levels = names(pal_km))),
              aes(x, y, color = grupo, group = grupo)) + geom_line(linewidth = 1) +
    scale_color_manual(values = pal_km, name = NULL) +
    guides(color = guide_legend(ncol = 1)) +
    theme_void() + theme(legend.position = "bottom",
                         legend.text = element_text(size = 11))
  patchwork::wrap_elements(full = cowplot::get_legend(g))
}

figura_km_2x2 <- function(flag_evento, fecha_evento, ylab, ylim = c(0, 1),
                          covars = c("edad_tx", "risk_f")) {
  b_gb   <- km_extendida("gbm",   flag_evento, fecha_evento, "gB",    ylab, ylim, TRUE,  covars = covars)
  b_gh   <- km_extendida("ghm",   flag_evento, fecha_evento, "gH",    ylab, ylim, TRUE,  covars = covars)
  b_gbgh <- km_extendida("gbghm", flag_evento, fecha_evento, "gB/gH", ylab, ylim, TRUE,  covars = covars)
  fig <- (b_gb$bloque | b_gh$bloque) / (b_gbgh$bloque | leyenda_km()) &
    theme(plot.background = element_rect(fill = "white", color = NA))
  fig
}
## Para guardar una figura concreta (no se guarda nada automaticamente):
## ggsave("~/Desktop/Figuras_tesis/Figura_KMext_neoplasia_2x2.png", fig_km_onco,
##        width = 24, height = 20, units = "cm", dpi = 300, bg = "white")

    ## Sustituyen a km_extendida() y figura_km_2x2() de [7.12] [b]. Requieren [7.12] [a].

pal_km <- c("Genotipos únicos" = "#2E8B57", "Genotipos mixtos" = "firebrick3")

km_extendida <- function(clasif, flag_evento, fecha_evento, titulo,
                         ylab = "Probabilidad libre de evento",
                         ylim = c(0, 1), mostrar_y = TRUE, tiempos_riesgo = NULL,
                         covars = c("edad_tx", "risk_f"), base = 12) {
  
  capture.output(a <- suppressWarnings(analisis_td(clasif, flag_evento, fecha_evento, titulo,
                                                   covars = covars, solo_principal = TRUE)))
  d9 <- a$d_principal %>%
    mutate(grupo = factor(ifelse(mix_td == 1, "Genotipos mixtos", "Genotipos únicos"),
                          levels = names(pal_km)))
  sf <- survfit(Surv(tstart, tstop, ev_td) ~ grupo, data = d9)
  
  cur <- data.frame(t = sf$time, s = sf$surv,
                    grupo = rep(names(sf$strata), sf$strata)) %>%
    mutate(grupo = factor(sub("grupo=", "", grupo), levels = names(pal_km)))
  ini <- cur %>% group_by(grupo) %>%
    summarise(t = min(d9$tstart[d9$grupo == first(grupo)]), s = 1, .groups = "drop")
  cur <- bind_rows(ini, cur) %>% arrange(grupo, t)
  
  if (is.null(tiempos_riesgo)) {
    xmax_dat <- max(d9$tstop, na.rm = TRUE)
    tiempos_riesgo <- seq(0, 2 * ceiling(xmax_dat / 2), 2)
  }
  r <- a$res[grepl("PRINCIPAL", a$res$Modelo), ]
  etiqueta_hr <- sprintf("p = %s", r$p)
  
  xmax    <- max(tiempos_riesgo)
  xbreaks <- tiempos_riesgo
  ybreaks <- seq(ylim[1], ylim[2], by = if (diff(ylim) > 0.5) 0.2 else 0.1)
  
  p <- ggplot(cur, aes(t, s, color = grupo)) +
    geom_step(linewidth = 1.1, lineend = "round") +
    scale_color_manual(values = pal_km, name = NULL) +
    scale_x_continuous(breaks = xbreaks, limits = c(0, xmax),
                       expand = expansion(mult = c(0.04, 0.02))) +
    scale_y_continuous(limits = ylim, breaks = ybreaks,
                       labels = function(x) sub("\\.", ",", x),
                       expand = expansion(mult = c(0, 0.02))) +
    annotate("text", x = xmax, y = ylim[2] - 0.02 * diff(ylim), label = etiqueta_hr,
             hjust = 1, vjust = 1, size = base * 0.3) +
    labs(x = "Años desde el trasplante", y = if (mostrar_y) ylab else NULL,
         title = titulo) +
    theme_classic(base_size = base) +
    theme(plot.title         = element_text(size = base, face = "bold", hjust = 0.5),
          legend.position    = "none",
          panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
          axis.line          = element_line(color = "black", linewidth = 0.5),
          axis.ticks         = element_line(color = "black", linewidth = 0.5),
          axis.ticks.length  = unit(3, "pt"),
          axis.text          = element_text(color = "black", size = base * 0.95),
          axis.title         = element_text(size = base),
          plot.margin        = margin(4, 8, 2, 4))
  
  riesgo <- lapply(tiempos_riesgo, function(tt) {
    d9 %>% filter(tstart <= tt, tstop > tt) %>% count(grupo, .drop = FALSE) %>% mutate(t = tt)
  }) %>% bind_rows() %>%
    tidyr::complete(grupo, t = tiempos_riesgo, fill = list(n = 0))
  tab <- ggplot(riesgo, aes(t, grupo, label = n, color = grupo)) +
    geom_text(size = base * 0.32) +
    scale_color_manual(values = pal_km) +
    scale_x_continuous(limits = c(0, xmax), breaks = xbreaks,
                       expand = expansion(mult = c(0.04, 0.02))) +
    scale_y_discrete(limits = rev(names(pal_km)), labels = c("Mixtos", "Únicos")) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL, title = "Pacientes en riesgo") +
    theme_void(base_size = base) +
    theme(plot.title  = element_text(size = base * 0.8, hjust = 0),
          legend.position = "none",
          axis.text.y = element_text(size = base * 0.8, hjust = 1, color = "black",
                                     margin = margin(r = 4)),
          plot.margin = margin(0, 8, 4, 4))
  
  list(plot = p, tabla = tab,
       bloque = (p / tab) + plot_layout(heights = c(4, 1)), res = a$res)
}

leyenda_km <- function(base = 12) {
  g <- ggplot(data.frame(x = c(1, 2, 1, 2), y = c(1, 1, 2, 2),
                         grupo = factor(rep(names(pal_km), each = 2), levels = names(pal_km))),
              aes(x, y, color = grupo, group = grupo)) + geom_line(linewidth = 1.1) +
    scale_color_manual(values = pal_km, name = NULL) +
    guides(color = guide_legend(ncol = 1)) +
    theme_void() + theme(legend.position = "bottom",
                         legend.text = element_text(size = base * 0.95),
                         legend.key.width = unit(24, "pt"))
  patchwork::wrap_elements(full = cowplot::get_legend(g))
}

figura_km_2x2 <- function(flag_evento, fecha_evento, ylab, ylim = c(0, 1),
                          covars = c("edad_tx", "risk_f"), base = 12) {
  b_gb   <- km_extendida("gbm",   flag_evento, fecha_evento, "gB",    ylab, ylim, TRUE, covars = covars, base = base)
  b_gh   <- km_extendida("ghm",   flag_evento, fecha_evento, "gH",    ylab, ylim, TRUE, covars = covars, base = base)
  b_gbgh <- km_extendida("gbghm", flag_evento, fecha_evento, "gB/gH", ylab, ylim, TRUE, covars = covars, base = base)
  # Los cuatro cuadrantes con el mismo ancho: el panel gB/gH no se estira sobre la celda de la leyenda
  (b_gb$bloque | b_gh$bloque) / (b_gbgh$bloque | leyenda_km(base)) +
    plot_layout(widths = c(1, 1), heights = c(1, 1)) &
    theme(plot.background = element_rect(fill = "white", color = NA))
}

## ---- Figura 38. Infecciones oportunistas no relacionadas con CMV ----
# ── 6. Preparar datos ─────────────────────────────────────────────────────────
df_opo <- pacientes %>%
  mutate(
    Ftx_dt           = as.Date(Ftx),
    F_opo_dt         = as.Date(F_opo),
    Ffinseg_dt       = as.Date(Ffinseg),
    fecha_fin_real   = if_else(inf_opo == 1, F_opo_dt, Ffinseg_dt),
    t_analisis_anios = as.numeric(fecha_fin_real - Ftx_dt) / 365.25,
    evento_opo       = if_else(inf_opo == 1, 1L, 0L),
    risk_f           = factor(risk, levels = c("MEDIO", "ALTO")),
    gbghm_grupo      = factor(gbghm, levels = c(0,1),
                              labels = c("Genotipos únicos",
                                         "Genotipos mixtos"))
  ) %>%
  filter(t_analisis_anios > 0)

# ── 7. Descriptivo ────────────────────────────────────────────────────────────
cat("\n── Infecciones oportunistas por grupo ──\n")
df_opo %>%
  group_by(gbghm_grupo) %>%
  summarise(
    n     = n(),
    n_opo = sum(evento_opo),
    pct   = round(n_opo / n * 100, 1),
    .groups = "drop"
  ) %>% print()

cat("\n── Total pacientes con IO ──\n")
cat(sprintf("n = %d (%.1f%%)\n",
            sum(pacientes$inf_opo == 1, na.rm = TRUE),
            100 * mean(pacientes$inf_opo == 1, na.rm = TRUE)))

cat("\n── Número de episodios por paciente ──\n")
table(pacientes$n_opo[pacientes$inf_opo == 1], useNA = "always") %>% print()

cat("\n── Etiología IO ──\n")
table(pacientes$det_opo, useNA = "always") %>%
  sort(decreasing = TRUE) %>% print()

cat("\n── Test chi-cuadrado IO vs grupo ──\n")
print(chisq.test(table(df_opo$gbghm_grupo, df_opo$evento_opo)))

# ── 8. KM + log-rank ─────────────────────────────────────────────────────────
surv_opo <- Surv(df_opo$t_analisis_anios, df_opo$evento_opo)
km_opo   <- survfit(surv_opo ~ gbghm_grupo, data = df_opo)
lr_opo   <- survdiff(surv_opo ~ gbghm_grupo, data = df_opo)
p_lr     <- 1 - pchisq(lr_opo$chisq, df = 1)
cat(sprintf("\n── Log-rank IO: p = %.3f ──\n", p_lr))

cat("\n── Supervivencia libre IO a 5 años ──\n")
print(summary(km_opo, times = 5))

# ── 9. Cox univariante y ajustado ─────────────────────────────────────────────
cox_opo_crudo <- coxph(surv_opo ~ gbghm_grupo, data = df_opo)
cox_opo_ajust <- coxph(surv_opo ~ gbghm_grupo + edad_tx + risk_f,
                       data = df_opo)

cat("\n── Cox IO crudo ──\n")
hr_oc <- exp(coef(cox_opo_crudo)[1])
ci_oc <- exp(confint(cox_opo_crudo)[1,])
p_oc  <- summary(cox_opo_crudo)$coefficients[1,5]
cat(sprintf("HR: %.2f (%.2f–%.2f)  p = %.3f\n",
            hr_oc, ci_oc[1], ci_oc[2], p_oc))

cat("\n── Cox IO ajustado ──\n")
data.frame(
  HR    = round(exp(coef(cox_opo_ajust)), 2),
  IC_lo = round(exp(confint(cox_opo_ajust))[,1], 2),
  IC_hi = round(exp(confint(cox_opo_ajust))[,2], 2),
  p     = round(summary(cox_opo_ajust)$coefficients[,5], 3)
) %>% print()

# ── IO segun genotipo mixto: gB, gH y gB/gH ──────────────────────────────────

# Funcion: descriptivo + test para una variable de mezcla
io_por_genotipo <- function(data, var_mix, etiqueta) {
  d <- data %>%
    mutate(grupo = factor(.data[[var_mix]], levels = c(0, 1),
                          labels = c("\u00danicos", "Mixtos"))) %>%
    filter(!is.na(grupo))
  
  tab <- table(d$grupo, d$evento_opo)
  
  resumen <- d %>%
    group_by(grupo) %>%
    summarise(
      n     = n(),
      n_opo = sum(evento_opo),
      pct   = round(100 * n_opo / n, 1),
      .groups = "drop"
    )
  
  # Fisher (robusto con n pequenos) y chi-cuadrado de Yates
  p_fisher <- fisher.test(tab)$p.value
  p_chisq  <- suppressWarnings(chisq.test(tab)$p.value)
  
  cat(sprintf("\n══════ %s ══════\n", etiqueta))
  print(resumen)
  cat(sprintf("p (Fisher)        = %s\n",
              sub("\\.", ",", sprintf("%.3f", p_fisher))))
  cat(sprintf("p (chi-cuadrado)  = %s\n",
              sub("\\.", ",", sprintf("%.3f", p_chisq))))
  
  invisible(resumen)
}

# Aplicar a las tres clasificaciones (usa df_opo, que ya tiene evento_opo)
io_por_genotipo(df_opo, "gbm",   "gB mixtos")
io_por_genotipo(df_opo, "ghm",   "gH mixtos")
io_por_genotipo(df_opo, "gbghm", "gB/gH mixtos")

# [Nuevo] Datos de la Figura 38: episodios por etiología a partir de det_opo.
# En el script antiguo df_counts y colores no se definían (venían de otra sesión).
# Cada paciente aporta n_opo episodios; las etiologías se leen de det_opo.
# Si det_opo menciona más etiologías que n_opo, se conservan las últimas n_opo
# (afecta a un paciente: "VEB (PTLD) + IFI pulmonar", con n_opo = 1, que cuenta como IFI).
clasificar_opo <- function(txt) {
  partes <- trimws(unlist(strsplit(txt, ",|\\+|;")))
  partes <- partes[partes != ""]
  case_when(
    grepl("VVZ", partes, ignore.case = TRUE)            ~ "VVZ",
    grepl("candidiasis", partes, ignore.case = TRUE)    ~ "Candidiasis",
    grepl("IFI", partes)                                ~ "IFI",
    grepl("VEB|PTLD|PLTD", partes)                      ~ "VEB/PTLD",
    grepl("TBC", partes)                                ~ "TBC",
    grepl("Cryptosporidium", partes)                    ~ "Cryptosporidium",
    grepl("VHH8|Kaposi", partes)                        ~ "VHH8/Sarcoma de Kaposi",
    grepl("JC", partes)                                 ~ "Virus JC",
    TRUE                                                ~ NA_character_
  )
}

episodios_opo <- pacientes %>%
  filter(inf_opo == 1, !is.na(det_opo)) %>%
  select(NHC, det_opo, n_opo) %>%
  mutate(categoria = lapply(det_opo, clasificar_opo)) %>%
  mutate(categoria = map2(categoria, n_opo, ~ tail(.x[!is.na(.x)], .y))) %>%
  tidyr::unnest(categoria)

df_counts <- episodios_opo %>%
  count(categoria) %>%
  mutate(pct = round(100 * n / sum(n), 1))

stopifnot(sum(df_counts$n) == sum(pacientes$n_opo, na.rm = TRUE))
print(df_counts %>% arrange(desc(n)))

colores <- c(
  "VVZ"                    = "#4E79A7",
  "Candidiasis"            = "#F28E2B",
  "IFI"                    = "#8CD17D",
  "VEB/PTLD"               = "#E15759",
  "TBC"                    = "#7B3F9E",
  "VHH8/Sarcoma de Kaposi" = "#2CB5A0",
  "Virus JC"               = "#8C8C8C",
  "Cryptosporidium"        = "#EDC948"
)
# ── Etiquetas personalizadas (edita aquí a tu gusto) ─────────────────────────
etiquetas_custom <- c(
  "VVZ"                    = "Herpes zóster",
  "Candidiasis"            = "Candidiasis mucocutánea",
  "IFI"                    = "Infección fúngica invasiva",
  "VEB/PTLD"               = "PTLD asociado a VEB",
  "TBC"                    = "Tuberculosis diseminada",
  "Cryptosporidium"        = "Criptosporidiosis",
  "VHH8/Sarcoma de Kaposi" = "Sarcoma de Kaposi (VHH-8)",
  "Virus JC"               = "LMP por virus JC"
)

figura_opo_barras <- ggplot(df_counts %>% arrange(n),
                            aes(x = reorder(categoria, n), y = n,
                                fill = categoria)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.2) +
  geom_text(
    aes(label = paste0("n=", n, " (",
                       gsub("\\.", ",", as.character(pct)), "%)")),
    hjust = -0.1, size = 3.8, color = "grey20"
  )  +
  scale_fill_manual(values = colores, labels = etiquetas_custom) +
  scale_x_discrete(labels = etiquetas_custom) +  # ← etiquetas eje Y
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.45))
  ) +
  coord_flip() +
  labs(x = NULL, y = "Número de episodios") +
  theme_classic(base_size = 13) +
  theme(
    legend.position    = "none",
    axis.text.y        = element_text(size = 11, color = "grey20"),
    axis.text.x        = element_text(size = 10, color = "grey40"),
    axis.title.x       = element_text(size = 11, color = "grey20"),
    axis.line          = element_line(color = "grey60", linewidth = 0.4),
    axis.ticks         = element_line(color = "grey60", linewidth = 0.3),
    plot.margin        = margin(10, 20, 10, 10)
  )

figura_opo_barras

guardar_fig(figura_opo_barras, "Figura_38", ancho = 16, alto = 9)

## ---- Modelo dependiente del tiempo: infecciones oportunistas (para la tabla resumen) ----
res_inf_opo <- list(
  gB   = analisis_td("gbm",   "inf_opo", "F_opo", "Infecciones oportunistas no CMV | gB mixtos"),
  gH   = analisis_td("ghm",   "inf_opo", "F_opo", "Infecciones oportunistas no CMV | gH mixtos"),
  gBgH = analisis_td("gbghm", "inf_opo", "F_opo", "Infecciones oportunistas no CMV | gB/gH mixtos")
)

## Descomentar para ver los nueve modelos de cada clasificación
# analisis_td("gbm",   "inf_opo", "F_opo", "Infecciones oportunistas no CMV | gB mixtos", solo_principal = FALSE)
# analisis_td("ghm",   "inf_opo", "F_opo", "Infecciones oportunistas no CMV | gH mixtos", solo_principal = FALSE)
# analisis_td("gbghm", "inf_opo", "F_opo", "Infecciones oportunistas no CMV | gB/gH mixtos", solo_principal = FALSE)

## ---- Figura 39. Episodios de infección por virus respiratorios comunitarios por año ----
VR <- read_excel(f_bd, 
                 sheet = "VR")
head(VR)


VR <- read_excel(f_bd, sheet = "VR") %>%
  select(-gbghm, -genotipo) %>%
  left_join(pacientes %>% select(NHC, gbghm), by = "NHC") %>%
  mutate(año = as.integer(año),
         mes = as.integer(mes),
         genotipo_label = ifelse(gbghm == 0, "Genotipos únicos", "Genotipos mixtos"),
         virus = factor(virus, levels = c("SARS-CoV-2", "Influenza", "VRS",
                                          "Rhinovirus", "Metapneumovirus",
                                          "Parainfluenza", "Adenovirus",
                                          "Coronavirus estacional")))

        
cat("=== PREVALENCIA INFECCIÓN VR ===\n")
cat("Total pacientes con inf_VR:", sum(pacientes$inf_VR, na.rm = TRUE), "\n")
cat("Prevalencia global:", round(mean(pacientes$inf_VR, na.rm = TRUE)*100, 1), "%\n\n")

prev_geno <- pacientes %>%
  group_by(gbghm) %>%
  summarise(n     = n(),
            n_inf = sum(inf_VR, na.rm = TRUE),
            prev  = round(n_inf / n * 100, 1), .groups = "drop")
print(prev_geno)

cat("\nTest Fisher prevalencia VR por genotipo:\n")
print(fisher.test(table(pacientes$gbghm, pacientes$inf_VR)))

cat("\n=== EPISODIOS POR PACIENTE ===\n")

n_ep_pac <- VR %>% count(NHC, gbghm, genotipo_label, name = "n_ep")

n_ep_all <- pacientes %>%
  select(NHC, gbghm) %>%
  left_join(n_ep_pac %>% select(NHC, n_ep), by = "NHC") %>%
  replace_na(list(n_ep = 0)) %>%
  mutate(genotipo_label = ifelse(gbghm == 0, "Genotipos únicos", "Genotipos mixtos"))

cat("Total episodios - únicos:", sum(n_ep_all$n_ep[n_ep_all$gbghm == 0]), "\n")
cat("Total episodios - mixtos:", sum(n_ep_all$n_ep[n_ep_all$gbghm == 1]), "\n")
cat("Total episodios global:",   sum(n_ep_all$n_ep), "\n\n")

cat("Mediana episodios (pacientes infectados):\n")
n_ep_pac %>%
  group_by(genotipo_label) %>%
  summarise(mediana = median(n_ep),
            p25     = quantile(n_ep, 0.25),
            p75     = quantile(n_ep, 0.75),
            .groups = "drop") %>%
  print()

cat("\nTest Wilcoxon episodios (todos los pacientes, incluyendo 0):\n")
print(wilcox.test(n_ep ~ gbghm, data = n_ep_all))

cat("\n=== DISTRIBUCIÓN POR VIRUS ===\n")
dist_virus <- VR %>%
  count(genotipo_label, virus, .drop = FALSE) %>%
  group_by(genotipo_label) %>%
  mutate(total = sum(n),
         pct   = round(n / total * 100, 1)) %>%
  ungroup()
print(dist_virus, n = 30)

cat("\n=== ANÁLISIS TEMPORAL ===\n")
VR_con_año <- VR %>% filter(!is.na(año))
n_post  <- sum(VR_con_año$año >= 2020)
n_total <- nrow(VR_con_año)
cat("Episodios 2020-2022:", n_post, "/", n_total,
    "(", round(n_post / n_total * 100, 1), "%)\n")

cat("\n=== COVID-19 ===\n")
covid_nhc <- VR %>% filter(virus == "SARS-CoV-2") %>% distinct(NHC)
pacientes <- pacientes %>%
  mutate(covid = ifelse(NHC %in% covid_nhc$NHC, 1, 0))

pacientes %>%
  group_by(gbghm) %>%
  summarise(n       = n(),
            n_covid = sum(covid),
            pct     = round(n_covid / n * 100, 1), .groups = "drop") %>%
  print()

cat("\nTest Fisher COVID por genotipo:\n")
print(fisher.test(table(pacientes$gbghm, pacientes$covid)))

cat("\nNeumonía COVID por genotipo:\n")
pacientes %>%
  group_by(gbghm) %>%
  summarise(n          = n(),
            n_neumonia = sum(neumoniaCOVID, na.rm = TRUE),
            pct        = round(n_neumonia / n * 100, 1),
            .groups    = "drop") %>%
  print()

cat("\n=== VR Y RECHAZO AGUDO ===\n")

cat("Mediana episodios VR según RA:\n")
n_ep_all %>%
  left_join(pacientes %>% select(NHC, RA), by = "NHC") %>%
  group_by(RA) %>%
  summarise(n       = n(),
            mediana = median(n_ep),
            p25     = quantile(n_ep, 0.25),
            p75     = quantile(n_ep, 0.75),
            .groups = "drop") %>%
  print()

print(wilcox.test(n_ep ~ RA,
                  data = n_ep_all %>%
                    left_join(pacientes %>% select(NHC, RA), by = "NHC")))

cat("\nOR infección VR ~ RA en genotipos únicos:\n")
pac_unicos <- pacientes %>% filter(gbghm == 0)
print(fisher.test(table(pac_unicos$inf_VR, pac_unicos$RA)))

cat("\nOR infección VR ~ RA en genotipos mixtos:\n")
pac_mixtos <- pacientes %>% filter(gbghm == 1)
print(fisher.test(table(pac_mixtos$inf_VR, pac_mixtos$RA)))

cat("\n=== VR + RA → DCIP ===\n")
pacientes %>%
  mutate(grupo = case_when(
    inf_VR == 1 & RA == 1 ~ "VR + RA",
    inf_VR == 0 & RA == 0 ~ "Sin VR ni RA",
    inf_VR == 1 & RA == 0 ~ "Solo VR",
    inf_VR == 0 & RA == 1 ~ "Solo RA"
  )) %>%
  group_by(grupo) %>%
  summarise(n      = n(),
            n_DCIP = sum(DCIP, na.rm = TRUE),
            pct    = round(n_DCIP / n * 100, 1),
            .groups = "drop") %>%
  print()

virus_colors <- c(
  "SARS-CoV-2"             = "#EE6363",
  "Influenza"              = "dodgerblue2",
  "VRS"                    = "seagreen3",
  "Rhinovirus"             = "#9B5DE5",
  "Metapneumovirus"        = "#F77F00",
  "Parainfluenza"          = "#00E5EE",
  "Adenovirus"             = "#FF69B4",
  "Coronavirus estacional" = "#ADB5BD"
)

años_rng <- 2010:2022

df_todos <- VR_con_año %>%
  mutate(grupo_label = "Todos los pacientes (n = 146)") %>%
  count(año, virus, grupo_label)

df_geno <- VR_con_año %>%
  mutate(grupo_label = ifelse(gbghm == 0,
                              "Genotipos gB+gH únicos (n = 95)",
                              "Genotipos gB+gH mixtos (n = 51)")) %>%
  count(año, virus, grupo_label)

df_fig <- bind_rows(df_todos, df_geno) %>%
  mutate(grupo_label = factor(grupo_label,
                              levels = c("Todos los pacientes (n = 146)",
                                         "Genotipos gB+gH únicos (n = 95)",
                                         "Genotipos gB+gH mixtos (n = 51)")))

df_completo <- expand.grid(
  año         = años_rng,
  virus       = factor(names(virus_colors), levels = names(virus_colors)),
  grupo_label = factor(levels(df_fig$grupo_label), levels = levels(df_fig$grupo_label))
) %>%
  left_join(df_fig, by = c("año", "virus", "grupo_label")) %>%
  replace_na(list(n = 0))

fig_VR <- ggplot(df_completo,
                 aes(x = año, y = n, fill = virus)) +
  geom_col(width = 0.75, position = "stack") +
  geom_vline(xintercept = 2019.5, linetype = "dashed",
             color = "gray40", linewidth = 0.45) +
  annotate("text", x = 2019.65, y = Inf,
           label = "Pandemia\nCOVID-19",
           hjust = 0, vjust = 1.3,
           size = 3.8, color = "gray40") +
  facet_wrap(~ grupo_label, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = virus_colors, breaks = names(virus_colors)) +
  scale_x_continuous(breaks = años_rng, labels = as.character(años_rng)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Año",
       y = "Núm. infecciones",
       fill = NULL) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 13),
    axis.text.x      = element_text(size = 11),
    axis.text.y      = element_text(size = 11),
    axis.title       = element_text(size = 12),
    legend.position  = "bottom",
    legend.text      = element_text(size = 11),
    panel.grid       = element_blank()
  ) +
  guides(fill = guide_legend(
    nrow = 2,
    byrow = TRUE,
    keywidth  = 0.8,
    keyheight = 0.8
  ))

print(fig_VR)

# [guardado antiguo] ggsave("~/Desktop/Figuras_tesis/Figura_R37.png", fig_VR,
# [guardado antiguo]        width = 27, height = 16, units = "cm",
# [guardado antiguo]        dpi = 300)


# Para cada paciente con RA y VR, ver qué ocurrió antes
pacientes %>%
  filter(RA == 1, inf_VR == 1) %>%
  select(NHC, F_RA1, gbghm) %>%
  left_join(VR %>% select(NHC, mes, año) %>%
              mutate(F_VR_primera = as.Date(paste(año, mes, "01", sep="-"))) %>%
              group_by(NHC) %>%
              summarise(F_VR_primera = min(F_VR_primera, na.rm = TRUE)),
            by = "NHC") %>%
  mutate(F_RA1 = as.Date(F_RA1),
         orden = case_when(
           F_RA1 < F_VR_primera ~ "RA antes de VR",
           F_VR_primera < F_RA1 ~ "VR antes de RA",
           TRUE ~ "simultáneo"
         )) %>%
  count(orden)

guardar_fig(fig_VR, "Figura_39", ancho = 27, alto = 16)

## ---- Rechazo agudo: modelo dependiente del tiempo ----
## 21 de 26 rechazos ocurren antes del primer genotipado: quedan 5 eventos evaluables.
## El modelo se ejecuta para documentarlo, pero no es interpretable y no se dibuja KM.
res_RA <- list(
  gB   = analisis_td("gbm",   "RA", "F_RA1", "Rechazo agudo | gB mixtos"),
  gH   = analisis_td("ghm",   "RA", "F_RA1", "Rechazo agudo | gH mixtos"),
  gBgH = analisis_td("gbghm", "RA", "F_RA1", "Rechazo agudo | gB/gH mixtos")
)

## Descomentar para ver los nueve modelos de cada clasificación
# analisis_td("gbm",   "RA", "F_RA1", "Rechazo agudo | gB mixtos", solo_principal = FALSE)
# analisis_td("ghm",   "RA", "F_RA1", "Rechazo agudo | gH mixtos", solo_principal = FALSE)
# analisis_td("gbghm", "RA", "F_RA1", "Rechazo agudo | gB/gH mixtos", solo_principal = FALSE)

## ---- Figura 40. Genotipados de CMV respecto al primer episodio de rechazo agudo ----
##
## Autocontenido: requiere `pacientes` y `dnaemias` cargados (no depende de [7.12]).
## Una fila por paciente con rechazo agudo. El cero es el primer rechazo agudo;
## cada punto es un episodio de DNAemia genotipado, situado en días desde el rechazo
## y coloreado según su patrón gB/gH (único o mixto).


unidad <- "meses"  # "dias" o "meses"
div    <- if (unidad == "dias") 1 else 30.44
etq_x  <- if (unidad == "dias") "Días desde el primer rechazo agudo" else
  "Meses desde el primer rechazo agudo"

# ── 1. Episodios genotipados ─────────────────────────────────────────────────
epis_gt <- dnaemias %>%
  mutate(Fini = as.Date(Fini)) %>%
  filter(!is.na(Fini), (!is.na(gbsum) & gbsum >= 1) | (!is.na(ghsum) & ghsum >= 1)) %>%
  transmute(NHC, Fini,
            patron = factor(ifelse(gbghmix == 1, "Genotipos mixtos gB/gH", "Genotipos únicos gB/gH"),
                            levels = c("Genotipos únicos gB/gH", "Genotipos mixtos gB/gH")))

# ── 2. Pacientes con rechazo agudo (excluido el hiperagudo) ──────────────────
ra <- pacientes %>%
  transmute(NHC, Ftx = as.Date(Ftx), F_RA1 = as.Date(F_RA1),
            F_RA2 = as.Date(F_RA2), F_RA3 = as.Date(F_RA3), RA) %>%
  filter(RA == 1, !is.na(F_RA1), as.numeric(F_RA1 - Ftx) > 0)

df_fig <- ra %>%
  inner_join(epis_gt, by = "NHC") %>%
  mutate(t = as.numeric(Fini - F_RA1) / div)

# Orden de filas: por el primer episodio genotipado respecto al rechazo
orden <- df_fig %>% group_by(NHC) %>%
  summarise(t_primero = min(t), .groups = "drop") %>%
  arrange(t_primero) %>%
  mutate(id = sprintf("P%02d", row_number()),
         # Categoría temporal: concomitante si el primer genotipado cae dentro
         # del mes anterior o posterior al rechazo (±30,44 días)
         grupo_pac = case_when(
           abs(t_primero * div) <= 30.44 ~ "Proximidad temporal (±30 días)",
           t_primero > 0               ~ "DNAemia >30 días después del rechazo",
           TRUE                        ~ "DNAemia >30 días antes del rechazo"))
df_fig <- df_fig %>% inner_join(orden, by = "NHC") %>%
  mutate(id = factor(id, levels = orden$id))

# Segundos y terceros rechazos (marcas adicionales)
df_ra_rec <- ra %>%
  select(NHC, F_RA1, F_RA2, F_RA3) %>%
  tidyr::pivot_longer(c(F_RA2, F_RA3), names_to = "episodio", values_to = "F_rec") %>%
  filter(!is.na(F_rec)) %>%
  mutate(t = as.numeric(F_rec - F_RA1) / div,
         episodio = "Episodios adicionales de rechazo agudo") %>%
  inner_join(orden %>% select(NHC, id), by = "NHC") %>%
  mutate(id = factor(id, levels = orden$id))
cat(sprintf("Pacientes con rechazo recurrente: %d (episodios adicionales: %d)\n",
            n_distinct(df_ra_rec$NHC), nrow(df_ra_rec)))

# ── 3. Resumen numérico ──────────────────────────────────────────────────────
n_tot   <- nrow(orden)
n_ant   <- sum(orden$grupo_pac == "DNAemia >30 días después del rechazo")
n_conc  <- sum(orden$grupo_pac == "Proximidad temporal (±30 días)")
n_post  <- sum(orden$grupo_pac == "DNAemia >30 días antes del rechazo")
pct     <- function(k) sub("\\.", ",", sprintf("%.1f", 100 * k / n_tot))
cat(sprintf("Pacientes con rechazo: %d\n  DNAemia >30 dias despues del rechazo: %d (%s %%)\n  proximidad temporal (+/-30 dias): %d (%s %%)\n  DNAemia >30 dias antes del rechazo: %d (%s %%)\n",
            n_tot, n_ant, pct(n_ant), n_conc, pct(n_conc), n_post, pct(n_post)))
cat(sprintf("%s hasta el primer genotipado (DNAemia >30 dias despues):\n", etq_x))
print(summary(orden$t_primero[orden$grupo_pac == "DNAemia >30 días después del rechazo"]))
cat("Intervalo en los de proximidad temporal (dias):\n")
print(round(orden$t_primero[orden$grupo_pac == "Proximidad temporal (±30 días)"] * div))

# ── 4. Figura (media página A4 en Word: 16 x 15 cm) ──────────────────────────
niv <- c("Genotipos únicos gB/gH", "Genotipos mixtos gB/gH", "Episodios adicionales de rechazo agudo")
pal <- c("Genotipos únicos gB/gH" = "#2E8B57",
         "Genotipos mixtos gB/gH" = "firebrick3",
         "Episodios adicionales de rechazo agudo" = "black")
df_fig    <- df_fig    %>% mutate(tipo = factor(as.character(patron), levels = niv))
df_ra_rec <- df_ra_rec %>% mutate(tipo = factor(episodio, levels = niv))
xr <- range(c(df_fig$t, df_ra_rec$t))

# Trazo por paciente que une episodios genotipados y rechazos posteriores
df_linea <- bind_rows(df_fig %>% select(id, t),
                      df_ra_rec %>% select(id, t)) %>%
  arrange(id, t)

banda   <- 30.44 / div   # ±1 mes en las unidades del eje
x_llave <- 12 * banda    # posición de las llaves (12 meses a la derecha del cero)
tick    <- 1.2 * banda   # longitud de los extremos de la llave
fig_ra_tiempo <- ggplot(df_fig, aes(x = t, y = id)) +
  annotate("rect", xmin = -banda, xmax = banda, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_line(data = df_linea, aes(x = t, y = id, group = id),
            color = "grey75", linewidth = 0.35) +
  geom_point(aes(color = tipo), size = 2.3) +
  geom_point(data = df_ra_rec, aes(x = t, y = id, color = tipo),
             shape = 24, fill = "white", size = 2.1, stroke = 0.7) +
  annotate("text", x = 0, y = n_tot + 0.8, label = "Primer rechazo agudo",
           size = 3.8, angle = 90, hjust = 1, vjust = -0.4) +
  annotate("text", x = xr[2], y = n_post + n_conc + 7, hjust = 1, vjust = 0, size = 3.7,
           lineheight = 0.9,
           label = sprintf("DNAemia iniciada >30 días\ndespués del rechazo: %d/%d (%s %%)",
                           n_ant, n_tot, pct(n_ant))) +
  # Llaves para los grupos concomitante y posterior (a la izquierda)
  annotate("segment", x = x_llave, xend = x_llave,
           y = n_post + 0.7, yend = n_post + n_conc + 0.3, linewidth = 0.5) +
  annotate("segment", x = x_llave, xend = x_llave - tick,
           y = c(n_post + 0.7, n_post + n_conc + 0.3),
           yend = c(n_post + 0.7, n_post + n_conc + 0.3), linewidth = 0.5) +
  annotate("text", x = x_llave + tick, y = n_post + n_conc / 2 + 0.5, hjust = 0, vjust = 0.5,
           size = 3.7, lineheight = 0.9,
           label = sprintf("Proximidad temporal (±30 días):\n%d/%d (%s %%)",
                           n_conc, n_tot, pct(n_conc))) +
  annotate("segment", x = x_llave, xend = x_llave,
           y = 0.7, yend = n_post + 0.3, linewidth = 0.5) +
  annotate("segment", x = x_llave, xend = x_llave - tick,
           y = c(0.7, n_post + 0.3), yend = c(0.7, n_post + 0.3), linewidth = 0.5) +
  annotate("text", x = x_llave + tick, y = n_post / 2 + 0.5, hjust = 0, vjust = 0.5,
           size = 3.7, lineheight = 0.9,
           label = sprintf("DNAemia iniciada >30 días\nantes del rechazo: %d/%d (%s %%)",
                           n_post, n_tot, pct(n_post))) +
  scale_color_manual(values = pal, name = NULL, drop = FALSE) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE,
                              override.aes = list(shape = c(16, 16, 24),
                                                  fill  = c(NA, NA, "white"),
                                                  size  = c(2.3, 2.3, 2.1)))) +
  scale_x_continuous(breaks = if (unidad == "meses") seq(-36, 132, 12) else waiver(),
                     labels = function(x) sub("-", "\u2212",
                                              format(x, big.mark = ".", decimal.mark = ",",
                                                     scientific = FALSE, trim = TRUE)),
                     expand = expansion(mult = c(0.04, 0.04))) +
  scale_y_discrete(expand = expansion(add = c(0.8, 1.2))) +
  labs(x = etq_x, y = "Pacientes con rechazo agudo") +
  theme_classic(base_size = 14) +
  theme(legend.position   = "bottom",
        legend.text       = element_text(size = 10.5, margin = margin(l = 3, r = 14)),
        legend.key.width  = unit(8, "pt"),
        legend.margin     = margin(t = 0),
        legend.justification = "center",
        legend.spacing.y  = unit(1, "pt"),
        axis.text.y       = element_text(size = 10, color = "black"),
        axis.text.x       = element_text(size = 11.5, color = "black"),
        axis.title        = element_text(size = 12.5),
        axis.ticks.y      = element_blank())

print(fig_ra_tiempo)

# ── 5. Exportar (media página A4 en Word, 16 x 15 cm) ────────────────────────
## PNG (600 ppp) para Word, TIFF (LZW, 600 ppp) para revistas y SVG (vectorial).
## SVG con svglite y fix_text_size = FALSE: sin ese argumento svglite añade
## atributos textLength que Word interpreta mal y muestra la figura encogida en una esquina.
## Requiere install.packages("svglite").
# [guardado antiguo] ggsave("~/Desktop/Figuras_tesis/Figura_41_rechazo_vs_genotipado.png", fig_ra_tiempo,
# [guardado antiguo]        width = 16, height = 15, units = "cm", dpi = 600, bg = "white")
# [guardado antiguo] ggsave("~/Desktop/Figuras_tesis/Figura_41_rechazo_vs_genotipado.tiff", fig_ra_tiempo,
# [guardado antiguo]        device = "tiff", compression = "lzw",
# [guardado antiguo]        width = 16, height = 15, units = "cm", dpi = 600, bg = "white")
# [guardado antiguo] ggsave("~/Desktop/Figuras_tesis/Figura_41_rechazo_vs_genotipado.svg", fig_ra_tiempo,
# [guardado antiguo]        device = svglite::svglite, fix_text_size = FALSE,
# [guardado antiguo]        width = 16, height = 15, units = "cm", bg = "white")

guardar_fig(fig_ra_tiempo, "Figura_40", ancho = 16, alto = 15)

## ---- DCIP: modelo dependiente del tiempo ----
res_DCIP <- list(
  gB   = analisis_td("gbm",   "DCIP", "F_DCIP", "Disfunción crónica del injerto pulmonar | gB mixtos"),
  gH   = analisis_td("ghm",   "DCIP", "F_DCIP", "Disfunción crónica del injerto pulmonar | gH mixtos"),
  gBgH = analisis_td("gbghm", "DCIP", "F_DCIP", "Disfunción crónica del injerto pulmonar | gB/gH mixtos")
)

## Descomentar para ver los nueve modelos de cada clasificación
# analisis_td("gbm",   "DCIP", "F_DCIP", "Disfunción crónica del injerto pulmonar | gB mixtos", solo_principal = FALSE)
# analisis_td("ghm",   "DCIP", "F_DCIP", "Disfunción crónica del injerto pulmonar | gH mixtos", solo_principal = FALSE)
# analisis_td("gbghm", "DCIP", "F_DCIP", "Disfunción crónica del injerto pulmonar | gB/gH mixtos", solo_principal = FALSE)

## ---- Figura 41. Supervivencia libre de DCIP según genotipos mixtos ----
fig_km_DCIP <- figura_km_2x2("DCIP", "F_DCIP", "Libre de DCIP", ylim = c(0.4, 1))
print(fig_km_DCIP)

## Exportar a SVG para una página A4 en Word (16 x 22 cm)
# [guardado antiguo] ggsave("~/Desktop/Figuras_tesis/Figura_42_KMext_DCIP_2x2.svg", fig_km_DCIP,
# [guardado antiguo]        device = svglite::svglite, fix_text_size = FALSE,
# [guardado antiguo]        width = 16, height = 22, units = "cm", bg = "white")

guardar_fig(fig_km_DCIP, "Figura_41", ancho = 16, alto = 22)

## ---- Neoplasia: exclusión de neoplasias del explante y modelo dependiente del tiempo ----
explante <- !is.na(pacientes$F_onco) & as.Date(pacientes$F_onco) == as.Date(pacientes$Ftx)
pacientes <- pacientes %>%
  mutate(onco       = ifelse(explante, 0L, onco),
         F_onco     = if_else(explante, as.Date(NA), as.Date(F_onco)),
         onco_tipo  = ifelse(explante, NA, onco_tipo),
         onco_tipo2 = ifelse(explante, NA, onco_tipo2),
         onco_det   = ifelse(explante, NA, onco_det))

res_onco <- list(
  gB   = analisis_td("gbm",   "onco", "F_onco", "Neoplasia global | gB mixtos"),
  gH   = analisis_td("ghm",   "onco", "F_onco", "Neoplasia global | gH mixtos"),
  gBgH = analisis_td("gbghm", "onco", "F_onco", "Neoplasia global | gB/gH mixtos")
)

## Descomentar para ver los nueve modelos de cada clasificación
# analisis_td("gbm",   "onco", "F_onco", "Neoplasia global | gB mixtos", solo_principal = FALSE)
# analisis_td("ghm",   "onco", "F_onco", "Neoplasia global | gH mixtos", solo_principal = FALSE)
# analisis_td("gbghm", "onco", "F_onco", "Neoplasia global | gB/gH mixtos", solo_principal = FALSE)

## ---- Figura 42. Supervivencia libre de neoplasia según genotipos mixtos ----
## Mismo estilo que la Figura 42 (funciones km_extendida / figura_km_2x2 de [7.12] [b2])
fig_km_onco <- figura_km_2x2("onco", "F_onco", "Supervivencia libre de neoplasia", ylim = c(0, 1))
print(fig_km_onco)

## Exportar a SVG para una página A4 en Word (16 x 22 cm)
# [guardado antiguo] ggsave("~/Desktop/Figuras_tesis/Figura_KMext_neoplasia_2x2.svg", fig_km_onco,
# [guardado antiguo]        device = svglite::svglite, fix_text_size = FALSE,
# [guardado antiguo]        width = 16, height = 22, units = "cm", bg = "white")

guardar_fig(fig_km_onco, "Figura_42", ancho = 16, alto = 22)

## ---- Tumor sólido no cutáneo: modelo dependiente del tiempo ----
res_solido <- list(
  gB   = analisis_td("gbm",   "solido", "F_onco", "Tumor sólido no cutáneo | gB mixtos"),
  gH   = analisis_td("ghm",   "solido", "F_onco", "Tumor sólido no cutáneo | gH mixtos"),
  gBgH = analisis_td("gbghm", "solido", "F_onco", "Tumor sólido no cutáneo | gB/gH mixtos")
)

## Descomentar para ver los nueve modelos de cada clasificación
# analisis_td("gbm",   "solido", "F_onco", "Tumor sólido no cutáneo | gB mixtos", solo_principal = FALSE)
# analisis_td("ghm",   "solido", "F_onco", "Tumor sólido no cutáneo | gH mixtos", solo_principal = FALSE)
# analisis_td("gbghm", "solido", "F_onco", "Tumor sólido no cutáneo | gB/gH mixtos", solo_principal = FALSE)

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

## ---- Mortalidad: modelo dependiente del tiempo ----
## Ajuste homogéneo con el resto de desenlaces: edad al trasplante + perfil serológico (risk_f)
res_exitus <- list(
  gB   = analisis_td("gbm",   "exitus", "Fexitus", "Mortalidad | gB mixtos"),
  gH   = analisis_td("ghm",   "exitus", "Fexitus", "Mortalidad | gH mixtos"),
  gBgH = analisis_td("gbghm", "exitus", "Fexitus", "Mortalidad | gB/gH mixtos")
)

## Descomentar para ver los nueve modelos de cada clasificación
# analisis_td("gbm",   "exitus", "Fexitus", "Mortalidad | gB mixtos", solo_principal = FALSE)
# analisis_td("ghm",   "exitus", "Fexitus", "Mortalidad | gH mixtos", solo_principal = FALSE)
# analisis_td("gbghm", "exitus", "Fexitus", "Mortalidad | gB/gH mixtos", solo_principal = FALSE)

## ---- Figura 44. Supervivencia global según genotipos mixtos ----
## Mismo estilo que la Figura 42 (funciones km_extendida / figura_km_2x2 de [7.12] [b2])
fig_km_exitus <- figura_km_2x2("exitus", "Fexitus", "Supervivencia global", ylim = c(0, 1))
print(fig_km_exitus)

## Exportar a SVG para una página A4 en Word (16 x 22 cm)
# [guardado antiguo] ggsave("~/Desktop/Figuras_tesis/Figura_KMext_mortalidad_2x2.svg", fig_km_exitus,
# [guardado antiguo]        device = svglite::svglite, fix_text_size = FALSE,
# [guardado antiguo]        width = 16, height = 22, units = "cm", bg = "white")

guardar_fig(fig_km_exitus, "Figura_44", ancho = 16, alto = 22)

## ---- Tabla 21. Riesgo de neoplasia y tumor sólido no cutáneo (tabla resumen de todos los desenlaces) ----
# Va al final porque reúne los modelos de todos los desenlaces (rechazo, DCIP,
# neoplasia, tumor sólido, infecciones oportunistas y mortalidad).
tabla_principal <- bind_rows(
  lapply(list(RA = res_RA, DCIP = res_DCIP, onco = res_onco, solido = res_solido,
              inf_opo = res_inf_opo, exitus = res_exitus), function(r)
                bind_rows(lapply(names(r), function(k)
                  cbind(Clasificacion = k, r[[k]]$res[grepl("PRINCIPAL", r[[k]]$res$Modelo), ])))),
  .id = "Desenlace") %>%
  select(Desenlace, Clasificacion, n, eventos, eventos_expuesto, HR, IC95, p)

etiq_des <- c(RA      = "Rechazo agudo",
              DCIP    = "Disfunción crónica del injerto pulmonar",
              onco    = "Neoplasia global",
              solido  = "Tumor sólido no cutáneo",
              inf_opo = "Infección oportunista no CMV",
              exitus  = "Mortalidad")
etiq_cla <- c(gB = "gB", gH = "gH", gBgH = "gB/gH")

tab_ft <- tabla_principal %>%
  mutate(
    Desenlace     = etiq_des[Desenlace],
    Clasificacion = etiq_cla[Clasificacion],
    `HR ajustada (IC 95 %)` = paste0(HR, " (", IC95, ")"),
    p = ifelse(grepl("^<", p), p, sub("^0,", "0,", p))
  ) %>%
  select(Desenlace, Clasificacion, n, eventos, eventos_expuesto,
         `HR ajustada (IC 95 %)`, p)

ft_principal <- flextable(tab_ft) %>%
  set_header_labels(
    Desenlace        = "Desenlace",
    Clasificacion    = "Clasificación",
    n                = "n",
    eventos          = "Eventos",
    eventos_expuesto = "Eventos en expuestos",
    p                = "p"
  ) %>%
  merge_v(j = "Desenlace") %>%
  valign(j = "Desenlace", valign = "top") %>%
  align(j = c("n", "eventos", "eventos_expuesto", "p"), align = "center", part = "all") %>%
  align(j = "HR ajustada (IC 95 %)", align = "center", part = "all") %>%
  bold(part = "header") %>%
  bold(i = ~ grepl("^<", p) | suppressWarnings(as.numeric(sub(",", ".", p))) < 0.05,
       j = c("HR ajustada (IC 95 %)", "p")) %>%
  flextable::font(fontname = "Source Sans Pro", part = "all") %>%
  border_remove() %>%
  hline_top(part = "header", border = fp_border(width = 1)) %>%
  hline_bottom(part = "header", border = fp_border(width = 0.5)) %>%
  hline_bottom(part = "body", border = fp_border(width = 1)) %>%
  hline(i = ~ !duplicated(Desenlace, fromLast = TRUE),
        border = fp_border(width = 0.25, color = "grey60")) %>%
  add_footer_lines(paste0(
    "Modelo de Cox con entrada diferida en el primer episodio con genotipado válido ",
    "del gen estudiado y exposición dependiente del tiempo (paso a mixto en la primera ",
    "detección de genotipos mixtos). Ajustado por edad al trasplante y perfil serológico ",
    "D/R (mortalidad: serología D/R en tres categorías). n: pacientes libres del desenlace ",
    "en la entrada; HR: hazard ratio; IC: intervalo de confianza. En negrita, p < 0,05.")) %>%
  fontsize(size = 8, part = "footer") %>%
  italic(part = "footer") %>%
  autofit()

ft_principal

## Exportar a Word (opcional):
# save_as_docx(ft_principal, path = "~/Desktop/Tabla_HR_principales.docx")

guardar_tabla(ft_principal, "Tabla_21")


# ---- Fin ----
cat("\nScript completado. Figuras y tablas en:", normalizePath(dir_salida), "\n")
