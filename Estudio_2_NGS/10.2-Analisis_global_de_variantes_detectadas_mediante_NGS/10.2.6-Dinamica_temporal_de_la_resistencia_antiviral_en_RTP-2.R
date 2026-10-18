# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.6. Dinámica temporal de la resistencia antiviral en RTP-2
# Extracto del script completo (ver README.md del estudio)

# ---- 10.2.6. Dinámica temporal de la resistencia antiviral en RTP-2 ----

## ---- Preparación común: cargas virales, tratamientos y variantes de RTP-2 ----
# Guarda las figuras con guardar_rtp2() (png + svg) en la carpeta de salidas.
## Cambios respecto a la versión anterior:
##  - Las MRA y las VND se toman de all_data_txp (misma tabla que los anexos), no de listas
##    escritas a mano (vsi_extra). Desaparecen dia_remap, F396L, D711Y, F300L, P704S, P712H,
##    Q323R y D605E (ahora PN). T700A no forma parte de all_data_txp (no superó el filtro
##    de cobertura) y no se representa.
##  - Las muestras sin resultado de UL54 o UL56 (cruces) se detectan desde los datos.
##  - Las VND representadas se seleccionan con un criterio explícito: localización en las
##    regiones funcionales asociadas a resistencia (regiones_resist); el resto queda en el anexo.
##  - Los tratamientos y las cargas virales (df_rtp2, trat_df) siguen siendo datos clínicos
##    introducidos manualmente.


fecha_txp_rtp2 <- as.Date("2022-06-17")
nhc_rtp2       <- "NHC_RTP-2"

# ── Colores y tamaños ─────────────────────────────────────────────────────────
col_GCV <- "#3D7A5A"; col_VGCV <- "#2E86AB"; col_MBV <- "#B8860B"; col_LEF <- "#C0652B"
col_FOS <- "#9B2335"; col_LMV  <- "#7B5EA7"; col_CMVIg <- "#A8C8E0"
colores_fill  <- c("UL97" = alpha("#4472C4", 0.80), "UL54" = alpha("#E6A817", 0.80), "UL56" = alpha("#E8735A", 0.80))
colores_color <- c("UL97" = "#4472C4", "UL54" = "#C8880A", "UL56" = "#C94A30")
tam_5  <- c("\u22655%" = 2, "\u226510%" = 3, "\u226515%" = 5, "\u226520%" = 7)
tam_10 <- c("\u226510%" = 3, "\u226515%" = 5, "\u226520%" = 7)

# ── Eje temporal (índice 1..51 sobre los días con determinación de CV) ───────
niveles_dia <- c(
  "1","26","42","49","55","62","69","75","84","94","108","122",
  "159","181","187","201","214","227","244","271","279","291","304",
  "361","371","382","388","390","395","402","409","418","423","433",
  "445","458","466","479","497","531","571","630","668","697","709",
  "742","774","801","850","872","909")
n_niv          <- length(niveles_dia)
pos_map        <- tibble(dia = as.numeric(niveles_dia), pos = seq_along(niveles_dia))
interpolar_pos <- function(d) approx(pos_map$dia, pos_map$pos, xout = d, rule = 2)$y

# ── Datos clínicos RTP-2 (carga viral) ───────────────────────────────────────
df_rtp2 <- tibble(
  dia = as.numeric(niveles_dia),
  log_cv = c(0.000, 0.000,
             4.278, 4.842, 5.042, 4.389, 3.891, 3.461, 3.234, 2.876, 0.000, 2.241,
             2.740, 3.557, 3.389, 3.461, 3.042, 2.885, 3.346,
             2.699, 3.279, 3.653, 4.008,
             3.741, 4.276, 4.610, 4.521, 4.422, 4.157, 4.312,
             0.000, 3.828, 3.872, 3.638,
             3.752, 4.381, 4.362, 4.364, 4.173,
             4.203, 3.712, 3.328, 2.965, 2.845, 2.823, 2.823,
             3.215, 3.034, 2.645, 3.002, 0.000)
)

# ── Muestras NGS de RTP-2 desde los datos: día y número de orden ─────────────
muestras_rtp2 <- all_data_txp %>%
  filter(NHC == nhc_rtp2) %>%
  mutate(Freg = as.Date(Freg), dia = as.numeric(Freg - fecha_txp_rtp2)) %>%
  distinct(Npet, dia) %>% arrange(dia) %>%
  mutate(muestra = row_number())
stopifnot(all(muestras_rtp2$dia %in% as.numeric(niveles_dia)))   # los días deben estar en el eje

df_rtp2 <- df_rtp2 %>%
  left_join(muestras_rtp2 %>% select(dia, muestra), by = "dia") %>%
  mutate(secuenciada = !is.na(muestra),
         x_pos       = interpolar_pos(dia),
         tipo_punto  = ifelse(secuenciada, "Muestra secuenciada", "Muestra no secuenciada"))

# Muestras sin resultado por gen (cruces): analizadas para RTP-2 pero sin ninguna fila del gen
genes_por_muestra <- all_data_txp %>%
  filter(NHC == nhc_rtp2) %>%
  mutate(Freg = as.Date(Freg), dia = as.numeric(Freg - fecha_txp_rtp2)) %>%
  distinct(dia, Protein)
dias_sin_gen <- function(gen) setdiff(muestras_rtp2$dia, genes_por_muestra$dia[genes_por_muestra$Protein == gen])

# ── Variantes de RTP-2 desde los datos ───────────────────────────────────────
variantes_rtp2 <- all_data_txp %>%
  filter(NHC == nhc_rtp2, class %in% c("MRA", "VSI", "NEW"), Prevalence >= 5) %>%
  mutate(
    Freg  = as.Date(Freg),
    dia   = as.numeric(Freg - fecha_txp_rtp2),
    x_pos = interpolar_pos(dia),
    gen_color = factor(as.character(Protein), levels = c("UL97", "UL54", "UL56")),
    clase_var = factor(ifelse(class == "NEW", "VND", class), levels = c("MRA", "VSI", "VND")),
    Prevalence_group = factor(case_when(
      Prevalence >= 20 ~ "\u226520%", Prevalence >= 15 ~ "\u226515%",
      Prevalence >= 10 ~ "\u226510%", TRUE ~ "\u22655%"),
      levels = c("\u22655%", "\u226510%", "\u226515%", "\u226520%")),
    Mutation = as.character(Mutation)
  ) %>%
  select(dia, x_pos, Protein, gen_color, Mutation, clase_var, Prevalence, Prevalence_group)

cat("RTP-2: muestras", nrow(muestras_rtp2), "| variantes (filas)", nrow(variantes_rtp2),
    "| únicas", n_distinct(variantes_rtp2$Mutation), "\n")
print(variantes_rtp2 %>% count(gen_color, clase_var))

# ── Regiones funcionales asociadas a resistencia (definidas a priori) ─────────
## Ajustar a los intervalos que se declaran en Material y métodos.
regiones_resist <- tibble(
  Protein = c("UL97", "UL54", "UL56"),
  ini     = c(335, 300, 230),
  fin     = c(665, 1000, 370)
)
variantes_rtp2 <- variantes_rtp2 %>%
  mutate(Position = as.numeric(gsub("[^0-9]", "", Mutation))) %>%
  left_join(regiones_resist, by = "Protein") %>%
  mutate(en_region = Position >= ini & Position <= fin) %>%
  select(-ini, -fin)

# Subconjunto para las figuras: todas las MRA + VND situadas en regiones de resistencia
variantes_rtp2_fig <- variantes_rtp2 %>% filter(clase_var == "MRA" | en_region)
cat("Variantes representadas:", n_distinct(variantes_rtp2_fig$Mutation), "de",
    n_distinct(variantes_rtp2$Mutation), "\n")
print(variantes_rtp2 %>% distinct(Protein, Mutation, clase_var, en_region) %>%
        count(Protein, clase_var, en_region))

# ── Tratamientos ─────────────────────────────────────────────────────────────
trat_df <- tibble(
  label  = c("GCV", "VGCV", "GCV", "MBV", "LEF", "FOS", "LMV", "CMV-Ig"),
  xmin_d = c( 42,    69,    271,   291,   361,   382,   409,    62),
  xmax_d = c( 68,   270,    290,   381,   381,   408,   432,   909),
  fila   = c("sup", "sup",  "sup", "sup", "top", "sup", "sup", "inf"),
  fill_col = c(col_GCV, col_VGCV, col_GCV, col_MBV, col_LEF, col_FOS, col_LMV, col_CMVIg),
  txt_col  = c(rep("white", 7), "#1A3A52")
)

make_rounded_bar <- function(xmin, xmax, ymid, h, fill_col, id, n_arc = 40) {
  r <- min(h / 2, (xmax - xmin) / 2)
  theta_l <- seq(pi / 2, 3 * pi / 2, length.out = n_arc)
  theta_r <- seq(-pi / 2, pi / 2, length.out = n_arc)
  tibble(x = c(xmin + r + r * cos(theta_l), xmax - r + r * cos(theta_r)),
         y = c(ymid + r * sin(theta_l), ymid + r * sin(theta_r)),
         group = id, fill_col = fill_col)
}

# ── Función de figura ────────────────────────────────────────────────────────
figura_rtp2 <- function(datos_mut, tamaños, y_mut_base = -0.4, y_mut_step = -0.55,
                        y_inf_trat = 5.60, y_sup_trat = 6.00, y_top_trat = 6.38,
                        bar_h = 0.35, y_lim_sup = 7, base_size = 13,
                        tam_etq = 3.8, tam_num = 3.5) {
  
  datos_mut <- datos_mut %>%
    mutate(gen_orden = as.integer(gen_color)) %>%
    arrange(gen_orden, clase_var, Mutation) %>%      # MRA antes que VND dentro de cada gen
    mutate(Mutation = factor(Mutation, levels = unique(Mutation)))
  mut_levels <- levels(datos_mut$Mutation)
  mut_y_map  <- setNames(y_mut_base + (seq_along(mut_levels) - 1) * y_mut_step, mut_levels)
  datos_mut  <- datos_mut %>% mutate(y_pos = mut_y_map[as.character(Mutation)])
  
  mut_tick_df <- datos_mut %>% distinct(Mutation, gen_color) %>%
    mutate(label = as.character(Mutation), y_pos = mut_y_map[label],
           col_text = colores_color[as.character(gen_color)])
  
  # Cruces: muestras sin resultado del gen, en las filas de ese gen
  cruces <- bind_rows(lapply(c("UL54", "UL56"), function(g) {
    muts <- mut_tick_df$label[mut_tick_df$gen_color == g]
    if (length(muts) == 0) return(NULL)
    expand_grid(dia = dias_sin_gen(g), Mutation = muts) %>%
      mutate(x_pos = interpolar_pos(dia), y_pos = mut_y_map[Mutation])
  }))
  
  y_trat <- c(sup = y_sup_trat, top = y_top_trat, inf = y_inf_trat)
  trat <- trat_df %>%
    mutate(y_row = y_trat[fila], xmin_pos = interpolar_pos(xmin_d), xmax_pos = interpolar_pos(xmax_d),
           x_mid = (xmin_pos + xmax_pos) / 2, ancho = xmax_pos - xmin_pos,
           txt_size = case_when(ancho > 3.5 ~ 3.4, ancho > 1.8 ~ 3.0, TRUE ~ 2.6))
  polys_trat <- pmap_dfr(list(trat$xmin_pos, trat$xmax_pos, trat$y_row, rep(bar_h, nrow(trat)),
                              trat$fill_col, seq_len(nrow(trat))), make_rounded_bar)
  
  mra_data <- datos_mut %>% filter(clase_var == "MRA")
  vnd_data <- datos_mut %>% filter(clase_var != "MRA")
  genes_presentes <- levels(droplevels(datos_mut$gen_color))
  y_lim_inf <- min(mut_y_map) - 0.4
  y_corte   <- y_inf_trat - 0.25
  
  ggplot() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = y_corte, ymax = y_lim_sup, fill = "white", color = NA) +
    geom_vline(xintercept = seq_along(niveles_dia), color = "#DDDDDD", linewidth = 0.25) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = y_corte, ymax = y_lim_sup, fill = "white", color = NA) +
    geom_polygon(data = polys_trat, aes(x, y, group = group), fill = polys_trat$fill_col, color = NA) +
    geom_text(data = trat, aes(x_mid, y_row, label = label), color = trat$txt_col,
              size = trat$txt_size, fontface = "bold") +
    annotate("text", x = rep(0.55, 5), y = 1:5, label = as.character(1:5), hjust = 1, size = 3.5, color = "grey35") +
    annotate("segment", x = rep(0.75, 5), xend = rep(0.97, 5), y = 1:5, yend = 1:5, color = "grey70", linewidth = 0.25) +
    geom_line(data = df_rtp2, aes(x_pos, log_cv, group = 1), linewidth = 0.6, color = "darkred") +
    geom_point(data = df_rtp2, aes(x_pos, log_cv, color = tipo_punto), size = 2.5) +
    geom_text(data = df_rtp2 %>% filter(!is.na(muestra)), aes(x_pos, log_cv, label = muestra),
              vjust = -0.8, size = tam_num, color = "darkred", fontface = "bold") +
    geom_segment(data = mut_tick_df, aes(x = -0.3, xend = n_niv + 0.5, y = y_pos, yend = y_pos),
                 color = "#DDDDDD", linewidth = 0.25) +
    geom_text(data = mut_tick_df, aes(x = -0.4, y = y_pos, label = label), color = mut_tick_df$col_text,
              hjust = 1, size = tam_etq, fontface = "bold") +
    annotate("text", x = -1.2, y = 3.0, label = "Log CV (UI/mL)", angle = 90, size = 4.0,
             color = "grey30", fontface = "bold", hjust = 0.5, vjust = 0) +
    geom_point(data = mra_data, aes(x_pos, y_pos, size = Prevalence_group, color = gen_color, fill = gen_color),
               shape = 21, stroke = 1, alpha = 0.80) +
    geom_point(data = vnd_data, aes(x_pos, y_pos, size = Prevalence_group, color = gen_color),
               shape = 21, fill = "white", stroke = 1) +
    geom_point(data = cruces, aes(x_pos, y_pos), shape = 4, size = 2.4, stroke = 0.8, color = "grey35") +
    scale_color_manual(values = c(colores_color, "Muestra secuenciada" = "darkred", "Muestra no secuenciada" = "black"),
                       breaks = c(genes_presentes, "Muestra secuenciada", "Muestra no secuenciada"), name = NULL) +
    scale_fill_manual(values = colores_fill, guide = "none") +
    scale_size_manual(values = tamaños, name = "Umbral") +
    scale_x_continuous(limits = c(-3.5, n_niv + 0.5), breaks = seq_along(niveles_dia),
                       labels = niveles_dia, expand = c(0, 0)) +
    scale_y_continuous(limits = c(y_lim_inf, y_lim_sup), breaks = 0:6, labels = NULL, expand = c(0, 0)) +
    guides(size = guide_legend(order = 1, nrow = 1,
                               override.aes = list(shape = 21, color = "grey50", fill = "grey70", stroke = 1)),
           color = guide_legend(order = 2, nrow = 1,
                                override.aes = list(shape = c(rep(21, length(genes_presentes)), 16, 16),
                                                    size  = c(rep(5, length(genes_presentes)), 4, 4),
                                                    fill  = c(colores_fill[genes_presentes], "darkred", "black"),
                                                    stroke = c(rep(1, length(genes_presentes)), 0, 0)))) +
    labs(x = "D\u00edas postrasplante", y = NULL) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(),
          axis.title.x = element_text(face = "bold", margin = margin(t = 12)),
          legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal",
          legend.title = element_text(size = 11, face = "bold"), legend.text = element_text(size = 11),
          legend.spacing.x = unit(0.25, "cm"),
          plot.margin = margin(5, 15, 5, 10),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))
}

guardar_rtp2 <- function(p, nombre, ancho, alto) {
  ggsave(file.path(dir_salida, paste0(nombre, ".png")), p,
         width = ancho, height = alto, units = "cm", dpi = 300, bg = "white")
  ggsave(file.path(dir_salida, paste0(nombre, ".svg")), p,
         device = svglite::svglite, fix_text_size = FALSE,
         width = ancho, height = alto, units = "cm", bg = "white")
}

## ---- Figura 49. Evolución de la carga viral y detección de MRA ----
mra_rtp2 <- variantes_rtp2 %>% filter(clase_var == "MRA")
cat("MRA en RTP-2:\n"); mra_rtp2 %>% select(dia, Protein, Mutation, Prevalence) %>% arrange(Protein, Mutation, dia) %>% print(n = Inf)

p_mra <- figura_rtp2(mra_rtp2, tam_5, y_mut_base = -0.4,
                     y_inf_trat = 6.28, y_sup_trat = 6.58, y_top_trat = 6.84, bar_h = 0.24, y_lim_sup = 7.25)
grid::grid.newpage(); print(p_mra)
guardar_rtp2(p_mra, "Figura_49", 32, 15)

## ---- Figura 50. Evolución de la carga viral y detección de MRA y VND ----
p_5 <- figura_rtp2(variantes_rtp2_fig, tam_5,
                   y_mut_base = -0.1, y_mut_step = -0.5,
                   base_size = 16, tam_etq = 5.2, tam_num = 4.5) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
        axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 12)),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text  = element_text(size = 14),
        legend.key.size = unit(0.9, "cm"))
grid::grid.newpage(); print(p_5)
guardar_rtp2(p_5, "Figura_50", 42, 24)

## ---- Complementaria. MRA y VND en RTP-2 al umbral ≥ 10% ----
variantes_rtp2_10 <- variantes_rtp2_fig %>% filter(Prevalence >= 10)
p_10 <- figura_rtp2(variantes_rtp2_10, tam_10, y_mut_base = -0.8)
grid::grid.newpage(); print(p_10)
guardar_rtp2(p_10, "Complementaria_RTP2_MRA_VND_10pct", 35, 21)


