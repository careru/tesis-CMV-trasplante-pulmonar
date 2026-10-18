# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: ANEXO.3. RTP-2 evolutivo, MRA, VSI y VND ≥ 10% (versión anterior)
# Extracto del script completo (ver README.md del estudio)

## ---- RTP-2 evolutivo, MRA, VSI y VND ≥ 10% (versión anterior) ----
# ══════════════════════════════════════════════════════════════════════════════
#  COLORES
# ══════════════════════════════════════════════════════════════════════════════
col_GCV   <- "#3D7A5A"
col_VGCV  <- "#2E86AB"
col_MBV   <- "#B8860B"
col_LEF   <- "#C0652B"
col_FOS   <- "#9B2335"
col_LMV   <- "#7B5EA7"
col_CMVIg <- "#A8C8E0"

colores_fill  <- c("UL97" = alpha("#4472C4", 0.80),
                   "UL54" = alpha("#E6A817", 0.80),
                   "UL56" = alpha("#E8735A", 0.80))
colores_color <- c("UL97" = "#4472C4",
                   "UL54" = "#C8880A",
                   "UL56" = "#C94A30")

# ── Solo a partir de ≥10% ────────────────────────────────────────────────────
tamaños <- c("\u226510%" = 3, "\u226515%" = 5, "\u226520%" = 7)

# ══════════════════════════════════════════════════════════════════════════════
#  EJE TEMPORAL — índice numérico 1..51
# ══════════════════════════════════════════════════════════════════════════════
niveles_dia <- c(
  "1","26","42","49","55","62","69","75","84","94","108","122",
  "159","181","187","201","214","227","244","271","279","291","304",
  "361","371","382","388","390","395","402","409","418","423","433",
  "445","458","466","479","497","531","571","630","668","697","709",
  "742","774","801","850","872","909"
)
n_niv   <- length(niveles_dia)

pos_map        <- tibble(dia = as.numeric(niveles_dia), pos = seq_along(niveles_dia))
interpolar_pos <- function(d) approx(pos_map$dia, pos_map$pos, xout = d, rule = 2)$y

pos_evr <- interpolar_pos(279)

# ══════════════════════════════════════════════════════════════════════════════
#  DATOS CLÍNICOS RTP-2
# ══════════════════════════════════════════════════════════════════════════════
df_rtp2 <- tibble(
  dia = c(
    1,   26,  42,  49,  55,  62,  69,  75,  84,  94,  108, 122,
    159, 181, 187, 201, 214, 227, 244, 271, 279, 291, 304,
    361, 371, 382, 388, 390, 395, 402, 409, 418, 423, 433,
    445, 458, 466, 479, 497, 531, 571, 630, 668, 697, 709,
    742, 774, 801, 850, 872, 909
  ),
  log_cv = c(
    0.000, 0.000,
    4.278, 4.842, 5.042, 4.389, 3.891, 3.461, 3.234, 2.876, 0.000, 2.241,
    2.740, 3.557, 3.389, 3.461, 3.042, 2.885, 3.346,
    2.699, 3.279, 3.653, 4.008,
    3.741, 4.276, 4.610, 4.521, 4.422, 4.157, 4.312,
    0.000, 3.828, 3.872, 3.638,
    3.752, 4.381, 4.362, 4.364, 4.173,
    4.203, 3.712, 3.328, 2.965, 2.845, 2.823, 2.823,
    3.215, 3.034, 2.645, 3.002, 0.000
  ),
  muestra = c(
    NA, NA,
    1, 2, 3, 4, NA, NA, NA, NA, NA, NA,
    NA, NA, NA, NA, NA, NA, NA,
    NA, NA, NA, 5,
    6, 7, 8, 9, 10, 11, NA,
    NA, NA, NA, NA,
    NA, 12, 13, 14, 15,
    NA, NA, NA, NA, NA, NA, NA,
    NA, NA, NA, NA, NA
  ),
  secuenciada = c(
    F, F,
    T, T, T, T, F, F, F, F, F, F,
    F, F, F, F, F, F, F,
    F, F, F, T,
    T, T, T, T, T, T, F,
    F, F, F, F,
    F, T, T, T, T,
    F, F, F, F, F, F, F,
    F, F, F, F, F
  )
) %>%
  mutate(
    x_pos      = interpolar_pos(dia),
    tipo_punto = ifelse(secuenciada, "Muestra secuenciada", "Muestra no secuenciada")
  )

# ══════════════════════════════════════════════════════════════════════════════
#  DATOS DE MUTACIONES
# ══════════════════════════════════════════════════════════════════════════════
dia_remap <- c(
  "41"="42","48"="49","54"="55","61"="62",
  "303"="304","360"="361","370"="371",
  "381"="382","387"="388","389"="390","394"="395",
  "457"="402","465"="466","478"="479","496"="497"
)
remap_dia <- function(d) {
  s <- as.character(d)
  ifelse(s %in% names(dia_remap), dia_remap[s], s)
}

vsi_extra <- tibble(
  Mutation = c(
    "F396L","D605E",
    rep("A497P", 12),
    "W424*","G340C","D481N",
    "Y281H","M697I",
    "F300L","Q323R","M491T","P704S","D711Y","P712H","D854G","S913G",
    "M491T","F702L","F817L"
  ),
  dia = c(
    54, 465,
    48,54,303,360,370,381,387,389,394,465,478,496,
    61,478,478,
    478,478,
    61,61,61,61,61,61,61,61,
    303,303,303
  ),
  Prevalence = c(
    12.9, 5.7,
    5.13,5.24,6.39,5.51,6.89,5.77,7.71,6.0,5.42,5.86,5.16,NA,
    23.9,6.2,5.1,
    8.0,8.4,
    7.9,25.0,10.2,7.3,8.6,8.2,19.6,15.1,
    6.9,5.2,5.2
  ),
  Protein = c(
    "UL54-VSI","UL97-VSI",
    rep("UL97-VND",12),
    "UL97-VND","UL97-VND","UL97-VND",
    "UL56-VND","UL56-VND",
    rep("UL54-VND",8),
    rep("UL54-VND",3)
  )
) %>%
  filter(!is.na(Prevalence)) %>%
  mutate(
    Prevalence_group = factor(case_when(
      Prevalence>=20~"\u226520%",Prevalence>=15~"\u226515%",
      Prevalence>=10~"\u226510%",TRUE~"\u22655%"),
      levels=c("\u22655%","\u226510%","\u226515%","\u226520%")),
    Protein  = factor(Protein,
                      levels=c("UL97","UL54","UL56",
                               "UL54-VSI","UL97-VSI",
                               "UL97-VND","UL54-VND","UL56-VND")),
    dia_str  = remap_dia(dia),
    x_pos    = interpolar_pos(as.numeric(dia_str))
  )

mra_rtp2_base <- mra_rtp2 %>%
  filter(!is.na(Mutation), !is.na(Prevalence_group)) %>%
  mutate(
    dia_str = remap_dia(dia),
    x_pos   = interpolar_pos(as.numeric(dia_str))
  )

mra_rtp2_all <- bind_rows(mra_rtp2_base, vsi_extra)

mra_rtp2_plot <- mra_rtp2_all %>%
  filter(!is.na(Mutation), !is.na(Prevalence_group)) %>%
  mutate(
    Mutation  = as.character(Mutation),
    gen_color = case_when(
      Protein %in% c("UL97","UL97-VSI","UL97-VND") ~ "UL97",
      Protein %in% c("UL54","UL54-VSI","UL54-VND") ~ "UL54",
      Protein %in% c("UL56","UL56-VND")             ~ "UL56"
    ),
    gen_color = factor(gen_color, levels=c("UL97","UL54","UL56")),
    clase_var = case_when(
      Protein %in% c("UL97","UL54","UL56")             ~ "MRA",
      Protein %in% c("UL97-VSI","UL54-VSI")            ~ "VSI",
      Protein %in% c("UL97-VND","UL54-VND","UL56-VND") ~ "VND"
    ),
    clase_var = factor(clase_var, levels=c("MRA","VSI","VND")),
    gen_orden = case_when(
      gen_color=="UL97" ~ 1L,
      gen_color=="UL54" ~ 2L,
      gen_color=="UL56" ~ 3L
    )
  ) %>%
  # ── Filtrar: solo prevalencia ≥10% ────────────────────────────────────────
  filter(Prevalence_group %in% c("\u226510%","\u226515%","\u226520%")) %>%
  arrange(gen_orden, Mutation) %>%
  mutate(Mutation = factor(Mutation, levels=unique(Mutation)))

# ── Posiciones Y para mutaciones ─────────────────────────────────────────────
mut_levels <- levels(mra_rtp2_plot$Mutation)
y_mut_base <- -0.8
y_mut_step <- -0.55
mut_y_map  <- setNames(
  y_mut_base + (seq_along(mut_levels)-1) * y_mut_step,
  mut_levels
)
mra_rtp2_plot <- mra_rtp2_plot %>%
  mutate(y_pos = mut_y_map[as.character(Mutation)])

mut_gen_map <- mra_rtp2_plot %>%
  select(Mutation, gen_color) %>%
  distinct() %>%
  mutate(Mutation = as.character(Mutation))

x_mut_label <- -0.3

# Cruces UL54 no analizado (muestras 6 y 14)
cruces_ul54_muts <- mra_rtp2_plot %>%
  filter(Protein %in% c("UL54","UL54-VSI","UL54-VND")) %>%
  pull(Mutation) %>% as.character() %>% unique()
cruces_m6_m14 <- bind_rows(
  tibble(x_pos=interpolar_pos(361), y_pos=mut_y_map[cruces_ul54_muts], Mutation=cruces_ul54_muts),
  tibble(x_pos=interpolar_pos(479), y_pos=mut_y_map[cruces_ul54_muts], Mutation=cruces_ul54_muts)
)

# Cruces UL56 para muestras 1–7
ul56_muts <- mra_rtp2_plot %>%
  filter(gen_color=="UL56") %>%
  pull(Mutation) %>% as.character() %>% unique()
dias_m1_7 <- c(42, 49, 55, 62, 304, 361, 371)
cruces_ul56_m1_7 <- bind_rows(lapply(dias_m1_7, function(d) {
  tibble(x_pos=interpolar_pos(d),
         y_pos=mut_y_map[ul56_muts],
         Mutation=ul56_muts)
}))

mra_data <- mra_rtp2_plot %>% filter(clase_var=="MRA")
vsi_data <- mra_rtp2_plot %>% filter(clase_var %in% c("VSI","VND"))

mut_tick_df <- tibble(
  label = mut_levels,
  y_pos = mut_y_map[mut_levels]
) %>%
  left_join(mut_gen_map, by=c("label"="Mutation")) %>%
  mutate(col_text = case_when(
    gen_color=="UL97" ~ "#4472C4",
    gen_color=="UL54" ~ "#C8880A",
    gen_color=="UL56" ~ "#C94A30"
  ))

# ══════════════════════════════════════════════════════════════════════════════
#  BARRAS DE ANTIVIRAL
# ══════════════════════════════════════════════════════════════════════════════
y_inf_trat <- 5.60; y_sup_trat <- 6.00; y_top_trat <- 6.38; bar_h <- 0.35

trat_df <- tibble(
  label    = c("GCV",  "VGCV", "GCV",  "MBV",  "LEF",  "FOS",  "LMV",  "CMV-Ig"),
  xmin_d   = c( 42,     69,     271,    291,    361,    382,    409,    62),
  xmax_d   = c( 68,    270,     290,    381,    381,    408,    432,   909),
  y_row    = c(y_sup_trat,y_sup_trat,y_sup_trat,y_sup_trat,
               y_top_trat,y_sup_trat,y_sup_trat,y_inf_trat),
  fill_col = c(col_GCV,col_VGCV,col_GCV,col_MBV,col_LEF,col_FOS,col_LMV,col_CMVIg),
  txt_col  = c("white","white","white","white","white","white","white","#1A3A52")
) %>%
  mutate(
    xmin_pos = interpolar_pos(xmin_d),
    xmax_pos = interpolar_pos(xmax_d),
    x_mid    = (xmin_pos + xmax_pos) / 2,
    ancho    = xmax_pos - xmin_pos,
    txt_size = case_when(ancho>3.5~3.4, ancho>1.8~3.0, TRUE~2.6)
  )

make_rounded_bar <- function(xmin, xmax, ymid, h, fill_col, id, n_arc=40) {
  r <- min(h/2, (xmax-xmin)/2)
  theta_l <- seq(pi/2, 3*pi/2, length.out=n_arc)
  theta_r <- seq(-pi/2, pi/2,  length.out=n_arc)
  tibble(
    x        = c(xmin+r+r*cos(theta_l), xmax-r+r*cos(theta_r)),
    y        = c(ymid+r*sin(theta_l),   ymid+r*sin(theta_r)),
    group    = id,
    fill_col = fill_col
  )
}

polys_trat <- pmap_dfr(
  list(trat_df$xmin_pos, trat_df$xmax_pos, trat_df$y_row,
       rep(bar_h, nrow(trat_df)), trat_df$fill_col, seq_len(nrow(trat_df))),
  function(xmin,xmax,ymid,h,fill_col,id) make_rounded_bar(xmin,xmax,ymid,h,fill_col,id)
)

y_lim_inf <- min(mut_y_map) - 0.1
y_lim_sup <- 7

# ══════════════════════════════════════════════════════════════════════════════
#  FIGURA
# ══════════════════════════════════════════════════════════════════════════════
p <- ggplot() +
  
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=5.50, ymax=y_lim_sup,
           fill="white", color=NA) +
  
  geom_vline(xintercept=seq_along(niveles_dia),
             color="#DDDDDD", linewidth=0.25, ymin=0, ymax=1) +
  
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=5.35, ymax=y_lim_sup,
           fill="white", color=NA) +
  
  geom_polygon(data=polys_trat,
               aes(x=x, y=y, group=group, fill=I(fill_col)),
               color=NA, inherit.aes=FALSE) +
  geom_text(data=trat_df %>% mutate(r=row_number()),
            aes(x=x_mid, y=y_row, label=label),
            color=trat_df$txt_col,
            size=trat_df$txt_size, fontface="bold",
            hjust=0.5, vjust=0.5, inherit.aes=FALSE) +
  
  
  # ── Etiquetas numéricas eje Y (Log CV) ───────────────────────────────────
  annotate("text",
           x     = rep(0.55, 5),
           y     = 1:5,
           label = as.character(1:5),
           hjust = 1, size = 3.5,
           color = "grey35") +
  
  # ── Marcas de tick ────────────────────────────────────────────────────────
  annotate("segment",
           x    = rep(0.75, 5), xend = rep(0.97, 5),
           y    = 1:5,          yend = 1:5,
           color = "grey70", linewidth = 0.25) +
  
  geom_line(data=df_rtp2, aes(x=x_pos, y=log_cv, group=1),
            linewidth=0.6, color="darkred") +
  geom_point(data=df_rtp2,
             aes(x=x_pos, y=log_cv, color=tipo_punto, shape=tipo_punto),
             size=2.5) +
  geom_text(data=df_rtp2 %>% filter(!is.na(muestra)),
            aes(x=x_pos, y=log_cv, label=muestra),
            vjust=-0.8, size=3.5, color="darkred", fontface="bold") +
  
  geom_segment(data=mut_tick_df,
               aes(x=x_mut_label, xend=n_niv+0.5, y=y_pos, yend=y_pos),
               color="#DDDDDD", linewidth=0.25, inherit.aes=FALSE) +
  
  geom_text(data=mut_tick_df,
            aes(x=x_mut_label - 0.1, y=y_pos, label=label),
            color=mut_tick_df$col_text,
            hjust=1, size=3.8, fontface="bold",
            inherit.aes=FALSE) +
  
  annotate("text", x=-1.2, y=2.5,
           label="Log CV (UI/mL)", angle=90, size=4.0,
           color="grey30", fontface="bold", hjust=0.5, vjust=0) +
  
  geom_point(data=mra_data,
             aes(x=x_pos, y=y_pos, size=Prevalence_group,
                 color=gen_color, fill=gen_color),
             shape=21, stroke=1, alpha=0.80) +
  geom_point(data=vsi_data,
             aes(x=x_pos, y=y_pos, size=Prevalence_group, color=gen_color),
             shape=21, fill="white", stroke=1) +
  geom_point(data=cruces_m6_m14,
             aes(x=x_pos, y=y_pos),
             shape=4, size=2.4, stroke=0.8,
             color="grey35", inherit.aes=FALSE) +
  geom_point(data=cruces_ul56_m1_7,
             aes(x=x_pos, y=y_pos),
             shape=4, size=2.4, stroke=0.8,
             color="grey35", inherit.aes=FALSE) +
  
  scale_color_manual(
    values = c(colores_color,
               "Muestra secuenciada"    = "darkred",
               "Muestra no secuenciada" = "black"),
    breaks = c("UL97","UL54",
               "Muestra secuenciada","Muestra no secuenciada"),
    labels = c("UL97","UL54          ",
               "Muestra secuenciada","Muestra no secuenciada"),
    name   = NULL
  ) +
  scale_fill_manual(
    values = colores_fill,
    breaks = c("UL97","UL54"),
    name   = NULL, guide="none"
  ) +
  scale_shape_manual(
    values = c("Muestra secuenciada"=16, "Muestra no secuenciada"=16),
    guide  = "none"
  ) +
  scale_size_manual(values=tamaños, name="Umbral") +
  
  scale_x_continuous(
    limits = c(-3.5, n_niv + 0.5),
    breaks = seq_along(niveles_dia),
    labels = niveles_dia,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    limits = c(y_lim_inf, y_lim_sup),
    breaks = seq(0, 6, by=1),
    labels = NULL,
    expand = c(0,0)
  ) +
  
  guides(
    size = guide_legend(
      order=1, nrow=1,
      override.aes=list(shape=21, color=rep("grey50",3),
                        fill=rep("grey70",3), stroke=1)
    ),
    color = guide_legend(
      order=2, nrow=1,
      override.aes=list(
        shape  = c(21, 21, 16, 16),
        size   = c( 5,  5,  4,  4),
        alpha  = c(0.8,0.8, 1,  1),
        fill   = c(alpha("#4472C4",0.8), alpha("#E6A817",0.8),
                   "darkred", "black"),
        stroke = c(1, 1, 0, 0),
        color  = c("#4472C4","#C8880A","darkred","black")
      ),
      breaks=c("UL97","UL54",
               "Muestra secuenciada","Muestra no secuenciada")
    )
  ) +
  
  labs(x="D\u00edas postrasplante", y=NULL) +
  coord_cartesian(clip="off") +
  
  theme_minimal(base_size=13) +
  theme(
    axis.text.x       = element_text(angle=45, hjust=1, size=9),
    axis.ticks.x      = element_line(color="grey90", linewidth=0.3),
    axis.text.y       = element_blank(),
    axis.ticks.y      = element_blank(),
    panel.grid        = element_blank(),
    axis.title.x      = element_text(face="bold", margin=margin(t=12)),
    axis.title.y      = element_blank(),
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.box        = "horizontal",
    legend.title      = element_text(size=11, face="bold"),
    legend.text       = element_text(size=11),
    legend.spacing.x  = unit(0.25,"cm"),
    plot.margin       = margin(5, 15, 5, 10),
    plot.background   = element_rect(fill="white", color=NA),
    panel.background  = element_rect(fill="white", color=NA)
  )

grid::grid.newpage()
print(p)

# [guardado antiguo] png(filename="~/Desktop/Figuras_tesis/figura_RTP2_evolucion_VSI_10pct.png",
# [guardado antiguo]     width=35, height=21, units="cm", res=300, bg="white")
# [guardado antiguo] print(p)
# [guardado antiguo] dev.off()


# [guardado antiguo] ggsave("~/Desktop/Figuras_tesis/figura_RTP2_evolucion_VSI_10pct.png", p,
# [guardado antiguo]        width = 35, height = 21, units = "cm", dpi = 300, bg = "white")
# [guardado antiguo] ggsave("~/Desktop/Figuras_tesis/figura_RTP2_evolucion_VSI_10pct.svg", p,
# [guardado antiguo]        device = svglite::svglite, fix_text_size = FALSE,
# [guardado antiguo]        width = 38, height = 22, units = "cm", bg = "white")

