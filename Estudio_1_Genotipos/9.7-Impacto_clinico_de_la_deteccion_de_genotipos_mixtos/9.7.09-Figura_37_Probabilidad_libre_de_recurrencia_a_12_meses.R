# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.7.09. Figura 37. Probabilidad libre de recurrencia a 12 meses según el patrón genotípico
# Extracto del script completo (ver README.md del estudio)

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

