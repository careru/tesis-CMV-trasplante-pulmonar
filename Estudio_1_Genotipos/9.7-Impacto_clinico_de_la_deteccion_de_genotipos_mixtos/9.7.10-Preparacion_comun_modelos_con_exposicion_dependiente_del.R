# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.7.10. Preparación común: modelos con exposición dependiente del tiempo y curvas KM extendidas
# Extracto del script completo (ver README.md del estudio)

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

