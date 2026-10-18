# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.7.15. Figura 40. Genotipados de CMV respecto al primer episodio de rechazo agudo
# Extracto del script completo (ver README.md del estudio)

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

