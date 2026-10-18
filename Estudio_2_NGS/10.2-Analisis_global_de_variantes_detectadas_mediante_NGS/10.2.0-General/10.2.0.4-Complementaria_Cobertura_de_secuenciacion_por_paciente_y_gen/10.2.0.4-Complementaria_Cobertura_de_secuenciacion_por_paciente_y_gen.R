# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.2.0.4. Complementaria. Cobertura de secuenciación por paciente y gen
# Extracto del script completo (ver README.md del estudio)

## ---- Complementaria. Cobertura de secuenciación por paciente y gen ----
# ── Pacientes y número de muestras por paciente ───────────────────────────────
pacientes_cob <- tibble::tribble(
  ~RTP,      ~n,
  "RTP-1",    1,
  "RTP-2",   15,
  "RTP-3",    1,
  "RTP-4",    1,
  "RTP-5",    1,
  "RTP-6",    1,
  "RTP-7",    1,
  "RTP-8",    1,
  "RTP-9",    1,
  "RTP-10",   1,
  "RTP-11",   2
)
genes <- c("UL97", "UL54", "UL56")

# ── Estado por muestra (válido / cobertura insuficiente / no analizado) ───────
estado_muestra <- function(rtp, gen, n) {
  if (gen == "UL97") return(rep("valido", n))
  if (gen == "UL54") {
    v <- rep("valido", n)
    if (rtp == "RTP-2") v[c(6, 14)] <- "insuf"   # muestras 6 y 14 (días +361, +479)
    return(v)
  }
  # UL56: analizado solo en RTP-2 (8 de 15), RTP-7 (1) y RTP-11 (2)
  if (rtp == "RTP-2")  return(c(rep("valido", 8), rep("no", 7)))
  if (rtp == "RTP-7")  return(rep("valido", n))
  if (rtp == "RTP-11") return(rep("valido", n))
  rep("no", n)
}

# ── Expandir a una fila por muestra ───────────────────────────────────────────
datos <- pacientes_cob %>%
  tidyr::crossing(Gen = factor(genes, levels = genes)) %>%
  mutate(estado = pmap(list(RTP, as.character(Gen), n), estado_muestra)) %>%
  unnest(estado) %>%
  group_by(RTP, Gen) %>%
  mutate(idx = row_number()) %>%
  ungroup()

# ── Maquetación de la cuadrícula (5 casillas por fila dentro de cada celda) ───
ncol_cell <- 5
gap_pac   <- 0.7
gene_base <- c(UL97 = 0, UL54 = 7, UL56 = 14)   # x de inicio de cada gen

pac_layout <- pacientes_cob %>%
  mutate(
    RTP   = factor(RTP, levels = pacientes_cob$RTP),
    n_sub = ceiling(n / ncol_cell)
  ) %>%
  arrange(RTP) %>%
  mutate(
    y_top    = c(0, -head(cumsum(n_sub + gap_pac), -1)),
    y_centro = y_top - (n_sub - 1) / 2
  )

datos <- datos %>%
  left_join(pac_layout %>% select(RTP, y_top), by = "RTP") %>%
  mutate(
    col    = (idx - 1) %% ncol_cell,
    sub    = (idx - 1) %/% ncol_cell,
    x      = gene_base[as.character(Gen)] + col,
    y      = y_top - sub,
    estado = factor(estado, levels = c("valido", "insuf", "no"))
  )

# ── Etiquetas y resumen (resumen calculado desde los datos) ──────────────────
lab_genes <- tibble(Gen = names(gene_base),
                    x = gene_base + (ncol_cell - 1) / 2,
                    y = max(datos$y) + 1.4)

resumen <- datos %>%
  group_by(Gen) %>%
  summarise(validos = sum(estado == "valido"), total = n(), .groups = "drop") %>%
  mutate(
    pct     = validos / total * 100,
    pct_txt = ifelse(abs(pct - 100) < 1e-9, "100",
                     formatC(pct, format = "f", digits = 1, decimal.mark = ",")),
    label   = paste0(validos, "/", total, " \u00b7 ", pct_txt, "%"),
    x       = gene_base[as.character(Gen)] + (ncol_cell - 1) / 2,
    y       = min(datos$y) - 1.6
  )

# ── Colores por estado ────────────────────────────────────────────────────────
col_fill <- c(valido = "#CDEBC5", insuf = "#F6CDA0", no = "#E6E6E6")
col_line <- c(valido = "#6FB16A", insuf = "#DC8A3C", no = "#C2C2C2")

# ── Figura ────────────────────────────────────────────────────────────────────
p_cobertura <- ggplot(datos, aes(x, y)) +
  geom_tile(aes(fill = estado, color = estado),
            width = 0.82, height = 0.82, linewidth = 0.5) +
  geom_text(data = pac_layout,
            aes(x = -4.6, y = y_centro, label = RTP),
            hjust = 0, fontface = "bold", size = 3.6,
            color = "grey20", inherit.aes = FALSE) +
  geom_text(data = pac_layout,
            aes(x = -0.5, y = y_centro, label = paste0("n=", n)),
            hjust = 1, size = 3.0, color = "grey50", inherit.aes = FALSE) +
  geom_text(data = lab_genes,
            aes(x = x, y = y, label = Gen),
            fontface = "bold", size = 4.3, color = "grey20", inherit.aes = FALSE) +
  geom_text(data = resumen,
            aes(x = x, y = y, label = label),
            fontface = "bold", size = 3.4, color = "grey25", inherit.aes = FALSE) +
  scale_fill_manual(
    values = col_fill,
    labels = c("Resultado válido", "Cobertura insuficiente", "No analizado"),
    name   = NULL
  ) +
  scale_color_manual(values = col_line, guide = "none") +
  scale_x_continuous(limits = c(-5, 19), expand = c(0, 0)) +
  coord_equal(clip = "off") +
  theme_void(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size = 11),
    legend.key.size = unit(0.5, "cm"),
    plot.margin     = margin(10, 15, 10, 15)
  )

p_cobertura

# [guardado antiguo] ggsave(
# [guardado antiguo]   filename = "~/Desktop/Figuras_tesis/figura_cobertura_NGS_paciente_gen.png",
# [guardado antiguo]   plot     = p_cobertura,
# [guardado antiguo]   width    = 20, height = 22, units = "cm", dpi = 600, bg = "white"
# [guardado antiguo] )

guardar_fig(p_cobertura, "Complementaria_cobertura_paciente_gen", ancho = 20, alto = 22, dpi = 600)


