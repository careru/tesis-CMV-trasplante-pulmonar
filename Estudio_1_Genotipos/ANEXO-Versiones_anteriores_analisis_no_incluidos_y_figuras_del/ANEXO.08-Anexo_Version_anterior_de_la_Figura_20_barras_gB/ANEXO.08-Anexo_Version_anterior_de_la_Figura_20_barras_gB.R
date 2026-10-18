# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.08. Anexo: Versión anterior de la Figura 20 (barras gB)
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Versión anterior de la Figura 20 (barras gB) ----
# --- 1) Formatear nombres de combinaciones gB ---
clean_combo <- function(x) {
  x <- gsub("GB", "gB", x)
  x <- gsub("\\+", " + ", x)      # espacio alrededor del +
  x <- gsub("\\s+", " ", x)       # limpiar dobles espacios
  trimws(x)
}

desglose_1$combo <- clean_combo(desglose_1$combo)
desglose_2$combo <- clean_combo(desglose_2$combo)
desglose_3$combo <- clean_combo(desglose_3$combo)

# --- 2) Añadir grupo correspondiente ---
desglose_1$grupo <- "Genotipos gB únicos"
desglose_2$grupo <- "Genotipos gB mixtos"
desglose_3$grupo <- "Genotipos gB mixtos"

# --- 3) Combinar ---
desglose_all <- bind_rows(desglose_1, desglose_2, desglose_3)

# --- 4) Crear una variable auxiliar para ordenar combinaciones ---
desglose_all <- desglose_all %>%
  mutate(
    elementos = strsplit(combo, " \\+ "),
    n_elem    = lengths(elementos),            # 1, 2 o 3
    combo_ord = sapply(elementos, function(v) paste(sort(v), collapse = " + "))
  )

# --- 5) Orden deseado ---
desglose_all <- desglose_all %>%
  arrange(
    n_elem,           # primero 1 genotipo, luego 2, luego 3
    combo_ord         # orden alfabético dentro del grupo
  )

# Reconvertir combo a factor con el nuevo orden
desglose_all$combo <- factor(desglose_all$combo, levels = unique(desglose_all$combo))

# --- 6) Calcular porcentajes dentro del grupo y crear etiquetas ---
desglose_all <- desglose_all %>%
  group_by(grupo) %>%
  mutate(
    pct_grupo = round(100 * n / sum(n), 1),
    pct_label = if_else(
      pct_grupo %% 1 == 0,
      as.character(as.integer(pct_grupo)),
      formatC(pct_grupo, digits = 1, format = "f")
    ),
    etiqueta = sprintf("%d\n(%s%%)", n, pct_label),
    total_grupo = sum(n)
  ) %>%
  ungroup() %>%
  mutate(
    strip_lab = case_when(
      grupo == "Genotipos gB únicos" ~ sprintf("Genotipos gB únicos\n(n=%d)",
                                               unique(total_grupo[grupo == "Genotipos gB únicos"])),
      grupo == "Genotipos gB mixtos" ~ sprintf("Genotipos gB mixtos\n(n=%d)",
                                               unique(total_grupo[grupo == "Genotipos gB mixtos"]))
    ),
    strip_lab = factor(strip_lab, levels = unique(strip_lab))
  )

# --- 7) Plot final ---
distribuciones_gB <- ggplot(desglose_all, aes(x = combo, y = n, fill = grupo)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = etiqueta),
            vjust = -0.4, size = 3.3, lineheight = 0.9) +
  facet_grid(~ strip_lab, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("Genotipos gB únicos" = "#4F94CD",
                               "Genotipos gB mixtos" = "tan2")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(x = "", y = "Número de muestras") +
  theme_classic(base_size = 13) +
  theme(
    plot.margin      = margin(t = 18, r = 15, b = 10, l = 10),
    legend.position  = "none",
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.line        = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks       = element_line(linewidth = 0.3, colour = "grey50"),
    axis.title.y     = element_text(size = 11),
    strip.text       = element_text(face = "bold"),
    panel.spacing.x  = unit(2.5, "lines")
  )

distribuciones_gB

