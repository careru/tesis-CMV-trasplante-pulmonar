# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: ANEXO.4. Rejilla genotipo × episodio (artículo de diversidad; usa datos del Estudio 1)
# Extracto del script completo (ver README.md del estudio)

## ---- Rejilla genotipo × episodio (artículo de diversidad; usa datos del Estudio 1) ----
# Requiere el objeto dnaemias del script de genotipos; aquí no está disponible.
## ===ARTICULO_diveerisdad_ Rejilla genotipo × episodio (JID, B/N, presencia) con * en estables ====
library(dplyr); library(tidyr); library(ggplot2)
gt_up   <- c("GB1","GB2","GB3","GB4","GH1","GH2")
gt_nm   <- c("gB1","gB2","gB3","gB4","gH1","gH2")
orden_td<- c("gB1","gB2","gB3","gB4","gH1","gH2")   # arriba → abajo

# 1) Pacientes con >1 episodio
ep_por_paciente <- dnaemias %>% filter(!is.na(gbghmix)) %>%
  group_by(NHC) %>% summarise(n = n(), .groups = "drop") %>% filter(n > 1)

# 2) Episodios ordenados cronológicamente
ep <- dnaemias %>% filter(!is.na(gbghmix)) %>%
  mutate(Fini = as.Date(Fini)) %>%
  semi_join(ep_por_paciente, by = "NHC") %>%
  arrange(NHC, Fini) %>%
  group_by(NHC) %>% mutate(ep_idx = row_number(), n_ep = n()) %>% ungroup()

# 3) Pacientes que cambian de conjunto genotípico entre episodios
changers <- ep %>% group_by(NHC) %>%
  summarise(np = n_distinct(paste(gb1,gb2,gb3,gb4,gh1,gh2)), .groups = "drop") %>%
  filter(np > 1) %>% pull(NHC)

# 4) Orden de paneles por nº de episodios ASCENDENTE (primero los de 2)
ord <- ep %>% filter(NHC %in% changers) %>%
  distinct(NHC, n_ep) %>%
  left_join(pacientes %>% select(NHC, serDR), by = "NHC") %>%
  arrange(n_ep, NHC) %>%
  mutate(facet = paste0("P", row_number(), " (", serDR, ")"),
         facet = factor(facet, levels = facet))

# 5) Presencia de cada genotipo por episodio
ep_pres <- muestras_ep %>% filter(NHC %in% changers) %>%
  select(NHC, Fini, Freg, all_of(gt_up)) %>% mutate(Fini = as.Date(Fini)) %>%
  pivot_longer(all_of(gt_up), names_to = "g", values_to = "pres") %>%
  group_by(NHC, Fini, g) %>%
  summarise(present = as.integer(any(pres == 1, na.rm = TRUE)), .groups = "drop")

# 5b) Genotipos ESTABLES por paciente (positivos en TODOS los episodios)
stable_gt <- ep %>% filter(NHC %in% changers) %>%
  group_by(NHC) %>%
  summarise(across(c(gb1,gb2,gb3,gb4,gh1,gh2),
                   ~ as.integer(all(. == 1, na.rm = TRUE))), .groups = "drop") %>%
  pivot_longer(-NHC, names_to = "g_low", values_to = "stable") %>%
  mutate(g = toupper(g_low))

# 6) Rejilla completa + estado + marca de estable
grid <- ep %>% filter(NHC %in% changers) %>%
  distinct(NHC, Fini, ep_idx, n_ep) %>%
  tidyr::crossing(g = gt_up) %>%
  left_join(ep_pres, by = c("NHC","Fini","g")) %>%
  left_join(stable_gt %>% select(NHC, g, stable), by = c("NHC","g")) %>%
  left_join(ord %>% select(NHC, facet), by = "NHC") %>%
  mutate(gt = gt_nm[match(g, gt_up)],
         present = ifelse(is.na(present), 0L, present),
         status = case_when(present == 0      ~ "absent",
                            grepl("^gB", gt)  ~ "gB genotype",
                            TRUE              ~ "gH genotype"),
         status = factor(status, levels = c("gB genotype","gH genotype","absent")),
         star = ifelse(!is.na(stable) & stable == 1 & present == 1, "*", ""))

# 7) Figura
fig_grid <- ggplot(grid, aes(x = factor(ep_idx), y = gt, fill = status)) +
  geom_tile(color = "grey75", width = 0.95, height = 0.95, linewidth = 0.3) +
  geom_text(aes(label = star), color = "white", size = 3.4,
            fontface = "bold", vjust = 0.72) +
  scale_fill_manual(values = c("gB genotype" = "grey15",
                               "gH genotype" = "grey55",
                               "absent"      = "white"),
                    breaks = c("gB genotype","gH genotype"), name = NULL) +
  scale_y_discrete(limits = rev(orden_td)) +          # gB1 arriba … gH2 abajo
  facet_wrap(~ facet, scales = "free_x", ncol = 4) +
  labs(x = "Genotyped DNAemia episode (chronological order)", y = NULL) +
  theme_minimal(base_size = 9, base_family = "sans") +
  theme(legend.position = "none",
        legend.text   = element_text(size = 8),
        panel.grid    = element_blank(),
        strip.text    = element_text(face = "bold", size = 7.5),
        axis.text     = element_text(color = "black", size = 7),
        axis.title.x  = element_text(margin = margin(t = 9)),
        axis.text.x   = element_text(color = "black", size = 7, margin = margin(t = 4)),
        axis.ticks    = element_blank(),
        panel.spacing = unit(0.5, "lines"))

fig_grid
# [guardado antiguo] ggsave("Figure_genotype_tracks_BW.tiff", fig_grid,
# [guardado antiguo]        width = 8.5, height = 9.3, dpi = 300, compression = "lzw")

}  # fin del anexo
