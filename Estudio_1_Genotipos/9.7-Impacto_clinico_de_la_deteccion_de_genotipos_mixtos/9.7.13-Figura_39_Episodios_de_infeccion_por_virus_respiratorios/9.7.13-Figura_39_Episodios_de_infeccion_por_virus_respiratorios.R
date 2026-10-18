# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.7.13. Figura 39. Episodios de infección por virus respiratorios comunitarios por año
# Extracto del script completo (ver README.md del estudio)

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

