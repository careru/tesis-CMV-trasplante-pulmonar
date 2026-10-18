# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: ANEXO.33. Anexo: Artículo: donuts
# Extracto del script completo (ver README.md del estudio)

## ---- Anexo: Artículo: donuts ----
# Tres donuts (n = 234 episodios): A. gB  B. gH  C. gB/gH. Autocontenida.
library(dplyr); library(tidyr); library(ggplot2); library(patchwork)

COL_MIXED <- "#E8A07C"   # coral claro para categorias "Mixed ..."
TOP_N    <- 4            # combinaciones gB/gH unicas mostradas individualmente
LAB_SIZE <- 3.6          # mismo tamanio de etiqueta en los tres donuts

df_ep <- dnaemias %>%
  mutate(across(c(gb1:gb4, gh1:gh2), ~ replace_na(., 0))) %>%
  filter(gbsum >= 1, ghsum >= 1) %>%
  mutate(
    mixed  = gbsum > 1 | ghsum > 1,
    gb_cat = case_when(gbsum > 1 ~ "Mixed gB",
                       gb1 == 1 ~ "gB1", gb2 == 1 ~ "gB2",
                       gb3 == 1 ~ "gB3", gb4 == 1 ~ "gB4"),
    gh_cat = case_when(ghsum > 1 ~ "Mixed gH",
                       gh1 == 1 ~ "gH1", gh2 == 1 ~ "gH2"),
    gB_one = case_when(gb1 == 1 ~ "gB1", gb2 == 1 ~ "gB2",
                       gb3 == 1 ~ "gB3", gb4 == 1 ~ "gB4"),
    gH_one = case_when(gh1 == 1 ~ "gH1", gh2 == 1 ~ "gH2"),
    combo_one = paste0(gB_one, "/", gH_one)
  )
n_ep    <- nrow(df_ep)
fmt_pct <- function(x) formatC(x, format = "f", digits = 1)   # decimal con punto

tab_gb <- df_ep %>% count(gb_cat, name = "n") %>%
  mutate(cat = factor(gb_cat, levels = c("gB1","gB2","gB3","gB4","Mixed gB"))) %>%
  arrange(cat) %>% mutate(pct = 100 * n / sum(n))
tab_gh <- df_ep %>% count(gh_cat, name = "n") %>%
  mutate(cat = factor(gh_cat, levels = c("gH1","gH2","Mixed gH"))) %>%
  arrange(cat) %>% mutate(pct = 100 * n / sum(n))
top_combos <- df_ep %>% filter(!mixed) %>% count(combo_one, sort = TRUE) %>%
  slice_head(n = TOP_N) %>% pull(combo_one)
orden_combo <- c(top_combos, "Other single", "Mixed gB/gH")
tab_combo <- df_ep %>%
  mutate(combo_cat = case_when(
    mixed                     ~ "Mixed gB/gH",
    combo_one %in% top_combos ~ combo_one,
    TRUE                      ~ "Other single")) %>%
  count(combo_cat, name = "n") %>%
  mutate(cat = factor(combo_cat, levels = orden_combo)) %>%
  arrange(cat) %>% mutate(pct = 100 * n / sum(n))

cat(sprintf("n = %d episodes\n", n_ep))
cat("Top single combos (most frequent first):\n"); print(top_combos)
cat("\n-- gB --\n");    print(tab_gb    %>% select(cat, n, pct))
cat("\n-- gH --\n");    print(tab_gh    %>% select(cat, n, pct))
cat("\n-- combo --\n"); print(tab_combo %>% select(cat, n, pct))

col_gb <- c("gB1" = "#D6EAF8", "gB2" = "#AED6F1", "gB3" = "#85C1E9",
            "gB4" = "#5DADE2", "Mixed gB" = COL_MIXED)
col_gh <- c("gH1" = "#E8DAEF", "gH2" = "#BB8FCE", "Mixed gH" = COL_MIXED)
green_ramp <- c("#4CAF6E", "#76C68F", "#A1DBB3", "#C8EBD3", "#E2F4E8")
greens <- green_ramp[seq_len(length(top_combos))]; names(greens) <- top_combos
col_combo <- c(greens, "Other single" = "#AEC2B5", "Mixed gB/gH" = COL_MIXED)

donut <- function(tab, cols, centro, titulo) {
  tab <- tab %>% mutate(lab = paste0(cat, "\n", fmt_pct(pct), "%"))
  ggplot(tab, aes(x = 2, y = pct, fill = cat)) +
    geom_col(width = 1, color = "white", linewidth = 0.6) +
    geom_text(aes(label = lab), position = position_stack(vjust = 0.5),
              size = LAB_SIZE, color = "black", lineheight = 0.9) +
    coord_polar(theta = "y") + xlim(0.5, 2.5) +
    scale_fill_manual(values = cols, guide = "none") +
    annotate("text", x = 0.5, y = 0, label = centro, fontface = "bold", size = 4.5) +
    ggtitle(titulo) + theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.15, size = 14,
                                    margin = margin(t = 6, b = -20)),
          plot.title.position = "plot")
}

gA <- donut(tab_gb,    col_gb,    sprintf("gB\nn = %d", n_ep),    "A")
gB <- donut(tab_gh,    col_gh,    sprintf("gH\nn = %d", n_ep),    "B")
gC <- donut(tab_combo, col_combo, sprintf("gB/gH\nn = %d", n_ep), "C")
print(gA + gB + gC)

}  # fin del anexo


