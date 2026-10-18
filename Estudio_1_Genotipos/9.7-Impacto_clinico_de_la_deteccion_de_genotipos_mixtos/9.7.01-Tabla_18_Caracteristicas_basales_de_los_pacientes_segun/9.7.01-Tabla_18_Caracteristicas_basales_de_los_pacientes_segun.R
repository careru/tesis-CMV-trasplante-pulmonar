# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.7.01. Tabla 18. Características basales de los pacientes según genotipos mixtos gB/gH
# Extracto del script completo (ver README.md del estudio)

## ---- Tabla 18. Características basales de los pacientes según genotipos mixtos gB/gH ----
# Este bloque vuelve a leer pacientes_sel: es la base de pacientes del resto de 9.7.
pacientes <- read_excel(f_bd,
                        sheet = "pacientes_sel")

# ── 1. Preparar datos ─────────────────────────────────────────────────────────
pacientes_tab <- pacientes %>%
  mutate(
    Ftx         = as.Date(Ftx),
    Ffinseg     = as.Date(Ffinseg),
    t_seg_anios = as.numeric(Ffinseg - Ftx) / 365.25,
    grupo       = factor(gbghm, levels = c(0,1),
                         labels = c("Genotipos únicos", "Genotipos mixtos"))
  )

g0  <- pacientes_tab %>% filter(grupo == "Genotipos únicos")
g1  <- pacientes_tab %>% filter(grupo == "Genotipos mixtos")
tot <- pacientes_tab
n0  <- nrow(g0); n1 <- nrow(g1); nt <- nrow(tot)

# ── 2. Funciones auxiliares ───────────────────────────────────────────────────
med_ric <- function(x) {
  x <- x[!is.na(x)]
  sprintf("%.1f (%.1f\u2013%.1f)",
          median(x), quantile(x, 0.25), quantile(x, 0.75))
}

n_pct <- function(x, n_total) {
  sprintf("%d (%.1f%%)", sum(x, na.rm = TRUE),
          100 * sum(x, na.rm = TRUE) / n_total)
}

fmt_p <- function(p) {
  ifelse(p < 0.001, "<0,001",
         sub("\\.", ",", sprintf("%.3f", p)))
}

p_cont <- function(x, g) {
  fmt_p(wilcox.test(x ~ g, exact = FALSE)$p.value)
}

p_cat <- function(x, g) {
  tab <- table(x, g)
  test <- tryCatch(
    fisher.test(tab, simulate.p.value = TRUE, B = 10000),
    error = function(e) chisq.test(tab)
  )
  fmt_p(test$p.value)
}

# ── 3. Construir filas ────────────────────────────────────────────────────────
filas <- bind_rows(
  
  tibble(variable = "n",
         total = as.character(nt),
         unicos = as.character(n0),
         mixtos = as.character(n1),
         p = ""),
  
  # ── Demográficas ──────────────────────────────────────────────────────────
  tibble(variable = "CARACTERÍSTICAS DEMOGRÁFICAS",
         total = "", unicos = "", mixtos = "", p = ""),
  
  tibble(variable = "Edad receptor (años), mediana (RIC)",
         total  = med_ric(tot$edad_tx),
         unicos = med_ric(g0$edad_tx),
         mixtos = med_ric(g1$edad_tx),
         p = p_cont(pacientes_tab$edad_tx, pacientes_tab$grupo)),
  
  tibble(variable = "Edad donante (años), mediana (RIC)",
         total  = med_ric(tot$edad_don),
         unicos = med_ric(g0$edad_don),
         mixtos = med_ric(g1$edad_don),
         p = p_cont(pacientes_tab$edad_don, pacientes_tab$grupo)),
  
  tibble(variable = "Sexo femenino receptor, n (%)",
         total  = n_pct(tot$SexF == 1, nt),
         unicos = n_pct(g0$SexF  == 1, n0),
         mixtos = n_pct(g1$SexF  == 1, n1),
         p = p_cat(pacientes_tab$SexF, pacientes_tab$grupo)),
  
  tibble(variable = "Sexo femenino donante, n (%)",
         total  = n_pct(tot$SexF_don == 1, nt),
         unicos = n_pct(g0$SexF_don  == 1, n0),
         mixtos = n_pct(g1$SexF_don  == 1, n1),
         p = p_cat(pacientes_tab$SexF_don, pacientes_tab$grupo)),
  
  # ── Seguimiento ──────────────────────────────────────────────────────────
  tibble(variable = "CARACTERÍSTICAS DE SEGUIMIENTO",
         total = "", unicos = "", mixtos = "", p = ""),
  
  tibble(variable = "Tiempo de seguimiento (años), mediana (RIC)",
         total  = med_ric(tot$t_seg_anios),
         unicos = med_ric(g0$t_seg_anios),
         mixtos = med_ric(g1$t_seg_anios),
         p = p_cont(pacientes_tab$t_seg_anios, pacientes_tab$grupo)),
  
  tibble(variable = "Episodios CMV totales, mediana (RIC)",
         total  = med_ric(tot$nepi_tot),
         unicos = med_ric(g0$nepi_tot),
         mixtos = med_ric(g1$nepi_tot),
         p = p_cont(pacientes_tab$nepi_tot, pacientes_tab$grupo)),
  
  tibble(variable = "Episodios CMV genotipados, mediana (RIC)",
         total  = med_ric(tot$nepi_gtc),
         unicos = med_ric(g0$nepi_gtc),
         mixtos = med_ric(g1$nepi_gtc),
         p = p_cont(pacientes_tab$nepi_gtc, pacientes_tab$grupo)),
  
  tibble(variable = "Determinaciones CV totales, mediana (RIC)",
         total  = med_ric(tot$Npcrs),
         unicos = med_ric(g0$Npcrs),
         mixtos = med_ric(g1$Npcrs),
         p = p_cont(pacientes_tab$Npcrs, pacientes_tab$grupo)),
  
  tibble(variable = "Antigenemias totales, mediana (RIC)",
         total  = med_ric(tot$Nag),
         unicos = med_ric(g0$Nag),
         mixtos = med_ric(g1$Nag),
         p = p_cont(pacientes_tab$Nag, pacientes_tab$grupo)),
  
  # ── Tipo de trasplante ────────────────────────────────────────────────────
  tibble(variable = "Tipo de trasplante, n (%)",
         total = "", unicos = "", mixtos = "",
         p = p_cat(pacientes_tab$tipotx, pacientes_tab$grupo)),
  
  tibble(variable = "   Unipulmonar",
         total  = n_pct(tot$tipotx == "UNI", nt),
         unicos = n_pct(g0$tipotx  == "UNI", n0),
         mixtos = n_pct(g1$tipotx  == "UNI", n1), p = ""),
  
  tibble(variable = "   Bipulmonar",
         total  = n_pct(tot$tipotx == "BI", nt),
         unicos = n_pct(g0$tipotx  == "BI", n0),
         mixtos = n_pct(g1$tipotx  == "BI", n1), p = ""),
  
  # ── Enfermedad de base ────────────────────────────────────────────────────
  tibble(variable = "Enfermedad de base, n (%)",
         total = "", unicos = "", mixtos = "",
         p = p_cat(pacientes_tab$ebase, pacientes_tab$grupo)),
  
  tibble(variable = "   EPID",
         total  = n_pct(tot$ebase == "EPID", nt),
         unicos = n_pct(g0$ebase  == "EPID", n0),
         mixtos = n_pct(g1$ebase  == "EPID", n1), p = ""),
  
  tibble(variable = "   EPOC",
         total  = n_pct(tot$ebase == "EPOC", nt),
         unicos = n_pct(g0$ebase  == "EPOC", n0),
         mixtos = n_pct(g1$ebase  == "EPOC", n1), p = ""),
  
  tibble(variable = "   FPI",
         total  = n_pct(tot$ebase == "FPI", nt),
         unicos = n_pct(g0$ebase  == "FPI", n0),
         mixtos = n_pct(g1$ebase  == "FPI", n1), p = ""),
  
  tibble(variable = "   FQ",
         total  = n_pct(tot$ebase == "FQ", nt),
         unicos = n_pct(g0$ebase  == "FQ", n0),
         mixtos = n_pct(g1$ebase  == "FQ", n1), p = ""),
  
  tibble(variable = "   HAP",
         total  = n_pct(tot$ebase == "HAP", nt),
         unicos = n_pct(g0$ebase  == "HAP", n0),
         mixtos = n_pct(g1$ebase  == "HAP", n1), p = ""),
  
  tibble(variable = "   Otros",
         total  = n_pct(!tot$ebase %in% c("EPID","EPOC","FPI","FQ","HAP"), nt),
         unicos = n_pct(!g0$ebase  %in% c("EPID","EPOC","FPI","FQ","HAP"), n0),
         mixtos = n_pct(!g1$ebase  %in% c("EPID","EPOC","FPI","FQ","HAP"), n1),
         p = ""),
  
  # ── Perfil serológico ─────────────────────────────────────────────────────
  tibble(variable = "Perfil serológico D/R, n (%)",
         total = "", unicos = "", mixtos = "",
         p = p_cat(pacientes_tab$serDR, pacientes_tab$grupo)),
  
  tibble(variable = "   D-/R+",
         total  = n_pct(tot$serDR == "D-/R+", nt),
         unicos = n_pct(g0$serDR  == "D-/R+", n0),
         mixtos = n_pct(g1$serDR  == "D-/R+", n1), p = ""),
  
  tibble(variable = "   D+/R+",
         total  = n_pct(tot$serDR == "D+/R+", nt),
         unicos = n_pct(g0$serDR  == "D+/R+", n0),
         mixtos = n_pct(g1$serDR  == "D+/R+", n1), p = ""),
  
  tibble(variable = "   D+/R-",
         total  = n_pct(tot$serDR == "D+/R-", nt),
         unicos = n_pct(g0$serDR  == "D+/R-", n0),
         mixtos = n_pct(g1$serDR  == "D+/R-", n1), p = ""),
  
  # ── Inmunosupresión ───────────────────────────────────────────────────────
  tibble(variable = "INMUNOSUPRESIÓN BASAL, n (%)",
         total = "", unicos = "", mixtos = "", p = ""),
  
  tibble(variable = "   Tacrolimus",
         total  = n_pct(grepl("TACRO", tot$IS_base), nt),
         unicos = n_pct(grepl("TACRO", g0$IS_base),  n0),
         mixtos = n_pct(grepl("TACRO", g1$IS_base),  n1),
         p = p_cat(grepl("TACRO", pacientes_tab$IS_base),
                   pacientes_tab$grupo)),
  
  tibble(variable = "   Ciclosporina",
         total  = n_pct(grepl("CICLO", tot$IS_base), nt),
         unicos = n_pct(grepl("CICLO", g0$IS_base),  n0),
         mixtos = n_pct(grepl("CICLO", g1$IS_base),  n1),
         p = p_cat(grepl("CICLO", pacientes_tab$IS_base),
                   pacientes_tab$grupo)),
  
  tibble(variable = "   Micofenolato",
         total  = n_pct(grepl("MICO", tot$IS_base), nt),
         unicos = n_pct(grepl("MICO", g0$IS_base),  n0),
         mixtos = n_pct(grepl("MICO", g1$IS_base),  n1),
         p = p_cat(grepl("MICO", pacientes_tab$IS_base),
                   pacientes_tab$grupo)),
  
  tibble(variable = "   Azatioprina",
         total  = n_pct(grepl("AZA", tot$IS_base), nt),
         unicos = n_pct(grepl("AZA", g0$IS_base),  n0),
         mixtos = n_pct(grepl("AZA", g1$IS_base),  n1),
         p = p_cat(grepl("AZA", pacientes_tab$IS_base),
                   pacientes_tab$grupo)),
  
  tibble(variable = "   Esteroides",
         total  = "146 (100%)",
         unicos = "95 (100%)",
         mixtos = "51 (100%)",
         p = "\u2014")
)

# ── 4. Flextable ──────────────────────────────────────────────────────────────
filas_seccion <- which(filas$variable %in% c(
  "CARACTERÍSTICAS DEMOGRÁFICAS",
  "CARACTERÍSTICAS DE SEGUIMIENTO",
  "INMUNOSUPRESIÓN BASAL, n (%)"
))

filas_padre <- which(filas$variable %in% c(
  "Tipo de trasplante, n (%)",
  "Enfermedad de base, n (%)",
  "Perfil serológico D/R, n (%)"
))

colnames(filas) <- c("Variable",
                     "Total\n(n = 146)",
                     "Genotipos únicos\n(n = 95)",
                     "Genotipos mixtos\n(n = 51)",
                     "p")

# p significativas para resaltar
p_sig <- c("<0,001")

ft <- flextable(filas) %>%
  width(j = "Variable",                   width = 4.5) %>%
  width(j = "Total\n(n = 146)",           width = 2.0) %>%
  width(j = "Genotipos únicos\n(n = 95)", width = 2.0) %>%
  width(j = "Genotipos mixtos\n(n = 51)", width = 2.0) %>%
  width(j = "p",                          width = 0.9) %>%
  align(align = "left",   part = "all") %>%
  align(j = c(2,3,4,5), align = "center", part = "all") %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header") %>%
  bg(bg = "white", part = "all") %>%
  bold(i = filas_seccion, part = "body") %>%
  bg(i = filas_seccion, bg = "#F0F0F0", part = "body") %>%
  italic(i = filas_padre, part = "body") %>%
  color(i = ~ p %in% p_sig, j = "p",
        color = "firebrick3", part = "body") %>%
  bold(i  = ~ p %in% p_sig, j = "p", part = "body") %>%
  border_remove() %>%
  hline_top(part = "header",
            border = fp_border(color = "black", width = 1.2)) %>%
  hline(part = "header",
        border = fp_border(color = "black", width = 0.8)) %>%
  hline_bottom(part = "body",
               border = fp_border(color = "black", width = 1.2)) %>%
  hline(i = filas_seccion - 1,
        border = fp_border(color = "grey70", width = 0.4),
        part = "body") %>%
  padding(padding.top = 2, padding.bottom = 2,
          padding.left = 4, padding.right = 4, part = "all") %>%
  add_footer_lines(values = c(
    "RIC: rango intercuartílico; EPID: enfermedad pulmonar intersticial difusa; EPOC: enfermedad pulmonar obstructiva crónica; FPI: fibrosis pulmonar idiopática; FQ: fibrosis quística; HAP: hipertensión arterial pulmonar.",
    "Variables continuas: test de Wilcoxon. Variables categóricas: test exacto de Fisher o chi-cuadrado.",
    "Valores de p en rojo indican diferencia estadísticamente significativa (p < 0,05)."
  )) %>%
  fontsize(size = 8, part = "footer") %>%
  color(color = "grey30", part = "footer") %>%
  set_caption(
    caption = "Tabla R11. Características basales de los pacientes según presencia de genotipos mixtos gB/gH."
  ) %>%
  set_table_properties(layout = "autofit")

ft

guardar_tabla(ft, "Tabla_18")

