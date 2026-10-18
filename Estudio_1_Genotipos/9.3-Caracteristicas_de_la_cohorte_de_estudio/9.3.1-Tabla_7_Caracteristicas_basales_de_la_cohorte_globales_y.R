# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.3.1. Tabla 7. Características basales de la cohorte globales y por perfiles de riesgo
# Extracto del script completo (ver README.md del estudio)

## ---- Tabla 7. Características basales de la cohorte globales y por perfiles de riesgo ----
pacientes <- read_excel(f_bd,
                        sheet = "pacientes_sel")

pacientes <- pacientes %>%
  mutate(
    edad_receptor           = as.numeric(difftime(Ftx, Fnac, units = "days")) / 365.25,
    # Seguimiento hasta Ffinseg (fallecimiento, traslado o 31/12/2022), como en Métodos
    tiempo_seguimiento_anos = as.numeric(difftime(as.Date(Ffinseg), as.Date(Ftx),
                                                  units = "days")) / 365.25,
    serDR = factor(serDR, levels = c("D-/R+", "D+/R+", "D+/R-")),
    tac   = grepl("TACRO", IS_base),
    cicl  = grepl("CICLO", IS_base),
    mmf   = grepl("MICO",  IS_base),
    aza   = grepl("AZA",   IS_base)
  )

gl <- pacientes
dm <- pacientes %>% filter(serDR == "D-/R+")
pp <- pacientes %>% filter(serDR == "D+/R+")
pm <- pacientes %>% filter(serDR == "D+/R-")

# ── Helpers ────────────────────────────────────────────────────────────────────
med <- function(x) sprintf("%.1f ± %.1f", mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))
mdn <- function(x) sprintf("%.1f (%.1f–%.1f)",
                           median(x, na.rm=TRUE),
                           quantile(x, 0.25, na.rm=TRUE),
                           quantile(x, 0.75, na.rm=TRUE))
pct <- function(x, tot=length(x)) sprintf("%d (%.1f%%)",
                                          sum(x, na.rm=TRUE),
                                          100*sum(x, na.rm=TRUE)/tot)
fp  <- function(p) if(is.na(p)) "-" else if(p<0.001) "<0,001" else gsub("\\.", ",", sprintf("%.3f", p))

eb  <- function(g, cat) {
  n <- sum(g$ebase == cat, na.rm=TRUE)
  sprintf("%d (%.1f%%)", n, 100*n/nrow(g))
}

# ── p-valores ──────────────────────────────────────────────────────────────────
p_edad_rx  <- fp(summary(aov(edad_receptor           ~ serDR, data=pacientes))[[1]][["Pr(>F)"]][1])
p_edad_don <- fp(summary(aov(edad_don                ~ serDR, data=pacientes))[[1]][["Pr(>F)"]][1])
p_sex_rx   <- fp(chisq.test(table(pacientes$SexF,     pacientes$serDR))$p.value)
p_sex_don  <- fp(chisq.test(table(pacientes$SexF_don, pacientes$serDR))$p.value)
p_tipotx   <- fp(chisq.test(table(pacientes$tipotx,   pacientes$serDR))$p.value)
p_ebase    <- fp(fisher.test(table(pacientes$ebase,   pacientes$serDR),
                             simulate.p.value=TRUE, B=10000)$p.value)
p_seg      <- fp(kruskal.test(tiempo_seguimiento_anos ~ serDR, data=pacientes)$p.value)

p_tac  <- tryCatch(fp(chisq.test(table(pacientes$tac,  pacientes$serDR), simulate.p.value=TRUE)$p.value), error=function(e) "-")
p_cicl <- tryCatch(fp(fisher.test(table(pacientes$cicl, pacientes$serDR))$p.value), error=function(e) "-")
p_mmf  <- tryCatch(fp(chisq.test(table(pacientes$mmf,  pacientes$serDR), simulate.p.value=TRUE)$p.value), error=function(e) "-")
p_aza  <- tryCatch(fp(fisher.test(table(pacientes$aza,  pacientes$serDR))$p.value), error=function(e) "-")

# ── Tabla ──────────────────────────────────────────────────────────────────────
df <- data.frame(
  Caracteristica = c(
    "Edad del receptor (años), media ± DE",
    "Edad del donante (años), media ± DE",
    "Género masculino del receptor, n (%)",
    "Género masculino del donante, n (%)",
    "Tipo de trasplante, n (%)",
    "   Unilateral",
    "   Bilateral",
    "Enfermedad basal, n (%)",
    "   EPOC / enfisema",
    "   EPID",
    "   Fibrosis quística",
    "   Hipertensión pulmonar",
    "   Otros",
    "Inmunosupresión inicial, n (%)",
    "   Esteroides",
    "   Tacrolimus",
    "   Ciclosporina",
    "   Micofenolato mofetilo / ácido micofenólico",
    "   Azatioprina",
    "Tiempo de seguimiento (años), mediana (RIC)"
  ),
  Global = c(
    med(gl$edad_receptor), med(gl$edad_don),
    pct(gl$SexF==0), pct(gl$SexF_don==0),
    "", pct(gl$tipotx=="UNI"), pct(gl$tipotx=="BI"),
    "", eb(gl,"EPOC"), eb(gl,"EPID"), eb(gl,"FQ"), eb(gl,"HAP"), eb(gl,"OTROS"),
    "", pct(rep(TRUE,nrow(gl))), pct(gl$tac), pct(gl$cicl), pct(gl$mmf), pct(gl$aza),
    mdn(gl$tiempo_seguimiento_anos)
  ),
  DmRp = c(
    med(dm$edad_receptor), med(dm$edad_don),
    pct(dm$SexF==0, nrow(dm)), pct(dm$SexF_don==0, nrow(dm)),
    "", pct(dm$tipotx=="UNI", nrow(dm)), pct(dm$tipotx=="BI", nrow(dm)),
    "", eb(dm,"EPOC"), eb(dm,"EPID"), eb(dm,"FQ"), eb(dm,"HAP"), eb(dm,"OTROS"),
    "", pct(rep(TRUE,nrow(dm)),nrow(dm)), pct(dm$tac,nrow(dm)),
    pct(dm$cicl,nrow(dm)), pct(dm$mmf,nrow(dm)), pct(dm$aza,nrow(dm)),
    mdn(dm$tiempo_seguimiento_anos)
  ),
  DpRp = c(
    med(pp$edad_receptor), med(pp$edad_don),
    pct(pp$SexF==0, nrow(pp)), pct(pp$SexF_don==0, nrow(pp)),
    "", pct(pp$tipotx=="UNI", nrow(pp)), pct(pp$tipotx=="BI", nrow(pp)),
    "", eb(pp,"EPOC"), eb(pp,"EPID"), eb(pp,"FQ"), eb(pp,"HAP"), eb(pp,"OTROS"),
    "", pct(rep(TRUE,nrow(pp)),nrow(pp)), pct(pp$tac,nrow(pp)),
    pct(pp$cicl,nrow(pp)), pct(pp$mmf,nrow(pp)), pct(pp$aza,nrow(pp)),
    mdn(pp$tiempo_seguimiento_anos)
  ),
  DpRm = c(
    med(pm$edad_receptor), med(pm$edad_don),
    pct(pm$SexF==0, nrow(pm)), pct(pm$SexF_don==0, nrow(pm)),
    "", pct(pm$tipotx=="UNI", nrow(pm)), pct(pm$tipotx=="BI", nrow(pm)),
    "", eb(pm,"EPOC"), eb(pm,"EPID"), eb(pm,"FQ"), eb(pm,"HAP"), eb(pm,"OTROS"),
    "", pct(rep(TRUE,nrow(pm)),nrow(pm)), pct(pm$tac,nrow(pm)),
    pct(pm$cicl,nrow(pm)), pct(pm$mmf,nrow(pm)), pct(pm$aza,nrow(pm)),
    mdn(pm$tiempo_seguimiento_anos)
  ),
  p = c(
    p_edad_rx, p_edad_don, p_sex_rx, p_sex_don,
    p_tipotx, "", "",
    p_ebase, "", "", "", "", "",
    "", "-", p_tac, p_cicl, p_mmf, p_aza,
    p_seg
  ),
  stringsAsFactors = FALSE
)

colnames(df) <- c("Característica", "Global\n(n=146)", "D-/R+\n(n=15)",
                  "D+/R+\n(n=110)", "D+/R-\n(n=21)", "p")

# ── Flextable ──────────────────────────────────────────────────────────────────
apartados <- c("Tipo de trasplante, n (%)", "Enfermedad basal, n (%)",
               "Inmunosupresión inicial, n (%)")

ft_tabla_7 <- flextable(df) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  bold(i = ~ Característica %in% apartados, part = "body") %>%
  italic(i = ~ grepl("^   ", Característica), j = "Característica", part = "body") %>%
  align(j = "Característica", align = "left",   part = "all") %>%
  align(j = -1,               align = "center", part = "all") %>%
  fontsize(part = "all", size = 10) %>%
  padding(part = "all", padding = 3) %>%
  add_footer_lines("*p < 0,05. EPID: enfermedad pulmonar intersticial difusa; EPOC: enfermedad pulmonar obstructiva crónica. Variables categóricas comparadas mediante test exacto de Fisher con simulación de Monte Carlo (B=10.000). Cada fármaco inmunosupresor se analizó como variable dicotómica independiente.") %>%
  fontsize(part = "footer", size = 9) %>%
  autofit()

ft_tabla_7

guardar_tabla(ft_tabla_7, "Tabla_07")


