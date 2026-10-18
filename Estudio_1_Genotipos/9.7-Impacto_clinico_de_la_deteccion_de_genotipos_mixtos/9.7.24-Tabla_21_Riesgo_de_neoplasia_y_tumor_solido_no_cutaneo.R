# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.7.24. Tabla 21. Riesgo de neoplasia y tumor sólido no cutáneo (tabla resumen de todos los desenlaces)
# Extracto del script completo (ver README.md del estudio)

## ---- Tabla 21. Riesgo de neoplasia y tumor sólido no cutáneo (tabla resumen de todos los desenlaces) ----
# Va al final porque reúne los modelos de todos los desenlaces (rechazo, DCIP,
# neoplasia, tumor sólido, infecciones oportunistas y mortalidad).
tabla_principal <- bind_rows(
  lapply(list(RA = res_RA, DCIP = res_DCIP, onco = res_onco, solido = res_solido,
              inf_opo = res_inf_opo, exitus = res_exitus), function(r)
                bind_rows(lapply(names(r), function(k)
                  cbind(Clasificacion = k, r[[k]]$res[grepl("PRINCIPAL", r[[k]]$res$Modelo), ])))),
  .id = "Desenlace") %>%
  select(Desenlace, Clasificacion, n, eventos, eventos_expuesto, HR, IC95, p)

etiq_des <- c(RA      = "Rechazo agudo",
              DCIP    = "Disfunción crónica del injerto pulmonar",
              onco    = "Neoplasia global",
              solido  = "Tumor sólido no cutáneo",
              inf_opo = "Infección oportunista no CMV",
              exitus  = "Mortalidad")
etiq_cla <- c(gB = "gB", gH = "gH", gBgH = "gB/gH")

tab_ft <- tabla_principal %>%
  mutate(
    Desenlace     = etiq_des[Desenlace],
    Clasificacion = etiq_cla[Clasificacion],
    `HR ajustada (IC 95 %)` = paste0(HR, " (", IC95, ")"),
    p = ifelse(grepl("^<", p), p, sub("^0,", "0,", p))
  ) %>%
  select(Desenlace, Clasificacion, n, eventos, eventos_expuesto,
         `HR ajustada (IC 95 %)`, p)

ft_principal <- flextable(tab_ft) %>%
  set_header_labels(
    Desenlace        = "Desenlace",
    Clasificacion    = "Clasificación",
    n                = "n",
    eventos          = "Eventos",
    eventos_expuesto = "Eventos en expuestos",
    p                = "p"
  ) %>%
  merge_v(j = "Desenlace") %>%
  valign(j = "Desenlace", valign = "top") %>%
  align(j = c("n", "eventos", "eventos_expuesto", "p"), align = "center", part = "all") %>%
  align(j = "HR ajustada (IC 95 %)", align = "center", part = "all") %>%
  bold(part = "header") %>%
  bold(i = ~ grepl("^<", p) | suppressWarnings(as.numeric(sub(",", ".", p))) < 0.05,
       j = c("HR ajustada (IC 95 %)", "p")) %>%
  flextable::font(fontname = "Source Sans Pro", part = "all") %>%
  border_remove() %>%
  hline_top(part = "header", border = fp_border(width = 1)) %>%
  hline_bottom(part = "header", border = fp_border(width = 0.5)) %>%
  hline_bottom(part = "body", border = fp_border(width = 1)) %>%
  hline(i = ~ !duplicated(Desenlace, fromLast = TRUE),
        border = fp_border(width = 0.25, color = "grey60")) %>%
  add_footer_lines(paste0(
    "Modelo de Cox con entrada diferida en el primer episodio con genotipado válido ",
    "del gen estudiado y exposición dependiente del tiempo (paso a mixto en la primera ",
    "detección de genotipos mixtos). Ajustado por edad al trasplante y perfil serológico ",
    "D/R (mortalidad: serología D/R en tres categorías). n: pacientes libres del desenlace ",
    "en la entrada; HR: hazard ratio; IC: intervalo de confianza. En negrita, p < 0,05.")) %>%
  fontsize(size = 8, part = "footer") %>%
  italic(part = "footer") %>%
  autofit()

ft_principal

## Exportar a Word (opcional):
# save_as_docx(ft_principal, path = "~/Desktop/Tabla_HR_principales.docx")

guardar_tabla(ft_principal, "Tabla_21")


