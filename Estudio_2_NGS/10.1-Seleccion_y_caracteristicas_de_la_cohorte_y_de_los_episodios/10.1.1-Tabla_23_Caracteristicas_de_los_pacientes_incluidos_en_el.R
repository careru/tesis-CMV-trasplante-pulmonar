# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 10.1.1. Tabla 23. Características de los pacientes incluidos en el análisis NGS (cifras del texto)
# Extracto del script completo (ver README.md del estudio)

## ---- Tabla 23. Características de los pacientes incluidos en el análisis NGS (cifras del texto) ----
# [Nuevo] La Tabla 23 se maquetó a mano en Word. Este bloque reproduce las cifras del
# texto de 10.1 a partir de CMV_NGS_pacientes_DEF.xlsx (hojas "TXP" y "txp_episodios")
# y compara el perfil serológico D/R entre las dos hojas.
rtp_tx <- read_excel(f_clinica, sheet = "TXP")
rtp_ep <- suppressWarnings(read_excel(f_clinica, sheet = "txp_episodios")) %>%
  filter(!is.na(Fnac))                     # una fila por paciente: el episodio analizado
stopifnot(nrow(rtp_tx) == 11, nrow(rtp_ep) == 11)

med_ric <- function(x, dec = 1) {
  q <- quantile(x, c(.5, .25, .75), na.rm = TRUE)
  sprintf("%s (RIC: %s–%s)", round(q[1], dec), round(q[2], dec), round(q[3], dec))
}

edad_tx <- round(rtp_tx$edad)              # la tabla de la tesis usa la edad redondeada
cat("Sexo: ", sum(rtp_tx$sexF == 0), "hombres y", sum(rtp_tx$sexF == 1), "mujeres\n")
cat("Edad al trasplante:", sprintf("%.1f ± %.1f años", mean(edad_tx), sd(edad_tx)), "\n")
cat("Tipo de trasplante:\n"); print(table(rtp_tx$tipotx))
cat("Meses desde el trasplante hasta el episodio:", med_ric(rtp_tx$tx_ep),
    sprintf("; rango %.1f–%.1f", min(rtp_tx$tx_ep), max(rtp_tx$tx_ep)), "\n")
cat("Carga viral pico (UI/mL):", med_ric(rtp_ep$cvpico_ev, 0), "\n")
cat("Síntomas:\n"); print(table(rtp_tx$clinica))
cat("Días hasta el aclaramiento (10 pacientes sin resistencia):",
    med_ric(rtp_ep$dur[rtp_ep$TXP != 2], 1), "\n")

# Perfil serológico D/R: las dos hojas no coinciden en RTP-3 y RTP-5, y ninguna
# coincide del todo con la Tabla 23 de la tesis (RTP-1 D−/R+, RTP-3 D+/R−).
# Conviene confirmarlo en la historia clínica antes de la defensa.
serologia <- tibble(ID = rtp_tx$ID, hoja_TXP = rtp_tx$serDR, hoja_episodios = rtp_ep$SeroD,
                    tesis_tabla23 = c("D-/R+", "D+/R-", "D+/R-", rep("D+/R+", 1), "D-/R+",
                                      rep("D+/R+", 6)))
print(serologia %>% filter(hoja_TXP != hoja_episodios | hoja_TXP != tesis_tabla23))


