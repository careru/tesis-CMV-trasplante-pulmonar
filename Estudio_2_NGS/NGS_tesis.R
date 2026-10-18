################################################################################
#
#   TESIS DOCTORAL · CMV EN RECEPTORES DE TRASPLANTE PULMONAR
#   Estudio 2: caracterización de variantes de UL97, UL54 y UL56 mediante NGS
#
#   Carmen Alhena Reyes Ruiz · Hospital Universitario 12 de Octubre · UAM
#
#   Organización: una sección por apartado de Resultados (10.1 y 10.2.1 a 10.2.6)
#   y, dentro de cada una, un bloque por Figura, Tabla o Anexo con la numeración
#   de la tesis. Todas las secciones se pliegan y despliegan (Alt+O / Alt+Shift+O).
#   El script corre de principio a fin; las salidas se guardan en 'dir_salida' con
#   su número de la tesis (Figura_45.png, Tabla_24.docx, Anexo_I.docx...).
#   Las tablas y figuras de trabajo que no aparecen numeradas en la tesis se
#   guardan con el prefijo "Complementaria_".
#
#   Reorganizado a partir de "Estudio2_NGS_CMV.R". El código de cada análisis es
#   el del script original; los cambios están marcados con [Nuevo] o [Corregido].
#
################################################################################


# ---- 0. Configuración ----

## ---- Rutas y archivos de datos ----
dir_datos  <- "~/Documents/CMV/ESTUDIO VARIANTES NGS"   # carpeta de trabajo
dir_aa     <- file.path(dir_datos, "NGS_CMV/AA")       # CSV de DeepChek (uno por muestra)
dir_salida <- file.path(dir_datos, "salidas_tesis")

f_muestras <- file.path(dir_datos, "allsamples.dta")             # muestras con NHC
f_clinica  <- file.path(dir_datos, "CMV_NGS_pacientes_DEF.xlsx")  # datos clínicos (Tabla 23)

ejecutar_anexo <- FALSE   # TRUE para ejecutar también los bloques del anexo

dir.create(dir_salida, showWarnings = FALSE)

## ---- Librerías ----
# Mismo orden de carga que el script original (data.table después de tidyverse).
library(tidyverse)
library(openxlsx)
library(readxl)
library(haven)
library(data.table)
library(flextable)
library(officer)
library(grid)
library(gridExtra)
library(patchwork)

## ---- Funciones para guardar figuras y tablas ----
guardar_fig <- function(p, nombre, ancho = 16, alto = 12, dpi = 400) {
  ggsave(file.path(dir_salida, paste0(nombre, ".png")), p,
         width = ancho, height = alto, units = "cm", dpi = dpi, bg = "white")
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(file.path(dir_salida, paste0(nombre, ".svg")), p,
           device = svglite::svglite, fix_text_size = FALSE,
           width = ancho, height = alto, units = "cm", bg = "white")
  }
  invisible(p)
}

guardar_tabla <- function(ft, nombre, apaisada = FALSE) {
  if (apaisada) {
    doc <- read_docx() %>% body_add_flextable(ft) %>% body_end_section_landscape()
    print(doc, target = file.path(dir_salida, paste0(nombre, ".docx")))
  } else {
    flextable::save_as_docx(ft, path = file.path(dir_salida, paste0(nombre, ".docx")))
  }
  invisible(ft)
}


# ---- 1. Preparación de los datos ----

## ---- Unión de los CSV de DeepChek (una fila por variante y muestra) ----
# ── Configuración ─────────────────────────────────────────────────────────────
path        <- dir_aa
output_file <- file.path(dir_salida, "1_out_mutaciones.xlsx")
extra_mra   <- c("L545F")

# ── Carga y fusión de archivos ────────────────────────────────────────────────
files <- list.files(path, pattern = "\\.csv$", full.names = TRUE)
cat("Archivos encontrados:", length(files), "\n")

all_data <- map(files, ~ {
  df <- read_delim(.x, delim = ";", col_types = cols(.default = "c")) %>%
    rename(Npet = 1) %>%
    mutate(Npet = trimws(Npet))
  
  names(df) <- str_replace(
    names(df),
    "Is mutation interest \\(CMV variant database .*\\)",
    "Resistencia"
  )
  names(df) <- str_replace(names(df), "Passed Expert System", "Filtro")
  
  return(df)
}) %>%
  bind_rows()

cat("Muestras únicas (pre-filtro):", n_distinct(all_data$Npet), "\n")

# ── Selección de columnas ─────────────────────────────────────────────────────
all_data <- all_data %>%
  dplyr::select(
    Npet, Protein, Position, Mutation, Prevalence, Resistencia, Filtro,
    matches("Mutation count FW|Mutation\\.count\\.FW"),
    matches("Mutation count RV|Mutation\\.count\\.RV"),
    matches("Position coverage FW|Position\\.coverage\\.FW"),
    matches("Position coverage RV|Position\\.coverage\\.RV"),
    matches("^Q-Score$|^Q\\.Score$")
  )

# ── Renombrar a nombres cortos ────────────────────────────────────────────────
old_names <- names(all_data)
names(all_data)[str_detect(old_names, "Mutation.*count.*FW|Mutation\\.count\\.FW")]       <- "mut_fw"
names(all_data)[str_detect(old_names, "Mutation.*count.*RV|Mutation\\.count\\.RV")]       <- "mut_rv"
names(all_data)[str_detect(old_names, "Position.*coverage.*FW|Position\\.coverage\\.FW")] <- "cov_fw"
names(all_data)[str_detect(old_names, "Position.*coverage.*RV|Position\\.coverage\\.RV")] <- "cov_rv"
names(all_data)[str_detect(old_names, "Q.Score|Q-Score")]                                 <- "q_score"

# ── Convertir a numérico y calcular métricas derivadas ───────────────────────
all_data <- all_data %>%
  mutate(
    across(c(mut_fw, mut_rv, cov_fw, cov_rv, q_score), as.numeric),
    cobertura_total = cov_fw + cov_rv,
    balance_ok      = mut_fw > 0 & mut_rv > 0
  )

cat("Columnas de calidad disponibles:",
    paste(c("mut_fw","mut_rv","cov_fw","cov_rv","q_score",
            "cobertura_total","balance_ok"), collapse = ", "), "\n")

# ── Filtros ───────────────────────────────────────────────────────────────────
n_dash <- sum(all_data$Mutation == "-", na.rm = TRUE)
cat("Filas con Mutation == '-' eliminadas:", n_dash, "\n")

all_data <- all_data %>%
  filter(Mutation != "-") %>%
  filter(Filtro == "Yes")

cat("Muestras únicas (post-filtro Expert System):", n_distinct(all_data$Npet), "\n")
cat("Variantes totales (post-filtro Expert System):", nrow(all_data), "\n")

# ── Corrección de Resistencia ─────────────────────────────────────────────────
all_data <- all_data %>%
  mutate(Resistencia = if_else(str_trim(Mutation) %in% extra_mra, "Yes", Resistencia))
cat("Variantes marcadas como MRA adicional:",
    sum(str_trim(all_data$Mutation) %in% extra_mra), "\n")

nas <- sum(is.na(all_data$Resistencia))
cat("NAs en Resistencia:", nas, "\n")
if (nas > 0) {
  all_data <- all_data %>% mutate(Resistencia = replace_na(Resistencia, "No"))
  cat("NAs corregidos a 'No'. Comprobación:", sum(is.na(all_data$Resistencia)), "\n")
}

# ── Exportar ──────────────────────────────────────────────────────────────────
write.xlsx(all_data, file = output_file)
cat("Archivo exportado:", output_file, "\n")

## ---- Fusión con los datos de las muestras (NHC) ----
# ── Exclusiones ───────────────────────────────────────────────────────────────
npet_excluidos <- unique(c(
  "102564799", "102344108", "102498654", "102500899", "213553115",
  "213645499", "334497189", "666063012", "102578305", "102473733",
  "102523258", "610121903", "610104837", "888328133", "102527866",
  "213797155", "666062877", "102775723", "102554490", "666075028",
  "111636931", "111631934", "213640305", "210744884", "210744746",
  "210744882", "210707909", "210734985", "210740947", "210740946",
  "666047983", "666050147", "666072527", "130384461", "210740916",
  "210743105", "120395862", "102856101", "102379827", "610099233",
  "102432523", "102344508", "610126260"
))
cat("Exclusiones únicas definidas:", length(npet_excluidos), "\n")

# ── Importar datos clínicos ───────────────────────────────────────────────────
pacientes <- read_dta(f_muestras)
setDT(pacientes)
pacientes <- pacientes %>% mutate(Npet = as.character(Npet))

# ── Aplicar exclusiones antes del join ───────────────────────────────────────
n_antes        <- n_distinct(all_data$Npet)
all_data_clean <- all_data %>% filter(!Npet %in% npet_excluidos)
n_despues      <- n_distinct(all_data_clean$Npet)
cat("Muestras antes de exclusiones:", n_antes,             "\n")
cat("Muestras excluidas:            ", n_antes - n_despues, "\n")
cat("Muestras tras exclusiones:     ", n_despues,           "\n")

# ── Fusión con datos clínicos ─────────────────────────────────────────────────
all_data_nhc <- all_data_clean %>%
  left_join(pacientes %>% dplyr::select(Npet, NHC, CV, Freg), by = "Npet") %>%
  filter(!is.na(NHC))

# ── Diagnóstico ───────────────────────────────────────────────────────────────
cat("\n── Pacientes (allsamples.dta) ──────────────────\n")
cat("Npet únicos:", n_distinct(pacientes$Npet), "\n")
cat("NHC únicos: ", n_distinct(pacientes$NHC),  "\n")
cat("\n── all_data_nhc (tras exclusiones + join) ───────\n")
cat("Npet únicos:", n_distinct(all_data_nhc$Npet), "\n")
cat("NHC únicos: ", n_distinct(all_data_nhc$NHC),  "\n")

# ── Comprobación: no deben quedar muestras sin match ─────────────────────────
sin_match_restantes <- setdiff(unique(all_data_clean$Npet), unique(pacientes$Npet))
if (length(sin_match_restantes) > 0) {
  cat("\n⚠️  Muestras aún sin match en allsamples.dta:",
      paste(sin_match_restantes, collapse = ", "), "\n")
} else {
  cat("\nTodas las muestras restantes tienen match en allsamples.dta.\n")
}

## ---- Filtro analítico: Expert System y frecuencia de variante ≥ 5% ----
all_data_nhc <- all_data_nhc %>%
  mutate(Prevalence = as.numeric(Prevalence)) %>%
  filter(Filtro == "Yes", Prevalence >= 5, cobertura_total >= 100)

cat("── Tras filtro Prevalence ≥5% ──────────────────\n")
cat("Variantes totales:", nrow(all_data_nhc),          "\n")
cat("Npet únicos:      ", n_distinct(all_data_nhc$Npet), "\n")
cat("NHC únicos:       ", n_distinct(all_data_nhc$NHC),  "\n")

cat("\nDistribución por gen:\n")
all_data_nhc %>%
  count(Protein, name = "n_variantes") %>%
  arrange(Protein) %>%
  print()

## ---- Selección de los receptores de trasplante pulmonar (RTP) ----
pacientes_txp <- c(
  "NHC_RTP-2", "NHC_RTP-11", "NHC_RTP-1", "NHC_RTP-9", "NHC_RTP-4",
  "NHC_RTP-10", "NHC_RTP-6", "NHC_RTP-5", "NHC_RTP-3", "NHC_RTP-8", "NHC_RTP-7"
)

all_data_txp <- all_data_nhc %>%
  filter(NHC %in% pacientes_txp)

cat("── Subcohorte TXP ──────────────────────────────\n")
cat("Pacientes TXP definidos:", length(pacientes_txp),         "\n")
cat("Pacientes TXP con datos:", n_distinct(all_data_txp$NHC),  "\n")
cat("Npet únicos TXP:        ", n_distinct(all_data_txp$Npet), "\n")
cat("Variantes totales TXP:  ", nrow(all_data_txp),            "\n")

cat("\nDistribución por gen (TXP):\n")
all_data_txp %>%
  count(Protein, name = "n_variantes") %>%
  arrange(Protein) %>%
  print()

sin_datos <- setdiff(pacientes_txp, unique(all_data_txp$NHC))
if (length(sin_datos) > 0) {
  cat("\n⚠️  Pacientes TXP sin variantes tras filtro:",
      paste(sin_datos, collapse = ", "), "\n")
} else {
  cat("\nTodos los pacientes TXP tienen datos tras el filtro.\n")
}

## ---- Clasificación CHARMD: MRA, VSI, PN y VND ----
# ══════════════════════════════════════════════════════════════════════════════
# VECTORES CHARMD — extraídos de CHARMD_UL54/UL56/UL97.xlsx (Tilloy et al. 2024)
# ══════════════════════════════════════════════════════════════════════════════

# ── UL54 ──────────────────────────────────────────────────────────────────────
UL54_MRA <- c(
  "D301N", "E303G", "E303D", "N408D", "N408H", "N408K", "N408S", "N410K",
  "F412C", "F412L", "F412S", "F412V", "D413A", "D413E", "D413N", "D413del",
  "D413Y", "P488R", "N495K", "K500N", "L501F", "L501I", "T503A", "T503I",
  "A505V", "K513E", "K513N", "K513Q", "K513R", "D515E", "D515Y", "L516R",
  "L516W", "I521T", "P522A", "P522S", "P522T", "C524del", "V526L", "C539G",
  "C539R", "D542E", "A543V", "L545F", "L545S", "L545W", "T552S", "T552N",
  "Q578H", "Q578L", "S585A", "D588E", "D588N", "C590F", "F595I", "T700A",
  "V715A", "V715M", "I726T", "E756D", "E756K", "E756Q", "L773V", "L776M",
  "V781I", "V787A", "V787E", "V787L", "L802M", "K805Q", "A809V", "V812L",
  "T813S", "T821I", "P829S", "A834P", "T838A", "G841A", "G841S", "M844T",
  "M844V", "A928T", "V946L", "L957F", "D981del", "981-982del", "A987G"
)
UL54_VSI <- c(
  "Q229K", "E235G", "D247N", "D262N", "D271A", "D284E", "D288N", "N345S",
  "Y380C", "F396L", "L424V", "F460L", "V476G", "V482G", "A492D", "A543S",
  "R581H", "V654G", "S660G", "G667N", "A692G", "F718L", "F718S", "I726V",
  "E793V", "Q795P", "Q795R", "L802V", "G822D", "M827I", "M828V", "P859A",
  "L862F", "V902G", "K947E", "M959T"
)
UL54_POL <- c(
  "N4H", "N4Y", "S8N", "S8T", "A15S", "A17P", "R21C", "S24L", "S32P",
  "V95E", "H142Y", "G143S", "D163H", "T242K", "C256S", "C256L", "N267Y",
  "A269V", "A336T", "I341T", "G347D", "T351A", "V355A", "F357I", "P375L",
  "L394F", "K426R", "G441S", "V450G", "F460S", "S464F", "H465Y", "N467S",
  "A469T", "A473V", "V483A", "R512C", "D515G", "Y518C", "P522L", "Q541P",
  "G604S", "S612N", "P617S", "I622L", "A626V", "P628A", "P628L", "G629S",
  "A631G", "S633F", "V634A", "M640R", "A647S", "S651E", "G653S", "S655L",
  "S660N", "S663N", "F669L", "S676G", "G678S", "V681del", "N685S", "H686Y",
  "A688V", "G689R", "T691S", "T691A", "A692V", "A692S", "A692F", "A693T",
  "S695T", "Q697H", "I722V", "L737M", "V759M", "A786P", "R792S", "R800C",
  "R847H", "N855D", "H863R", "Q868R", "D870H", "V873L", "G874R", "in884T",
  "A885S", "A885T", "in885S", "in885T", "P887S", "L890F", "L890S", "T892I",
  "S897L", "N898D", "E899K", "E903G", "G920S", "V927M", "S932N", "V953A",
  "R954H", "S956Y", "A972V", "E989Q", "G993C", "T995K", "N998D", "S1000L",
  "D1005N", "R1006C", "R1006G", "A1012V", "L1020I", "R1052C", "D1055L",
  "R1070G", "Q1071R", "L1074M", "V1079I", "T1108A", "N1116H", "A1122T",
  "P1129Q", "G1133S", "A1138T", "S1146N", "N1147S", "R1149T", "G1151del",
  "P1153S", "L1156del", "S1162L", "S1235T", "C1241G", "A1108T", "T1122A",
  "L897S", "T885A"
)

# ── UL56 ──────────────────────────────────────────────────────────────────────
UL56_MRA <- c(
  "C25F", "V231A", "V231L", "N232Y", "V236A", "V236L", "V236M", "E237D",
  "L241P", "T244K", "L254F", "L257F", "L257I", "K258E", "F261C", "F261L",
  "Y321C", "C325F", "C325R", "C325Y", "C325W", "M329T", "A365S", "N368D",
  "R369G", "R369M", "R369S", "R369T", "R396S"
)
UL56_VSI <- c(
  "S229F", "L328V", "F345L"
)
UL56_POL <- c(
  "E2K", "V12L", "V33A", "A36V", "E75G", "L122P", "R129H", "A221V",
  "D242G", "A327V", "E393K", "A425V", "G436E", "M442T", "A444G", "S445N",
  "N446del", "N446S", "449-451del", "450-453del", "T452I", "S454N", "G460V",
  "A464V", "A464T", "V471A", "V476A", "E485G", "R507C", "H509N", "Q577R",
  "V582M", "D586N", "A648R", "A648S", "G649D", "L658F", "H698Q", "C721Y",
  "S749N", "V778A", "A779T", "A779V", "A788T", "V793A", "V793P", "P800L",
  "P803A", "Y806T", "T811P", "A812E", "A812G", "G813A", "A820Q", "L831M",
  "G840A", "S848N", "V425A", "N586D"
)

# ── UL97 ──────────────────────────────────────────────────────────────────────
UL97_MRA <- c(
  "L337M", "F342S", "F342Y", "G343A", "V353A", "K355del", "V356G", "K359E",
  "K359Q", "L397R", "L405P", "T409M", "H411L", "H411N", "H411Y", "D456N",
  "M460I", "M460T", "M460V", "V466G", "C480F", "C480R", "C518Y", "H520Q",
  "P521L", "590-600del", "591-594del", "591-607del", "C592F", "C592G",
  "A594E", "A594G", "A594P", "A594S", "A594T", "A594V", "L595F", "L595S",
  "L595W", "L595del", "595-603del", "E596G", "E596Y", "E596Q", "599-603del",
  "K599T", "601-602del", "601-603del", "C603W", "C603Y", "C607Y", "I610T",
  "A613V", "Y617del"
)
UL97_VSI <- c(
  "A427V", "V466M", "H469Y", "R476C", "A478V", "N510S", "M550I", "A582V",
  "A588V", "590-603del", "A591D", "A591V", "E596del", "N597del", "597-598del",
  "597-599del", "G598S", "K599E", "K599R", "L600I", "L600del", "T601M",
  "C603R", "C603S", "C607F", "Y617H", "G623S", "E655K", "T659I",
  "V665I", "A674T"
)
UL97_POL <- c(
  "S7F", "L14F", "Q19E", "Q54E", "N68D", "N108S", "S108N", "R112C", "R112H",
  "D118N", "L126Q", "P132L", "S133P", "V171I", "D184G", "G215A", "A238T",
  "I244V", "P247S", "S249C", "F294C", "K316T", "E317A", "D329H", "V356W",
  "V356L", "K356V", "T429I", "A440V", "A442G", "Q449K", "Q449R", "M460L",
  "G561A", "H587Y", "M615V", "A639T", "A669V", "T75A", "Q126L", "V244I",
  "D68N", "S249N",
  "D605E"   # Reclasificada de VSI a PN: sin efecto fenotipico sobre ganciclovir (Sun et al. 2012)
)

cat("\n── Vectores CHARMD (Tilloy et al. 2024) ────────\n")
cat(sprintf("UL54 → MRA: %d | VSI: %d | POL: %d\n",
            length(UL54_MRA), length(UL54_VSI), length(UL54_POL)))
cat(sprintf("UL56 → MRA: %d | VSI: %d | POL: %d\n",
            length(UL56_MRA), length(UL56_VSI), length(UL56_POL)))
cat(sprintf("UL97 → MRA: %d | VSI: %d | POL: %d\n",
            length(UL97_MRA), length(UL97_VSI), length(UL97_POL)))

# ── Asignación de clase (prioridad: MRA > VSI > POL > NEW) ───────────────────
all_data_nhc <- all_data_nhc %>%
  mutate(class = case_when(
    Protein == "UL54" & Mutation %in% UL54_MRA ~ "MRA",
    Protein == "UL54" & Mutation %in% UL54_VSI ~ "VSI",
    Protein == "UL54" & Mutation %in% UL54_POL ~ "POL",
    Protein == "UL56" & Mutation %in% UL56_MRA ~ "MRA",
    Protein == "UL56" & Mutation %in% UL56_VSI ~ "VSI",
    Protein == "UL56" & Mutation %in% UL56_POL ~ "POL",
    Protein == "UL97" & Mutation %in% UL97_MRA ~ "MRA",
    Protein == "UL97" & Mutation %in% UL97_VSI ~ "VSI",
    Protein == "UL97" & Mutation %in% UL97_POL ~ "POL",
    TRUE ~ "NEW"
  ))

# Regenerar all_data_txp con la columna class ya asignada
all_data_txp <- all_data_nhc %>%
  filter(NHC %in% pacientes_txp)

cat("\nall_data_txp actualizado con class:",
    n_distinct(all_data_txp$NHC), "pacientes |",
    nrow(all_data_txp), "variantes\n")

cat("\n── Distribución de class por gen (todos) ───────\n")
all_data_nhc %>%
  count(Protein, class) %>%
  arrange(Protein, class) %>%
  print()

cat("\n── Variantes NEW por gen (no en CHARMD) ────────\n")
all_data_nhc %>%
  filter(class == "NEW") %>%
  count(Protein, Mutation, sort = TRUE) %>%
  print(n = Inf)

## ---- Identificadores RTP-1 a RTP-11 y tablas de trabajo ----
rtp_map <- tibble::tribble(
  ~NHC,      ~RTP,
  "NHC_RTP-1",   "RTP-1",
  "NHC_RTP-2",  "RTP-2",
  "NHC_RTP-3",    "RTP-3",
  "NHC_RTP-4",    "RTP-4",
  "NHC_RTP-5",   "RTP-5",
  "NHC_RTP-6",  "RTP-6",
  "NHC_RTP-7",   "RTP-7",
  "NHC_RTP-8",  "RTP-8",
  "NHC_RTP-9",  "RTP-9",
  "NHC_RTP-10",  "RTP-10",
  "NHC_RTP-11",  "RTP-11"
)

all_data_nhc <- all_data_nhc %>% left_join(rtp_map, by = "NHC")
all_data_txp <- all_data_txp %>% left_join(rtp_map, by = "NHC")

UL97_data <- all_data_nhc %>% filter(Protein == "UL97")
UL54_data <- all_data_nhc %>% filter(Protein == "UL54")
UL56_data <- all_data_nhc %>% filter(Protein == "UL56")


# ---- 10.1. Selección y características de la cohorte y de los episodios ----

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


# ---- 10.2. Análisis global de variantes detectadas mediante NGS ----

## ---- Muestras analizadas por gen (26 muestras; UL97 26, UL54 24, UL56 11) ----
npet_txp_total <- all_data_nhc %>%
  filter(NHC %in% pacientes_txp) %>%
  distinct(Npet) %>%
  nrow()
cat("Total peticiones únicas (Npet) en subcohorte TXP:", npet_txp_total, "\n\n")

cat("── Npet con resultado válido por gen ────────────────────────────────────\n")
all_data_nhc %>%
  filter(NHC %in% pacientes_txp) %>%
  group_by(Protein) %>%
  summarise(
    n_npet      = n_distinct(Npet),
    pct         = round(n_distinct(Npet) / npet_txp_total * 100, 1),
    n_pacientes = n_distinct(NHC),
    .groups     = "drop"
  ) %>%
  filter(Protein %in% c("UL97", "UL54", "UL56")) %>%
  arrange(Protein) %>%
  print()

npet_ul97 <- all_data_nhc %>%
  filter(NHC %in% pacientes_txp, Protein == "UL97") %>% distinct(Npet)
npet_ul54 <- all_data_nhc %>%
  filter(NHC %in% pacientes_txp, Protein == "UL54") %>% distinct(Npet)
npet_ul56 <- all_data_nhc %>%
  filter(NHC %in% pacientes_txp, Protein == "UL56") %>% distinct(Npet)

cat("\n── Npet en UL97 pero no en UL54 (cobertura insuficiente u omitidas) ─────\n")
npet_solo_ul97 <- setdiff(npet_ul97$Npet, npet_ul54$Npet)
cat("N =", length(npet_solo_ul97), "\n")
cat(npet_solo_ul97, sep = "\n")

cat("\n── Npet en UL54 pero no en UL97 ─────────────────────────────────────────\n")
npet_solo_ul54 <- setdiff(npet_ul54$Npet, npet_ul97$Npet)
cat("N =", length(npet_solo_ul54), "\n")
cat(npet_solo_ul54, sep = "\n")

cat("\n── Npet con UL56 analizado ───────────────────────────────────────────────\n")
cat("N =", nrow(npet_ul56), "\n")
cat("Porcentaje sobre total TXP:",
    round(nrow(npet_ul56) / npet_txp_total * 100, 1), "%\n")

cat("\n── Desglose por Npet: genes analizados (TRUE/FALSE) ─────────────────────\n")
all_data_nhc %>%
  filter(NHC %in% pacientes_txp) %>%
  distinct(Npet, NHC) %>%
  left_join(npet_ul97 %>% mutate(UL97 = TRUE), by = "Npet") %>%
  left_join(npet_ul54 %>% mutate(UL54 = TRUE), by = "Npet") %>%
  left_join(npet_ul56 %>% mutate(UL56 = TRUE), by = "Npet") %>%
  mutate(across(c(UL97, UL54, UL56), ~replace_na(., FALSE))) %>%
  left_join(
    all_data_nhc %>% filter(NHC %in% pacientes_txp) %>% distinct(Npet, RTP),
    by = "Npet"
  ) %>%
  arrange(RTP, Npet) %>%
  print(n = Inf)

# ── Diagnóstico calidad: resumen de cobertura y Q-Score en subcohorte TXP ────
cat("\n── Resumen de métricas de calidad (subcohorte TXP, Prevalence ≥5%) ──────\n")
all_data_txp %>%
  filter(Protein %in% c("UL97", "UL54", "UL56")) %>%
  group_by(Protein) %>%
  summarise(
    cob_media   = round(mean(cobertura_total, na.rm = TRUE), 0),
    cob_mediana = round(median(cobertura_total, na.rm = TRUE), 0),
    cob_min     = round(min(cobertura_total, na.rm = TRUE), 0),
    q_mediana   = round(median(q_score, na.rm = TRUE), 1),
    pct_balance = round(mean(balance_ok, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) %>%
  rename(
    "Gen"             = Protein,
    "Cob. media (×)"  = cob_media,
    "Cob. mediana (×)"= cob_mediana,
    "Cob. mínima (×)" = cob_min,
    "Q-Score mediana" = q_mediana,
    "% balance FW/RV" = pct_balance
  ) %>%
  print()

## ---- Complementaria. Variantes por umbral, gen y clase (RTP) ----
thresholds <- c(5, 10, 15, 20)
genes      <- c("UL97", "UL54", "UL56")

calcular_celda <- function(data, gene, clase, thr) {
  d <- data %>%
    filter(Protein == gene, class == clase, Prevalence >= thr)
  n_var <- n_distinct(d$Mutation)
  n_pac <- n_distinct(d$NHC)
  if (n_var == 0) return("0 / 0")
  paste0(n_var, " / ", n_pac)
}

# ── Construir dataframe ───────────────────────────────────────────────────────
filas <- map_dfr(genes, function(gene) {
  total_pac <- n_distinct(filter(all_data_txp, Protein == gene)$NHC)
  map_dfr(thresholds, function(thr) {
    tibble(
      Gen    = if (thr == thresholds[1]) paste0(gene, "  (n = ", total_pac, " pac.)") else "",
      Umbral = paste0("\u2265 ", thr, "%"),
      MRA    = calcular_celda(all_data_txp, gene, "MRA",   thr),
      POL    = calcular_celda(all_data_txp, gene, "POL",   thr),
      VSI    = calcular_celda(all_data_txp, gene, "VSI",   thr),
      Nuevas = calcular_celda(all_data_txp, gene, "NEW",   thr)
    )
  })
})

# ── Identificar filas con MRA > 0 para colorear ───────────────────────────────
filas_mra <- which(filas$MRA != "0 / 0")

# ── Bordes ────────────────────────────────────────────────────────────────────
b_azul <- officer::fp_border(color = "#2E5FA3", width = 1.5)
b_sep  <- officer::fp_border(color = "#2E5FA3", width = 1.2)
b_gris <- officer::fp_border(color = "#CCCCCC", width = 0.5)
b_head <- officer::fp_border(color = "#FFFFFF", width = 0.8)

# ── Flextable ─────────────────────────────────────────────────────────────────
ft <- flextable::flextable(filas) %>%
  
  # Cabeceras
  flextable::set_header_labels(
    Gen = "Gen", Umbral = "Umbral",
    MRA = "MRA", POL = "POL", VSI = "VSI", Nuevas = "Nuevas"
  ) %>%
  flextable::add_header_row(
    values    = c("", "", "Variantes \u00fanicas / Pacientes afectados", "", "", ""),
    colwidths = c(1, 1, 1, 1, 1, 1),
    top       = TRUE
  ) %>%
  flextable::merge_at(i = 1, j = 3:6, part = "header") %>%
  
  # Fusión Gen
  flextable::merge_at(i = 1:4,  j = 1, part = "body") %>%
  flextable::merge_at(i = 5:8,  j = 1, part = "body") %>%
  flextable::merge_at(i = 9:12, j = 1, part = "body") %>%
  
  # Alineación
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = 1:2, align = "left", part = "body") %>%
  flextable::valign(valign = "center", part = "all") %>%
  
  # Tipografía base
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1, part = "body") %>%
  
  # Cabecera: fondo azul oscuro, texto blanco
  flextable::bg(part = "header", bg = "#2E5FA3") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::border_inner_v(border = b_head, part = "header") %>%
  
  # Cuerpo: sombreado alternado por bloque gen
  flextable::bg(i = 1:4,  bg = "#EEF4FB", part = "body") %>%
  flextable::bg(i = 5:8,  bg = "#FFFFFF", part = "body") %>%
  flextable::bg(i = 9:12, bg = "#EEF4FB", part = "body") %>%
  
  # Columna Umbral: fondo gris neutro para guía visual
  flextable::bg(j = 2, bg = "#F5F5F5", part = "body") %>%
  flextable::color(j = 2, color = "#555555", part = "body") %>%
  flextable::italic(j = 2, part = "body") %>%
  
  # MRA con color si > 0
  flextable::bold(i = filas_mra, j = 3, part = "body") %>%
  flextable::color(i = filas_mra, j = 3, color = "#B22222", part = "body") %>%
  
  # Bordes
  flextable::border_outer(border = b_azul, part = "all") %>%
  flextable::border_inner_h(border = b_gris, part = "body") %>%
  flextable::border_inner_v(border = b_gris, part = "body") %>%
  flextable::hline(i = c(4, 8), border = b_sep, part = "body") %>%
  flextable::hline(i = 1, border = b_sep, part = "header") %>%
  
  # Anchos
  flextable::width(j = 1,   width = 4.0, unit = "cm") %>%
  flextable::width(j = 2,   width = 1.8, unit = "cm") %>%
  flextable::width(j = 3:6, width = 2.8, unit = "cm") %>%
  
  # Padding
  flextable::padding(padding = 6, part = "all") %>%
  
  # Nota al pie
  flextable::add_footer_lines(
    "Valores expresados como: variantes \u00fanicas detectadas / pacientes con al menos una variante de esa clase al umbral indicado. MRA: mutaci\u00f3n de resistencia a antivirales; POL: polimorfismo natural; VSI: variante de significado incierto; Nuevas: variantes no clasificadas en CHARMD (Tilloy et al., 2024). An\u00e1lisis restringido a la subcohorte de trasplante pulmonar (n\u00famero de pacientes con datos por gen indicado entre par\u00e9ntesis)."
  ) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft

guardar_tabla(ft, "Complementaria_variantes_umbral_RTP")

## ---- Complementaria. Variantes por umbral, gen y clase (pacientes no RTP) ----
# Complemento: todos los pacientes con datos que NO son TXP
all_data_notxp <- all_data_nhc %>%
  filter(!NHC %in% pacientes_txp)

cat("── Subcohorte NO TXP ───────────────────────────\n")
cat("Pacientes NO TXP con datos:", n_distinct(all_data_notxp$NHC), "\n")
cat("Npet únicos NO TXP:        ", n_distinct(all_data_notxp$Npet), "\n")
cat("Variantes totales NO TXP:  ", nrow(all_data_notxp), "\n")

# ── Construir dataframe ───────────────────────────────────────────────────────
filas_notxp <- map_dfr(genes, function(gene) {
  total_pac <- n_distinct(filter(all_data_notxp, Protein == gene)$NHC)
  map_dfr(thresholds, function(thr) {
    tibble(
      Gen    = if (thr == thresholds[1]) paste0(gene, "  (n = ", total_pac, " pac.)") else "",
      Umbral = paste0("\u2265 ", thr, "%"),
      MRA    = calcular_celda(all_data_notxp, gene, "MRA", thr),
      POL    = calcular_celda(all_data_notxp, gene, "POL", thr),
      VSI    = calcular_celda(all_data_notxp, gene, "VSI", thr),
      Nuevas = calcular_celda(all_data_notxp, gene, "NEW", thr)
    )
  })
})

filas_mra_notxp <- which(filas_notxp$MRA != "0 / 0")

ft_notxp <- flextable::flextable(filas_notxp) %>%
  flextable::set_header_labels(
    Gen = "Gen", Umbral = "Umbral",
    MRA = "MRA", POL = "POL", VSI = "VSI", Nuevas = "Nuevas"
  ) %>%
  flextable::add_header_row(
    values    = c("", "", "Variantes \u00fanicas / Pacientes afectados", "", "", ""),
    colwidths = c(1, 1, 1, 1, 1, 1),
    top       = TRUE
  ) %>%
  flextable::merge_at(i = 1, j = 3:6, part = "header") %>%
  flextable::merge_at(i = 1:4,  j = 1, part = "body") %>%
  flextable::merge_at(i = 5:8,  j = 1, part = "body") %>%
  flextable::merge_at(i = 9:12, j = 1, part = "body") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = 1:2, align = "left", part = "body") %>%
  flextable::valign(valign = "center", part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1, part = "body") %>%
  flextable::bg(part = "header", bg = "#2E5FA3") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::border_inner_v(border = b_head, part = "header") %>%
  flextable::bg(i = 1:4,  bg = "#EEF4FB", part = "body") %>%
  flextable::bg(i = 5:8,  bg = "#FFFFFF", part = "body") %>%
  flextable::bg(i = 9:12, bg = "#EEF4FB", part = "body") %>%
  flextable::bg(j = 2, bg = "#F5F5F5", part = "body") %>%
  flextable::color(j = 2, color = "#555555", part = "body") %>%
  flextable::italic(j = 2, part = "body") %>%
  flextable::bold(i = filas_mra_notxp, j = 3, part = "body") %>%
  flextable::color(i = filas_mra_notxp, j = 3, color = "#B22222", part = "body") %>%
  flextable::border_outer(border = b_azul, part = "all") %>%
  flextable::border_inner_h(border = b_gris, part = "body") %>%
  flextable::border_inner_v(border = b_gris, part = "body") %>%
  flextable::hline(i = c(4, 8), border = b_sep, part = "body") %>%
  flextable::hline(i = 1, border = b_sep, part = "header") %>%
  flextable::width(j = 1,   width = 4.0, unit = "cm") %>%
  flextable::width(j = 2,   width = 1.8, unit = "cm") %>%
  flextable::width(j = 3:6, width = 2.8, unit = "cm") %>%
  flextable::padding(padding = 6, part = "all") %>%
  flextable::add_footer_lines(
    "Valores expresados como: variantes \u00fanicas detectadas / pacientes con al menos una variante de esa clase al umbral indicado. MRA: mutaci\u00f3n de resistencia a antivirales; POL: polimorfismo natural; VSI: variante de significado incierto; Nuevas: variantes no clasificadas en CHARMD (Tilloy et al., 2024). An\u00e1lisis restringido a pacientes no pertenecientes a la subcohorte de trasplante pulmonar."
  ) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft_notxp

guardar_tabla(ft_notxp, "Complementaria_variantes_umbral_noRTP")

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


# ---- 10.2.1. Variantes asociadas a resistencia antiviral (MRA) ----

## ---- Tabla 24. Características de las MRA detectadas en la paciente RTP-2 ----
nhc_rtp2 <- "NHC_RTP-2"

# ── Datos fijos (CHARMD) ──────────────────────────────────────────────────────
charmd <- tibble(
  Protein    = c("UL97", "UL97", "UL97", "UL54"),
  Mutation   = c("A594V", "T409M", "H411Y", "T503I"),
  Resistencia_fenotipica = c(
    "Ganciclovir",
    "Maribavir",
    "Maribavir",
    "Ganciclovir y Cidofovir"
  )
)

# ── Calcular metricas desde all_data_txp ──────────────────────────────────────
mra_raw <- all_data_txp %>%
  filter(class == "MRA", Prevalence >= 5)

mra_stats <- mra_raw %>%
  group_by(Protein, Mutation) %>%
  summarise(
    n_muestras   = n_distinct(Npet),
    prev_media   = round(mean(Prevalence, na.rm = TRUE), 1),
    prev_de      = round(sd(Prevalence,   na.rm = TRUE), 1),
    cv_mediana   = round(median(as.numeric(CV), na.rm = TRUE), 0),
    cv_p25       = round(quantile(as.numeric(CV), 0.25, na.rm = TRUE), 0),
    cv_p75       = round(quantile(as.numeric(CV), 0.75, na.rm = TRUE), 0),
    cob_mediana  = round(median(cobertura_total, na.rm = TRUE), 0),
    cob_p25      = round(quantile(cobertura_total, 0.25, na.rm = TRUE), 0),
    cob_p75      = round(quantile(cobertura_total, 0.75, na.rm = TRUE), 0),
    cob_min      = round(min(cobertura_total, na.rm = TRUE), 0),
    q_med        = round(median(q_score, na.rm = TRUE), 0),
    q_p25        = round(quantile(q_score, 0.25, na.rm = TRUE), 0),
    q_p75        = round(quantile(q_score, 0.75, na.rm = TRUE), 0),
    n_balance_ok = sum(balance_ok, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(charmd, by = c("Protein", "Mutation"))

# ── Umbral minimo de deteccion ────────────────────────────────────────────────
umbral_min_mra <- map_dfr(c(5, 10, 15, 20), function(thr) {
  all_data_txp %>%
    filter(class == "MRA", Prevalence >= thr) %>%
    distinct(Protein, Mutation) %>%
    mutate(umbral = thr)
}) %>%
  group_by(Protein, Mutation) %>%
  summarise(umbral_min = paste0("\u2265", min(umbral), "%"), .groups = "drop")

mra_stats <- mra_stats %>%
  left_join(umbral_min_mra, by = c("Protein", "Mutation"))

# ── Total muestras validas por gen en RTP-2 (calculado) ───────────────────────
n_muestras_gen <- all_data_nhc %>%
  filter(NHC == nhc_rtp2, Protein %in% c("UL97", "UL54")) %>%
  group_by(Protein) %>%
  summarise(n_total = n_distinct(Npet), .groups = "drop")

cat("── Denominadores de la Tabla 25 ──\n")
print(n_muestras_gen)

n_ul97 <- n_muestras_gen$n_total[n_muestras_gen$Protein == "UL97"]
n_ul54 <- n_muestras_gen$n_total[n_muestras_gen$Protein == "UL54"]

# ── Helper de formato numerico sin espacios sobrantes ─────────────────────────
fnum <- function(x) format(x, big.mark = ".", trim = TRUE)

mra_stats <- mra_stats %>%
  left_join(n_muestras_gen, by = "Protein") %>%
  mutate(
    pct   = round(n_muestras / n_total * 100, 1),
    col_n = paste0(n_muestras, " (",
                   format(pct, decimal.mark = ",", trim = TRUE), "%)"),
    col_prev = if_else(
      is.na(prev_de),
      paste0(format(prev_media, decimal.mark = ",", trim = TRUE), " (n = 1)"),
      paste0(format(prev_media, decimal.mark = ",", trim = TRUE),
             " \u00b1 ", format(prev_de, decimal.mark = ",", trim = TRUE))
    ),
    col_cv = if_else(
      is.na(cv_mediana),
      "\u2014",
      paste0(fnum(cv_mediana), " (", fnum(cv_p25), "\u2013", fnum(cv_p75), ")")
    ),
    col_cob = if_else(
      n_muestras == 1,
      paste0(fnum(cob_mediana), "\u00d7"),
      paste0(fnum(cob_mediana), "\u00d7 (", fnum(cob_p25), "\u2013",
             fnum(cob_p75), "; m\u00edn. ", fnum(cob_min), "\u00d7)")
    ),
    col_q = if_else(
      n_muestras == 1,
      as.character(q_med),
      paste0(q_med, " (", q_p25, "\u2013", q_p75, ")")
    ),
    col_balance = paste0(n_balance_ok, "/", n_muestras)
  )

# ── Orden explicito ───────────────────────────────────────────────────────────
orden_mra <- tibble(
  Protein  = c("UL97", "UL97", "UL97", "UL54"),
  Mutation = c("A594V", "T409M", "H411Y", "T503I")
)

filas_mra <- orden_mra %>%
  left_join(mra_stats, by = c("Protein", "Mutation")) %>%
  mutate(
    Gen = case_when(
      Mutation == "A594V" ~ paste0("UL97\n(n = ", n_ul97, " muestras)"),
      Mutation == "T503I" ~ paste0("UL54\n(n = ", n_ul54, " muestras)"),
      TRUE ~ ""
    )
  ) %>%
  dplyr::select(
    Gen, Mutation, col_n, col_prev, col_cv,
    col_cob, col_q, col_balance, umbral_min, Resistencia_fenotipica
  )

names(filas_mra) <- c(
  "Gen", "MRA", "Muestras, n (%)",
  "Prevalencia (media \u00b1 DE, %)",
  "CV en muestras con MRA\nmediana (RIC), UI/mL",
  "Cobertura\n(mediana, RIC; m\u00edn.)",
  "Q-score\n(mediana, RIC)",
  "Balance FW/RV\n(n muestras)",
  "Umbral m\u00ednimo de\ndetecci\u00f3n",
  "Resistencia fenot\u00edpica*"
)

# ── Bordes ────────────────────────────────────────────────────────────────────
b_azul <- officer::fp_border(color = "#2E5FA3", width = 1.5)
b_sep  <- officer::fp_border(color = "#2E5FA3", width = 1.2)
b_gris <- officer::fp_border(color = "#CCCCCC", width = 0.5)
b_head <- officer::fp_border(color = "#FFFFFF", width = 0.8)

# ── Flextable ─────────────────────────────────────────────────────────────────
ft_mra <- flextable::flextable(filas_mra) %>%
  flextable::merge_at(i = 1:3, j = 1, part = "body") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = 1:2, align = "left", part = "body") %>%
  flextable::valign(valign = "center", part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1:2, part = "body") %>%
  flextable::bg(part = "header", bg = "#2E5FA3") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::border_inner_v(border = b_head, part = "header") %>%
  flextable::bg(i = 1:3, bg = "#EEF4FB", part = "body") %>%
  flextable::bg(i = 4,   bg = "#FFFFFF", part = "body") %>%
  flextable::border_outer(border = b_azul, part = "all") %>%
  flextable::border_inner_h(border = b_gris, part = "body") %>%
  flextable::border_inner_v(border = b_gris, part = "body") %>%
  flextable::hline(i = 3, border = b_sep, part = "body") %>%
  flextable::width(j = 1,  width = 2.4, unit = "cm") %>%
  flextable::width(j = 2,  width = 1.8, unit = "cm") %>%
  flextable::width(j = 3,  width = 2.0, unit = "cm") %>%
  flextable::width(j = 4,  width = 2.5, unit = "cm") %>%
  flextable::width(j = 5,  width = 2.8, unit = "cm") %>%
  flextable::width(j = 6,  width = 3.0, unit = "cm") %>%
  flextable::width(j = 7,  width = 1.8, unit = "cm") %>%
  flextable::width(j = 8,  width = 2.1, unit = "cm") %>%
  flextable::width(j = 9,  width = 2.0, unit = "cm") %>%
  flextable::width(j = 10, width = 3.0, unit = "cm") %>%
  flextable::padding(padding = 5, part = "all") %>%
  flextable::add_footer_lines(
    "*Resistencia fenot\u00edpica seg\u00fan la base de datos CHARMD (Tilloy et al., 2024). CV: carga viral; DE: desviaci\u00f3n est\u00e1ndar; RIC: rango intercuart\u00edlico. El umbral m\u00ednimo de detecci\u00f3n indica el umbral de frecuencia de variante m\u00e1s bajo al que la variante es detectable. Cobertura: n\u00famero total de lecturas forward + reverse que cubren la posici\u00f3n, expresado como mediana, rango intercuart\u00edlico y valor m\u00ednimo. Q-score: puntuaci\u00f3n de calidad de llamada. Balance FW/RV: n\u00famero de muestras con soporte en ambas cadenas sobre el total de muestras con la variante."
  ) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft_mra

guardar_tabla(ft_mra, "Tabla_24", apaisada = TRUE)

## ---- Complementaria. Muestras con MRA según umbral ----
# ── Calcular mra_por_umbral ───────────────────────────────────────────────────
thresholds <- c(5, 10, 15, 20)

mra_por_umbral <- map_dfr(thresholds, function(thr) {
  all_data_txp %>%
    filter(class == "MRA", Prevalence >= thr) %>%
    group_by(Protein) %>%
    summarise(n_muestras = n_distinct(Npet), .groups = "drop") %>%
    mutate(umbral = paste0("≥", thr, "%"))
}) %>%
  filter(Protein %in% c("UL97", "UL54")) %>%
  mutate(
    umbral  = factor(umbral, levels = paste0("≥", thresholds, "%")),
    Protein = factor(Protein, levels = c("UL97", "UL54"))
  )

# ── Parámetros ────────────────────────────────────────────────────────────────
colores <- c("UL97" = "#4472C4", "UL54" = "#FFC000")

ref_sanger <- tibble(
  Protein     = c("UL97", "UL54"),
  ref_value   = c(8, 6),
  color_linea = c("#104E8B", "darkgoldenrod3")
)

mra_por_umbral <- mra_por_umbral %>%
  mutate(vjust_label = if_else(Protein == "UL54", 1.8, -1))

# ── Figura ────────────────────────────────────────────────────────────────────
fig_mra_umbrales <- ggplot(mra_por_umbral, aes(x = umbral, y = n_muestras,
                           color = Protein, group = Protein)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  geom_text(
    aes(label = n_muestras, vjust = vjust_label),
    size     = 13 / .pt,
    fontface = "bold",
    color    = "black",
    show.legend = FALSE
  ) +
  geom_hline(
    data = ref_sanger,
    aes(yintercept = ref_value),
    color     = alpha(ref_sanger$color_linea, 0.5),
    linetype  = "dashed",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  scale_color_manual(values = colores, name = NULL) +
  scale_y_continuous(
    limits = c(0, 13),
    breaks = seq(0, 12, by = 3),
    expand = c(0, 0)
  ) +
  labs(
    x = "Umbral de prevalencia",
    y = "N muestras ≥1 MRA"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position    = "top",
    legend.text        = element_text(face = "bold", size = 13),
    axis.title         = element_text(face = "bold", size = 13),
    axis.text          = element_text(color = "black", size = 13),
    axis.title.x       = element_text(face = "bold", size = 13, margin = margin(t = 10)),
    axis.title.y       = element_text(face = "bold", size = 13),
    panel.grid.major.y = element_line(color = "#EEEEEE", linewidth = 0.4)
  )

# [guardado antiguo] ggsave(
# [guardado antiguo]   filename = "~/Desktop/Figuras_tesis/figura_MRA_umbrales.png",
# [guardado antiguo]   plot     = last_plot(),
# [guardado antiguo]   width    = 10,
# [guardado antiguo]   height   = 10,
# [guardado antiguo]   units    = "cm",
# [guardado antiguo]   dpi      = 300
# [guardado antiguo] )

guardar_fig(fig_mra_umbrales, "Complementaria_MRA_umbrales", ancho = 10, alto = 10, dpi = 300)


# ---- 10.2.2. Polimorfismos naturales (PN) ----

## ---- Datos de PN ----
all_data_txp %>%
  filter(class == "POL", Prevalence >= 5) %>%
  group_by(Protein) %>%
  summarise(
    n_variantes = n_distinct(Mutation),
    n_pacientes = n_distinct(NHC),
    .groups = "drop"
  )

all_data_txp %>%
  filter(class == "POL", Prevalence >= 5) %>%
  group_by(Protein, Mutation) %>%
  summarise(
    NHC_list = paste(sort(unique(NHC)), collapse = ", "),
    n_pac    = n_distinct(NHC),
    .groups  = "drop"
  ) %>%
  arrange(Protein, desc(n_pac)) %>%
  print(n = Inf)

## ---- Anexo I. Polimorfismos naturales identificados mediante NGS ----
# El bloque guarda el Word apaisado (Anexo_I.docx) en la carpeta de salidas.
## Requiere all_data_txp (con class y columnas DeepChek). Ejecutar antes el bloque del
## Anexo II para disponer de fmt_n y fmt_1, o dejar las dos funciones de abajo.


fmt_n <- function(x) format(round(x), big.mark = ".", decimal.mark = ",", trim = TRUE)
fmt_1 <- function(x) sub("\\.", ",", sprintf("%.1f", x))

denom <- all_data_txp %>% group_by(Protein) %>%
  summarise(N_m = n_distinct(Npet), N_p = n_distinct(RTP), .groups = "drop")
etq_gen <- setNames(sprintf("p%s (%d muestras, %d pacientes)", denom$Protein, denom$N_m, denom$N_p), denom$Protein)

pn <- all_data_txp %>%
  filter(class == "POL", Prevalence >= 5) %>%
  mutate(mut_tot = mut_fw + mut_rv, balance = pmin(mut_fw, mut_rv) / mut_tot) %>%
  group_by(Protein, Mutation) %>%
  summarise(
    n_m = n_distinct(Npet), n_p = n_distinct(RTP),
    pacientes = paste(unique(RTP[order(as.numeric(sub("RTP-", "", RTP)))]), collapse = ", "),
    prev_mean = mean(Prevalence), prev_sd = sd(Prevalence),
    cob_med = median(cobertura_total), cob_q1 = quantile(cobertura_total, .25), cob_q3 = quantile(cobertura_total, .75),
    lect_min = min(mut_tot), fw_rv_min = paste0(mut_fw[which.min(mut_tot)], "/", mut_rv[which.min(mut_tot)]),
    bal_min = min(balance), q_min = min(q_score), .groups = "drop") %>%
  left_join(denom, by = "Protein") %>%
  mutate(
    Posicion   = as.numeric(gsub("[^0-9]", "", Mutation)),
    Gen        = factor(etq_gen[Protein], levels = etq_gen[c("UL97", "UL54", "UL56")]),
    Muestras   = sprintf("%d (%s)", n_m, fmt_1(100 * n_m / N_m)),
    Pacientes  = sprintf("%d (%s)", n_p, fmt_1(100 * n_p / N_p)),
    Frecuencia = ifelse(n_m == 1, fmt_1(prev_mean), sprintf("%s \u00b1 %s", fmt_1(prev_mean), fmt_1(prev_sd))),
    Cobertura  = ifelse(n_m == 1, paste0(fmt_n(cob_med), " \u00d7"),
                        sprintf("%s \u00d7 (%s\u2013%s)", fmt_n(cob_med), fmt_n(cob_q1), fmt_n(cob_q3))),
    Lecturas   = sprintf("%d (%s)", lect_min, fw_rv_min),
    Hebra      = fmt_1(100 * bal_min),
    Q          = as.character(q_min)
  ) %>%
  arrange(Gen, Posicion) %>%
  select(Gen, Mutation, Muestras, Pacientes, pacientes, Frecuencia, Cobertura, Lecturas, Hebra, Q)

cat("PN por gen:\n"); print(table(pn$Gen))

# ── Flextable con fila combinada por gen ─────────────────────────────────────
pn_grp <- as_grouped_data(pn, groups = "Gen")

ft_pn <- as_flextable(pn_grp, hide_grouplabel = TRUE) %>%
  set_header_labels(
    Mutation = "Variante", Muestras = "Muestras, n (%)", Pacientes = "Pacientes, n (%)",
    pacientes = "Pacientes (Tabla X)",
    Frecuencia = "Frecuencia de variante (%), media \u00b1 DE",
    Cobertura = "Cobertura por posici\u00f3n, mediana (RIC)",
    Lecturas = "Lecturas mutadas, m\u00ednimo (fw/rv)",
    Hebra = "Hebra minoritaria, m\u00ednimo (%)",
    Q = "Q-score, m\u00ednimo") %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "body") %>%
  bold(part = "header") %>%
  bold(i = ~ !is.na(Gen), bold = TRUE) %>%                       # filas de gen en negrita
  bg(i = ~ !is.na(Gen), bg = "#1F3B73") %>%                      # fondo azul oscuro en las filas de gen
  color(i = ~ !is.na(Gen), color = "white") %>%
  bg(part = "header", bg = "#1F3B73") %>% color(part = "header", color = "white") %>%
  bold(j = "Mutation", bold = TRUE) %>%
  flextable::font(fontname = "Source Sans Pro", part = "all") %>%
  fontsize(size = 7.5, part = "all") %>%
  border_inner(border = fp_border(color = "grey70", width = 0.5)) %>%
  border_outer(border = fp_border(color = "grey40", width = 1)) %>%
  add_footer_lines(paste0(
    "PN: polimorfismos naturales (variantes recogidas en CHARMD, versi\u00f3n actualizada a 20 de abril de 2026, sin efecto ",
    "sobre la susceptibilidad antiviral; incluye D605E, reclasificada como PN). Muestras y pacientes sobre el total con ",
    "resultado v\u00e1lido para cada gen. Frecuencia y cobertura resumen las muestras en que se detect\u00f3 la variante. ",
    "Lecturas mutadas: n\u00famero m\u00ednimo de lecturas que soportan la variante, con su reparto entre hebras directa (fw) e ",
    "inversa (rv). Hebra minoritaria: proporci\u00f3n m\u00ednima de lecturas mutadas en la hebra menos representada. ",
    "El filtro por Q-score del sistema experto estaba desactivado.")) %>%
  fontsize(size = 6.5, part = "footer") %>%
  set_table_properties(layout = "autofit", width = 1)

ft_pn

# ── Exportar a Word en horizontal (carpeta Descargas) ─────────────────────────
doc <- read_docx() %>%
  body_add_par("Anexo. Polimorfismos naturales identificados mediante NGS", style = "heading 2") %>%
  body_add_flextable(ft_pn) %>%
  body_end_section_landscape()
print(doc, target = file.path(dir_salida, "Anexo_I.docx"))

## ---- Complementaria. Tabla de PN con cobertura y Q-score ----
# ── Datos ─────────────────────────────────────────────────────────────────────
pol_stats <- all_data_txp %>%
  filter(class == "POL", Prevalence >= 5) %>%
  group_by(Protein, Mutation) %>%
  summarise(
    n_muestras   = n_distinct(Npet),
    n_pacientes  = n_distinct(NHC),
    prev_media   = round(mean(Prevalence, na.rm = TRUE), 1),
    prev_de      = round(sd(Prevalence, na.rm = TRUE), 1),
    cob_mediana  = round(median(cobertura_total, na.rm = TRUE), 0),
    cob_p25      = round(quantile(cobertura_total, 0.25, na.rm = TRUE), 0),
    cob_p75      = round(quantile(cobertura_total, 0.75, na.rm = TRUE), 0),
    cob_min      = round(min(cobertura_total, na.rm = TRUE), 0),
    q_med        = round(median(q_score, na.rm = TRUE), 0),
    q_p25        = round(quantile(q_score, 0.25, na.rm = TRUE), 0),
    q_p75        = round(quantile(q_score, 0.75, na.rm = TRUE), 0),
    .groups = "drop"
  )

# ── Umbral máximo de detección ────────────────────────────────────────────────
umbral_max_pol <- map_dfr(c(5, 10, 15, 20), function(thr) {
  all_data_txp %>%
    filter(class == "POL", Prevalence >= thr) %>%
    distinct(Protein, Mutation) %>%
    mutate(umbral = thr)
}) %>%
  group_by(Protein, Mutation) %>%
  summarise(
    umbral_max = paste0("\u2265", max(umbral), "%"),
    .groups = "drop"
  )

# ── Total muestras y pacientes por gen ────────────────────────────────────────
n_total_gen <- tibble(
  Protein         = c("UL97", "UL54", "UL56"),
  n_tot_muestras  = c(26, 24, 11),
  n_tot_pacientes = c(11, 11, 3)
)

# ── Unir todo ─────────────────────────────────────────────────────────────────
pol_stats <- pol_stats %>%
  left_join(umbral_max_pol, by = c("Protein", "Mutation")) %>%
  left_join(n_total_gen,    by = "Protein") %>%
  mutate(
    pct_muestras  = round(n_muestras  / n_tot_muestras  * 100, 1),
    pct_pacientes = round(n_pacientes / n_tot_pacientes * 100, 1),
    col_prev = if_else(
      is.na(prev_de),
      paste0(format(prev_media, decimal.mark = ","), " (n=1)"),
      paste0(format(prev_media, decimal.mark = ","),
             " \u00b1 ",
             format(prev_de, decimal.mark = ","))
    ),
    col_muestras  = paste0(n_muestras, " (",
                           format(pct_muestras,  decimal.mark = ","), "%)"),
    col_pacientes = paste0(n_pacientes, " (",
                           format(pct_pacientes, decimal.mark = ","), "%)"),
    col_cobertura = paste0(
      format(cob_mediana, big.mark = "."), "\u00d7",
      " (", format(cob_p25, big.mark = "."),
      "\u2013", format(cob_p75, big.mark = "."), "\u00d7;",
      " m\u00edn. ", format(cob_min, big.mark = "."), "\u00d7)"
    ),
    col_qscore = paste0(
      q_med, " (", q_p25, "\u2013", q_p75, ")"
    )
  )

# ── Orden ─────────────────────────────────────────────────────────────────────
pol_stats <- pol_stats %>%
  mutate(Protein = factor(Protein, levels = c("UL97", "UL54", "UL56"))) %>%
  arrange(Protein, desc(n_muestras))

# ── Dataframe para flextable ──────────────────────────────────────────────────
pol_ft <- pol_stats %>%
  group_by(Protein) %>%
  mutate(
    Gen = if_else(row_number() == 1,
                  paste0(as.character(Protein),
                         "\n(n = ", n_tot_muestras, " muestras)"),
                  "")
  ) %>%
  ungroup() %>%
  dplyr::select(Gen, Mutation, col_muestras, col_pacientes,
                col_prev, umbral_max, col_cobertura, col_qscore)

names(pol_ft) <- c(
  "Gen", "PN",
  "Muestras, n (%)",
  "Pacientes, n (%)",
  "Prevalencia\n(media \u00b1 DE, %)",
  "Umbral m\u00e1x.\ndetecci\u00f3n",
  "Cobertura\n(mediana, RIC; m\u00edn.)",
  "Q-Score\n(mediana, RIC)"
)

# ── Filas por bloque ──────────────────────────────────────────────────────────
n_ul97 <- sum(pol_stats$Protein == "UL97")
n_ul54 <- sum(pol_stats$Protein == "UL54")
n_ul56 <- sum(pol_stats$Protein == "UL56")

filas_ul97 <- 1:n_ul97
filas_ul54 <- (n_ul97 + 1):(n_ul97 + n_ul54)
filas_ul56 <- (n_ul97 + n_ul54 + 1):(n_ul97 + n_ul54 + n_ul56)
filas_sep  <- c(n_ul97, n_ul97 + n_ul54)

# ── Bordes ────────────────────────────────────────────────────────────────────
b_azul <- officer::fp_border(color = "#2E5FA3", width = 1.5)
b_sep  <- officer::fp_border(color = "#2E5FA3", width = 1.2)
b_gris <- officer::fp_border(color = "#CCCCCC", width = 0.5)
b_head <- officer::fp_border(color = "#FFFFFF", width = 0.8)

# ── Flextable ─────────────────────────────────────────────────────────────────
ft_pol <- flextable::flextable(pol_ft) %>%
  flextable::merge_at(i = filas_ul97, j = 1, part = "body") %>%
  flextable::merge_at(i = filas_ul54, j = 1, part = "body") %>%
  flextable::merge_at(i = filas_ul56, j = 1, part = "body") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = 1:2, align = "left", part = "body") %>%
  flextable::valign(valign = "center", part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1:2, part = "body") %>%
  flextable::bg(part = "header", bg = "#2E5FA3") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::border_inner_v(border = b_head, part = "header") %>%
  flextable::bg(i = filas_ul97, bg = "#EEF4FB", part = "body") %>%
  flextable::bg(i = filas_ul54, bg = "#FFFFFF", part = "body") %>%
  flextable::bg(i = filas_ul56, bg = "#EEF4FB", part = "body") %>%
  flextable::border_outer(border = b_azul, part = "all") %>%
  flextable::border_inner_h(border = b_gris, part = "body") %>%
  flextable::border_inner_v(border = b_gris, part = "body") %>%
  flextable::hline(i = filas_sep, border = b_sep, part = "body") %>%
  flextable::width(j = 1, width = 3.0, unit = "cm") %>%
  flextable::width(j = 2, width = 2.0, unit = "cm") %>%
  flextable::width(j = 3, width = 2.5, unit = "cm") %>%
  flextable::width(j = 4, width = 2.5, unit = "cm") %>%
  flextable::width(j = 5, width = 3.0, unit = "cm") %>%
  flextable::width(j = 6, width = 2.0, unit = "cm") %>%
  flextable::width(j = 7, width = 3.5, unit = "cm") %>%
  flextable::width(j = 8, width = 2.5, unit = "cm") %>%
  flextable::padding(padding = 6, part = "all") %>%
  flextable::add_footer_lines(paste0(
    "PN: polimorfismo natural. DE: desviaci\u00f3n est\u00e1ndar. RIC: rango intercuart\u00edlico. ",
    "Los porcentajes de muestras y pacientes se calculan sobre el total con datos v\u00e1lidos por gen ",
    "(UL97: 26 muestras, 11 pacientes; UL54: 24 muestras, 11 pacientes; UL56: 11 muestras, 3 pacientes). ",
    "El umbral m\u00e1ximo de detecci\u00f3n indica el umbral de frecuencia de variante m\u00e1s alto al que la variante es detectable. ",
    "Cobertura: suma de lecturas forward + reverse en la posici\u00f3n (mediana y RIC entre todas las muestras con la variante; m\u00ednimo observado). ",
    "Q-Score: mediana y RIC de la puntuaci\u00f3n de calidad de llamada de variante."
  )) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft_pol

guardar_tabla(ft_pol, "Complementaria_PN", apaisada = TRUE)

## ---- Complementaria. Cobertura de los PN ----
# ── Datos ─────────────────────────────────────────────────────────────────────
pol_cob <- all_data_txp %>%
  filter(class == "POL", Prevalence >= 5) %>%
  mutate(
    Position = as.integer(str_extract(Mutation, "\\d+")),
    Protein  = factor(Protein, levels = c("UL97", "UL54", "UL56"))
  ) %>%
  group_by(Protein, Mutation, Position) %>%
  summarise(
    mediana  = median(cobertura_total, na.rm = TRUE),
    p25      = quantile(cobertura_total, 0.25, na.rm = TRUE),
    p75      = quantile(cobertura_total, 0.75, na.rm = TRUE),
    n        = n_distinct(Npet),
    prev_max = max(Prevalence, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  arrange(Protein, Position) %>%
  mutate(
    Mutation  = factor(Mutation, levels = unique(Mutation)),
    unica     = n == 1,
    baja_prev = prev_max < 20
  )

# ── Colores por gen ───────────────────────────────────────────────────────────
col_gen <- c("UL97" = "#4472C4", "UL54" = "#C8880A", "UL56" = "#C0392B")

# ── Shapes: rombo (18) si prev_max < 20%, círculo (16/21) si no ──────────────
# Combinamos unica y baja_prev en un indicador de shape
pol_cob <- pol_cob %>%
  mutate(
    shape_tipo = case_when(
      baja_prev & unica  ~ "baja_unica",
      baja_prev & !unica ~ "baja_multi",
      !baja_prev & unica ~ "norm_unica",
      TRUE               ~ "norm_multi"
    )
  )

shape_map <- c(
  "norm_multi" = 16,   # círculo relleno
  "norm_unica" = 21,   # círculo hueco
  "baja_multi" = 18,   # rombo relleno
  "baja_unica" = 23    # rombo hueco
)

# ── Plot ──────────────────────────────────────────────────────────────────────
p_cob_pol <- ggplot(pol_cob,
                    aes(x = Mutation, y = mediana, color = Protein)) +
  
  # RIC (solo si n > 1)
  geom_errorbar(
    data = pol_cob %>% filter(!unica),
    aes(ymin = p25, ymax = p75),
    width = 0.35, linewidth = 0.6, alpha = 0.7
  ) +
  
  # Puntos con shape según tipo
  geom_point(
    aes(shape = shape_tipo,
        fill  = ifelse(unica | shape_tipo == "baja_unica", "white", Protein)),
    size = 3.2, stroke = 0.8
  ) +
  scale_shape_manual(
    values = shape_map,
    labels = c(
      "norm_multi" = "Frecuencia variante \u226520%",
      "norm_unica" = "Frecuencia variante \u226520%, muestra \u00fanica",
      "baja_multi" = "Frecuencia variante <20%",
      "baja_unica" = "Frecuencia variante <20%, muestra \u00fanica"
    ),
    name = NULL
  ) +
  scale_fill_manual(values = c(col_gen, "white" = "white"), guide = "none") +
  scale_color_manual(values = col_gen, guide = "none") +
  
  scale_y_continuous(
    labels  = function(x) format(x, big.mark = ".", decimal.mark = ","),
    expand  = expansion(mult = c(0.05, 0.12))
  ) +
  
  facet_wrap(~ Protein, scales = "free_x", nrow = 1) +
  
  labs(
    x       = NULL,
    y       = "Cobertura de posici\u00f3n (x)",
    caption = paste0(
    )
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1,
                                    size = 12, face = "bold"),
    axis.text.y      = element_text(size = 10),
    strip.text       = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "#EEF4FB", color = NA),
    panel.spacing    = unit(1.2, "cm"),
    legend.position  = "bottom",
    legend.text      = element_text(size = 10),
    plot.caption     = element_text(size = 10, color = "grey40",
                                    hjust = 0, margin = margin(t = 8)),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin      = margin(t = 10, r = 15, b = 10, l = 10)
  )

p_cob_pol

# [guardado antiguo] ggsave(
# [guardado antiguo]   "~/Desktop/Figuras_tesis/figura_cobertura_PN.png",
# [guardado antiguo]   p_cob_pol,
# [guardado antiguo]   width = 32, height = 14, units = "cm", dpi = 400,
# [guardado antiguo]   bg = "white"
# [guardado antiguo] )

guardar_fig(p_cob_pol, "Complementaria_cobertura_PN", ancho = 32, alto = 14, dpi = 400)


# ---- 10.2.3. Variantes de significado incierto (VSI) ----

## ---- Datos y tabla de VSI (vacía: D605E se reclasificó como PN) ----
all_data_txp %>%
  filter(class == "VSI", Prevalence >= 5) %>%
  group_by(Protein, Mutation) %>%
  summarise(
    n_muestras  = n_distinct(Npet),
    n_pacientes = n_distinct(NHC),
    NHC_list    = paste(sort(unique(NHC)), collapse = ", "),
    prev_media  = round(mean(Prevalence, na.rm = TRUE), 1),
    prev_de     = round(sd(Prevalence, na.rm = TRUE), 1),
    prev_min    = round(min(Prevalence, na.rm = TRUE), 1),
    prev_max    = round(max(Prevalence, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(Protein, desc(n_muestras)) %>%
  print(n = Inf)


all_data_txp %>%
  filter(class == "VSI", Prevalence >= 5) %>%
  select(Npet, Freg, Protein, Mutation, Prevalence, CV) %>%
  arrange(Protein)

# ── Datos VSI ─────────────────────────────────────────────────────────────────
vsi_stats <- all_data_txp %>%
  filter(class == "VSI", Prevalence >= 5) %>%
  group_by(Protein, Mutation) %>%
  summarise(
    n_muestras  = n_distinct(Npet),
    n_pacientes = n_distinct(NHC),
    cob_mediana  = round(median(cobertura_total, na.rm = TRUE), 0),
    cob_p25      = round(quantile(cobertura_total, 0.25, na.rm = TRUE), 0),
    cob_p75      = round(quantile(cobertura_total, 0.75, na.rm = TRUE), 0),
    cob_min      = round(min(cobertura_total, na.rm = TRUE), 0),
    q_med        = round(median(q_score, na.rm = TRUE), 0),
    q_p25        = round(quantile(q_score, 0.25, na.rm = TRUE), 0),
    q_p75        = round(quantile(q_score, 0.75, na.rm = TRUE), 0),
    n_balance_ok = sum(balance_ok, na.rm = TRUE),
    prev_media  = round(mean(Prevalence, na.rm = TRUE), 1),
    prev_de     = round(sd(Prevalence, na.rm = TRUE), 1),
    cv_mediana  = round(median(as.numeric(CV), na.rm = TRUE), 0),
    cv_p25      = round(quantile(as.numeric(CV), 0.25, na.rm = TRUE), 0),
    cv_p75      = round(quantile(as.numeric(CV), 0.75, na.rm = TRUE), 0),
    .groups = "drop"
  )

# ── Umbral máximo de detección ────────────────────────────────────────────────
umbral_max_vsi <- map_dfr(c(5, 10, 15, 20), function(thr) {
  all_data_txp %>%
    filter(class == "VSI", Prevalence >= thr) %>%
    distinct(Protein, Mutation) %>%
    mutate(umbral = thr)
}) %>%
  group_by(Protein, Mutation) %>%
  summarise(
    umbral_max = paste0("\u2265", max(umbral), "%"),
    .groups = "drop"
  )

# ── Total muestras por gen (denominador RTP-2) ────────────────────────────────
n_muestras_gen <- tibble(
  Protein = c("UL97", "UL54"),
  n_total = c(14, 12)
)

# ── Unir todo ─────────────────────────────────────────────────────────────────
vsi_stats <- vsi_stats %>%
  left_join(umbral_max_vsi, by = c("Protein", "Mutation")) %>%
  left_join(n_muestras_gen, by = "Protein") %>%
  mutate(
    pct      = round(n_muestras / n_total * 100, 1),
    col_n    = paste0(n_muestras, " (",
                      format(pct, decimal.mark = ","), "%)"),
    col_prev = if_else(
      is.na(prev_de),
      paste0(format(prev_media, decimal.mark = ","), " (n=1)"),
      paste0(format(prev_media, decimal.mark = ","),
             " \u00b1 ",
             format(prev_de, decimal.mark = ","))
    ),
    col_cv = if_else(
      is.na(cv_mediana),
      "\u2014",
      paste0(format(cv_mediana, big.mark = "."),
             " (", format(cv_p25, big.mark = "."),
             "\u2013", format(cv_p75, big.mark = "."), ")")
    )
  )

# ── Orden: UL97 primero, luego UL54 ──────────────────────────────────────────
vsi_stats <- vsi_stats %>%
  mutate(Protein = factor(Protein, levels = c("UL97", "UL54"))) %>%
  arrange(Protein)

# ── Dataframe para flextable ──────────────────────────────────────────────────
vsi_ft <- vsi_stats %>%
  group_by(Protein) %>%
  mutate(
    col_cob = if_else(
      n_muestras == 1,
      paste0(format(cob_mediana, big.mark = "."), "×"),
      paste0(
        format(cob_mediana, big.mark = "."), "× (",
        format(cob_p25, big.mark = "."), "–",
        format(cob_p75, big.mark = "."), "; mín. ",
        format(cob_min, big.mark = "."), "×)"
      )
    ),
    col_q = if_else(
      n_muestras == 1,
      as.character(q_med),
      paste0(q_med, " (", q_p25, "–", q_p75, ")")
    ),
    col_balance = paste0(n_balance_ok, "/", n_muestras),
    Gen = if_else(row_number() == 1,
                  paste0(as.character(Protein),
                         "\n(n = ", n_total, " muestras)"),
                  "")
  ) %>%
  ungroup() %>%
  dplyr::select(
    Gen, Mutation, col_n, col_prev, col_cv,
    col_cob, col_q, col_balance, umbral_max
  )
names(vsi_ft) <- c(
  "Gen", "VSI",
  "Muestras, n (%)",
  "Prevalencia (media ± DE, %)",
  "CV en muestra con VSI\n(UI/mL)",
  "Cobertura\n(mediana, RIC; mín.)",
  "Q-score\n(mediana, RIC)",
  "Balance FW/RV\n(n muestras)",
  "Umbral máximo de\ndetección"
)

# ── Bordes ────────────────────────────────────────────────────────────────────
b_azul <- officer::fp_border(color = "#2E5FA3", width = 1.5)
b_sep  <- officer::fp_border(color = "#2E5FA3", width = 1.2)
b_gris <- officer::fp_border(color = "#CCCCCC", width = 0.5)
b_head <- officer::fp_border(color = "#FFFFFF", width = 0.8)

n_ul97 <- sum(vsi_stats$Protein == "UL97")
n_ul54 <- sum(vsi_stats$Protein == "UL54")
n_total_filas <- nrow(vsi_ft)

filas_ul97 <- if (n_ul97 > 0) 1:n_ul97 else integer(0)
filas_ul54 <- if (n_ul54 > 0) (n_ul97 + 1):n_total_filas else integer(0)

ft_vsi <- flextable::flextable(vsi_ft) %>%
  { if (length(filas_ul97) > 1) flextable::merge_at(., i = filas_ul97, j = 1, part = "body") else . } %>%
  { if (length(filas_ul54) > 1) flextable::merge_at(., i = filas_ul54, j = 1, part = "body") else . } %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = 1:2, align = "left", part = "body") %>%
  flextable::valign(valign = "center", part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1:2, part = "body") %>%
  flextable::bg(part = "header", bg = "#2E5FA3") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::border_inner_v(border = b_head, part = "header") %>%
  { if (length(filas_ul97) > 0) flextable::bg(., i = filas_ul97, bg = "#EEF4FB", part = "body") else . } %>%
  { if (length(filas_ul54) > 0) flextable::bg(., i = filas_ul54, bg = "#FFFFFF",  part = "body") else . } %>%
  flextable::border_outer(border = b_azul, part = "all") %>%
  flextable::border_inner_h(border = b_gris, part = "body") %>%
  flextable::border_inner_v(border = b_gris, part = "body") %>%
  { if (n_ul97 > 0 & n_ul54 > 0) flextable::hline(., i = n_ul97, border = b_sep, part = "body") else . } %>%
  flextable::width(j = 1, width = 2.6, unit = "cm") %>%
  flextable::width(j = 2, width = 1.8, unit = "cm") %>%
  flextable::width(j = 3, width = 2.2, unit = "cm") %>%
  flextable::width(j = 4, width = 2.7, unit = "cm") %>%
  flextable::width(j = 5, width = 2.6, unit = "cm") %>%
  flextable::width(j = 6, width = 3.0, unit = "cm") %>%
  flextable::width(j = 7, width = 2.0, unit = "cm") %>%
  flextable::width(j = 8, width = 2.2, unit = "cm") %>%
  flextable::width(j = 9, width = 2.1, unit = "cm") %>%
  flextable::padding(padding = 6, part = "all") %>%
  flextable::add_footer_lines(
    "VSI: variante de significado incierto. DE: desviaci\u00f3n est\u00e1ndar. Los porcentajes se calculan sobre el total de muestras v\u00e1lidas del paciente RTP-2 por gen (UL97: 14 muestras; UL54: 12 muestras). El umbral m\u00e1ximo de detecci\u00f3n indica el umbral de prevalencia al\u00e9lica m\u00e1s alto al que la variante es detectable."
  ) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft_vsi


# ---- 10.2.4. Variantes no descritas previamente (VND) ----

## ---- Datos de VND ----
# ── Resumen actualizado de variantes NEW ──────────────────────────────────────
new_vars <- all_data_txp %>%
  filter(class == "NEW", Prevalence >= 5) %>%
  mutate(Position = as.numeric(Position)) %>%
  group_by(Protein, Mutation, Position) %>%
  summarise(
    n_muestras  = n_distinct(Npet),
    n_pacientes = n_distinct(NHC),
    prev_media  = round(mean(Prevalence, na.rm = TRUE), 1),
    prev_max    = round(max(Prevalence, na.rm = TRUE), 1),
    .groups = "drop"
  )

cat("Total NEW por gen:\n")
new_vars %>% count(Protein) %>% print()

cat("\nNEW por gen y umbral:\n")
map_dfr(c(5, 10, 15, 20), function(thr) {
  all_data_txp %>%
    filter(class == "NEW", Prevalence >= thr) %>%
    group_by(Protein) %>%
    summarise(n_variantes = n_distinct(Mutation),
              n_pacientes = n_distinct(NHC), .groups = "drop") %>%
    mutate(umbral = paste0("≥", thr, "%"))
}) %>%
  arrange(Protein, umbral) %>%
  print(n = Inf)


new_vars_dom <- new_vars %>%
  mutate(
    dominio = case_when(
      Protein == "UL97" & Position >= 335 & Position <= 370 ~ "P-loop",
      Protein == "UL97" & Position >= 409 & Position <= 411 ~ "Gatekeeper",
      Protein == "UL97" & Position >= 440 & Position <= 480 ~ "Loop catalítico",
      Protein == "UL97" & Position >= 520 & Position <= 607 ~ "Región canónica GCV",
      Protein == "UL54" & Position >= 408 & Position <= 412 ~ "Región IV",
      Protein == "UL54" & Position >= 501 & Position <= 545 ~ "Delta-C",
      Protein == "UL54" & Position >= 700 & Position <= 715 ~ "Región II",
      Protein == "UL54" & Position >= 725 & Position <= 964 ~ "Unión ATP/dNTP",
      Protein == "UL54" & Position >= 802 & Position <= 821 ~ "Región III",
      Protein == "UL54" & Position >= 910 & Position <= 913 ~ "Región I",
      Protein == "UL54" & Position >= 962 & Position <= 987 ~ "Región V",
      Protein == "UL54" & Position >= 1153 & Position <= 1159 ~ "NLS-A",
      Protein == "UL54" & Position >= 1213 & Position <= 1242 ~ "UL44 binding",
      Protein == "UL56" & Position == 25 ~ "N-terminal",
      Protein == "UL56" & Position >= 202 & Position <= 220 ~ "LAGLIDADG",
      Protein == "UL56" & Position >= 231 & Position <= 369 ~ "Región LMV",
      Protein == "UL56" & Position >= 514 & Position <= 522 ~ "Interacción UL89",
      Protein == "UL56" & Position >= 671 & Position <= 680 ~ "Interacción UL89",
      Protein == "UL56" & Position >= 757 & Position <= 764 ~ "Interacción UL89",
      TRUE ~ "Fuera de dominio funcional"
    ),
    en_dominio = dominio != "Fuera de dominio funcional"
  ) %>%
  arrange(Protein, desc(en_dominio), desc(n_pacientes), desc(prev_max))

cat("── En dominio funcional ─────────────────────────\n")
new_vars_dom %>% filter(en_dominio) %>% print(n = Inf)

cat("\n── Fuera de dominio pero relevantes ─────────────\n")
new_vars_dom %>%
  filter(!en_dominio, n_pacientes >= 2 | prev_max >= 15) %>%
  print(n = Inf)

## ---- Complementaria. VND en dominios funcionales (todos los genes) ----
# ── Datos: variantes nuevas en dominio funcional ──────────────────────────────
new_dom_resumen <- new_vars_dom %>%
  filter(en_dominio) %>%
  mutate(
    rtp = case_when(
      Protein == "UL54" & Mutation %in% c(
        "D711Y","D854G","E530D","E846*","F702L","F763L","F817L",
        "G967C","G971C","I833V","I840V","M1219I","M959V","P704S",
        "P712H","P744H","P797H","R1157S","R876L","R901L","S913G",
        "G340C","Y281H") ~ "RTP-2",
      Protein == "UL54" & Mutation %in% c("A714V","L736M","R1157M","T760A") ~ "RTP-1",
      Protein == "UL54" & Mutation %in% c("A835T","V772A") ~ "RTP-10",
      Protein == "UL54" & Mutation == "Q511R" ~ "RTP-3",
      Protein == "UL56" & Mutation == "N330S" ~ "RTP-11",
      Protein == "UL56" & Mutation == "Y281H" ~ "RTP-2",
      Protein == "UL97" & Mutation == "D574G" ~ "RTP-1",
      Protein == "UL97" & Mutation == "S604P" ~ "RTP-10",
      Protein == "UL97" & Mutation == "G340C" ~ "RTP-2",
      TRUE ~ "otro"
    )
  ) %>%
  arrange(Protein, dominio, desc(prev_max)) %>%
  mutate(
    Protein = factor(Protein, levels = c("UL54", "UL56", "UL97")),
    col_prev = paste0(format(prev_max, decimal.mark = ","), "%")
  ) %>%
  select(Protein, Mutation, dominio, prev_max, rtp)

# ── Construir dataframe para flextable ────────────────────────────────────────
ft_new_dom <- new_dom_resumen %>%
  group_by(Protein) %>%
  mutate(
    Gen = if_else(row_number() == 1, as.character(Protein), "")
  ) %>%
  ungroup() %>%
  select(Gen, Mutation, dominio, prev_max, rtp)

names(ft_new_dom) <- c(
  "Gen", "Variante",
  "Dominio funcional",
  "Prevalencia máx. (%)",
  "Paciente"
)

# ── Filas por bloque gen ──────────────────────────────────────────────────────
n_ul54 <- sum(new_dom_resumen$Protein == "UL54")
n_ul56 <- sum(new_dom_resumen$Protein == "UL56")
n_ul97 <- sum(new_dom_resumen$Protein == "UL97")

filas_ul54 <- 1:n_ul54
filas_ul56 <- (n_ul54 + 1):(n_ul54 + n_ul56)
filas_ul97 <- (n_ul54 + n_ul56 + 1):(n_ul54 + n_ul56 + n_ul97)
filas_sep  <- c(n_ul54, n_ul54 + n_ul56)

# ── Bordes ────────────────────────────────────────────────────────────────────
b_azul <- officer::fp_border(color = "#2E5FA3", width = 1.5)
b_sep  <- officer::fp_border(color = "#2E5FA3", width = 1.2)
b_gris <- officer::fp_border(color = "#CCCCCC", width = 0.5)
b_head <- officer::fp_border(color = "#FFFFFF", width = 0.8)

# ── Flextable ─────────────────────────────────────────────────────────────────
ft_new <- flextable::flextable(ft_new_dom) %>%
  flextable::merge_at(i = filas_ul54, j = 1, part = "body") %>%
  flextable::merge_at(i = filas_ul56, j = 1, part = "body") %>%
  flextable::merge_at(i = filas_ul97, j = 1, part = "body") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = 1:3, align = "left", part = "body") %>%
  flextable::valign(valign = "center", part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1:2, part = "body") %>%
  flextable::bg(part = "header", bg = "#2E5FA3") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::border_inner_v(border = b_head, part = "header") %>%
  flextable::bg(i = filas_ul54, bg = "#EEF4FB", part = "body") %>%
  flextable::bg(i = filas_ul56, bg = "#FFFFFF", part = "body") %>%
  flextable::bg(i = filas_ul97, bg = "#EEF4FB", part = "body") %>%
  flextable::border_outer(border = b_azul, part = "all") %>%
  flextable::border_inner_h(border = b_gris, part = "body") %>%
  flextable::border_inner_v(border = b_gris, part = "body") %>%
  flextable::hline(i = filas_sep, border = b_sep, part = "body") %>%
  flextable::width(j = 1, width = 2.0, unit = "cm") %>%
  flextable::width(j = 2, width = 2.5, unit = "cm") %>%
  flextable::width(j = 3, width = 4.5, unit = "cm") %>%
  flextable::width(j = 4, width = 2.5, unit = "cm") %>%
  flextable::width(j = 5, width = 2.0, unit = "cm") %>%
  flextable::padding(padding = 5, part = "all") %>%
  flextable::add_footer_lines(
    "Tabla restringida a variantes no clasificadas en CHARMD (Tilloy et al., 2024) detectadas en al menos un dominio funcional conocido, al umbral de prevalencia \u2265 5%. Prevalencia m\u00e1xima expresada como el valor m\u00e1ximo observado entre todas las muestras en que se detect\u00f3 la variante."
  ) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft_new

guardar_tabla(ft_new, "Complementaria_VND_dominios")

## ---- Figura 45. VND en UL97 ----
dominios_ul97 <- data.frame(
  start = c(230),
  end   = c(620),
  label = c("Dominio quinasa"),
  fill  = c("lightskyblue")
)

hotspots_ul97 <- data.frame(
  start = c(335, 397, 519, 590),
  end   = c(370, 480, 521, 607),
  label = c("", "", "", ""),
  fill  = c("#B03060", "#8E44AD", "#E67E22", "#E67E22")
)

ul97_umbrales <- all_data_txp %>%
  filter(class == "NEW", Protein == "UL97") %>%
  group_by(Mutation, Position) %>%
  summarise(
    prev_max   = max(Prevalence, na.rm = TRUE),
    umbral_max = case_when(
      max(Prevalence) >= 20 ~ 20,
      max(Prevalence) >= 15 ~ 15,
      max(Prevalence) >= 10 ~ 10,
      TRUE                  ~ 5
    ),
    .groups = "drop"
  ) %>%
  mutate(Position = as.numeric(Position))

ul97 <- new_vars_dom %>%
  filter(Protein == "UL97") %>%
  left_join(ul97_umbrales %>% select(Mutation, umbral_max), by = "Mutation")

nivel <- c("5" = 0.6, "10" = 1.1, "15" = 1.6, "20" = 2.1)
paso  <- 0.4

df <- ul97 %>%
  mutate(
    y_nivel    = nivel[as.character(umbral_max)],
    en_quinasa = Position >= 230 & Position <= 620,
    y_point    = ifelse(en_quinasa, -y_nivel, y_nivel),
    color      = ifelse(en_quinasa, "#C0392B", "#2980B9"),
    label_mut  = case_when(
      en_quinasa                        ~ Mutation,
      prev_max >= 10 | n_pacientes >= 4 ~ Mutation,
      TRUE                              ~ NA_character_
    )
  )

df_lab_azul <- df %>%
  filter(!is.na(label_mut), !en_quinasa) %>%
  arrange(Position) %>%
  mutate(prox_group = cumsum(c(TRUE, diff(Position) > 30))) %>%
  group_by(prox_group) %>%
  mutate(
    within_rank = row_number(),
    y_lab       = 3.0 + (within_rank - 1) * 0.60
  ) %>%
  ungroup()

df_lab_rojo <- df %>%
  filter(!is.na(label_mut), en_quinasa) %>%
  arrange(Position) %>%
  mutate(prox_group = cumsum(c(TRUE, diff(Position) > 30))) %>%
  group_by(prox_group) %>%
  mutate(
    within_rank = row_number(),
    y_lab       = -(3.0 + (within_rank - 1) * 0.60)
  ) %>%
  ungroup()

ymax_plot <- max(
  max(df_lab_azul$y_lab, na.rm = TRUE),
  max(abs(df_lab_rojo$y_lab), na.rm = TRUE)
) + 0.5

niveles_ref <- data.frame(
  y     = c(0.6, 1.1, 1.6, 2.1, -0.6, -1.1, -1.6, -2.1),
  label = rep(c("\u22655%", "\u226510%", "\u226515%", "\u226520%"), 2),
  lado  = rep(c("arriba", "abajo"), each = 4)
)

p_ul97 <- ggplot() +
  
  # Líneas de referencia — empiezan después de las etiquetas, acaban en aa 707
  geom_segment(data = niveles_ref,
               aes(x = -42, xend = 707, y = y, yend = y),
               color = "grey85", linewidth = 0.3, linetype = "dashed") +
  
  # Etiquetas de prevalencia
  geom_text(data = niveles_ref %>% filter(lado == "arriba"),
            aes(x = -55, y = y, label = label),
            hjust = 1, size = 4.0, color = "grey50", fontface = "bold") +
  geom_text(data = niveles_ref %>% filter(lado == "abajo"),
            aes(x = -55, y = y, label = label),
            hjust = 1, size = 4.0, color = "grey50", fontface = "bold") +
  
  # Dominios
  geom_rect(data = dominios_ul97,
            aes(xmin = start, xmax = end, ymin = 0.0, ymax = 0.40, fill = fill),
            alpha = 0.9, color = "white", linewidth = 0.25) +
  geom_text(data = dominios_ul97 %>% mutate(mid = (start + end) / 2),
            aes(x = mid, y = 0.20, label = label),
            size = 4.0, color = "black", fontface = "bold",
            vjust = 0.5, hjust = 0.5) +
  
  # Hotspots
  geom_rect(data = hotspots_ul97,
            aes(xmin = start, xmax = end, ymin = -0.40, ymax = 0.0, fill = fill),
            alpha = 0.9, color = "white", linewidth = 0.25) +
  geom_text(data = hotspots_ul97 %>% mutate(mid = (start + end) / 2),
            aes(x = mid, y = -0.20, label = label),
            size = 3, color = "white", fontface = "bold",
            vjust = 0.5, hjust = 0.5) +
  
  scale_fill_identity() +
  
  # Backbone aa 1–707 con tapas
  geom_segment(aes(x = 1,   xend = 707, y = 0,     yend = 0),
               color = "grey10", linewidth = 1.0) +
  geom_segment(aes(x = 1,   xend = 1,   y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  geom_segment(aes(x = 707, xend = 707, y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  
  # Etiquetas extremos backbone
  annotate("text", x = 1,   y = -0.28, label = "1",   size = 4.0,
           color = "grey10", vjust = 1) +
  annotate("text", x = 707, y = -0.28, label = "707", size = 4.0,
           color = "grey10", vjust = 1) +
  
  # Stems
  geom_segment(data = df,
               aes(x = Position, xend = Position,
                   y = ifelse(en_quinasa, -0.42, 0.42),
                   yend = y_point, color = color),
               linewidth = 0.45, alpha = 0.8) +
  
  # Puntos
  geom_point(data = df,
             aes(x = Position, y = y_point, color = color, size = prev_max),
             alpha = 0.88) +
  scale_color_identity() +
  scale_size_continuous(range = c(1.2, 4.5), guide = "none") +
  
  # Guías azules
  geom_segment(data = df_lab_azul,
               aes(x = Position, xend = Position,
                   y = y_point + 0.08, yend = y_lab - 0.08),
               color = "#2980B9", linewidth = 0.3, alpha = 0.4,
               linetype = "dotted") +
  
  # Guías rojas
  geom_segment(data = df_lab_rojo,
               aes(x = Position, xend = Position,
                   y = y_point - 0.08, yend = y_lab + 0.08),
               color = "#C0392B", linewidth = 0.3, alpha = 0.4,
               linetype = "dotted") +
  
  # Etiquetas azules
  geom_text(data = df_lab_azul,
            aes(x = Position, y = y_lab, label = label_mut),
            color = "#2980B9", size = 4.0, fontface = "bold",
            hjust = 0.5, vjust = 0) +
  
  # Etiquetas rojas
  geom_text(data = df_lab_rojo,
            aes(x = Position, y = y_lab, label = label_mut),
            color = "#C0392B", size = 4.0, fontface = "bold",
            hjust = 0.5, vjust = 1) +
  
  scale_x_continuous(
    limits = c(-100, 730),
    breaks = NULL,
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(limits = c(-ymax_plot, ymax_plot)) +
  labs(title = NULL, x = "", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(t = 8, r = 12, b = 8, l = 35)
  )

leyenda_hotspots <- ggplot(data.frame(
  x     = c(1.5, 3.7, 6),
  y     = c(1, 1, 1),
  label = c("GCV/MBV-R (335-370)",
            "MBV-R (397-480)",
            "GCV-R cod\u00f3n 520 + (590-607)"),
  color = c("#B03060", "#8E44AD", "#E67E22")
)) +
  geom_point(aes(x, y, color = color), size = 4, shape = 15) +
  geom_text(aes(x + 0.15, y, label = label),
            hjust = 0, size = 4.0, color = "grey20") +
  scale_color_identity() +
  scale_x_continuous(limits = c(0.8, 9.5)) +
  scale_y_continuous(limits = c(0.5, 1.5)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(t = 0, b = 6, l = 35)
  )

figura_ul97 <- p_ul97 / leyenda_hotspots +
  plot_layout(heights = c(10, 1))

figura_ul97
# [guardado antiguo] ggsave(
# [guardado antiguo]   "~/Desktop/Figuras_tesis/figura_lollipop_UL97.png",
# [guardado antiguo]   figura_ul97,
# [guardado antiguo]   width = 30, height = 15, units = "cm", dpi = 400,
# [guardado antiguo]   bg = "white"
# [guardado antiguo] )

guardar_fig(figura_ul97, "Figura_45", ancho = 30, alto = 15, dpi = 400)

## ---- Complementaria. VND del dominio quinasa de UL97 ----
ul97_quinasa_datos <- all_data_txp %>%
  filter(class == "NEW", Protein == "UL97", Prevalence >= 5) %>%
  mutate(Position = as.numeric(Position)) %>%
  filter(Position >= 230, Position <= 620) %>%
  group_by(Mutation, Position) %>%
  summarise(
    Paciente       = paste(sort(unique(RTP)), collapse = ", "),
    N_muestras     = n_distinct(Npet),
    Prev_min       = round(min(Prevalence), 1),
    Prev_max       = round(max(Prevalence), 1),
    # ── Métricas de calidad agregadas ────────────────────────────────────
    cob_mediana = round(median(cobertura_total, na.rm = TRUE), 0),
    cob_p25     = round(quantile(cobertura_total, 0.25, na.rm = TRUE), 0),
    cob_p75     = round(quantile(cobertura_total, 0.75, na.rm = TRUE), 0),
    cob_min     = round(min(cobertura_total, na.rm = TRUE), 0),
    q_med = round(median(q_score, na.rm = TRUE), 0),
    q_p25 = round(quantile(q_score, 0.25, na.rm = TRUE), 0),
    q_p75 = round(quantile(q_score, 0.75, na.rm = TRUE), 0),
    n_balance_ok   = sum(balance_ok, na.rm = TRUE),   # muestras con FW y RV > 0
    .groups = "drop"
  ) %>%
  mutate(
    Prevalencia = ifelse(
      Prev_min == Prev_max,
      paste0(format(Prev_max, decimal.mark = ","), "%"),
      paste0(format(Prev_min, decimal.mark = ","), "–",
             format(Prev_max, decimal.mark = ","), "%")
    ),
    Cobertura = paste0(
      format(cob_mediana, big.mark = "."), "× (",
      format(cob_p25, big.mark = "."), "–",
      format(cob_p75, big.mark = "."), "×; mín. ",
      format(cob_min, big.mark = "."), "×)"
    ),
    QScore = paste0(
      q_med, " (",
      q_p25, "–",
      q_p75, ")"
    ),
    Balance = paste0(n_balance_ok, "/", N_muestras)
  ) %>%
  arrange(Position)

hotspot_map <- tibble::tribble(
  ~Mutation, ~Hotspot,
  "K316R",   "No",
  "G340C",   "GCV/MBV-R (aa 335–370)",
  "R412Q",   "MBV-R (aa 397–480)",
  "W424*",   "MBV-R (aa 397–480)",
  "D481N",   "No",
  "A497P",   "No",
  "R516C",   "No",
  "D574G",   "No",
  "S604P",   "GCV-R (aa 590–607)",
  "L608P",   "No"
)

ul97_quinasa_datos <- ul97_quinasa_datos %>%
  left_join(hotspot_map, by = "Mutation") %>%
  select(Mutation, Paciente, N_muestras, Prevalencia,
         Cobertura, QScore, Balance, Hotspot)

names(ul97_quinasa_datos) <- c(
  "Mutación", "Paciente", "N\nmuestras", "Prevalencia",
  "Cobertura\n(mediana, RIC; mín.)", "Q-Score\n(mediana)",
  "Balance FW/RV\n(n muestras)", "Hotspot"
)

filas_hotspot <- which(ul97_quinasa_datos$Hotspot != "No")
filas_no      <- which(ul97_quinasa_datos$Hotspot == "No")
fila_w424     <- which(ul97_quinasa_datos[[1]] == "W424*")

b_azul <- officer::fp_border(color = "#2E5FA3", width = 1.5)
b_gris <- officer::fp_border(color = "#CCCCCC", width = 0.5)
b_head <- officer::fp_border(color = "#FFFFFF", width = 0.8)

ft_ul97_quinasa <- flextable::flextable(ul97_quinasa_datos) %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = c(1, 2, 8), align = "left", part = "body") %>%
  flextable::valign(valign = "center", part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1, part = "body") %>%
  flextable::bold(i = fila_w424, part = "body") %>%
  flextable::bg(i = fila_w424,   bg = "#FDECEA", part = "body") %>%
  flextable::bg(i = filas_hotspot[!filas_hotspot %in% fila_w424],
                bg = "#FEF9E7", part = "body") %>%
  flextable::bg(i = filas_no, bg = "#F0F4F8", part = "body") %>%
  flextable::bg(part = "header", bg = "#2E5FA3") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::border_inner_v(border = b_head, part = "header") %>%
  flextable::border_outer(border = b_azul, part = "all") %>%
  flextable::border_inner_h(border = b_gris, part = "body") %>%
  flextable::border_inner_v(border = b_gris, part = "body") %>%
  flextable::width(j = 1, width = 2.0, unit = "cm") %>%
  flextable::width(j = 2, width = 3.8, unit = "cm") %>%
  flextable::width(j = 3, width = 1.4, unit = "cm") %>%
  flextable::width(j = 4, width = 2.2, unit = "cm") %>%
  flextable::width(j = 5, width = 2.8, unit = "cm") %>%
  flextable::width(j = 6, width = 1.8, unit = "cm") %>%
  flextable::width(j = 7, width = 2.2, unit = "cm") %>%
  flextable::width(j = 8, width = 3.5, unit = "cm") %>%
  flextable::padding(padding = 5, part = "all") %>%
  flextable::add_footer_lines(paste0(
    "VND en dominio quinasa de UL97 (aa 230–620) al umbral \u22655%. ",
    "Cobertura: suma de lecturas forward + reverse en la posición (media entre muestras con la variante; mínimo observado). ",
    "Q-Score: puntuación de calidad de llamada de variante (mediana). ",
    "Balance FW/RV: número de muestras con lecturas en ambas cadenas sobre total de muestras con la variante. ",
    "W424*: codón de parada prematuro."
  )) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft_ul97_quinasa

guardar_tabla(ft_ul97_quinasa, "Complementaria_VND_UL97", apaisada = TRUE)

## ---- Figura 46. VND en UL54 ----
# ── Dominios UL54 ─────────────────────────────────────────────────────────────
dominios_ul54 <- data.frame(
  start = c(290, 696, 1150),
  end   = c(550, 990, 1243),
  label = c("Exonucleasa", "Polimerasa", "NLS + unión ppUL44"),
  fill  = c("#82E0AA", "#7FB3D3", "#E8A0D0")
)

# ── Hotspots UL54 ─────────────────────────────────────────────────────────────
hotspots_ul54 <- data.frame(
  start = c(300, 408, 498, 985,   598, 802,   700,   750, 781, 945),
  end   = c(325, 420, 545, 990,   600, 825,   715,   760, 790, 955),
  fill  = c(rep("#E67E22", 4),
            rep("#8E44AD", 2),
            "#E74C3C",
            rep("#C0392B", 3))
)

# ── Datos UL54 con umbral_max ─────────────────────────────────────────────────
ul54_umbrales <- all_data_txp %>%
  filter(class == "NEW", Protein == "UL54") %>%
  group_by(Mutation, Position) %>%
  summarise(
    prev_max   = max(Prevalence, na.rm = TRUE),
    umbral_max = case_when(
      max(Prevalence) >= 20 ~ 20,
      max(Prevalence) >= 15 ~ 15,
      max(Prevalence) >= 10 ~ 10,
      TRUE                  ~ 5
    ),
    .groups = "drop"
  ) %>%
  mutate(Position = as.numeric(Position))

ul54 <- new_vars_dom %>%
  filter(Protein == "UL54") %>%
  left_join(ul54_umbrales %>% select(Mutation, umbral_max), by = "Mutation")

# ── Regiones de resistencia ───────────────────────────────────────────────────
en_region_resistencia <- function(pos) {
  (pos >= 300 & pos <= 325) |
    (pos >= 408 & pos <= 420) |
    (pos >= 498 & pos <= 545) |
    (pos >= 598 & pos <= 600) |
    (pos >= 700 & pos <= 715) |
    (pos >= 750 & pos <= 760) |
    (pos >= 781 & pos <= 790) |
    (pos >= 802 & pos <= 825) |
    (pos >= 945 & pos <= 955) |
    (pos >= 985 & pos <= 990)
}

# ── Niveles fijos ─────────────────────────────────────────────────────────────
nivel <- c("5" = 0.6, "10" = 1.1, "15" = 1.6, "20" = 2.1)

df <- ul54 %>%
  mutate(
    en_resist = en_region_resistencia(as.numeric(Position)),
    y_nivel   = nivel[as.character(umbral_max)],
    y_point   = ifelse(en_resist, -y_nivel, y_nivel),
    color     = ifelse(en_resist, "#C0392B", "#2980B9"),
    label_mut = case_when(
      en_resist                         ~ Mutation,
      prev_max >= 10 | n_pacientes >= 4 ~ Mutation,
      TRUE                              ~ NA_character_
    )
  )

# ── Etiquetas con agrupación por proximidad ───────────────────────────────────
df_lab_azul <- df %>%
  filter(!is.na(label_mut), !en_resist) %>%
  arrange(Position) %>%
  mutate(prox_group = cumsum(c(TRUE, diff(Position) > 40))) %>%
  group_by(prox_group) %>%
  mutate(
    within_rank = row_number(),
    y_lab       = 3.0 + (within_rank - 1) * 0.60
  ) %>%
  ungroup()

df_lab_rojo <- df %>%
  filter(!is.na(label_mut), en_resist) %>%
  arrange(Position) %>%
  mutate(
    escalon = case_when(
      Position < 817  ~ (row_number() - 1) %% 3,
      TRUE            ~ (row_number() - 1 - sum(Position < 817)) %% 3
    ),
    y_lab = -(3.0 + escalon * 0.60)
  )%>%
  ungroup()

ymax_plot <- max(df_lab_azul$y_lab, na.rm = TRUE) + 0.5
ymin_plot <- min(df_lab_rojo$y_lab, na.rm = TRUE) - 0.5

# ── Niveles de referencia ─────────────────────────────────────────────────────
niveles_ref <- data.frame(
  y     = c(0.6, 1.1, 1.6, 2.1, -0.6, -1.1, -1.6, -2.1),
  label = rep(c("\u22655%", "\u226510%", "\u226515%", "\u226520%"), 2),
  lado  = rep(c("arriba", "abajo"), each = 4)
)

# ── Plot UL54 ─────────────────────────────────────────────────────────────────
p_ul54 <- ggplot() +
  
  # Líneas de referencia — de -55 a aa 1243
  geom_segment(data = niveles_ref,
               aes(x = -55, xend = 1243, y = y, yend = y),
               color = "grey85", linewidth = 0.3, linetype = "dashed") +
  
  # Etiquetas de prevalencia
  geom_text(data = niveles_ref %>% filter(lado == "arriba"),
            aes(x = -70, y = y, label = label),
            hjust = 1, size = 4.0, color = "grey50", fontface = "bold") +
  geom_text(data = niveles_ref %>% filter(lado == "abajo"),
            aes(x = -70, y = y, label = label),
            hjust = 1, size = 4.0, color = "grey50", fontface = "bold") +
  
  # Dominios
  geom_rect(data = dominios_ul54,
            aes(xmin = start, xmax = end, ymin = 0.0, ymax = 0.40, fill = fill),
            alpha = 0.85, color = "white", linewidth = 0.25) +
  geom_text(data = dominios_ul54 %>% mutate(mid = (start + end) / 2),
            aes(x = mid, y = 0.23, label = label),
            size = 4.0, color = "grey15", fontface = "bold",
            vjust = 0.5, hjust = 0.5) +
  
  # Hotspots
  geom_rect(data = hotspots_ul54,
            aes(xmin = start, xmax = end, ymin = -0.40, ymax = 0.0, fill = fill),
            alpha = 0.9, color = "white", linewidth = 0.25) +
  
  scale_fill_identity() +
  
  # Backbone aa 1–1243 con tapas
  geom_segment(aes(x = 1,    xend = 1243, y = 0,     yend = 0),
               color = "grey10", linewidth = 1.0) +
  geom_segment(aes(x = 1,    xend = 1,    y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  geom_segment(aes(x = 1243, xend = 1243, y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  
  # Etiquetas extremos backbone
  annotate("text", x = 1,    y = -0.28, label = "1",    size = 4.0,
           color = "grey10", vjust = 1) +
  annotate("text", x = 1243, y = -0.28, label = "1243", size = 4.0,
           color = "grey10", vjust = 1) +
  
  # Stems
  geom_segment(data = df,
               aes(x = Position, xend = Position,
                   y = ifelse(en_resist, -0.42, 0.42),
                   yend = y_point, color = color),
               linewidth = 0.45, alpha = 0.8) +
  
  # Puntos
  geom_point(data = df,
             aes(x = Position, y = y_point, color = color, size = prev_max),
             alpha = 0.88) +
  scale_color_identity() +
  scale_size_continuous(range = c(1.2, 4.5), guide = "none") +
  
  # Guías azules
  geom_segment(data = df_lab_azul,
               aes(x = Position, xend = Position,
                   y = y_point + 0.08, yend = y_lab - 0.08),
               color = "#2980B9", linewidth = 0.3, alpha = 0.4,
               linetype = "dotted") +
  
  # Guías rojas
  geom_segment(data = df_lab_rojo,
               aes(x = Position, xend = Position,
                   y = y_point - 0.08, yend = y_lab + 0.08),
               color = "#C0392B", linewidth = 0.3, alpha = 0.4,
               linetype = "dotted") +
  
  # Etiquetas azules
  geom_text(data = df_lab_azul,
            aes(x = Position, y = y_lab, label = label_mut),
            color = "#2980B9", size = 4.0, fontface = "bold",
            hjust = 0.5, vjust = 0) +
  
  # Etiquetas rojas
  geom_text(data = df_lab_rojo,
            aes(x = Position, y = y_lab, label = label_mut),
            color = "#C0392B", size = 4.0, fontface = "bold",
            hjust = 0.5, vjust = 1) +
  
  scale_x_continuous(
    limits = c(-130, 1380),
    breaks = NULL,
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(limits = c(ymin_plot, ymax_plot)) +
  labs(title = NULL, x = "", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(t = 8, r = 12, b = 8, l = 50)
  )

# ── Leyenda hotspots ──────────────────────────────────────────────────────────
leyenda_hotspots_ul54 <- ggplot(data.frame(
  x     = c(1.0, 4.2, 7.4, 10.6),
  y     = 1,
  label = c("GCV/CDV", "GCV/CDV/FOS", "FOS > GCV", "FOS"),
  color = c("#E67E22", "#8E44AD", "#E74C3C", "#C0392B")
)) +
  geom_point(aes(x, y, color = color), size = 4, shape = 15) +
  geom_text(aes(x + 0.18, y, label = label),
            hjust = 0, size = 4.0, color = "grey20") +
  scale_color_identity() +
  scale_x_continuous(limits = c(0.5, 13)) +
  scale_y_continuous(limits = c(0.5, 1.5)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(t = 0, b = 6, l = 50)
  )

# ── Composición final ─────────────────────────────────────────────────────────
figura_ul54 <- p_ul54 / leyenda_hotspots_ul54 +
  plot_layout(heights = c(10, 1)) +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  )

figura_ul54

# [guardado antiguo] ggsave(
# [guardado antiguo]   "~/Desktop/Figuras_tesis/figura_lollipop_UL54.png",
# [guardado antiguo]   figura_ul54,
# [guardado antiguo]   width = 34, height = 18, units = "cm", dpi = 400
# [guardado antiguo] )

guardar_fig(figura_ul54, "Figura_46", ancho = 34, alto = 18, dpi = 400)

## ---- Complementaria. VND de UL54 en regiones funcionales ----
# ── Funciones auxiliares UL54 (definir ANTES del bloque [5f]) ─────────────────

get_hotspot <- function(pos) {
  dplyr::case_when(
    (pos >= 300 & pos <= 325) |
      (pos >= 408 & pos <= 420) |
      (pos >= 498 & pos <= 545) |
      (pos >= 985 & pos <= 990)  ~ "GCV/CDV",
    (pos >= 598 & pos <= 600) |
      (pos >= 802 & pos <= 825)  ~ "GCV/CDV/FOS",
    (pos >= 700 & pos <= 715)  ~ "FOS > GCV",
    (pos >= 750 & pos <= 760) |
      (pos >= 781 & pos <= 790) |
      (pos >= 945 & pos <= 955)  ~ "FOS",
    TRUE                       ~ "No"
  )
}

get_dominio <- function(pos) {
  dplyr::case_when(
    pos >= 290  & pos <= 550  ~ "Exonucleasa",
    pos >= 696  & pos <= 990  ~ "Polimerasa",
    pos >= 1150 & pos <= 1243 ~ "NLS + unión ppUL44",
    TRUE                      ~ "Fuera de dominio"
  )
  
}

ul54_datos <- all_data_txp %>%
  filter(class == "NEW", Protein == "UL54", Prevalence >= 5) %>%
  mutate(Position = as.numeric(Position)) %>%
  group_by(Mutation, Position) %>%
  summarise(
    Paciente     = paste(sort(unique(RTP)), collapse = ", "),
    N_muestras   = n_distinct(Npet),
    N_pacientes  = n_distinct(NHC),
    Prev_min     = round(min(Prevalence), 1),
    Prev_max     = round(max(Prevalence), 1),
    cob_media    = round(mean(cobertura_total, na.rm = TRUE), 0),
    cob_min      = round(min(cobertura_total,  na.rm = TRUE), 0),
    q_med        = round(median(q_score, na.rm = TRUE), 0),
    n_balance_ok = sum(balance_ok, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Prevalencia = ifelse(
      Prev_min == Prev_max,
      paste0(format(Prev_max, decimal.mark = ","), "%"),
      paste0(format(Prev_min, decimal.mark = ","), "–",
             format(Prev_max, decimal.mark = ","), "%")
    ),
    Cobertura = paste0(format(cob_media, big.mark = "."), "×",
                       " (mín. ", format(cob_min, big.mark = "."), "×)"),
    QScore    = as.character(q_med),
    Balance   = paste0(n_balance_ok, "/", N_muestras),
    Hotspot   = get_hotspot(Position)
  ) %>%
  filter(
    Hotspot != "No" |
      (get_dominio(Position) %in% c("Exonucleasa", "Polimerasa", "NLS + unión ppUL44") &
         Prev_max >= 10)
  ) %>%
  arrange(Position) %>%
  select(Mutation, Paciente, N_muestras, Prevalencia,
         Cobertura, QScore, Balance, Hotspot)

names(ul54_datos) <- c(
  "Mutación", "Paciente", "N\nmuestras", "Prevalencia",
  "Cobertura\n(media, mín.)", "Q-Score\n(mediana)",
  "Balance FW/RV\n(n muestras)", "Hotspot"
)

filas_hotspot <- which(ul54_datos$Hotspot != "No")
filas_no      <- which(ul54_datos$Hotspot == "No")
filas_sanger  <- which(
  as.numeric(sub("%.*", "", sub(",", ".", ul54_datos$Prevalencia))) >= 20 |
    grepl("^[2-9][0-9]", ul54_datos$Prevalencia)
)

ft_ul54 <- flextable::flextable(ul54_datos) %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = c(1, 2, 8), align = "left", part = "body") %>%
  flextable::valign(valign = "center", part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1, part = "body") %>%
  flextable::bg(i = filas_hotspot, bg = "#FEF9E7", part = "body") %>%
  flextable::bg(i = filas_no,      bg = "#F0F4F8", part = "body") %>%
  flextable::bold(i = filas_sanger, part = "body") %>%
  flextable::bg(i = filas_sanger,   bg = "#FDECEA", part = "body") %>%
  flextable::bg(part = "header", bg = "#2E5FA3") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::border_inner_v(border = b_head, part = "header") %>%
  flextable::border_outer(border = b_azul, part = "all") %>%
  flextable::border_inner_h(border = b_gris, part = "body") %>%
  flextable::border_inner_v(border = b_gris, part = "body") %>%
  flextable::width(j = 1, width = 2.0, unit = "cm") %>%
  flextable::width(j = 2, width = 3.8, unit = "cm") %>%
  flextable::width(j = 3, width = 1.4, unit = "cm") %>%
  flextable::width(j = 4, width = 2.2, unit = "cm") %>%
  flextable::width(j = 5, width = 2.8, unit = "cm") %>%
  flextable::width(j = 6, width = 1.8, unit = "cm") %>%
  flextable::width(j = 7, width = 2.2, unit = "cm") %>%
  flextable::width(j = 8, width = 3.8, unit = "cm") %>%
  flextable::padding(padding = 5, part = "all") %>%
  flextable::add_footer_lines(paste0(
    "VND en UL54 seleccionadas al umbral \u22655% (en hotspot o en dominio funcional con prevalencia \u226515% o >3 pacientes). ",
    "Cobertura: lecturas totales FW+RV (media y mínimo entre muestras con la variante). ",
    "Q-Score: mediana de la puntuación de calidad de llamada. ",
    "Balance FW/RV: muestras con soporte en ambas cadenas. ",
    "Fondo amarillo: en hotspot de resistencia. Fondo azul: fuera de hotspot, incluida por dominio/prevalencia. ",
    "Fondo rojo: prevalencia \u226520% (equivalente Sanger). ",
    "Hotspot según regiones documentadas para GCV/CDV/FOS en UL54."
  )) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft_ul54

guardar_tabla(ft_ul54, "Complementaria_VND_UL54", apaisada = TRUE)

## ---- Figura 47. VND en UL56 ----
dominios_ul56 <- data.frame(
  start = c(4,   117, 280, 514, 671, 757),
  end   = c(20,  130, 310, 522, 680, 764),
  label = c("NLS", "endonucleasa", "Dedo de zinc", "Interacción pUL89", "Estab. 1", "Estab. 2"),
  fill  = c("cadetblue3", "cadetblue3", "cadetblue3", "#F0A500", "#F0A500", "#F0A500")
)

hotspots_ul56 <- data.frame(
  start = c(20,  230, 392),
  end   = c(30,  370, 398),
  label = c("",  "", ""),
  fill  = c("#B03060", "#B03060", "#B03060")
)

ul56_umbrales <- all_data_txp %>%
  filter(class == "NEW", Protein == "UL56") %>%
  group_by(Mutation, Position) %>%
  summarise(
    prev_max   = max(Prevalence, na.rm = TRUE),
    umbral_max = case_when(
      max(Prevalence) >= 20 ~ 20,
      max(Prevalence) >= 15 ~ 15,
      max(Prevalence) >= 10 ~ 10,
      TRUE                  ~ 5
    ),
    .groups = "drop"
  ) %>%
  mutate(Position = as.numeric(Position))

ul56 <- new_vars_dom %>%
  filter(Protein == "UL56") %>%
  left_join(ul56_umbrales %>% select(Mutation, umbral_max), by = "Mutation")

nivel <- c("5" = 0.6, "10" = 1.1, "15" = 1.6, "20" = 2.1)
paso  <- 0.4

df <- ul56 %>%
  mutate(
    y_nivel   = nivel[as.character(umbral_max)],
    en_lmv    = Position >= 231 & Position <= 412,
    y_point   = ifelse(en_lmv, -y_nivel, y_nivel),
    color     = ifelse(en_lmv, "#C0392B", "#2980B9"),
    label_mut = Mutation
  )

df_lab_azul_base <- df %>%
  filter(!is.na(label_mut), !en_lmv) %>%
  arrange(Position) %>%
  mutate(row_idx = row_number())

df_lab_azul <- df_lab_azul_base %>%
  mutate(
    y_lab_raw = 2.5 + (row_idx - 1) * paso,
    y_lab_raw = case_when(
      Position <= 30 ~ y_lab_raw + 0.8,
      TRUE           ~ y_lab_raw
    ),
    y_lab = case_when(
      Position >= 595 ~ y_lab_raw - 1.6,
      TRUE            ~ y_lab_raw
    )
  )

df_lab_rojo <- df %>%
  filter(!is.na(label_mut), en_lmv) %>%
  arrange(Position) %>%
  mutate(y_lab = -(0.9 + (row_number() - 1) * paso))

ymax_plot <- max(df_lab_azul$y_lab, na.rm = TRUE) + 0.5

niveles_ref <- data.frame(
  y     = c(0.6, 1.1, 1.6, 2.1, -0.6, -1.1, -1.6, -2.1),
  label = rep(c("\u22655%", "\u226510%", "\u226515%", "\u226520%"), 2),
  lado  = rep(c("arriba", "abajo"), each = 4)
)

p_ul56 <- ggplot() +
  geom_segment(data = niveles_ref,
               aes(x = -42, xend = 850, y = y, yend = y),
               color = "grey85", linewidth = 0.3, linetype = "dashed") +
  # Etiquetas prevalencia — más cerca del aa 1 (x = -55, límite -75)
  geom_text(data = niveles_ref %>% filter(lado == "arriba"),
            aes(x = -55, y = y, label = label),
            hjust = 1, size = 3.5, color = "grey50", fontface = "bold") +
  geom_text(data = niveles_ref %>% filter(lado == "abajo"),
            aes(x = -55, y = y, label = label),
            hjust = 1, size = 3.5, color = "grey50", fontface = "bold") +
  # Dominios
  geom_rect(data = dominios_ul56,
            aes(xmin = start, xmax = end, ymin = 0.0, ymax = 0.40, fill = fill),
            alpha = 0.9, color = "white", linewidth = 0.25) +
  geom_text(data = dominios_ul56 %>% mutate(mid = (start + end) / 2),
            aes(x = mid, y = 0.20, label = label),
            size = 3.5, color = "black", fontface = "bold",
            vjust = 0.5, hjust = 0.5, lineheight = 0.85) +
  # Hotspots
  geom_rect(data = hotspots_ul56,
            aes(xmin = start, xmax = end, ymin = -0.40, ymax = 0.0, fill = fill),
            alpha = 0.9, color = "white", linewidth = 0.25) +
  scale_fill_identity() +
  # Backbone aa 1–850
  geom_segment(aes(x = 1, xend = 850, y = 0, yend = 0),
               color = "grey10", linewidth = 1.0) +
  # Tapas extremos (palitos)
  geom_segment(aes(x = 1,   xend = 1,   y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  geom_segment(aes(x = 850, xend = 850, y = -0.12, yend = 0.12),
               color = "grey10", linewidth = 1.2) +
  # Stems
  geom_segment(data = df,
               aes(x = Position, xend = Position,
                   y = ifelse(en_lmv, -0.42, 0.42),
                   yend = y_point, color = color),
               linewidth = 0.45, alpha = 0.8) +
  # Puntos
  geom_point(data = df,
             aes(x = Position, y = y_point, color = color, size = prev_max),
             alpha = 0.88) +
  scale_color_identity() +
  scale_size_continuous(range = c(1.2, 4.5), guide = "none") +
  # Guías azules
  geom_segment(data = df_lab_azul,
               aes(x = Position, xend = Position,
                   y = y_point + 0.08, yend = y_lab - 0.08),
               color = "#2980B9", linewidth = 0.3, alpha = 0.4,
               linetype = "dotted") +
  # Guías rojas
  geom_segment(data = df_lab_rojo,
               aes(x = Position, xend = Position,
                   y = y_point - 0.08, yend = y_lab + 0.08),
               color = "#C0392B", linewidth = 0.3, alpha = 0.4,
               linetype = "dotted") +
  # Etiquetas azules
  geom_text(data = df_lab_azul,
            aes(x = Position, y = y_lab, label = label_mut),
            color = "#2980B9", size = 3.5, fontface = "bold",
            hjust = 0.5, vjust = 0) +
  # Etiquetas rojas
  geom_text(data = df_lab_rojo,
            aes(x = Position, y = y_lab, label = label_mut),
            color = "#C0392B", size = 3.5, fontface = "bold",
            hjust = 0.5, vjust = 1) +
  annotate("text", x = 1,   y = -0.28, label = "1",   size = 3.5, color = "grey10", vjust = 1) +
  annotate("text", x = 850, y = -0.28, label = "850", size = 3.5, color = "grey10", vjust = 1) +
  scale_x_continuous(
    limits = c(-100, 870),
    breaks = NULL,
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(limits = c(-2.5, ymax_plot)) +
  # Sin título de panel — el nombre "UL56" va en el pie de figura
  labs(title = NULL, x = "", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(t = 8, r = 12, b = 8, l = 35)
  )

leyenda_hotspots_ul56 <- ggplot(data.frame(
  x     = 1.5,
  y     = 1,
  label = "LMV-R",
  color = "#B03060"
)) +
  geom_point(aes(x, y, color = color), size = 4, shape = 15) +
  geom_text(aes(x + 0.15, y, label = label),
            hjust = 0, size = 3.5, color = "grey20") +
  scale_color_identity() +
  scale_x_continuous(limits = c(0.8, 4)) +
  scale_y_continuous(limits = c(0.5, 1.5)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(t = 0, b = 6, l = 35)
  )

figura_ul56 <- p_ul56 / leyenda_hotspots_ul56 +
  plot_layout(heights = c(10, 1))

figura_ul56

# [guardado antiguo] ggsave(
# [guardado antiguo]   "~/Desktop/Figuras_tesis/figura_lollipop_UL56.png",
# [guardado antiguo]   figura_ul56,
# [guardado antiguo]   width = 28, height = 15, units = "cm", dpi = 400,
# [guardado antiguo]   bg = "white"
# [guardado antiguo] )

guardar_fig(figura_ul56, "Figura_47", ancho = 28, alto = 15, dpi = 400)

## ---- Complementaria. VND de UL56 ----
ul56_dom_datos <- all_data_txp %>%
  filter(class == "NEW", Protein == "UL56", Prevalence >= 5) %>%
  mutate(Position = as.numeric(Position)) %>%
  group_by(Mutation, Position) %>%
  summarise(
    Paciente     = paste(sort(unique(RTP)), collapse = ", "),
    N_muestras   = n_distinct(Npet),
    Prev_min     = round(min(Prevalence), 1),
    Prev_max     = round(max(Prevalence), 1),
    cob_mediana = round(median(cobertura_total, na.rm = TRUE), 0),
    cob_p25     = round(quantile(cobertura_total, 0.25, na.rm = TRUE), 0),
    cob_p75     = round(quantile(cobertura_total, 0.75, na.rm = TRUE), 0),
    cob_min     = round(min(cobertura_total, na.rm = TRUE), 0),
    q_med       = round(median(q_score, na.rm = TRUE), 0),
    q_p25       = round(quantile(q_score, 0.25, na.rm = TRUE), 0),
    q_p75       = round(quantile(q_score, 0.75, na.rm = TRUE), 0),
    n_balance_ok = sum(balance_ok, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Prevalencia = ifelse(
      Prev_min == Prev_max,
      paste0(format(Prev_max, decimal.mark = ","), "%"),
      paste0(format(Prev_min, decimal.mark = ","), "–",
             format(Prev_max, decimal.mark = ","), "%")
    ),
    Cobertura = ifelse(
      N_muestras == 1,
      paste0(format(cob_mediana, big.mark = "."), "×"),
      paste0(
        format(cob_mediana, big.mark = "."), "× (",
        format(cob_p25, big.mark = "."), "–",
        format(cob_p75, big.mark = "."), "×; mín. ",
        format(cob_min, big.mark = "."), "×)"
      )
    ),
    QScore = ifelse(
      N_muestras == 1,
      as.character(q_med),
      paste0(q_med, " (", q_p25, "–", q_p75, ")")
    ),
    Balance = paste0(n_balance_ok, "/", N_muestras)
  ) %>%
  arrange(Position)

hotspot_map_ul56 <- tibble::tribble(
  ~Mutation, ~Hotspot,
  "Y281H",   "LMV-R (aa 230–370)",
  "N330S",   "LMV-R (aa 230–370)",
  "D480G",   "No",
  "I508V",   "No",
  "R513C",   "No",
  "Y595C",   "No",
  "N599S",   "No",
  "L612P",   "No",
  "M697I",   "No"
)

ul56_dom_datos <- ul56_dom_datos %>%
  left_join(hotspot_map_ul56, by = "Mutation") %>%
  mutate(Hotspot = replace_na(Hotspot, "No")) %>%
  select(Mutation, Paciente, N_muestras, Prevalencia,
         Cobertura, QScore, Balance, Hotspot)

names(ul56_dom_datos) <- c(
  "Mutación", "Paciente", "N\nmuestras", "Prevalencia",
  "Cobertura\n(mediana, RIC; mín.)", "Q-Score\n(mediana, RIC)",
  "Balance FW/RV\n(n muestras)", "Hotspot"
)

filas_hotspot <- which(ul56_dom_datos$Hotspot != "No")
filas_no      <- which(ul56_dom_datos$Hotspot == "No")
fila_r513c    <- which(ul56_dom_datos[[1]] == "R513C")

ft_ul56_dom <- flextable::flextable(ul56_dom_datos) %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = c(1, 2, 8), align = "left", part = "body") %>%
  flextable::valign(valign = "center", part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1, part = "body") %>%
  flextable::bold(i = fila_r513c, part = "body") %>%
  flextable::bg(i = fila_r513c, bg = "#FDECEA", part = "body") %>%
  flextable::bg(i = filas_hotspot[!filas_hotspot %in% fila_r513c],
                bg = "#FEF9E7", part = "body") %>%
  flextable::bg(i = filas_no[!filas_no %in% fila_r513c],
                bg = "#F0F4F8", part = "body") %>%
  flextable::bg(part = "header", bg = "#2E5FA3") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::border_inner_v(border = b_head, part = "header") %>%
  flextable::border_outer(border = b_azul, part = "all") %>%
  flextable::border_inner_h(border = b_gris, part = "body") %>%
  flextable::border_inner_v(border = b_gris, part = "body") %>%
  flextable::width(j = 1, width = 2.0, unit = "cm") %>%
  flextable::width(j = 2, width = 3.8, unit = "cm") %>%
  flextable::width(j = 3, width = 1.4, unit = "cm") %>%
  flextable::width(j = 4, width = 2.2, unit = "cm") %>%
  flextable::width(j = 5, width = 2.8, unit = "cm") %>%
  flextable::width(j = 6, width = 1.8, unit = "cm") %>%
  flextable::width(j = 7, width = 2.2, unit = "cm") %>%
  flextable::width(j = 8, width = 3.5, unit = "cm") %>%
  flextable::padding(padding = 5, part = "all") %>%
  flextable::add_footer_lines(paste0(
    "VND en UL56 al umbral \u22655%. ",
    "Cobertura: lecturas totales FW+RV en la posición (media y mínimo entre muestras con la variante). ",
    "Q-Score: mediana de la puntuación de calidad de llamada. ",
    "Balance FW/RV: muestras con soporte en ambas cadenas. ",
    "R513C: única variante con prevalencia \u226520% (equivalente Sanger). ",
    "Hotspot LMV-R según CHARMD (Tilloy et al., 2024)."
  )) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft_ul56_dom

guardar_tabla(ft_ul56_dom, "Complementaria_VND_UL56", apaisada = TRUE)

## ---- Complementaria. Tabla combinada de VND (todos los genes) ----
# ── Anotación de región de interés (prioridad: hotspot > dominio) ─────────────
annot_ul97 <- function(pos) dplyr::case_when(
  pos >= 335 & pos <= 370  ~ "Hotspot GCV/MBV-R (aa 335\u2013370)",
  pos >= 397 & pos <= 480  ~ "Hotspot MBV-R (aa 397\u2013480)",
  pos >= 519 & pos <= 521  ~ "Hotspot GCV-R (aa 519\u2013521)",
  pos >= 590 & pos <= 607  ~ "Hotspot GCV-R (aa 590\u2013607)",
  pos >= 230 & pos <= 620  ~ "Dominio quinasa",
  TRUE                     ~ ""
)
annot_ul54 <- function(pos) dplyr::case_when(
  pos >= 300 & pos <= 325  ~ "Hotspot GCV/CDV (aa 300\u2013325)",
  pos >= 408 & pos <= 420  ~ "Hotspot GCV/CDV (aa 408\u2013420)",
  pos >= 498 & pos <= 545  ~ "Hotspot GCV/CDV (aa 498\u2013545)",
  pos >= 985 & pos <= 990  ~ "Hotspot GCV/CDV (aa 985\u2013990)",
  pos >= 598 & pos <= 600  ~ "Hotspot GCV/CDV/FOS (aa 598\u2013600)",
  pos >= 802 & pos <= 825  ~ "Hotspot GCV/CDV/FOS (aa 802\u2013825)",
  pos >= 700 & pos <= 715  ~ "Hotspot FOS/GCV (aa 700\u2013715)",
  pos >= 750 & pos <= 760  ~ "Hotspot FOS (aa 750\u2013760)",
  pos >= 781 & pos <= 790  ~ "Hotspot FOS (aa 781\u2013790)",
  pos >= 945 & pos <= 955  ~ "Hotspot FOS (aa 945\u2013955)",
  pos >= 290 & pos <= 550  ~ "Dominio exonucleasa",
  pos >= 696 & pos <= 990  ~ "Dominio polimerasa",
  pos >= 1150 & pos <= 1243 ~ "Regi\u00f3n NLS / uni\u00f3n ppUL44",
  TRUE                     ~ ""
)
annot_ul56 <- function(pos) dplyr::case_when(
  pos >= 230 & pos <= 370  ~ "Hotspot LMV-R (aa 230\u2013370)",
  TRUE                     ~ ""
)

# ── Formateadores (coma decimal, punto de millar) ─────────────────────────────
pc  <- function(x) formatC(x, format = "f", digits = 1, decimal.mark = ",")
mil <- function(x) formatC(x, format = "f", digits = 0, big.mark = ".", decimal.mark = ",")

# ── Denominadores (de tu objeto n_total_gen, bloque [3b]) ─────────────────────
n_total_gen <- tibble(
  Protein         = c("UL97", "UL54", "UL56"),
  n_tot_muestras  = c(26, 24, 11),
  n_tot_pacientes = c(11, 11, 3)
)

# ── Construcción de la tabla ──────────────────────────────────────────────────
tabla_vnd <- all_data_txp %>%
  filter(class == "NEW", Prevalence >= 5) %>%
  mutate(Position = as.numeric(Position)) %>%
  group_by(Protein, Mutation, Position) %>%
  summarise(
    n_muestras  = n_distinct(Npet),
    n_pacientes = n_distinct(NHC),
    prev_media  = round(mean(Prevalence, na.rm = TRUE), 1),
    prev_de     = round(sd(Prevalence,   na.rm = TRUE), 1),
    cob_mediana = round(median(cobertura_total, na.rm = TRUE), 0),
    cob_p25     = round(quantile(cobertura_total, 0.25, na.rm = TRUE), 0),
    cob_p75     = round(quantile(cobertura_total, 0.75, na.rm = TRUE), 0),
    .groups = "drop"
  ) %>%
  left_join(n_total_gen, by = "Protein") %>%
  mutate(
    pct_muestras  = round(n_muestras  / n_tot_muestras  * 100, 1),
    pct_pacientes = round(n_pacientes / n_tot_pacientes * 100, 1),
    col_muestras  = paste0(n_muestras,  " (", pc(pct_muestras),  "%)"),
    col_pacientes = paste0(n_pacientes, " (", pc(pct_pacientes), "%)"),
    col_prev = if_else(
      is.na(prev_de),
      paste0(pc(prev_media), " (n=1)"),
      paste0(pc(prev_media), " \u00b1 ", pc(prev_de))
    ),
    col_cob = if_else(
      n_muestras == 1,
      paste0(mil(cob_mediana), "\u00d7"),
      paste0(mil(cob_mediana), "\u00d7 (", mil(cob_p25), "\u2013", mil(cob_p75), ")")
    ),
    dom_hot = case_when(
      Protein == "UL97" ~ annot_ul97(Position),
      Protein == "UL54" ~ annot_ul54(Position),
      Protein == "UL56" ~ annot_ul56(Position),
      TRUE              ~ ""
    ),
    Protein = factor(Protein, levels = c("UL97", "UL54", "UL56"))
  ) %>%
  arrange(Protein, Position)

# ── Diagnóstico (verifica que son 101) ────────────────────────────────────────
cat("Total VND (todas):", nrow(tabla_vnd), "\n")
tabla_vnd %>% count(Protein) %>% print()
cat("Con Dominio/Hotspot asignado:", sum(tabla_vnd$dom_hot != ""), "\n")

# ── Dataframe para flextable ──────────────────────────────────────────────────
vnd_ft <- tabla_vnd %>%
  transmute(
    Gen = as.character(Protein), Variante = Mutation,
    col_muestras, col_pacientes, col_prev, col_cob, dom_hot
  )
names(vnd_ft) <- c(
  "Gen", "Variante", "Muestras, n (%)", "Pacientes, n (%)",
  "Frecuencia al\u00e9lica\n(media \u00b1 DE)",
  "Cobertura de posici\u00f3n\nMediana (RIC)",
  "Dominio / Hotspot"
)

# ── Bloques por gen (para sombreado y separadores) ───────────────────────────
n_ul97 <- sum(tabla_vnd$Protein == "UL97")
n_ul54 <- sum(tabla_vnd$Protein == "UL54")
n_ul56 <- sum(tabla_vnd$Protein == "UL56")
filas_ul97 <- seq_len(n_ul97)
filas_ul54 <- n_ul97 + seq_len(n_ul54)
filas_ul56 <- n_ul97 + n_ul54 + seq_len(n_ul56)
filas_sep  <- c(n_ul97, n_ul97 + n_ul54)

b_azul <- officer::fp_border(color = "#2E5FA3", width = 1.5)
b_sep  <- officer::fp_border(color = "#2E5FA3", width = 1.2)
b_gris <- officer::fp_border(color = "#CCCCCC", width = 0.5)
b_head <- officer::fp_border(color = "#FFFFFF", width = 0.8)

# ── Flextable ─────────────────────────────────────────────────────────────────
ft_vnd <- flextable::flextable(vnd_ft) %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = c(1, 2, 7), align = "left", part = "body") %>%
  flextable::valign(valign = "center", part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1:2, part = "body") %>%
  flextable::bg(part = "header", bg = "#2E5FA3") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::border_inner_v(border = b_head, part = "header") %>%
  flextable::bg(i = filas_ul97, bg = "#EEF4FB", part = "body") %>%
  flextable::bg(i = filas_ul54, bg = "#FFFFFF", part = "body") %>%
  flextable::bg(i = filas_ul56, bg = "#EEF4FB", part = "body") %>%
  flextable::border_outer(border = b_azul, part = "all") %>%
  flextable::border_inner_h(border = b_gris, part = "body") %>%
  flextable::border_inner_v(border = b_gris, part = "body") %>%
  flextable::hline(i = filas_sep, border = b_sep, part = "body") %>%
  flextable::width(j = 1, width = 1.6, unit = "cm") %>%
  flextable::width(j = 2, width = 1.9, unit = "cm") %>%
  flextable::width(j = 3, width = 2.3, unit = "cm") %>%
  flextable::width(j = 4, width = 2.3, unit = "cm") %>%
  flextable::width(j = 5, width = 2.7, unit = "cm") %>%
  flextable::width(j = 6, width = 3.3, unit = "cm") %>%
  flextable::width(j = 7, width = 4.6, unit = "cm") %>%
  flextable::padding(padding = 4, part = "all") %>%
  flextable::add_footer_lines(paste0(
    "VND: variante no descrita en CHARMD (Tilloy et al., 2024), umbral \u22655%. ",
    "Porcentajes sobre muestras con secuenciaci\u00f3n v\u00e1lida por gen (UL97: 26; UL54: 24; UL56: 11) ",
    "y sobre pacientes con dato por gen (UL97: 11; UL54: 11; UL56: 3). ",
    "Frecuencia al\u00e9lica: media \u00b1 DE; en variantes con una sola muestra se indica el valor \u00fanico (n=1). ",
    "Cobertura: lecturas FW+RV en la posici\u00f3n, mediana (RIC). ",
    "Dominio / Hotspot: detalle solo en variantes localizadas en una regi\u00f3n de inter\u00e9s; en blanco si quedan fuera. ",
    "El asterisco (*) indica cod\u00f3n de parada."
  )) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::font(fontname = "Arial", part = "footer") %>%
  flextable::color(color = "#555555", part = "footer") %>%
  flextable::italic(part = "footer") %>%
  flextable::align(align = "left", part = "footer")

ft_vnd

guardar_tabla(ft_vnd, "Complementaria_VND_todas", apaisada = TRUE)

## ---- Anexo II. VND con región funcional y evidencia técnica ----
# El bloque guarda el Word apaisado (Anexo_II.docx) en la carpeta de salidas.
## Requiere all_data_txp (con class y columnas DeepChek: cobertura_total, mut_fw, mut_rv, q_score).
## Una fila por VND; frecuencia y cobertura resumen todas las muestras en que se detectó.


fmt_n <- function(x) format(round(x), big.mark = ".", decimal.mark = ",", trim = TRUE)
fmt_1 <- function(x) sub("\\.", ",", sprintf("%.1f", x))

# ── Regiones funcionales y hotspots (editar límites según Material y métodos) ──
## Se asigna la primera fila que contenga la posición: los hotspots van antes que los dominios.
regiones_funcionales <- tribble(
  ~Protein, ~ini,  ~fin,  ~etiqueta,
  # UL97
  "UL97",   335,   370,   "GCV/MBV-R (aa 335\u2013370)",
  "UL97",   397,   480,   "MBV-R (aa 397\u2013480)",
  "UL97",   590,   607,   "GCV-R (aa 590\u2013607)",
  "UL97",   335,   665,   "GCV-R (aa 335\u2013665)",
  "UL97",   337,   651,   "pUL97-KD",
  # UL54
  "UL54",   300,   325,   "GCV/CDV-R (aa 300\u2013325)",
  "UL54",   498,   545,   "GCV/CDV-R (aa 498\u2013545)",
  "UL54",   700,   715,   "FOS/GCV-R (aa 700\u2013715)",
  "UL54",   750,   760,   "FOS-R (aa 750\u2013760)",
  "UL54",   802,   825,   "GCV/CDV/FOS-R (aa 802\u2013825)",
  "UL54",   985,   990,   "GCV/CDV-R (aa 985\u2013990)",
  "UL54",   300,   545,   "pUL54-Exo",
  "UL54",   696,   988,   "pUL54-Pol",
  "UL54",  1153,  1242,   "pUL54-Cter/NLS",
  # UL56
  "UL56",   230,   370,   "LMV-R (aa 230\u2013370)"
)
asignar_region <- function(prot, pos) {
  r <- regiones_funcionales %>% filter(Protein == prot, pos >= ini, pos <= fin)
  if (nrow(r) == 0) "" else r$etiqueta[1]
}

# ── Denominadores por gen ────────────────────────────────────────────────────
denom <- all_data_txp %>% group_by(Protein) %>%
  summarise(N_m = n_distinct(Npet), N_p = n_distinct(RTP), .groups = "drop")

# ── Tabla ────────────────────────────────────────────────────────────────────
anexo2 <- all_data_txp %>%
  filter(class == "NEW", Prevalence >= 5) %>%
  mutate(mut_tot = mut_fw + mut_rv,
         balance = pmin(mut_fw, mut_rv) / mut_tot) %>%
  group_by(Protein, Mutation) %>%
  summarise(
    n_m       = n_distinct(Npet),
    n_p       = n_distinct(RTP),
    pacientes = paste(unique(RTP[order(as.numeric(sub("RTP-", "", RTP)))]), collapse = ", "),
    prev_mean = mean(Prevalence), prev_sd = sd(Prevalence),
    cob_med   = median(cobertura_total),
    cob_q1    = quantile(cobertura_total, .25), cob_q3 = quantile(cobertura_total, .75),
    lect_min  = min(mut_tot),
    fw_rv_min = paste0(mut_fw[which.min(mut_tot)], "/", mut_rv[which.min(mut_tot)]),
    bal_min   = min(balance),
    q_min     = min(q_score),
    .groups = "drop") %>%
  left_join(denom, by = "Protein") %>%
  mutate(
    Posicion  = as.numeric(gsub("[^0-9]", "", Mutation)),
    Region    = mapply(asignar_region, Protein, Posicion),
    Gen       = factor(Protein, levels = c("UL97", "UL54", "UL56"),
                       labels = c(sprintf("pUL97 (%d muestras, %d pacientes)", denom$N_m[denom$Protein == "UL97"], denom$N_p[denom$Protein == "UL97"]),
                                  sprintf("pUL54 (%d muestras, %d pacientes)", denom$N_m[denom$Protein == "UL54"], denom$N_p[denom$Protein == "UL54"]),
                                  sprintf("pUL56 (%d muestras, %d pacientes)", denom$N_m[denom$Protein == "UL56"], denom$N_p[denom$Protein == "UL56"]))),
    Muestras   = sprintf("%d (%s)", n_m, fmt_1(100 * n_m / N_m)),
    Pacientes  = sprintf("%d (%s)", n_p, fmt_1(100 * n_p / N_p)),
    Frecuencia = ifelse(n_m == 1, fmt_1(prev_mean),
                        sprintf("%s \u00b1 %s", fmt_1(prev_mean), fmt_1(prev_sd))),
    Cobertura  = ifelse(n_m == 1, paste0(fmt_n(cob_med), " \u00d7"),
                        sprintf("%s \u00d7 (%s\u2013%s)", fmt_n(cob_med), fmt_n(cob_q1), fmt_n(cob_q3))),
    Lecturas   = sprintf("%d (%s)", lect_min, fw_rv_min),
    Hebra      = fmt_1(100 * bal_min),
    Q          = as.character(q_min)
  ) %>%
  arrange(Gen, Posicion) %>%
  select(Gen, Mutation, Muestras, Pacientes, pacientes, Frecuencia, Cobertura, Lecturas, Hebra, Q, Region)

cat("VND por gen:\n"); print(table(anexo2$Gen))

# ── Flextable con fila combinada por gen ─────────────────────────────────────
anexo2_grp <- as_grouped_data(anexo2, groups = "Gen")

ft_anexo2 <- as_flextable(anexo2_grp, hide_grouplabel = TRUE) %>%
  set_header_labels(
    Mutation = "Variante",
    Muestras = "Muestras, n (%)", Pacientes = "Pacientes, n (%)",
    pacientes = "Pacientes (Tabla X)",
    Frecuencia = "Frecuencia de variante (%), media \u00b1 DE",
    Cobertura = "Cobertura por posici\u00f3n, mediana (RIC)",
    Lecturas = "Lecturas mutadas, m\u00ednimo (fw/rv)",
    Hebra = "Hebra minoritaria, m\u00ednimo (%)",
    Q = "Q-score, m\u00ednimo",
    Region = "Regi\u00f3n funcional / hotspot de resistencia"
  ) %>%
  align(align = "center", part = "all") %>%
  align(j = c(1, 10), align = "left", part = "body") %>%
  bold(part = "header") %>%
  bold(i = ~ !is.na(Gen), bold = TRUE) %>%
  bg(i = ~ !is.na(Gen), bg = "#1F3B73") %>%
  color(i = ~ !is.na(Gen), color = "white") %>%
  bg(part = "header", bg = "#1F3B73") %>% color(part = "header", color = "white") %>%
  bold(j = "Mutation", bold = TRUE) %>%
  flextable::font(fontname = "Source Sans Pro", part = "all") %>%
  fontsize(size = 7.5, part = "all") %>%
  border_inner(border = fp_border(color = "grey70", width = 0.5)) %>%
  border_outer(border = fp_border(color = "grey40", width = 1)) %>%
  add_footer_lines(paste0(
    "VND: variantes no descritas previamente (no recogidas en CHARMD, versi\u00f3n actualizada a 20 de abril de 2026). ",
    "Muestras y pacientes sobre el total con resultado v\u00e1lido para cada gen. Frecuencia y cobertura resumen las muestras en ",
    "que se detect\u00f3 la variante (valor \u00fanico cuando se detect\u00f3 en una sola). Lecturas mutadas: n\u00famero m\u00ednimo de lecturas que ",
    "soportan la variante en cualquiera de las muestras, con su reparto entre hebras directa (fw) e inversa (rv). Hebra minoritaria: ",
    "proporci\u00f3n m\u00ednima de lecturas mutadas en la hebra menos representada. El filtro por Q-score del sistema experto estaba ",
    "desactivado. pUL97-KD: dominio quinasa de pUL97; Exo: dominio exonucleasa; Pol: dominio polimerasa; Cter/NLS: extremo ",
    "C-terminal y se\u00f1al de localizaci\u00f3n nuclear; -R: regi\u00f3n asociada a resistencia; GCV: ganciclovir; CDV: cidofovir; ",
    "FOS: foscarnet; MBV: maribavir; LMV: letermovir. * Cod\u00f3n de terminaci\u00f3n prematuro.")) %>%
  fontsize(size = 6.5, part = "footer") %>%
  set_table_properties(layout = "autofit", width = 1)

ft_anexo2

# ── Exportar a Word en horizontal ────────────────────────────────────────────
doc <- read_docx() %>%
  body_add_par("Anexo II. Variantes no descritas previamente identificadas mediante NGS", style = "heading 2") %>%
  body_add_flextable(ft_anexo2) %>%
  body_end_section_landscape()
print(doc, target = file.path(dir_salida, "Anexo_II.docx"))


# ---- 10.2.5. Influencia del umbral de detección ----

## ---- Figura 48. Variantes únicas según el umbral de frecuencia de variante ----
# ── Paleta ──────────────────────────────────────────────────────────────────
pal_genes <- c("UL97" = "#4A90C4", "UL54" = "#E8A33D", "UL56" = "#D9534F")

# ── Datos ───────────────────────────────────────────────────────────────────
datos_umbral <- tribble(
  ~Gen,   ~Umbral,  ~Tipo,  ~Variantes,
  "UL97", "≥ 5%",   "MRA",   3,
  "UL97", "≥ 10%",  "MRA",   3,
  "UL97", "≥ 15%",  "MRA",   2,
  "UL97", "≥ 20%",  "MRA",   2,
  "UL54", "≥ 5%",   "MRA",   1,
  "UL54", "≥ 10%",  "MRA",   1,
  "UL54", "≥ 15%",  "MRA",   1,
  "UL54", "≥ 20%",  "MRA",   1,
  "UL56", "≥ 5%",   "MRA",   0,
  "UL56", "≥ 10%",  "MRA",   0,
  "UL56", "≥ 15%",  "MRA",   0,
  "UL56", "≥ 20%",  "MRA",   0,
  "UL97", "≥ 5%",   "PN",    9,      # incluye D605E (5,7%)
  "UL97", "≥ 10%",  "PN",    8,
  "UL97", "≥ 15%",  "PN",    8,
  "UL97", "≥ 20%",  "PN",    8,
  "UL54", "≥ 5%",   "PN",   12,
  "UL54", "≥ 10%",  "PN",   11,
  "UL54", "≥ 15%",  "PN",   11,
  "UL54", "≥ 20%",  "PN",   11,
  "UL56", "≥ 5%",   "PN",    3,
  "UL56", "≥ 10%",  "PN",    3,
  "UL56", "≥ 15%",  "PN",    3,
  "UL56", "≥ 20%",  "PN",    3,
  "UL97", "≥ 5%",   "VND",  25,
  "UL97", "≥ 10%",  "VND",   2,
  "UL97", "≥ 15%",  "VND",   2,
  "UL97", "≥ 20%",  "VND",   1,
  "UL54", "≥ 5%",   "VND",  67,
  "UL54", "≥ 10%",  "VND",  17,
  "UL54", "≥ 15%",  "VND",  13,
  "UL54", "≥ 20%",  "VND",   6,
  "UL56", "≥ 5%",   "VND",   9,
  "UL56", "≥ 10%",  "VND",   1,
  "UL56", "≥ 15%",  "VND",   1,
  "UL56", "≥ 20%",  "VND",   1
) %>%
  mutate(
    Umbral = factor(Umbral, levels = c("≥ 20%", "≥ 15%", "≥ 10%", "≥ 5%")),
    Tipo   = factor(Tipo,   levels = c("MRA", "PN", "VND")),
    Gen    = factor(Gen,    levels = c("UL97", "UL54", "UL56"))
  )

# ── Figura ──────────────────────────────────────────────────────────────────
p_umbral <- ggplot(datos_umbral, aes(x = Umbral, y = Variantes, fill = Gen)) +
  geom_col(
    position  = position_dodge(width = 0.8),
    width     = 0.7,
    color     = "white",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = if_else(Variantes > 0, as.character(Variantes), ""),
        color = Gen),
    position = position_dodge(width = 0.8),
    vjust    = -0.6,
    size     = 3.9,
    fontface = "bold",
    show.legend = FALSE
  ) +
  facet_wrap(~ Tipo, scales = "free_y", nrow = 1) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = pal_genes) +
  scale_color_manual(values = pal_genes) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.22))) +
  labs(
    x    = "Umbral de frecuencia de variante",
    y    = "Número de variantes únicas",
    fill = "Gen"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 11) +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = 12),
    axis.title.x       = element_text(face = "bold", size = 11,
                                      margin = margin(t = 8)),
    axis.title.y       = element_text(face = "bold", size = 11,
                                      margin = margin(r = 8)),
    axis.text          = element_text(color = "black", size = 11),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
    panel.spacing      = unit(0.7, "cm"),
    legend.position    = "bottom",
    legend.title       = element_text(face = "bold", size = 11),
    legend.text        = element_text(size = 11),
    legend.margin      = margin(t = 0),
    plot.margin        = margin(4, 6, 4, 4)
  )

p_umbral

# [guardado antiguo] ggsave(
# [guardado antiguo]   filename = "~/Desktop/Figuras_tesis/figura_47_variantes_umbral_NGS.png",
# [guardado antiguo]   plot     = p_umbral,
# [guardado antiguo]   width    = 16,
# [guardado antiguo]   height   = 8,
# [guardado antiguo]   units    = "cm",
# [guardado antiguo]   dpi      = 600
# [guardado antiguo] )

# [Nuevo] Los recuentos de la figura están escritos a mano. Se comprueba que coinciden
# con los datos (variantes únicas por gen, clase y umbral en all_data_txp).
recuento_datos <- expand_grid(Gen = c("UL97", "UL54", "UL56"), thr = c(5, 10, 15, 20),
                              Tipo = c("MRA", "PN", "VND")) %>%
  mutate(Variantes_datos = pmap_int(list(Gen, thr, Tipo), function(g, t, cl) {
    clase <- c(MRA = "MRA", PN = "POL", VND = "NEW")[[cl]]
    n_distinct(filter(all_data_txp, Protein == g, class == clase, Prevalence >= t)$Mutation)
  }),
  Umbral = paste0("≥ ", thr, "%"))
comprobacion_fig48 <- datos_umbral %>%
  mutate(across(c(Gen, Umbral, Tipo), as.character)) %>%
  left_join(recuento_datos, by = c("Gen", "Umbral", "Tipo"))
stopifnot(all(comprobacion_fig48$Variantes == comprobacion_fig48$Variantes_datos))
cat("Figura 48: los recuentos escritos a mano coinciden con los datos.\n")

guardar_fig(p_umbral, "Figura_48", ancho = 16, alto = 8, dpi = 600)


# ---- 10.2.6. Dinámica temporal de la resistencia antiviral en RTP-2 ----

## ---- Preparación común: cargas virales, tratamientos y variantes de RTP-2 ----
# Guarda las figuras con guardar_rtp2() (png + svg) en la carpeta de salidas.
## Cambios respecto a la versión anterior:
##  - Las MRA y las VND se toman de all_data_txp (misma tabla que los anexos), no de listas
##    escritas a mano (vsi_extra). Desaparecen dia_remap, F396L, D711Y, F300L, P704S, P712H,
##    Q323R y D605E (ahora PN). T700A no forma parte de all_data_txp (no superó el filtro
##    de cobertura) y no se representa.
##  - Las muestras sin resultado de UL54 o UL56 (cruces) se detectan desde los datos.
##  - Las VND representadas se seleccionan con un criterio explícito: localización en las
##    regiones funcionales asociadas a resistencia (regiones_resist); el resto queda en el anexo.
##  - Los tratamientos y las cargas virales (df_rtp2, trat_df) siguen siendo datos clínicos
##    introducidos manualmente.


fecha_txp_rtp2 <- as.Date("2022-06-17")
nhc_rtp2       <- "NHC_RTP-2"

# ── Colores y tamaños ─────────────────────────────────────────────────────────
col_GCV <- "#3D7A5A"; col_VGCV <- "#2E86AB"; col_MBV <- "#B8860B"; col_LEF <- "#C0652B"
col_FOS <- "#9B2335"; col_LMV  <- "#7B5EA7"; col_CMVIg <- "#A8C8E0"
colores_fill  <- c("UL97" = alpha("#4472C4", 0.80), "UL54" = alpha("#E6A817", 0.80), "UL56" = alpha("#E8735A", 0.80))
colores_color <- c("UL97" = "#4472C4", "UL54" = "#C8880A", "UL56" = "#C94A30")
tam_5  <- c("\u22655%" = 2, "\u226510%" = 3, "\u226515%" = 5, "\u226520%" = 7)
tam_10 <- c("\u226510%" = 3, "\u226515%" = 5, "\u226520%" = 7)

# ── Eje temporal (índice 1..51 sobre los días con determinación de CV) ───────
niveles_dia <- c(
  "1","26","42","49","55","62","69","75","84","94","108","122",
  "159","181","187","201","214","227","244","271","279","291","304",
  "361","371","382","388","390","395","402","409","418","423","433",
  "445","458","466","479","497","531","571","630","668","697","709",
  "742","774","801","850","872","909")
n_niv          <- length(niveles_dia)
pos_map        <- tibble(dia = as.numeric(niveles_dia), pos = seq_along(niveles_dia))
interpolar_pos <- function(d) approx(pos_map$dia, pos_map$pos, xout = d, rule = 2)$y

# ── Datos clínicos RTP-2 (carga viral) ───────────────────────────────────────
df_rtp2 <- tibble(
  dia = as.numeric(niveles_dia),
  log_cv = c(0.000, 0.000,
             4.278, 4.842, 5.042, 4.389, 3.891, 3.461, 3.234, 2.876, 0.000, 2.241,
             2.740, 3.557, 3.389, 3.461, 3.042, 2.885, 3.346,
             2.699, 3.279, 3.653, 4.008,
             3.741, 4.276, 4.610, 4.521, 4.422, 4.157, 4.312,
             0.000, 3.828, 3.872, 3.638,
             3.752, 4.381, 4.362, 4.364, 4.173,
             4.203, 3.712, 3.328, 2.965, 2.845, 2.823, 2.823,
             3.215, 3.034, 2.645, 3.002, 0.000)
)

# ── Muestras NGS de RTP-2 desde los datos: día y número de orden ─────────────
muestras_rtp2 <- all_data_txp %>%
  filter(NHC == nhc_rtp2) %>%
  mutate(Freg = as.Date(Freg), dia = as.numeric(Freg - fecha_txp_rtp2)) %>%
  distinct(Npet, dia) %>% arrange(dia) %>%
  mutate(muestra = row_number())
stopifnot(all(muestras_rtp2$dia %in% as.numeric(niveles_dia)))   # los días deben estar en el eje

df_rtp2 <- df_rtp2 %>%
  left_join(muestras_rtp2 %>% select(dia, muestra), by = "dia") %>%
  mutate(secuenciada = !is.na(muestra),
         x_pos       = interpolar_pos(dia),
         tipo_punto  = ifelse(secuenciada, "Muestra secuenciada", "Muestra no secuenciada"))

# Muestras sin resultado por gen (cruces): analizadas para RTP-2 pero sin ninguna fila del gen
genes_por_muestra <- all_data_txp %>%
  filter(NHC == nhc_rtp2) %>%
  mutate(Freg = as.Date(Freg), dia = as.numeric(Freg - fecha_txp_rtp2)) %>%
  distinct(dia, Protein)
dias_sin_gen <- function(gen) setdiff(muestras_rtp2$dia, genes_por_muestra$dia[genes_por_muestra$Protein == gen])

# ── Variantes de RTP-2 desde los datos ───────────────────────────────────────
variantes_rtp2 <- all_data_txp %>%
  filter(NHC == nhc_rtp2, class %in% c("MRA", "VSI", "NEW"), Prevalence >= 5) %>%
  mutate(
    Freg  = as.Date(Freg),
    dia   = as.numeric(Freg - fecha_txp_rtp2),
    x_pos = interpolar_pos(dia),
    gen_color = factor(as.character(Protein), levels = c("UL97", "UL54", "UL56")),
    clase_var = factor(ifelse(class == "NEW", "VND", class), levels = c("MRA", "VSI", "VND")),
    Prevalence_group = factor(case_when(
      Prevalence >= 20 ~ "\u226520%", Prevalence >= 15 ~ "\u226515%",
      Prevalence >= 10 ~ "\u226510%", TRUE ~ "\u22655%"),
      levels = c("\u22655%", "\u226510%", "\u226515%", "\u226520%")),
    Mutation = as.character(Mutation)
  ) %>%
  select(dia, x_pos, Protein, gen_color, Mutation, clase_var, Prevalence, Prevalence_group)

cat("RTP-2: muestras", nrow(muestras_rtp2), "| variantes (filas)", nrow(variantes_rtp2),
    "| únicas", n_distinct(variantes_rtp2$Mutation), "\n")
print(variantes_rtp2 %>% count(gen_color, clase_var))

# ── Regiones funcionales asociadas a resistencia (definidas a priori) ─────────
## Ajustar a los intervalos que se declaran en Material y métodos.
regiones_resist <- tibble(
  Protein = c("UL97", "UL54", "UL56"),
  ini     = c(335, 300, 230),
  fin     = c(665, 1000, 370)
)
variantes_rtp2 <- variantes_rtp2 %>%
  mutate(Position = as.numeric(gsub("[^0-9]", "", Mutation))) %>%
  left_join(regiones_resist, by = "Protein") %>%
  mutate(en_region = Position >= ini & Position <= fin) %>%
  select(-ini, -fin)

# Subconjunto para las figuras: todas las MRA + VND situadas en regiones de resistencia
variantes_rtp2_fig <- variantes_rtp2 %>% filter(clase_var == "MRA" | en_region)
cat("Variantes representadas:", n_distinct(variantes_rtp2_fig$Mutation), "de",
    n_distinct(variantes_rtp2$Mutation), "\n")
print(variantes_rtp2 %>% distinct(Protein, Mutation, clase_var, en_region) %>%
        count(Protein, clase_var, en_region))

# ── Tratamientos ─────────────────────────────────────────────────────────────
trat_df <- tibble(
  label  = c("GCV", "VGCV", "GCV", "MBV", "LEF", "FOS", "LMV", "CMV-Ig"),
  xmin_d = c( 42,    69,    271,   291,   361,   382,   409,    62),
  xmax_d = c( 68,   270,    290,   381,   381,   408,   432,   909),
  fila   = c("sup", "sup",  "sup", "sup", "top", "sup", "sup", "inf"),
  fill_col = c(col_GCV, col_VGCV, col_GCV, col_MBV, col_LEF, col_FOS, col_LMV, col_CMVIg),
  txt_col  = c(rep("white", 7), "#1A3A52")
)

make_rounded_bar <- function(xmin, xmax, ymid, h, fill_col, id, n_arc = 40) {
  r <- min(h / 2, (xmax - xmin) / 2)
  theta_l <- seq(pi / 2, 3 * pi / 2, length.out = n_arc)
  theta_r <- seq(-pi / 2, pi / 2, length.out = n_arc)
  tibble(x = c(xmin + r + r * cos(theta_l), xmax - r + r * cos(theta_r)),
         y = c(ymid + r * sin(theta_l), ymid + r * sin(theta_r)),
         group = id, fill_col = fill_col)
}

# ── Función de figura ────────────────────────────────────────────────────────
figura_rtp2 <- function(datos_mut, tamaños, y_mut_base = -0.4, y_mut_step = -0.55,
                        y_inf_trat = 5.60, y_sup_trat = 6.00, y_top_trat = 6.38,
                        bar_h = 0.35, y_lim_sup = 7, base_size = 13,
                        tam_etq = 3.8, tam_num = 3.5) {
  
  datos_mut <- datos_mut %>%
    mutate(gen_orden = as.integer(gen_color)) %>%
    arrange(gen_orden, clase_var, Mutation) %>%      # MRA antes que VND dentro de cada gen
    mutate(Mutation = factor(Mutation, levels = unique(Mutation)))
  mut_levels <- levels(datos_mut$Mutation)
  mut_y_map  <- setNames(y_mut_base + (seq_along(mut_levels) - 1) * y_mut_step, mut_levels)
  datos_mut  <- datos_mut %>% mutate(y_pos = mut_y_map[as.character(Mutation)])
  
  mut_tick_df <- datos_mut %>% distinct(Mutation, gen_color) %>%
    mutate(label = as.character(Mutation), y_pos = mut_y_map[label],
           col_text = colores_color[as.character(gen_color)])
  
  # Cruces: muestras sin resultado del gen, en las filas de ese gen
  cruces <- bind_rows(lapply(c("UL54", "UL56"), function(g) {
    muts <- mut_tick_df$label[mut_tick_df$gen_color == g]
    if (length(muts) == 0) return(NULL)
    expand_grid(dia = dias_sin_gen(g), Mutation = muts) %>%
      mutate(x_pos = interpolar_pos(dia), y_pos = mut_y_map[Mutation])
  }))
  
  y_trat <- c(sup = y_sup_trat, top = y_top_trat, inf = y_inf_trat)
  trat <- trat_df %>%
    mutate(y_row = y_trat[fila], xmin_pos = interpolar_pos(xmin_d), xmax_pos = interpolar_pos(xmax_d),
           x_mid = (xmin_pos + xmax_pos) / 2, ancho = xmax_pos - xmin_pos,
           txt_size = case_when(ancho > 3.5 ~ 3.4, ancho > 1.8 ~ 3.0, TRUE ~ 2.6))
  polys_trat <- pmap_dfr(list(trat$xmin_pos, trat$xmax_pos, trat$y_row, rep(bar_h, nrow(trat)),
                              trat$fill_col, seq_len(nrow(trat))), make_rounded_bar)
  
  mra_data <- datos_mut %>% filter(clase_var == "MRA")
  vnd_data <- datos_mut %>% filter(clase_var != "MRA")
  genes_presentes <- levels(droplevels(datos_mut$gen_color))
  y_lim_inf <- min(mut_y_map) - 0.4
  y_corte   <- y_inf_trat - 0.25
  
  ggplot() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = y_corte, ymax = y_lim_sup, fill = "white", color = NA) +
    geom_vline(xintercept = seq_along(niveles_dia), color = "#DDDDDD", linewidth = 0.25) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = y_corte, ymax = y_lim_sup, fill = "white", color = NA) +
    geom_polygon(data = polys_trat, aes(x, y, group = group), fill = polys_trat$fill_col, color = NA) +
    geom_text(data = trat, aes(x_mid, y_row, label = label), color = trat$txt_col,
              size = trat$txt_size, fontface = "bold") +
    annotate("text", x = rep(0.55, 5), y = 1:5, label = as.character(1:5), hjust = 1, size = 3.5, color = "grey35") +
    annotate("segment", x = rep(0.75, 5), xend = rep(0.97, 5), y = 1:5, yend = 1:5, color = "grey70", linewidth = 0.25) +
    geom_line(data = df_rtp2, aes(x_pos, log_cv, group = 1), linewidth = 0.6, color = "darkred") +
    geom_point(data = df_rtp2, aes(x_pos, log_cv, color = tipo_punto), size = 2.5) +
    geom_text(data = df_rtp2 %>% filter(!is.na(muestra)), aes(x_pos, log_cv, label = muestra),
              vjust = -0.8, size = tam_num, color = "darkred", fontface = "bold") +
    geom_segment(data = mut_tick_df, aes(x = -0.3, xend = n_niv + 0.5, y = y_pos, yend = y_pos),
                 color = "#DDDDDD", linewidth = 0.25) +
    geom_text(data = mut_tick_df, aes(x = -0.4, y = y_pos, label = label), color = mut_tick_df$col_text,
              hjust = 1, size = tam_etq, fontface = "bold") +
    annotate("text", x = -1.2, y = 3.0, label = "Log CV (UI/mL)", angle = 90, size = 4.0,
             color = "grey30", fontface = "bold", hjust = 0.5, vjust = 0) +
    geom_point(data = mra_data, aes(x_pos, y_pos, size = Prevalence_group, color = gen_color, fill = gen_color),
               shape = 21, stroke = 1, alpha = 0.80) +
    geom_point(data = vnd_data, aes(x_pos, y_pos, size = Prevalence_group, color = gen_color),
               shape = 21, fill = "white", stroke = 1) +
    geom_point(data = cruces, aes(x_pos, y_pos), shape = 4, size = 2.4, stroke = 0.8, color = "grey35") +
    scale_color_manual(values = c(colores_color, "Muestra secuenciada" = "darkred", "Muestra no secuenciada" = "black"),
                       breaks = c(genes_presentes, "Muestra secuenciada", "Muestra no secuenciada"), name = NULL) +
    scale_fill_manual(values = colores_fill, guide = "none") +
    scale_size_manual(values = tamaños, name = "Umbral") +
    scale_x_continuous(limits = c(-3.5, n_niv + 0.5), breaks = seq_along(niveles_dia),
                       labels = niveles_dia, expand = c(0, 0)) +
    scale_y_continuous(limits = c(y_lim_inf, y_lim_sup), breaks = 0:6, labels = NULL, expand = c(0, 0)) +
    guides(size = guide_legend(order = 1, nrow = 1,
                               override.aes = list(shape = 21, color = "grey50", fill = "grey70", stroke = 1)),
           color = guide_legend(order = 2, nrow = 1,
                                override.aes = list(shape = c(rep(21, length(genes_presentes)), 16, 16),
                                                    size  = c(rep(5, length(genes_presentes)), 4, 4),
                                                    fill  = c(colores_fill[genes_presentes], "darkred", "black"),
                                                    stroke = c(rep(1, length(genes_presentes)), 0, 0)))) +
    labs(x = "D\u00edas postrasplante", y = NULL) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(),
          axis.title.x = element_text(face = "bold", margin = margin(t = 12)),
          legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal",
          legend.title = element_text(size = 11, face = "bold"), legend.text = element_text(size = 11),
          legend.spacing.x = unit(0.25, "cm"),
          plot.margin = margin(5, 15, 5, 10),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))
}

guardar_rtp2 <- function(p, nombre, ancho, alto) {
  ggsave(file.path(dir_salida, paste0(nombre, ".png")), p,
         width = ancho, height = alto, units = "cm", dpi = 300, bg = "white")
  ggsave(file.path(dir_salida, paste0(nombre, ".svg")), p,
         device = svglite::svglite, fix_text_size = FALSE,
         width = ancho, height = alto, units = "cm", bg = "white")
}

## ---- Figura 49. Evolución de la carga viral y detección de MRA ----
mra_rtp2 <- variantes_rtp2 %>% filter(clase_var == "MRA")
cat("MRA en RTP-2:\n"); mra_rtp2 %>% select(dia, Protein, Mutation, Prevalence) %>% arrange(Protein, Mutation, dia) %>% print(n = Inf)

p_mra <- figura_rtp2(mra_rtp2, tam_5, y_mut_base = -0.4,
                     y_inf_trat = 6.28, y_sup_trat = 6.58, y_top_trat = 6.84, bar_h = 0.24, y_lim_sup = 7.25)
grid::grid.newpage(); print(p_mra)
guardar_rtp2(p_mra, "Figura_49", 32, 15)

## ---- Figura 50. Evolución de la carga viral y detección de MRA y VND ----
p_5 <- figura_rtp2(variantes_rtp2_fig, tam_5,
                   y_mut_base = -0.1, y_mut_step = -0.5,
                   base_size = 16, tam_etq = 5.2, tam_num = 4.5) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
        axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 12)),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text  = element_text(size = 14),
        legend.key.size = unit(0.9, "cm"))
grid::grid.newpage(); print(p_5)
guardar_rtp2(p_5, "Figura_50", 42, 24)

## ---- Complementaria. MRA y VND en RTP-2 al umbral ≥ 10% ----
variantes_rtp2_10 <- variantes_rtp2_fig %>% filter(Prevalence >= 10)
p_10 <- figura_rtp2(variantes_rtp2_10, tam_10, y_mut_base = -0.8)
grid::grid.newpage(); print(p_10)
guardar_rtp2(p_10, "Complementaria_RTP2_MRA_VND_10pct", 35, 21)


