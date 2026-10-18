# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 1.5. Clasificación CHARMD: MRA, VSI, PN y VND
# Extracto del script completo (ver README.md del estudio)

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

