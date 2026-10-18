# Estudio 2 · Variantes de UL97, UL54 y UL56 mediante NGS
# Apartado: 1.6. Identificadores RTP-1 a RTP-11 y tablas de trabajo
# Extracto del script completo (ver README.md del estudio)

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


