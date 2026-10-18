# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 0.2. Librerías
# Extracto del script completo (ver README.md del estudio)

## ---- Librerías ----
# MASS va primero para que dplyr::select no quede enmascarado.
library(MASS)
library(readxl)
library(haven)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(stringr)
library(lubridate)
library(scales)
library(ggplot2)
library(ggrepel)
library(ggtext)
library(ggpubr)
library(rstatix)
library(patchwork)
library(cowplot)
library(grid)
library(reshape2)
library(flextable)
library(officer)
library(survival)
library(survminer)
library(lme4)
library(lmerTest)
library(car)
library(dunn.test)
library(DiagrammeR)

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
recode <- dplyr::recode   # car::recode lo enmascara

