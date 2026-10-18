# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 1. Datos comunes: AUC recalculada desde las cargas seriadas
# Extracto del script completo (ver README.md del estudio)

# ---- 1. Datos comunes: AUC recalculada desde las cargas seriadas ----
## Trapecio sobre log10(CV) entre Fini y Ffin de cada episodio, con las fechas
## tal como están en dnaemias_sel. Negativas y "< LDC" = 0. Si no hay
## determinación en Ffin, se añade un 0 en esa fecha. Sustituye a la AUCCV de
## Stata, que sumaba el tramo entre la determinación previa y la primera positiva.
## Requiere DTA_CV_TXP.dta en el directorio de trabajo.

AUC_importada <- local({
  cargas <- read_dta(f_cargas) %>%
    zap_labels() %>%
    transmute(
      NHC   = as.numeric(NHC),
      fecha = as.Date(Freg),
      cv    = if_else(!is.na(CV) & CV > 0, as.numeric(CV), 0)
    ) %>%
    filter(!is.na(fecha)) %>%
    group_by(NHC, fecha) %>%
    summarise(cv = max(cv), .groups = "drop")
  
  episodios <- read_excel(f_bd,
                          sheet = "dnaemias_sel") %>%
    transmute(NHC  = as.numeric(NHC),
              Fini = as.Date(Fini),
              Ffin = as.Date(Ffin)) %>%
    filter(!is.na(Fini), !is.na(Ffin))
  
  auc_trapecio <- function(fecha, cv) {
    if (length(cv) < 2) return(NA_real_)
    y <- ifelse(cv > 0, pmax(log10(cv), 0), 0)
    t <- as.numeric(fecha)
    sum(diff(t) * (head(y, -1) + tail(y, -1)) / 2)
  }
  
  episodios %>%
    inner_join(cargas, by = "NHC", relationship = "many-to-many") %>%
    filter(fecha >= Fini, fecha <= Ffin) %>%
    group_by(NHC, Fini) %>%
    arrange(fecha, .by_group = TRUE) %>%
    summarise(
      AUC_importada = {
        f <- fecha; v <- cv
        if (!any(f == first(Ffin))) { f <- c(f, first(Ffin)); v <- c(v, 0) }
        auc_trapecio(f, v)
      },
      .groups = "drop"
    )
})

## Función para incorporar la AUC cada vez que se lee dnaemias_sel:
## AUCCV pasa a contener AUC_importada; la de Stata queda en AUCCV_stata.
incorporar_AUC <- function(df) {
  n_filas <- nrow(df)
  df <- df %>%
    select(-any_of(c("AUC_importada", "AUCCV_stata"))) %>%
    mutate(.NHC = as.numeric(NHC), .Fini = as.Date(Fini)) %>%
    left_join(AUC_importada %>% rename(.NHC = NHC, .Fini = Fini),
              by = c(".NHC", ".Fini"), relationship = "many-to-one") %>%
    select(-.NHC, -.Fini) %>%
    mutate(AUCCV_stata = AUCCV,
           AUCCV       = AUC_importada)
  stopifnot(nrow(df) == n_filas,
            sum(!is.na(df$AUCCV)) == sum(!is.na(df$AUCCV_stata)))
  df
}

## Comprobación (esperado: 220 episodios, mediana 73,0)
AUC_importada %>%
  summarise(episodios = n(),
            mediana   = median(AUC_importada),
            p25       = quantile(AUC_importada, 0.25),
            p75       = quantile(AUC_importada, 0.75))


