# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.7.04. Datos para el texto (9.7.2): tiempo hasta la primera DNAemia
# Extracto del script completo (ver README.md del estudio)

## ---- Datos para el texto (9.7.2): tiempo hasta la primera DNAemia ----
df_base <- dnaemias %>%
  dplyr::filter(primer_ev == 1) %>%
  dplyr::mutate(
    Ftx              = as.Date(Ftx),
    Fini             = as.Date(Fini),
    meses_tx_inicio  = as.numeric(Fini - Ftx) / 30.44,
    serDR = factor(serDR, levels = c("D+/R-", "D+/R+", "D-/R+"))
  ) %>%
  dplyr::filter(!is.na(serDR), !is.na(meses_tx_inicio))

df_A <- df_base %>%
  dplyr::select(-any_of(c("gbm", "ghm", "gbghm"))) %>%
  dplyr::left_join(pacientes %>% dplyr::select(NHC, gbm, ghm, gbghm),
                   by = "NHC")

df_B <- df_base

fmt_mediana_ric <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("NA")
  med <- median(x)
  ric <- quantile(x, probs = c(0.25, 0.75))
  sprintf("%.1f (RIC %.1f-%.1f)", med, ric[1], ric[2])
}

# >>> FUNCION CORREGIDA: extraer vectores directamente para wilcox.test <<<
comparar_grupos <- function(datos, var_mixto, serologia, etiqueta) {
  sub <- datos %>%
    dplyr::filter(serDR == serologia) %>%
    dplyr::filter(!is.na(.data[[var_mixto]])) %>%
    dplyr::filter(!is.na(meses_tx_inicio))
  
  if (nrow(sub) == 0) {
    cat(sprintf("  %s: sin datos\n", etiqueta))
    return(invisible(NULL))
  }
  
  # Extraer vectores directamente con [[ ]]
  vals     <- sub$meses_tx_inicio
  grupos   <- sub[[var_mixto]]
  
  mixtos <- vals[grupos == 1]
  unicos <- vals[grupos == 0]
  
  cat(sprintf("  %s:\n", etiqueta))
  cat(sprintf("    Mixtos  (n=%d): %s\n",
              length(mixtos), fmt_mediana_ric(mixtos)))
  cat(sprintf("    Unicos  (n=%d): %s\n",
              length(unicos), fmt_mediana_ric(unicos)))
  
  # >>> CORRECCION: pasar vectores directamente, no formula <<<
  if (length(mixtos) >= 1 && length(unicos) >= 1) {
    p <- tryCatch(
      stats::wilcox.test(x = mixtos, y = unicos, exact = FALSE)$p.value,
      error = function(e) {
        cat(sprintf("    ERROR en wilcox: %s\n", e$message))
        NA
      }
    )
    if (!is.na(p)) {
      if (p < 0.001) {
        cat("    p < 0.001\n")
      } else {
        cat(sprintf("    p = %.3f\n", p))
      }
    } else {
      cat("    p = NA\n")
    }
  }
}

for (grupo in c("D+/R+", "D-/R+", "D+/R-")) {
  cat("\n========================================\n")
  cat(sprintf("GRUPO %s\n", grupo))
  cat("========================================\n")
  
  n_grupo <- df_base %>%
    dplyr::filter(serDR == grupo) %>%
    dplyr::pull(NHC) %>%
    dplyr::n_distinct()
  cat(sprintf("N pacientes: %d\n\n", n_grupo))
  
  med_global <- df_base %>%
    dplyr::filter(serDR == grupo) %>%
    dplyr::pull(meses_tx_inicio)
  cat(sprintf("Mediana global: %s\n\n", fmt_mediana_ric(med_global)))
  
  cat("ABORDAJE A (nivel paciente):\n")
  comparar_grupos(df_A, "gbm",   grupo, "gB")
  comparar_grupos(df_A, "ghm",   grupo, "gH")
  comparar_grupos(df_A, "gbghm", grupo, "gB+gH")
  
  cat("\nABORDAJE B (nivel episodio):\n")
  comparar_grupos(df_B, "gbmix",   grupo, "gB")
  comparar_grupos(df_B, "ghmix",   grupo, "gH")
  comparar_grupos(df_B, "gbghmix", grupo, "gB+gH")
}

cat("\n========================================\n")
cat("DISTRIBUCION TEMPORAL D+/R+\n")
cat("========================================\n")

dist_dpr <- df_base %>%
  dplyr::filter(serDR == "D+/R+") %>%
  dplyr::mutate(tras_6m = meses_tx_inicio > 6) %>%
  dplyr::count(tras_6m) %>%
  dplyr::mutate(pct = round(100 * n / sum(n), 1))
print(dist_dpr)

