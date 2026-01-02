############################################################
# PharmaToxAI :: GEO PIPELINE
# Módulo de ingesta automática de datasets transcriptómicos
#
# Propósito dentro de PharmaToxAI:
#   - Buscar en GEO (NCBI) series (GSE) relevantes para un
#     término biológico/toxicológico (ej. "shiga toxin",
#     "cisplatin nephrotoxicity", "oxidative stress", etc.).
#   - Filtrar por organismo de interés (ej. Homo sapiens,
#     Mus musculus) para mantener coherencia biológica.
#   - Descargar matrices procesadas (SeriesMatrix) y, cuando
#     estén disponibles, archivos RAW (supplementary).
#   - Estandarizar la estructura de carpetas de entrada para:
#       * módulos in silico de PharmaToxAI
#         (EDA, DEG, enriquecimiento, QSAR sobre firmas)
#       * integración con resultados in vitro
#         (ensayos de toxicidad/actividad biológica).
#   - Generar una tabla resumen + log de descargas que sirva
#     como "catálogo de datasets" para priorizar qué GSE
#     entran al pipeline analítico de PharmaToxAI.
#
# Notas de integración:
#   - Este script SOLO hace ingesta/organización de datos
#     (data ingestion). No corre análisis estadístico.

############################################################


#############################
# CONFIGURACIÓN BÁSICA  #
#############################

# Directorio base donde se guardará todo
BASE_DIR    <- "/home/user/Cosmetics_genes/geo_datasets"

# Término de búsqueda principal (lo más general posible, evita que la búsqueda no encuentre nada)
SEARCH_TERM <- "cosmetics ingredients"

# Organismo (tal como lo usa NCBI)
ORGANISM    <- "Homo sapiens" # nombres comunes: Bos taurus, Mus musculus, Homo sapiens

# Tipo de entrada en GEO:
# - "gse[Filter]" = GEO Series (datasets)
ENTRY_TYPE  <- "gse[Filter]"

# lo bueno es que podemos agregar otros términos para construir un query más complejo

# Máximo de resultados (ajustar si esperás muchos datasets)
RETMAX      <- 10

# Nombre de archivos de salida
SUMMARY_CSV  <- file.path(BASE_DIR, "geo_search_results.csv")
DOWNLOAD_LOG <- file.path(BASE_DIR, "geo_download_log.csv")

# Flags para controlar qué se descarga
DOWNLOAD_MATRIX <- TRUE   # TRUE = baja SeriesMatrix (GEOquery::getGEO)
DOWNLOAD_RAW    <- TRUE   # TRUE = baja archivos suplementarios (getGEOSuppFiles)


################################
# CARGA DE LIBRERÍAS       #
################################

# Asumo que ya están instaladas.
# Si no: install.packages("rentrez"); BiocManager::install("GEOquery")

library(rentrez)
library(GEOquery)
library(dplyr)
library(readr)
library(purrr)
library(stringr)


################################
# FUNCIONES AUXILIARES      #
################################

################################
# FUNCIONES DE CÁLCULO Y ANÁLISIS
################################

# Calcular coeficiente de variación (CV) para cada gen
calculate_gene_cv <- function(expr_matrix) {
  # CV = sd / mean para cada fila (gen)
  gene_means <- rowMeans(expr_matrix, na.rm = TRUE)
  gene_sds   <- apply(expr_matrix, 1, sd, na.rm = TRUE)
  
  cv <- gene_sds / gene_means
  cv[is.infinite(cv) | is.na(cv)] <- NA_real_
  
  return(cv)
}

# Detectar muestras outlier basado en correlación
detect_outlier_samples <- function(expr_matrix, threshold = 0.85) {
  # Calcular matriz de correlación entre muestras
  cor_matrix <- cor(expr_matrix, use = "pairwise.complete.obs")
  
  # Para cada muestra, calcular correlación media con otras muestras
  mean_cors <- colMeans(cor_matrix, na.rm = TRUE)
  
  # Identificar outliers: correlación media < threshold
  outliers <- names(mean_cors)[mean_cors < threshold]
  
  tibble(
    Sample          = colnames(expr_matrix),
    Mean_correlation = mean_cors,
    Is_outlier      = colnames(expr_matrix) %in% outliers
  )
}

# Calcular estadísticas de datos faltantes
calculate_missing_stats <- function(expr_matrix) {
  total_values <- prod(dim(expr_matrix))
  missing_values <- sum(is.na(expr_matrix) | !is.finite(expr_matrix))
  missing_pct <- (missing_values / total_values) * 100
  
  # Por muestra
  missing_per_sample <- colSums(is.na(expr_matrix) | !is.finite(expr_matrix))
  
  # Por gen
  missing_per_gene <- rowSums(is.na(expr_matrix) | !is.finite(expr_matrix))
  
  list(
    total_missing_pct = missing_pct,
    missing_per_sample = missing_per_sample,
    missing_per_gene = missing_per_gene,
    samples_with_high_missing = names(missing_per_sample)[missing_per_sample > ncol(expr_matrix) * 0.5],
    genes_with_high_missing = names(missing_per_gene)[missing_per_gene > nrow(expr_matrix) * 0.5]
  )
}

# Identificar genes de baja varianza
identify_low_variance_genes <- function(expr_matrix, percentile = 0.1) {
  # Calcular varianza para cada gen
  gene_vars <- apply(expr_matrix, 1, var, na.rm = TRUE)
  
  # Identificar genes en el percentil más bajo
  threshold <- quantile(gene_vars, percentile, na.rm = TRUE)
  low_var_genes <- names(gene_vars)[gene_vars <= threshold & !is.na(gene_vars)]
  
  tibble(
    Gene = rownames(expr_matrix),
    Variance = gene_vars,
    Low_variance = rownames(expr_matrix) %in% low_var_genes
  )
}

# Calcular estadísticas de correlación entre muestras
calculate_sample_correlations <- function(expr_matrix) {
  cor_matrix <- cor(expr_matrix, use = "pairwise.complete.obs")
  
  # Extraer valores únicos (triángulo superior sin diagonal)
  cor_values <- cor_matrix[upper.tri(cor_matrix)]
  
  tibble(
    Min_correlation = min(cor_values, na.rm = TRUE),
    Q25_correlation = quantile(cor_values, 0.25, na.rm = TRUE),
    Median_correlation = median(cor_values, na.rm = TRUE),
    Q75_correlation = quantile(cor_values, 0.75, na.rm = TRUE),
    Max_correlation = max(cor_values, na.rm = TRUE),
    Mean_correlation = mean(cor_values, na.rm = TRUE)
  )
}

# Realizar análisis de componentes principales (PCA)
perform_pca_analysis <- function(expr_matrix) {
  # Transponer: muestras en filas, genes en columnas
  expr_t <- t(expr_matrix)
  
  # Remover columnas con varianza cero o NA
  valid_cols <- apply(expr_t, 2, function(x) var(x, na.rm = TRUE) > 0 & !all(is.na(x)))
  expr_t_clean <- expr_t[, valid_cols]
  
  # Imputar NAs con la mediana de cada columna (gen)
  for (i in 1:ncol(expr_t_clean)) {
    expr_t_clean[is.na(expr_t_clean[, i]), i] <- median(expr_t_clean[, i], na.rm = TRUE)
  }
  
  # Realizar PCA
  pca_result <- tryCatch({
    prcomp(expr_t_clean, center = TRUE, scale. = TRUE)
  }, error = function(e) {
    message("Error en PCA: ", conditionMessage(e))
    return(NULL)
  })
  
  return(pca_result)
}

################################
# FUNCIONES DE VISUALIZACIÓN
################################

# Crear gráfico de distribución de expresión
plot_expression_distribution <- function(expr_matrix, gse_id, output_dir) {
  pdf_file <- file.path(output_dir, paste0(gse_id, "_expression_distribution.pdf"))
  
  tryCatch({
    pdf(pdf_file, width = 12, height = 8)
    par(mfrow = c(1, 2))
    
    # Histograma
    vals <- as.numeric(expr_matrix)
    vals <- vals[is.finite(vals)]
    
    if (length(vals) > 0) {
      hist(vals, 
           breaks = 100, 
           main = paste(gse_id, "- Distribución de expresión"),
           xlab = "Valores de expresión",
           ylab = "Frecuencia",
           col = "steelblue",
           border = "white")
      
      # Gráfico de densidad
      plot(density(vals, na.rm = TRUE),
           main = paste(gse_id, "- Densidad de expresión"),
           xlab = "Valores de expresión",
           ylab = "Densidad",
           col = "darkred",
           lwd = 2)
      polygon(density(vals, na.rm = TRUE), col = rgb(1, 0, 0, 0.3), border = NA)
    }
    
    dev.off()
    message("   - Gráfico de distribución guardado: ", pdf_file)
    return(TRUE)
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    message("   ! Error creando gráfico de distribución: ", conditionMessage(e))
    return(FALSE)
  })
}

# Crear heatmap de correlación entre muestras
plot_sample_correlation_heatmap <- function(expr_matrix, gse_id, output_dir) {
  pdf_file <- file.path(output_dir, paste0(gse_id, "_sample_correlation.pdf"))
  
  tryCatch({
    cor_matrix <- cor(expr_matrix, use = "pairwise.complete.obs")
    
    pdf(pdf_file, width = 10, height = 10)
    
    # Crear heatmap personalizado
    image(1:ncol(cor_matrix), 1:nrow(cor_matrix), 
          t(cor_matrix),
          col = colorRampPalette(c("blue", "white", "red"))(100),
          xlab = "Muestras", ylab = "Muestras",
          main = paste(gse_id, "- Correlación entre muestras"),
          axes = FALSE)
    
    # Añadir escala de colores
    axis(1, at = 1:ncol(cor_matrix), labels = colnames(cor_matrix), las = 2, cex.axis = 0.7)
    axis(2, at = 1:nrow(cor_matrix), labels = rownames(cor_matrix), las = 2, cex.axis = 0.7)
    
    dev.off()
    message("   - Heatmap de correlación guardado: ", pdf_file)
    return(TRUE)
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    message("   ! Error creando heatmap: ", conditionMessage(e))
    return(FALSE)
  })
}

# Crear gráfico PCA
plot_pca <- function(pca_result, gse_id, output_dir) {
  if (is.null(pca_result)) {
    message("   ! No se puede crear gráfico PCA: resultado nulo")
    return(FALSE)
  }
  
  pdf_file <- file.path(output_dir, paste0(gse_id, "_pca.pdf"))
  
  tryCatch({
    pdf(pdf_file, width = 12, height = 8)
    par(mfrow = c(1, 2))
    
    # Gráfico de varianza explicada
    var_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
    barplot(var_explained[1:min(10, length(var_explained))],
            names.arg = paste0("PC", 1:min(10, length(var_explained))),
            main = paste(gse_id, "- Varianza explicada por PC"),
            xlab = "Componente Principal",
            ylab = "% Varianza explicada",
            col = "steelblue",
            las = 2)
    
    # Gráfico PC1 vs PC2
    plot(pca_result$x[, 1], pca_result$x[, 2],
         xlab = paste0("PC1 (", round(var_explained[1], 2), "%)"),
         ylab = paste0("PC2 (", round(var_explained[2], 2), "%)"),
         main = paste(gse_id, "- PCA Plot"),
         pch = 19,
         col = rgb(0, 0, 1, 0.5),
         cex = 1.5)
    
    # Añadir etiquetas de muestra
    text(pca_result$x[, 1], pca_result$x[, 2], 
         labels = rownames(pca_result$x),
         cex = 0.6, pos = 3)
    
    dev.off()
    message("   - Gráfico PCA guardado: ", pdf_file)
    return(TRUE)
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    message("   ! Error creando gráfico PCA: ", conditionMessage(e))
    return(FALSE)
  })
}

# Crear boxplot de expresión por muestra
plot_sample_boxplots <- function(expr_matrix, gse_id, output_dir) {
  pdf_file <- file.path(output_dir, paste0(gse_id, "_sample_boxplots.pdf"))
  
  tryCatch({
    pdf(pdf_file, width = max(12, ncol(expr_matrix) * 0.3), height = 8)
    
    boxplot(expr_matrix,
            main = paste(gse_id, "- Distribución de expresión por muestra"),
            xlab = "Muestras",
            ylab = "Valores de expresión",
            las = 2,
            col = rainbow(ncol(expr_matrix)),
            cex.axis = 0.7,
            outline = FALSE)
    
    dev.off()
    message("   - Boxplot de muestras guardado: ", pdf_file)
    return(TRUE)
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    message("   ! Error creando boxplot: ", conditionMessage(e))
    return(FALSE)
  })
}

# Crear gráfico de coeficiente de variación
plot_cv_distribution <- function(cv_values, gse_id, output_dir) {
  pdf_file <- file.path(output_dir, paste0(gse_id, "_cv_distribution.pdf"))
  
  tryCatch({
    cv_clean <- cv_values[is.finite(cv_values)]
    
    pdf(pdf_file, width = 10, height = 6)
    
    hist(cv_clean,
         breaks = 100,
         main = paste(gse_id, "- Distribución del Coeficiente de Variación"),
         xlab = "CV (SD/Mean)",
         ylab = "Frecuencia",
         col = "darkorange",
         border = "white")
    
    abline(v = median(cv_clean, na.rm = TRUE), col = "red", lwd = 2, lty = 2)
    legend("topright", 
           legend = paste("Mediana CV =", round(median(cv_clean, na.rm = TRUE), 3)),
           col = "red", lty = 2, lwd = 2)
    
    dev.off()
    message("   - Gráfico CV guardado: ", pdf_file)
    return(TRUE)
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    message("   ! Error creando gráfico CV: ", conditionMessage(e))
    return(FALSE)
  })
}

# Generar reporte de calidad completo
generate_quality_report <- function(gse_id, expr_matrix, output_dir) {
  message("   - Generando reporte de calidad para ", gse_id)
  
  # Calcular todas las métricas
  cv_values <- calculate_gene_cv(expr_matrix)
  outlier_info <- detect_outlier_samples(expr_matrix)
  missing_stats <- calculate_missing_stats(expr_matrix)
  low_var_genes <- identify_low_variance_genes(expr_matrix)
  cor_stats <- calculate_sample_correlations(expr_matrix)
  
  # Realizar PCA
  pca_result <- perform_pca_analysis(expr_matrix)
  
  # Crear todos los gráficos
  plot_expression_distribution(expr_matrix, gse_id, output_dir)
  plot_sample_correlation_heatmap(expr_matrix, gse_id, output_dir)
  plot_sample_boxplots(expr_matrix, gse_id, output_dir)
  plot_cv_distribution(cv_values, gse_id, output_dir)
  
  if (!is.null(pca_result)) {
    plot_pca(pca_result, gse_id, output_dir)
  }
  
  # Crear resumen de calidad en texto
  report_file <- file.path(output_dir, paste0(gse_id, "_quality_report.txt"))
  
  sink(report_file)
  cat("=" , rep("=", 60), "\n", sep = "")
  cat("REPORTE DE CALIDAD - ", gse_id, "\n")
  cat("=" , rep("=", 60), "\n\n", sep = "")
  
  cat("DIMENSIONES DEL DATASET:\n")
  cat("  - Número de genes/features: ", nrow(expr_matrix), "\n")
  cat("  - Número de muestras: ", ncol(expr_matrix), "\n\n")
  
  cat("ESTADÍSTICAS DE DATOS FALTANTES:\n")
  cat("  - Porcentaje total de datos faltantes: ", 
      round(missing_stats$total_missing_pct, 2), "%\n")
  cat("  - Muestras con >50% datos faltantes: ", 
      length(missing_stats$samples_with_high_missing), "\n")
  cat("  - Genes con >50% datos faltantes: ", 
      length(missing_stats$genes_with_high_missing), "\n\n")
  
  cat("ESTADÍSTICAS DE CORRELACIÓN ENTRE MUESTRAS:\n")
  cat("  - Correlación mínima: ", round(cor_stats$Min_correlation, 3), "\n")
  cat("  - Correlación mediana: ", round(cor_stats$Median_correlation, 3), "\n")
  cat("  - Correlación media: ", round(cor_stats$Mean_correlation, 3), "\n")
  cat("  - Correlación máxima: ", round(cor_stats$Max_correlation, 3), "\n\n")
  
  cat("MUESTRAS OUTLIER (correlación media < 0.85):\n")
  outliers <- outlier_info %>% filter(Is_outlier)
  if (nrow(outliers) > 0) {
    for (i in 1:nrow(outliers)) {
      cat("  - ", outliers$Sample[i], " (cor = ", 
          round(outliers$Mean_correlation[i], 3), ")\n", sep = "")
    }
  } else {
    cat("  - No se detectaron outliers\n")
  }
  cat("\n")
  
  cat("COEFICIENTE DE VARIACIÓN (CV):\n")
  cv_clean <- cv_values[is.finite(cv_values)]
  cat("  - CV mediano: ", round(median(cv_clean, na.rm = TRUE), 3), "\n")
  cat("  - CV medio: ", round(mean(cv_clean, na.rm = TRUE), 3), "\n")
  cat("  - Genes con CV alto (>1): ", sum(cv_clean > 1, na.rm = TRUE), "\n\n")
  
  cat("GENES DE BAJA VARIANZA (percentil <10%):\n")
  low_var_count <- sum(low_var_genes$Low_variance, na.rm = TRUE)
  cat("  - Número de genes de baja varianza: ", low_var_count, "\n")
  cat("  - Porcentaje: ", round(low_var_count / nrow(expr_matrix) * 100, 2), "%\n\n")
  
  if (!is.null(pca_result)) {
    var_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
    cat("ANÁLISIS DE COMPONENTES PRINCIPALES (PCA):\n")
    cat("  - Varianza explicada por PC1: ", round(var_explained[1], 2), "%\n")
    cat("  - Varianza explicada por PC2: ", round(var_explained[2], 2), "%\n")
    cat("  - Varianza acumulada PC1-PC2: ", 
        round(sum(var_explained[1:2]), 2), "%\n\n")
  }
  
  cat("=" , rep("=", 60), "\n", sep = "")
  cat("Fecha de generación: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("=" , rep("=", 60), "\n", sep = "")
  
  sink()
  
  message("   - Reporte de calidad guardado: ", report_file)
  
  # Retornar resumen en formato tibble
  return(tibble(
    GSE_ID = gse_id,
    N_genes = nrow(expr_matrix),
    N_samples = ncol(expr_matrix),
    Missing_pct = missing_stats$total_missing_pct,
    N_outlier_samples = sum(outlier_info$Is_outlier),
    Median_correlation = cor_stats$Median_correlation,
    Median_CV = median(cv_clean, na.rm = TRUE),
    N_low_var_genes = low_var_count,
    PC1_var = if (!is.null(pca_result)) round((pca_result$sdev[1]^2) / sum(pca_result$sdev^2) * 100, 2) else NA,
    PC2_var = if (!is.null(pca_result)) round((pca_result$sdev[2]^2) / sum(pca_result$sdev^2) * 100, 2) else NA
  ))
}

# Armar término de búsqueda para Entrez/GEO
build_geo_term <- function(search_term, organism, entry_type) {
  # Ejemplo: "shiga toxin AND Homo sapiens[ORGN] AND gse[Filter]"
  term <- paste0(search_term, " AND ", organism, "[ORGN] AND ", entry_type)
  return(term)
}

# Buscar GSE en GEO (db = "gds") vía rentrez
search_geo_series <- function(search_term,
                              organism,
                              entry_type,
                              retmax = 500) {
  term <- build_geo_term(search_term, organism, entry_type)
  message(">> Buscando en GEO con término: ", term)

  res <- entrez_search(
    db     = "gds",
    term   = term,
    retmax = retmax
  )

  if (length(res$ids) == 0) {
    message(">> No se encontraron GSE para el término dado.")
    return(NULL)
  }

  message(">> Se encontraron ",
          length(res$ids),
          " entradas (posibles GSE).")
  return(res$ids)  # ids numéricos del GDS
}

# Obtener metadatos de cada GSE
fetch_gse_metadata <- function(gse_ids,
                               search_term,
                               organism) {
  # gse_ids: vector de IDs numéricos (ej. "32710")
  message(">> Recuperando metadata de cada GSE...")

  summaries <- map(gse_ids, ~ entrez_summary(db = "gds", id = .x))

  df <- tibble(
    GSE_ID       = paste0("GSE", sapply(summaries, `[[`, "gse")),
    GDS_ID       = gse_ids,
    Title        = sapply(summaries, `[[`, "title"),
    Platform_org = sapply(summaries, `[[`, "platform_organism"),
    N_samples    = sapply(summaries, `[[`, "n_samples"),
    EntryType    = sapply(summaries, `[[`, "entrytype"),
    PubMed_ID    = sapply(summaries, function(x) {
      pmid <- x$pubmed_ids
      if (length(pmid) == 0)
        return(NA_character_)
      paste(pmid, collapse = ";")
    }),
    Summary      = sapply(summaries, `[[`, "summary"),
    Query_term   = search_term,
    Organism     = organism
  )

  return(df)
}

# Crear subcarpeta para un GSE
make_gse_dir <- function(base_dir, gse_id) {
  # gse_id tipo "GSE32710"
  gse_dir <- file.path(base_dir, gse_id)
  if (!dir.exists(gse_dir)) {
    dir.create(gse_dir, recursive = TRUE, showWarnings = FALSE)
  }
  return(gse_dir)
}

# Descargar matriz procesada (SeriesMatrix) de un GSE
#     Devuelve TRUE/FALSE según éxito
download_gse_matrix <- function(gse_id, dest_dir) {
  message("   - Descargando SeriesMatrix para ", gse_id, " ...")

  res <- tryCatch({
    gse_obj <- getGEO(
      GEO        = gse_id, # ID tipo "GSE32710"
      GSEMatrix  = TRUE, # TRUE para obtener ExpressionSet
      AnnotGPL   = DOWNLOAD_RAW # TRUE para anotar con plataforma
    )

    saveRDS(
      gse_obj,
      file = file.path(dest_dir, paste0(gse_id, "_GEOmatrix.rds"))
    )
    TRUE
  }, error = function(e) {
    message("   ! Error en SeriesMatrix de ", gse_id, ": ", conditionMessage(e))
    FALSE
  })

  return(res)
}

# descargar archivos suplementarios (RAW) de un GSE. Devuelve TRUE/FALSE según éxito
download_gse_raw <- function(gse_id, dest_dir) {
  message("   - Descargando RAW (supplementary files) para ", gse_id, " ...")

  res <- tryCatch({
    # baseDir evita cambiar el working directory
    getGEOSuppFiles(
      GSE        = gse_id,
      baseDir    = dest_dir,
      makeDirectory = TRUE
    )
    TRUE
  }, error = function(e) {
    message("   ! Error en RAW de ", gse_id, ": ", conditionMessage(e))
    FALSE
  })

  return(res)
}

#Pipeline completo para descargar TODO de un GSE
#     Devuelve un tibble con estado de descarga
process_single_gse <- function(gse_id,
                               base_dir,
                               download_matrix = TRUE,
                               download_raw    = TRUE) {
  message(">> Procesando ", gse_id, " ...")

  gse_dir <- make_gse_dir(base_dir, gse_id)

  ok_matrix <- NA
  ok_raw    <- NA

  if (download_matrix) {
    ok_matrix <- download_gse_matrix(gse_id, gse_dir)
  }
  if (download_raw) {
    ok_raw <- download_gse_raw(gse_id, gse_dir)
  }

  message(">> Listo ", gse_id)

  tibble(
    GSE_ID        = gse_id,
    Matrix_ok     = ok_matrix,
    RAW_ok        = ok_raw,
    Download_time = Sys.time()
  )
}


################################
# EJECUCIÓN DEL PIPELINE    #
################################

# Crear carpeta base si no existe
if (!dir.exists(BASE_DIR)) {
  dir.create(BASE_DIR, recursive = TRUE, showWarnings = FALSE)
}

# Buscar GSEs en GEO
gse_numeric_ids <- search_geo_series(
  search_term = SEARCH_TERM,
  organism    = ORGANISM,
  entry_type  = ENTRY_TYPE,
  retmax      = RETMAX
)

if (is.null(gse_numeric_ids)) {
  stop("No se encontraron GSE para el criterio de búsqueda. Fin del script.")
}

# Recuperar metadatos y armar tabla resumen
gse_metadata <- fetch_gse_metadata(
  gse_ids     = gse_numeric_ids,
  search_term = SEARCH_TERM,
  organism    = ORGANISM
)

# Guardar tabla resumen a CSV
write_csv(gse_metadata, SUMMARY_CSV)
message(">> Tabla resumen guardada en: ", SUMMARY_CSV)

# Extraer lista de IDs tipo "GSE32710"
gse_ids <- unique(gse_metadata$GSE_ID)

message(">> GSE a descargar:")
print(gse_ids)


################################
# DESCARGA MASIVA           #
################################

# Recorremos cada GSE y descargamos matrices + RAW, imap_dfr agrega índice (posición en el vector) para log de progreso
download_log <- purrr::imap_dfr(
  gse_ids,
  ~ {
    message("[", .y, "/", length(gse_ids), "]")
    process_single_gse(
      gse_id          = .x,
      base_dir        = BASE_DIR,
      download_matrix = DOWNLOAD_MATRIX,
      download_raw    = DOWNLOAD_RAW
    )
  }
)

# Guardar log de descargas
write_csv(download_log, DOWNLOAD_LOG)

message("======================================================")
message("PIPELINE GEO COMPLETADO")
message("Base dir: ", BASE_DIR)
message("Término de búsqueda: '", SEARCH_TERM, "' | Organismo: ", ORGANISM)
message("Resumen en: ", SUMMARY_CSV)
message("Log descargas: ", DOWNLOAD_LOG)
message("======================================================")

############################################################
# RESUMEN AUTOMÁTICO DE LOS GEOmatrix DESCARGADOS
#   es probable que esta parte no funcione bien, ya que depende de ls datos descargados y no tienen un estándar a seguir sino que depende de cómo lo subieron
############################################################

if (DOWNLOAD_MATRIX) {

  # Biobase solo se usa acá
  suppressPackageStartupMessages({
    library(Biobase)
  })

  message(">> Escaneando GEOmatrix descargados para resumen exploratorio...")

  # Buscar todos los GEOmatrix RDS en BASE_DIR
  rds_files <- list.files(
    path       = BASE_DIR,
    pattern    = "_GEOmatrix\\.rds$",
    recursive  = TRUE,
    full.names = TRUE
  )

  if (length(rds_files) == 0) {
    message(">> No se encontraron GEOmatrix RDS para resumir.")
  } else {

    geo_eset_summary <- purrr::map_dfr(rds_files, function(rds_path) {
      gse_id <- basename(dirname(rds_path))  # carpeta = GSE ID

      message("   - Resumiendo ", gse_id, " ...")

      # Intentar leer el RDS sin romper todo
      eset <- tryCatch({
        obj <- readRDS(rds_path)
        if (is.list(obj)) obj[[1]] else obj
      }, error = function(e) {
        message("     ! Error leyendo ", gse_id, ": ", conditionMessage(e))
        return(NULL)
      })

      if (is.null(eset) || !inherits(eset, "ExpressionSet")) {
        return(NULL)
      }

      # Extraer matrices básicas
      expr  <- Biobase::exprs(eset)
      pheno <- Biobase::pData(eset)

      n_features <- nrow(expr)
      n_samples  <- ncol(expr)

      # Estadísticos rápidos de expresión
      vals <- as.numeric(expr)
      vals <- vals[is.finite(vals)]

      if (length(vals) == 0) {
        min_expr    <- NA_real_
        median_expr <- NA_real_
        max_expr    <- NA_real_
        scale_guess <- NA_character_
      } else {
        min_expr    <- as.numeric(stats::quantile(vals, 0.01, na.rm = TRUE))
        median_expr <- stats::median(vals, na.rm = TRUE)
        max_expr    <- as.numeric(stats::quantile(vals, 0.99, na.rm = TRUE))

        # Heurística muy simple de escala: si max < 100  puede log2like, si no, crudo / no log
        scale_guess <- ifelse(max(vals, na.rm = TRUE) < 100,
                              "log2_like",
                              "raw_like")
      }

      # Contar columnas de pheno con pocos niveles (buenas candidatas a "grupo")
      if (nrow(pheno) > 0) {
        n_levels <- sapply(pheno, function(x) length(unique(as.character(x))))
        # columnas con entre 2 y 10 niveles
        group_like_cols <- names(n_levels)[n_levels >= 2 & n_levels <= 10]

        # columnas cuyo nombre sugiere tratamiento/grupo, poco probable que acierte
        group_keywords <- grepl(
          pattern = "treat|group|condit|time|dose|cell.?line",
          x       = tolower(colnames(pheno)),
          ignore.case = TRUE
        )
        group_keyword_cols <- colnames(pheno)[group_keywords]

        # ¿aparece el termino de busqueda en algún texto de pheno?
        pheno_text <- paste(apply(pheno, 2, as.character), collapse = " ")
        has_term  <- grepl(SEARCH_TERM, pheno_text, ignore.case = TRUE)

      } else {
        group_like_cols    <- character(0)
        group_keyword_cols <- character(0)
        has_term       <- FALSE
      }

      tibble::tibble(
        GSE_ID          = gse_id,
        RDS_path        = rds_path,
        N_features      = n_features,
        N_samples       = n_samples,
        Expr_min_01     = min_expr,
        Expr_median     = median_expr,
        Expr_max_99     = max_expr,
        Scale_guess     = scale_guess,                       # "log2_like" vs "raw_like"
        Group_like_cols = paste(group_like_cols, collapse = ";"),
        Group_key_cols  = paste(group_keyword_cols, collapse = ";"),
        Has_term_in_pheno = has_term
      )
    })

    # Guardar resumen si se obtuvo algo
    if (nrow(geo_eset_summary) > 0) {
      ESET_SCAN_CSV <- file.path(BASE_DIR, "geo_eset_scan_summary.csv")
      readr::write_csv(geo_eset_summary, ESET_SCAN_CSV)
      message(">> Resumen de GEOmatrix guardado en: ", ESET_SCAN_CSV)
    } else {
      message(">> No se pudo generar resumen de GEOmatrix (sin ExpressionSet válidos).")
    }
  }
}

############################################################
# ANÁLISIS DE CALIDAD Y GENERACIÓN DE GRÁFICOS
############################################################

if (DOWNLOAD_MATRIX) {
  message(">> Generando análisis de calidad y visualizaciones para datasets descargados...")
  
  # Buscar todos los GEOmatrix RDS en BASE_DIR
  rds_files <- list.files(
    path       = BASE_DIR,
    pattern    = "_GEOmatrix\\.rds$",
    recursive  = TRUE,
    full.names = TRUE
  )
  
  if (length(rds_files) == 0) {
    message(">> No se encontraron GEOmatrix RDS para analizar.")
  } else {
    
    # Crear directorio para reportes si no existe
    reports_dir <- file.path(BASE_DIR, "quality_reports")
    if (!dir.exists(reports_dir)) {
      dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    quality_summary <- purrr::map_dfr(rds_files, function(rds_path) {
      gse_id <- basename(dirname(rds_path))
      
      message("   - Analizando ", gse_id, " ...")
      
      # Intentar leer el RDS
      eset <- tryCatch({
        obj <- readRDS(rds_path)
        if (is.list(obj)) obj[[1]] else obj
      }, error = function(e) {
        message("     ! Error leyendo ", gse_id, ": ", conditionMessage(e))
        return(NULL)
      })
      
      if (is.null(eset) || !inherits(eset, "ExpressionSet")) {
        return(NULL)
      }
      
      # Extraer matriz de expresión
      expr <- tryCatch({
        Biobase::exprs(eset)
      }, error = function(e) {
        message("     ! Error extrayendo expresión de ", gse_id, ": ", conditionMessage(e))
        return(NULL)
      })
      
      if (is.null(expr) || nrow(expr) == 0 || ncol(expr) == 0) {
        message("     ! Matriz de expresión vacía para ", gse_id)
        return(NULL)
      }
      
      # Crear directorio para este GSE dentro de reports
      gse_report_dir <- file.path(reports_dir, gse_id)
      if (!dir.exists(gse_report_dir)) {
        dir.create(gse_report_dir, recursive = TRUE, showWarnings = FALSE)
      }
      
      # Generar reporte completo de calidad
      quality_metrics <- tryCatch({
        generate_quality_report(gse_id, expr, gse_report_dir)
      }, error = function(e) {
        message("     ! Error generando reporte para ", gse_id, ": ", conditionMessage(e))
        return(NULL)
      })
      
      return(quality_metrics)
    })
    
    # Guardar resumen consolidado si se obtuvo algo
    if (!is.null(quality_summary) && nrow(quality_summary) > 0) {
      QUALITY_SUMMARY_CSV <- file.path(BASE_DIR, "quality_analysis_summary.csv")
      readr::write_csv(quality_summary, QUALITY_SUMMARY_CSV)
      message(">> Resumen de calidad consolidado guardado en: ", QUALITY_SUMMARY_CSV)
      message(">> Reportes individuales y gráficos en: ", reports_dir)
    } else {
      message(">> No se pudo generar resumen de calidad.")
    }
  }
}

message("======================================================")
message("ANÁLISIS COMPLETO FINALIZADO")
message("======================================================")

