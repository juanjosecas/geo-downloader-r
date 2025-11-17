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
        Has_shiga_in_pheno = has_shiga
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
