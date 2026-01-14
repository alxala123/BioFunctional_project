# Especificar un mirror de CRAN
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Lista de paquetes requeridos
required_packages <- c(
  "shiny",              
  "shinydashboard",     
  "shinycssloaders",    
  "httr",              
  "readr",              
  "dplyr",              
  "tibble",            
  "future.apply",       
  "DT",                 
  "rvest",              
  "visNetwork",         
  "bslib",              
  "fastmap",            
  "shinyBS",            
  "shinyjs",            
  "plotly",
  "progressr",
  "rentrez",
  "tidyr",
  "stringr",
  "jsonlite",
  "rlang",
  "renv"
)

# Instalar y cargar los paquetes
installed_packages <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}



# Function to find the kegg structure
ancestors_kegg <- function(data) {
  relations <- unique(data$ONTOLOGY)
  processed_data <- data.frame()  
  
  for (i in seq_along(relations)) {
    relation <- relations[i] 
    
    # Obtain the reactions and interactions
    metabolic_domain <- ""   
    metabolic_subdomain <- ""    
    
    # read the website and extract the table
    url <- paste0("https://www.genome.jp/dbget-bin/www_bget?pathway:", relation)
    webpage <- read_html(url)
    
    # web scraping ont he table
    tables <- html_nodes(webpage, "table")
    table <- html_table(tables[[1]], fill = TRUE)
    # online
    # Find the line with the relevant info
    class_row_index <- which(table[, 1] == "Class")
    if (length(class_row_index) > 0) {
      class_row <- table[class_row_index, ]
      
      class_values <- unlist(strsplit(as.character(class_row), ";"))
      metabolic_domain <- class_values[2]
      metabolic_subdomain <- gsub("BRITE hierarchy", "", class_values[3])
    }
    
    # Create a dataframe with the processsed data
    subset_data <- subset(data, ONTOLOGY == relation)
    subset_data$metabolic_domain <- metabolic_domain
    subset_data$metabolic_subdomain <- metabolic_subdomain
    processed_data <- rbind(processed_data, subset_data)
  }
  processed_data$metabolic_domain <- ifelse(processed_data$metabolic_domain == "BRITE hierarchy", NA, processed_data$metabolic_domain)
  return(processed_data)
}


ancestors_gene_ontologies <- function(ontologies, groups) {
  # Función para obtener ancestros de QuickGO
  get_ancestors <- function(ontology) {
    url <- sprintf("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/%s/ancestors?relations=is_a%%2Cpart_of%%2Coccurs_in%%2Cregulates",
                   URLencode(ontology, reserved = TRUE))
    response <- tryCatch({
      content(GET(url, accept("application/json")), "parsed")
    }, error = function(e) {
      print(paste("Error al obtener datos para:", ontology))
      return(NULL)
    })
    
    # Depuración: Verificar respuesta de la API
    print(paste("Consultando:", ontology))
    print(response)
    
    # Validar respuesta y extraer ancestros
    if (!is.null(response$results) && length(response$results) > 0) {
      if (!response$results[[1]]$isObsolete) {
        ancestors <- response$results[[1]]$ancestors
        
        # Check if ancestors only contains the term itself
        if (length(ancestors) == 1 && ancestors[[1]] == ontology) {
          return(NA)  # Return NA if there are no true ancestors
        } else {
          return(ancestors)  # Return ancestors as they are
        }
      }
    }
    return(NA)  # If no ancestors found
  }
  
  # Obtener ancestros en paralelo
  plan(multisession)  # Usa varios núcleos (ajustar si hay problemas en Windows)
  ancestors <- future_lapply(ontologies, get_ancestors)
  
  # Convertir a caracteres y evitar valores NULL
  ancestors <- sapply(ancestors, function(x) ifelse(is.null(x), NA, toString(x)))
  
  # Crear dataframe asegurando que cada ontología mantenga su grupo
  data <- data.frame(GROUP = groups, ONTOLOGY = ontologies, ANCESTORS = ancestors, stringsAsFactors = FALSE)
  
  # Filtrar filas sin ancestros
  data <- data[!is.na(data$ANCESTORS) & data$ANCESTORS != "", ]
  
  return(data)
}





# Function to retrieve children of a GO term
get_children_quickgo <- function(go_id) {
  # Construct URL to retrieve children of a GO term
  base_url <- "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"
  url <- paste0(base_url, go_id, "/children")
  response <- httr::GET(url)
  if (httr::http_type(response) == "application/json") {
    children <- httr::content(response, "parsed")
    return(children)
  } else {
    stop("Error: The response is not in JSON format.")
  }               
}

# Function to retrieve information of children of a GO term
get_children_info <- function(go_term) {
  children_quickgo <- get_children_quickgo(go_term)
  children_list <- children_quickgo$results
  children_df <- data.frame(id = character(), name = character(), stringsAsFactors = FALSE)
  for (child_info in children_list[[1]]$children) {
    child_id <- child_info$id
    child_name <- child_info$name
    child_df <- data.frame(id = child_id, name = child_name, stringsAsFactors = FALSE)
    children_df <- rbind(children_df, child_df)
  }
  return(children_df)
}


# Función para obtener el nombre de la ontología desde QuickGO
get_ontology_name <- function(go_id) {
  if (is.na(go_id) || go_id == "") return(NA)
  
  url <- paste0("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/", go_id)
  response <- try(GET(url, accept_json()), silent = TRUE)
  
  if (inherits(response, "try-error") || status_code(response) != 200) return(NA)
  
  content <- content(response, as = "parsed", type = "application/json")
  
  if (!is.null(content$results) && length(content$results) > 0) {
    return(coalesce(content$results[[1]]$name, NA))
    
  }
  
  return(NA)
}

find_first_matching_ancestor <- function(ancestors) {
  if (is.null(ancestors) || length(ancestors) == 0) {
    return(NA)
  }
  
  for (ancestor in ancestors) {
    if (!is.null(children_info_mf) && ancestor %in% children_info_mf$id) {
      return(ancestor)
    }
    if (!is.null(children_info_bp) && ancestor %in% children_info_bp$id) {
      return(ancestor)
    }
    if (!is.null(children_info_cc) && ancestor %in% children_info_cc$id) {
      return(ancestor)
    }
  }
  return(NA)
}


# Define function to get the name of the first ancestor
get_first_ancestor_name <- function(ancestor_id) {
  if (ancestor_id %in% children_info_mf$id) {
    return(children_info_mf$name[children_info_mf$id == ancestor_id])
  } else if (ancestor_id %in% children_info_bp$id) {
    return(children_info_bp$name[children_info_bp$id == ancestor_id])
  } else if (ancestor_id %in% children_info_cc$id) {
    return(children_info_cc$name[children_info_cc$id == ancestor_id])
  } else {
    return(NA)
  }
}



# Define function remove_duplicate_observations before analyze_regulation
remove_duplicate_observations <- function(data) {
  # Encuentra las observaciones duplicadas por "ontology", "group" y "EA_VALUE"
  duplicated_rows <- duplicated(data[, c("ONTOLOGY", "GROUP", "EA_VALUE")]) | duplicated(data[, c("ONTOLOGY", "GROUP", "EA_VALUE")], fromLast = TRUE)
  
  # Cambia GOEA_up_DOWN a "NEUTRAL" en las observaciones duplicadas
  data[duplicated_rows, "UP_DOWN"] <- "NEUTRAL"
  
  # Elimina las filas duplicadas basadas en "ONTOLOGY", "GROUP_1", "GROUP_2", y "EA_VALUE"
  data <- data[!duplicated(data[, c("ONTOLOGY", "GROUP_1", "GROUP_2", "EA_VALUE")]), ]
  
  return(data)
}


# Define function analyze_regulation
analyze_regulation <- function(data) {
  
  # Split the data into groups based on the combination of "ONTOLOGY" and "GROUP" columns
  # Each unique combination of "ONTOLOGY" and "GROUP" will create a separate subset
  pairs <- split(data, paste(data$ONTOLOGY, data$GROUP))
  
  # Initialize a numeric vector to store the row names of the observations that will be kept
  observations_to_keep <- numeric()
  
  # Iterate over each subset (pair) of observations
  for (pair in pairs) {
    # Convert each pair (subset) into a data frame
    pair_df <- as.data.frame(pair)
    
    # Check if the subset contains exactly two observations
    if (nrow(pair_df) == 2) {
      # Extract the "EA_VALUE" for both observations
      EA_VALUE_1 <- pair_df[["EA_VALUE"]][1]
      EA_VALUE_2 <- pair_df[["EA_VALUE"]][2]
      
      # Compare the "EA_VALUE" of the two observations
      if (EA_VALUE_1 != EA_VALUE_2) {
        # If the values are different, keep the observation with the higher "EA_VALUE"
        observation_to_keep <- ifelse(EA_VALUE_1 > EA_VALUE_2, rownames(pair_df)[1], rownames(pair_df)[2])
        observations_to_keep <- c(observations_to_keep, observation_to_keep)
      } else {
        # If the "EA_VALUE" values are equal, keep both observations
        observations_to_keep <- c(observations_to_keep, rownames(pair_df))
      }
    } else if (nrow(pair_df) == 1) {
      # If there is only one observation in the subset, keep it
      observations_to_keep <- c(observations_to_keep, rownames(pair_df)[1])
    }
  }
  
  # Filter the original data to keep only the observations whose row names are in "observations_to_keep"
  data <- data[rownames(data) %in% observations_to_keep, ]
  
  # Return the filtered data after removing any duplicate observations (assumed functionality)
  return(remove_duplicate_observations(data))
}

# Función para limpiar el nombre del taxón o gen
clean_name <- function(name) {
  cleaned <- gsub(" ", "_", name)  # Reemplazar espacios por guiones bajos
  return(cleaned)
}

# Función para obtener datos de UniProt
get_uniprot_data <- function(taxon, gene, All_fields = FALSE) {
  # Limpiar los nombres de los parámetros
  cleaned_taxon <- clean_name(taxon)
  cleaned_gene <- clean_name(gene)
  
  # Construir la URL base para la consulta
  base_url <- "https://rest.uniprot.org/uniprotkb/stream?compressed=true"
  
  # Seleccionar los campos a obtener
  if (All_fields) {
    fields <- "&fields=accession,go_f,go_p,go_c,xref_kegg,other_fields_if_any"
  } else {
    fields <- "&fields=accession,go_f,go_p,go_c,xref_kegg"
  }
  
  # Crear la consulta combinada para taxón y gen
  query <- paste0("query=", URLencode(paste("(", cleaned_taxon, ")", "AND", "(", cleaned_gene, ")", "AND", "(reviewed:true)", sep = " ")))
  full_url <- paste0(base_url, "&", fields, "&format=tsv&", query)
  
  # Mostrar la URL de la solicitud HTTP para depuración
  cat("HTTP Request URL:", full_url, "\n")
  
  # Realizar la solicitud HTTP
  content <- httr::GET(full_url, write_disk(tf <- tempfile(fileext = ".tsv")))
  
  # Verificar si la solicitud fue exitosa
  if (content$status_code != 200) {
    cat("Error: HTTP request failed with status", content$status_code, "\n")
    return(NULL)
  } else {
    # Leer y devolver los datos
    data <- read.delim(tf, stringsAsFactors = FALSE)
    return(data)
  }
}

# Función principal GANGO
calculate_enrichment_scores_kegg <- function(kegg_data, gene_id_col = "KEGG_Pathway", alpha = 0.05) {
  count_cols <- grep("^count_", colnames(kegg_data), value = TRUE)
  if (length(count_cols) < 2) {
    stop("El dataset debe contener al menos dos columnas que comiencen con 'count_' para calcular el enriquecimiento.")
  }
  
  all_results <- data.frame()
  full_results <- data.frame()
  
  for (i in 1:(length(count_cols) - 1)) {
    for (j in (i + 1):length(count_cols)) {
      group1_col <- count_cols[i]
      group2_col <- count_cols[j]
      
      # Filtrar solo filas con conteo en alguno de los grupos
      kegg_data_filtered <- kegg_data %>% filter((.data[[group1_col]] > 0) | (.data[[group2_col]] > 0))
      
      gene_sets <- unique(kegg_data_filtered[[gene_id_col]])
      
      if (length(gene_sets) == 0) {
        cat("⚠️ No hay términos KEGG para analizar entre", group1_col, "y", group2_col, "\n")
        next
      }
      
      kegg_data_filtered <- kegg_data_filtered %>%
        mutate(logFC = log2(.data[[group2_col]] + 1) - log2(.data[[group1_col]] + 1))
      
      calculate_enrichment_score <- function(set_genes) {
        ranked_logFC_set <- kegg_data_filtered$logFC[kegg_data_filtered[[gene_id_col]] == set_genes]
        if (length(ranked_logFC_set) == 0) return(NA)
        NES <- sum(ranked_logFC_set) / sqrt(sum(kegg_data_filtered$logFC^2, na.rm = TRUE))
        return(NES)
      }
      
      create_contingency_table <- function(set_genes) {
        count_in_set_group1 <- sum(kegg_data_filtered[[group1_col]][kegg_data_filtered[[gene_id_col]] == set_genes], na.rm = TRUE)
        count_in_set_group2 <- sum(kegg_data_filtered[[group2_col]][kegg_data_filtered[[gene_id_col]] == set_genes], na.rm = TRUE)
        count_not_in_set_group1 <- sum(kegg_data_filtered[[group1_col]], na.rm = TRUE) - count_in_set_group1
        count_not_in_set_group2 <- sum(kegg_data_filtered[[group2_col]], na.rm = TRUE) - count_in_set_group2
        matrix(c(count_in_set_group1, count_not_in_set_group1, count_in_set_group2, count_not_in_set_group2), nrow = 2)
      }
      
      calculate_pvalue <- function(set_genes) {
        contingency_table <- create_contingency_table(set_genes)
        if (any(is.na(contingency_table))) return(NA)
        if (any(contingency_table < 5)) {
          return(fisher.test(contingency_table)$p.value)
        } else {
          return(suppressWarnings(chisq.test(contingency_table, correct = FALSE)$p.value))
        }
      }
      
      enrichment_scores <- sapply(gene_sets, calculate_enrichment_score)
      pvalues <- sapply(gene_sets, calculate_pvalue)
      
      results <- data.frame(
        KEGG_Pathway = gene_sets,
        EA_VALUE = enrichment_scores,
        pvalue = pvalues,
        GROUP_1 = group1_col,
        GROUP_2 = group2_col,
        GROUP = paste0(group1_col, "_vs_", group2_col),
        stringsAsFactors = FALSE
      )
      
      results$pvalue_corrected <- p.adjust(results$pvalue, method = "BH")
      results <- results %>%
        mutate(
          FDR = pvalue_corrected,
          UP_DOWN = case_when(
            EA_VALUE > 0 ~ "UP",
            EA_VALUE < 0 ~ "DOWN",
            TRUE ~ "NEUTRAL"
          )
        ) %>%
        arrange(FDR)
      
      full_results <- rbind(full_results, results)
      
      significant_results <- results %>%
        filter(!is.na(FDR) & FDR < alpha)
      
      all_results <- rbind(all_results, significant_results)
    }
  }
  
  if (nrow(all_results) == 0) {
    cat("No se encontraron términos KEGG significativos con el umbral de FDR:", alpha, "\n")
  }
  
  return(list(significant_results = all_results, full_results = full_results))
}


GANGO <- function(taxon, gene, group, All_fields = FALSE) {
  cat("=== INICIANDO GANGO ===\n")
  
  # Descargar lista de species KEGG válidas
  valid_species <- tryCatch({
    read.delim("https://rest.kegg.jp/list/organism", header = FALSE)
  }, error = function(e) {
    stop("❌ No se pudo descargar la lista de species KEGG.")
  })
  valid_codes <- valid_species$V2
  
  combined_data <- data.frame()
  
  # 1) Obtener datos de UniProt y combinar
  for (i in seq_along(taxon)) {
    cat(sprintf("[INFO] Procesando taxon %d/%d: %s - gen: %s\n",
                i, length(taxon), taxon[i], gene[i]))
    data <- get_uniprot_data(taxon[i], gene[i], All_fields)
    if (is.null(data) || nrow(data) == 0) {
      warning(sprintf("⚠️ No se encontraron datos para %s - %s", taxon[i], gene[i]))
    } else {
      data$TAXON <- taxon[i]
      data$GENE  <- gene[i]
      data$GROUP <- group[i]
      combined_data <- rbind(combined_data, data)
    }
  }
  
  if (nrow(combined_data) == 0) {
    cat("❌ No se encontraron datos en ningún caso.\n")
    return(NULL)
  }
  cat(sprintf("✅ Datos combinados: %d filas\n", nrow(combined_data)))
  
  # 2) Procesar GO
  cat("📦 Procesando anotaciones GO...\n")
  data_long <- combined_data %>%
    pivot_longer(cols = starts_with("Gene.Ontology"),
                 names_to = "Ontology_Type",
                 values_to = "Ontology") %>%
    separate_rows(Ontology, sep = "; ") %>%
    filter(!is.na(Ontology)) %>%
    mutate(
      Ontology_Type = case_when(
        str_detect(Ontology_Type, "molecular.function") ~ "molecular_function",
        str_detect(Ontology_Type, "biological.process")  ~ "biological_process",
        str_detect(Ontology_Type, "cellular.component")  ~ "cellular_component",
        TRUE ~ Ontology_Type
      )
    ) %>%
    select(Ontology, Ontology_Type, GROUP) %>%
    extract(Ontology,
            into = c("Ontology_Name", "Ontology_ID"),
            regex = "(.*)\\s*(\\[GO:\\d+\\])",
            remove = FALSE) %>%
    mutate(Ontology_ID = str_replace_all(Ontology_ID, "\\[|\\]", ""))
  
  data_grouped <- data_long %>%
    group_by(Ontology_Name, Ontology_ID, Ontology_Type, GROUP) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    pivot_wider(names_from = GROUP,
                values_from = Count,
                values_fill = 0) %>%
    rename_with(~ paste0("count_", .),
                -c(Ontology_Name, Ontology_ID, Ontology_Type))
  
  enrichment_results_GO <- calculate_enrichment_scores(data_grouped)
  enrichment_results_GO$Ontology_Name <- data_grouped$Ontology_Name
  enrichment_results_GO$Ontology_ID   <- data_grouped$Ontology_ID
  enrichment_results_GO$Ontology_Type <- data_grouped$Ontology_Type
  cat("✅ GO procesado correctamente\n")
  
  # 3) Procesar KEGG
  cat("🔗 Procesando rutas KEGG...\n")
  kegg_data <- combined_data %>%
    select(Entry, KEGG, TAXON, GENE, GROUP) %>%
    filter(!is.na(KEGG) & KEGG != "") %>%
    mutate(KEGG = str_replace_all(KEGG, ";$", "")) %>%
    separate_rows(KEGG, sep = ";") %>%
    mutate(
      KEGG    = str_trim(KEGG),
      SPECIES = sub(":.*", "", KEGG)
    ) %>%
    rename(GeneID = Entry, KEGG_ID = KEGG) %>%
    filter(SPECIES %in% valid_codes)
  
  if (nrow(kegg_data) == 0) {
    cat("⚠️ No se encontraron IDs KEGG válidos.\n")
    enrichment_results_KEGG <- NULL
  } else {
    species_needed <- unique(kegg_data$SPECIES)
    species_maps <- list()
    
    for (sp in species_needed) {
      cat(sprintf("📥 Descargando pathways para especie %s…\n", sp))
      link_url <- paste0("https://rest.kegg.jp/link/pathway/", sp)
      list_url <- paste0("https://rest.kegg.jp/list/pathway/", sp)
      
      species_maps[[sp]] <- tryCatch({
        link_df <- read.delim(link_url, header = FALSE, stringsAsFactors = FALSE)
        colnames(link_df) <- c("KEGG_ID", "KEGG_Pathway")
        link_df$KEGG_ID <- gsub(paste0(sp, ":"), "", link_df$KEGG_ID)
        link_df$KEGG_Pathway <- gsub("path:", "", link_df$KEGG_Pathway)
        
        list_df <- read.delim(list_url, header = FALSE, stringsAsFactors = FALSE)
        colnames(list_df) <- c("Pathway_Full", "KEGG_Name")
        list_df$KEGG_Pathway <- gsub("path:", "", list_df$Pathway_Full)
        
        merged <- merge(link_df, list_df[, c("KEGG_Pathway", "KEGG_Name")],
                        by = "KEGG_Pathway", all.x = TRUE)
        merged
      }, error = function(e) {
        warning(sprintf("❌ Falló descarga mapping o nombres para %s", sp))
        NULL
      })
    }
    
    # Mapear KEGG_ID a pathways y nombres
    kegg_data <- kegg_data %>%
      rowwise() %>%
      mutate(
        KEGG_Pathway = {
          df <- species_maps[[SPECIES]]
          if (!is.null(df)) {
            hits <- df$KEGG_Pathway[df$KEGG_ID == gsub(paste0(SPECIES, ":"), "", KEGG_ID)]
            if (length(hits) > 0) paste(unique(hits), collapse = ";") else NA_character_
          } else {
            NA_character_
          }
        },
        KEGG_Name = {
          df <- species_maps[[SPECIES]]
          if (!is.null(df)) {
            hits <- df$KEGG_Name[df$KEGG_ID == gsub(paste0(SPECIES, ":"), "", KEGG_ID)]
            if (length(hits) > 0) paste(unique(hits), collapse = ";") else NA_character_
          } else {
            NA_character_
          }
        }
      ) %>%
      ungroup()
    
    
    kegg_pathway_data <- kegg_data %>%
      separate_rows(KEGG_Pathway, KEGG_Name, sep = ";") %>%
      filter(!is.na(KEGG_Pathway) & KEGG_Pathway != "")
    
    if (nrow(kegg_pathway_data) == 0) {
      cat("⚠️ Ningún pathway KEGG mapeado.\n")
      enrichment_results_KEGG <- NULL
    } else {
      kegg_grouped <- kegg_pathway_data %>%
        group_by(KEGG_Pathway, KEGG_Name, GROUP) %>%
        summarise(Count = n(), .groups = 'drop') %>%
        pivot_wider(names_from = GROUP,
                    values_from = Count,
                    values_fill = 0) %>%
        rename_with(~ paste0("count_", .),
                    -c(KEGG_Pathway, KEGG_Name))
      
      enrichment_results_KEGG <- calculate_enrichment_scores_kegg(kegg_grouped)
      
      if (!is.null(enrichment_results_KEGG$significant_results)) {
        enrichment_results_KEGG$significant_results <- enrichment_results_KEGG$significant_results %>%
          left_join(kegg_grouped[, c("KEGG_Pathway", "KEGG_Name")], by = "KEGG_Pathway")
      }
      cat("✅ KEGG pathways procesados correctamente\n")
    }
  }
  
  cat("✅ GANGO finalizado con éxito\n")
  return(list(
    GO_results       = enrichment_results_GO,
    KEGG_results     = enrichment_results_KEGG,
    combined_data    = combined_data
  ))
}

# Calcular los puntajes de enriquecimiento
calculate_enrichment_scores <- function(go_data, gene_id_col = "Ontology_ID", alpha = 0.05) {
  # Verificar que existen columnas de recuento adecuadas
  count_cols <- grep("^count_", colnames(go_data), value = TRUE)
  if (length(count_cols) < 2) {
    stop("El dataset debe contener al menos dos columnas que comiencen con 'count_' para calcular el enriquecimiento.")
  }
  
  # Inicialización de resultados
  all_results <- data.frame()
  full_results <- data.frame()
  
  for (i in 1:(length(count_cols) - 1)) {
    for (j in (i + 1):length(count_cols)) {
      group1_col <- count_cols[i]
      group2_col <- count_cols[j]
      
      # Cálculo del log2 Fold Change (logFC) - para obtener una medición del cambio relativo
      go_data <- go_data %>%
        mutate(logFC = log2(.data[[group2_col]] + 1) - log2(.data[[group1_col]] + 1))
      
      # Asegúrate de que `ranked_logFC` esté correctamente definido
      rank_genes_logFC <- setNames(go_data$logFC, go_data[[gene_id_col]])
      
      # Verifica que `ranked_logFC` se haya definido correctamente
      if (length(rank_genes_logFC) == 0) {
        stop("ranked_logFC no se ha definido correctamente. Revisa los datos.")
      }
      
      gene_sets <- unique(go_data[[gene_id_col]])
      
      # Función para calcular el Enrichment Score (NES)
      calculate_enrichment_score <- function(set_genes) {
        ranked_logFC_set <- rank_genes_logFC[names(rank_genes_logFC) %in% set_genes]
        if (length(ranked_logFC_set) == 0) return(NA)
        NES <- sum(ranked_logFC_set) / sqrt(sum(rank_genes_logFC^2, na.rm = TRUE))
        return(NES)
      }
      
      # Función para crear la tabla de contingencia
      create_contingency_table <- function(set_genes) {
        count_in_set_group1 <- sum(go_data[[group1_col]][go_data[[gene_id_col]] %in% set_genes], na.rm = TRUE)
        count_in_set_group2 <- sum(go_data[[group2_col]][go_data[[gene_id_col]] %in% set_genes], na.rm = TRUE)
        count_not_in_set_group1 <- sum(go_data[[group1_col]], na.rm = TRUE) - count_in_set_group1
        count_not_in_set_group2 <- sum(go_data[[group2_col]], na.rm = TRUE) - count_in_set_group2
        matrix(c(count_in_set_group1, count_not_in_set_group1, count_in_set_group2, count_not_in_set_group2), nrow = 2)
      }
      
      # Función para calcular el p-valor
      calculate_pvalue <- function(set_genes) {
        contingency_table <- create_contingency_table(set_genes)
        
        # Verificar si la tabla tiene valores NA
        if (any(is.na(contingency_table))) {
          return(NA)
        }
        
        # Verificar si alguna celda tiene menos de 5, en ese caso usar Fisher
        if (any(contingency_table < 5)) {
          return(fisher.test(contingency_table)$p.value)
        } else {
          return(suppressWarnings(chisq.test(contingency_table, correct = FALSE)$p.value))
        }
      }
      
      # Calcular los scores de enriquecimiento y p-values
      enrichment_scores <- sapply(gene_sets, calculate_enrichment_score)
      pvalues <- sapply(gene_sets, calculate_pvalue)
      
      # Crear el dataframe de resultados
      results <- data.frame(
        ONTOLOGY = gene_sets,
        EA_VALUE = enrichment_scores,
        pvalue = pvalues,
        GROUP_1 = group1_col,
        GROUP_2 = group2_col,
        GROUP = paste0(group1_col, "_vs_", group2_col),
        row.names = NULL
      )
      
      # Corregir los p-valores con el método BH
      results$pvalue_corrected <- p.adjust(results$pvalue, method = "BH")
      
      # Añadir columnas adicionales y ordenar por FDR
      results <- results %>%
        mutate(
          FDR = pvalue_corrected,
          UP_DOWN = case_when(
            EA_VALUE > 0 ~ "UP",
            EA_VALUE < 0 ~ "DOWN",
            TRUE ~ "NEUTRAL"
          )
        ) %>%
        arrange(FDR)
      
      # Guardar todos los resultados completos
      full_results <- rbind(full_results, results)
      
      # Filtrar los resultados significativos
      significant_results <- results %>%
        filter(!is.na(FDR) & FDR < alpha)
      
      # Guardar los resultados significativos
      all_results <- rbind(all_results, significant_results)
    }
  }
  
  # Verificar si hay resultados significativos
  if (nrow(all_results) == 0) {
    cat("No se encontraron términos GO significativos con el umbral de FDR:", alpha, "\n")
  }
  
  # Ahora, realizar el merge para agregar Ontology_Name y Ontology_Type a significant_results
  significant_results <- merge(
    significant_results,                      # DataFrame con resultados significativos
    go_data %>% select(Ontology_ID, Ontology_Name, Ontology_Type),  # Seleccionar solo las columnas necesarias de go_data
    by.x = "ONTOLOGY",                        # Usamos ONTOLOGY como clave en significant_results
    by.y = "Ontology_ID",                     # Usamos Ontology_ID como clave en go_data
    all.x = TRUE                               # Mantener todas las filas de significant_results
  )
  # Ajuste de nombres de columnas en los resultados finales
  significant_results <- significant_results %>%
    rename(
      ONT_NAME = Ontology_Name,
      ONT_DESCRIPTION = Ontology_Type
    )
  
  
  # Retornar los resultados significativos y completos
  return(list(significant_results = significant_results, full_results = full_results))
}

# Función para leer archivo CSV en Shiny
read_taxon_file <- function(file) {
  data <- read.csv(file$datapath, stringsAsFactors = FALSE)
  return(data)
}

# Get children information for different ontology types
children_info_mf <- get_children_info("GO:0003674")
children_info_bp <- get_children_info("GO:0008150")
children_info_cc <- get_children_info("GO:0005575")

