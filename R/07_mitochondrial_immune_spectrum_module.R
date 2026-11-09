
# ============================================================
# MODULE 07
# MITOCHONDRIAL–IMMUNE INTEGRATED SPECTRUM
#
# Independent Shiny module.
# Does NOT call the existing MLP.
# Does NOT modify existing modules.
# ============================================================


# ------------------------------------------------------------
# 1. HELPER FUNCTIONS
# ------------------------------------------------------------

mito_read_csv <- function(file) {
  tryCatch(
    read.csv(
      file$datapath,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    error = function(e) {
      stop(
        paste("Cannot read CSV:", conditionMessage(e))
      )
    }
  )
}


mito_normalize_sample_id <- function(x) {
  
  x <- trimws(as.character(x))
  
  is_tcga <- grepl("^TCGA-", x)
  
  x[is_tcga] <- substr(
    x[is_tcga],
    1,
    16
  )
  
  x
}


mito_validate_sample_id <- function(data) {

  if (!"Sample_ID" %in% names(data)) {
    stop("The file must contain a Sample_ID column.")
  }

  ids <- trimws(as.character(data$Sample_ID))

  if (any(is.na(ids) | ids == "")) {
    stop("Sample_ID cannot contain missing or empty values.")
  }

  unique_ids <- unique(ids)

  if (length(unique_ids) != 1) {
    stop("Each uploaded profile must contain exactly one unique Sample_ID.")
  }

  unique_ids[[1]]
}

mito_validate_methylation <- function(data,
                                      require_sample_id = FALSE) {
  
  required <- c("CpG", "Beta")
  
  missing_columns <- setdiff(
    required,
    names(data)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      paste(
        "Missing columns:",
        paste(missing_columns, collapse = ", ")
      )
    )
  }
  
  if (require_sample_id) {
    mito_validate_sample_id(data)
  }
  
  data$CpG <- trimws(as.character(data$CpG))
  data$Beta <- suppressWarnings(
    as.numeric(data$Beta)
  )
  
  if (any(is.na(data$CpG) | data$CpG == "")) {
    stop("CpG identifiers cannot be empty.")
  }
  
  if (anyDuplicated(data$CpG)) {
    stop("Duplicate CpG identifiers were detected.")
  }
  
  if (any(!is.finite(data$Beta))) {
    stop("Beta values must be finite numeric values.")
  }
  
  if (any(data$Beta < 0 | data$Beta > 1)) {
    stop("Beta values must be between 0 and 1.")
  }
  
  data
}


mito_map_expression_genes <- function(data,
                                      annotation) {

  required_data <- c("Gene", "Expression")
  required_annotation <- c(
    "ENSEMBL", "Gene_Symbol", "Mapping_status",
    "Symbol_count", "Entrez_count", "ENTREZID_candidates"
  )

  missing_data <- setdiff(required_data, names(data))
  missing_annotation <- setdiff(required_annotation, names(annotation))

  if (length(missing_data) > 0) {
    stop("Missing expression columns: ", paste(missing_data, collapse = ", "))
  }

  if (length(missing_annotation) > 0) {
    stop("Missing annotation columns: ", paste(missing_annotation, collapse = ", "))
  }

  data$Gene <- trimws(as.character(data$Gene))
  data$Expression <- suppressWarnings(as.numeric(data$Expression))

  if (any(is.na(data$Gene) | data$Gene == "")) {
    stop("Gene identifiers cannot be empty.")
  }

  if (anyDuplicated(data$Gene)) {
    stop("Duplicate uploaded gene identifiers were detected.")
  }

  if (any(!is.finite(data$Expression))) {
    stop("Expression values must be finite numeric values.")
  }

  annotation$ENSEMBL <- sub(
    "\\.[0-9]+$",
    "",
    as.character(annotation$ENSEMBL)
  )

  annotation$Gene_Symbol <- trimws(as.character(annotation$Gene_Symbol))
  annotation$Mapping_status <- as.character(annotation$Mapping_status)
  annotation$Symbol_count <- suppressWarnings(as.numeric(annotation$Symbol_count))
  annotation$Entrez_count <- suppressWarnings(as.numeric(annotation$Entrez_count))
  annotation$ENTREZID_candidates <- trimws(as.character(annotation$ENTREZID_candidates))

  valid_symbols <- annotation[
    annotation$Mapping_status == "Unique gene symbol" &
      annotation$Symbol_count == 1 &
      annotation$Entrez_count == 1 &
      !is.na(annotation$ENTREZID_candidates) &
      annotation$ENTREZID_candidates != "" &
      !is.na(annotation$Gene_Symbol) &
      annotation$Gene_Symbol != "" &
      !is.na(annotation$ENSEMBL) &
      annotation$ENSEMBL != "",
    ,
    drop = FALSE
  ]

  if (nrow(valid_symbols) == 0) {
    stop("No genes passed the strict annotation criteria.")
  }

  if (anyDuplicated(valid_symbols$ENSEMBL)) {
    stop("The strict annotation contains duplicated Ensembl IDs.")
  }

  unique_symbol_rows <- !duplicated(valid_symbols$Gene_Symbol) &
    !duplicated(valid_symbols$Gene_Symbol, fromLast = TRUE)

  unambiguous_symbols <- valid_symbols[unique_symbol_rows, , drop = FALSE]

  symbol_map <- setNames(
    unambiguous_symbols$ENSEMBL,
    unambiguous_symbols$Gene_Symbol
  )

  symbol_by_ensembl <- setNames(
    unambiguous_symbols$Gene_Symbol,
    unambiguous_symbols$ENSEMBL
  )

  gene_ids <- sub("\\.[0-9]+$", "", data$Gene)

  is_ensembl <- grepl("^ENSG[0-9]+$", gene_ids)
  mapped_ids <- rep(NA_character_, length(gene_ids))

  # Ensembl IDs are retained only if they pass strict annotation.
  valid_ensembl <- gene_ids[is_ensembl] %in% valid_symbols$ENSEMBL
  mapped_ids[is_ensembl] <- ifelse(
    valid_ensembl,
    gene_ids[is_ensembl],
    NA_character_
  )

  # Gene symbols are mapped only through the strict one-to-one table.
  is_symbol <- !is_ensembl
  mapped_ids[is_symbol] <- unname(symbol_map[data$Gene[is_symbol]])

  data$ENSEMBL <- mapped_ids
  data <- data[!is.na(data$ENSEMBL), , drop = FALSE]

  if (nrow(data) == 0) {
    stop("No uploaded genes passed strict TCGA annotation mapping.")
  }

  if (anyDuplicated(data$ENSEMBL)) {
    stop("Multiple uploaded identifiers map to the same Ensembl gene.")
  }

  data$Gene_Symbol <- unname(symbol_by_ensembl[data$ENSEMBL])

  data
}
mito_calculate_ssgsea <- function(expression_data,
                                  gene_sets,
                                  component_map) {
  
  expression_matrix <- matrix(
    expression_data$Expression,
    ncol = 1,
    dimnames = list(
      expression_data$ENSEMBL,
      "Uploaded_sample"
    )
  )
  
  common_genes <- intersect(
    rownames(expression_matrix),
    unique(unlist(gene_sets))
  )
  
  if (length(common_genes) < 20) {
    stop(
      "Insufficient overlap with the pathway gene sets. At least 20 mapped genes are required."
    )
  }
  
  expression_matrix <- expression_matrix[
    common_genes,
    ,
    drop = FALSE
  ]
  
  parameters <- GSVA::ssgseaParam(
    exprData = expression_matrix,
    geneSets = gene_sets,
    minSize = 5,
    maxSize = Inf,
    alpha = 0.25,
    normalize = TRUE
  )
  
  scores <- GSVA::gsva(
    parameters,
    verbose = FALSE
  )
  
  result <- data.frame(
    pathway = rownames(scores),
    ssGSEA_score = as.numeric(scores[, 1]),
    stringsAsFactors = FALSE
  )
  
  result <- merge(
    result,
    component_map,
    by = "pathway",
    all.x = TRUE,
    sort = FALSE
  )

  if (anyNA(result$component)) {
    stop("Some pathways have no component classification.")
  }

  result <- result[
    match(rownames(scores), result$pathway),
    ,
    drop = FALSE
  ]
  
  result[order(
    result$ssGSEA_score,
    decreasing = TRUE
  ), ]
}


mito_plot_pathways <- function(data) {
  
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = reorder(pathway, ssGSEA_score),
      y = ssGSEA_score,
      fill = component
    )
  ) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = c(
        "Mitochondrial-associated" = "#4169E1",
        "Other metabolic-related" = "#D4AF37"
      )
    ) +
    ggplot2::labs(
      x = NULL,
      y = "ssGSEA score",
      fill = "Pathway category"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8),
      legend.position = "bottom"
    )
}


mito_find_cpg_associations <- function(cpg_data,
                                       robust_pairs) {
  
  result <- robust_pairs[
    robust_pairs$CpG_ID %in% cpg_data$CpG,
    ,
    drop = FALSE
  ]
  
  if (nrow(result) == 0) {
    return(result)
  }
  
  result <- result[
    order(
      abs(result$Pooled_rho),
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  head(result, 30)
}


# ------------------------------------------------------------
# 2. USER INTERFACE
# ------------------------------------------------------------

mitoImmuneSpectrumUI <- function(id) {
  
  ns <- shiny::NS(id)
  
  shiny::tagList(
    
    shiny::tags$h3(
      "Mitochondrial–Immune Integrated Spectrum"
    ),
    
    shiny::tags$p(
      "Explore methylation, gene expression and transcriptomic pathway enrichment using independent or matched molecular profiles."
    ),
    
    shiny::tabsetPanel(
      
      id = ns("profile_tabs"),
      
      # ------------------------------------------------------
      # TAB 1: EPIGENETIC
      # ------------------------------------------------------
      
      shiny::tabPanel(
        "Epigenetic profile",
        
        shiny::br(),
        
        shiny::fileInput(
          ns("cpg_file"),
          "Upload CpG methylation CSV",
          accept = ".csv"
        ),
        
        shiny::tags$p(
          "Required columns: CpG, Beta. Optional column: Sample_ID."
        ),
        
        shiny::verbatimTextOutput(
          ns("cpg_status")
        ),
        
        shiny::plotOutput(
          ns("beta_plot"),
          height = "300px"
        ),
        
        shiny::h4(
          "TCGA CpG–gene associations"
        ),
        
        shiny::tableOutput(
          ns("cpg_association_table")
        )
      ),
      
      # ------------------------------------------------------
      # TAB 2: TRANSCRIPTOMIC
      # ------------------------------------------------------
      
      shiny::tabPanel(
        "Transcriptomic profile",
        
        shiny::br(),
        
        shiny::fileInput(
          ns("expression_file"),
          "Upload gene expression CSV",
          accept = ".csv"
        ),
        
        shiny::tags$p(
          "Required columns: Gene, Expression. Gene can be an Ensembl ID or an unambiguous gene symbol. Optional column: Sample_ID."
        ),
        
        shiny::tags$p(
          "Use normalized/log-transformed expression values. The uploaded profile is scored independently using the stored TCGA pathway gene sets."
        ),
        
        shiny::verbatimTextOutput(
          ns("expression_status")
        ),
        
        shiny::plotOutput(
          ns("expression_pathway_plot"),
          height = "650px"
        ),
        
        shiny::tableOutput(
          ns("expression_pathway_table")
        )
      ),
      
      # ------------------------------------------------------
      # TAB 3: INTEGRATED
      # ------------------------------------------------------
      
      shiny::tabPanel(
        "Integrated profile",
        
        shiny::br(),
        
        shiny::fileInput(
          ns("integrated_cpg_file"),
          "Upload methylation CSV",
          accept = ".csv"
        ),
        
        shiny::tags$p(
          "Required columns: Sample_ID, CpG, Beta."
        ),
        
        shiny::fileInput(
          ns("integrated_expression_file"),
          "Upload expression CSV",
          accept = ".csv"
        ),
        
        shiny::tags$p(
          "Required columns: Sample_ID, Gene, Expression."
        ),
        
        shiny::tags$p(
          "Both files must represent exactly one sample. TCGA barcodes are matched using their first 16 characters."
        ),
        
        shiny::verbatimTextOutput(
          ns("integrated_status")
        ),
        
        shiny::h4(
          "Integrated pathway profile"
        ),
        
        shiny::plotOutput(
          ns("integrated_pathway_plot"),
          height = "650px"
        ),
        
        shiny::h4(
          "Methylation profile"
        ),
        
        shiny::plotOutput(
          ns("integrated_beta_plot"),
          height = "300px"
        ),
        
        shiny::h4(
          "Matched CpG–gene associations"
        ),
        
        shiny::tableOutput(
          ns("integrated_association_table")
        )
      )
    ),
    
    shiny::br(),
    
    shiny::tags$div(
      style = paste(
        "padding: 12px;",
        "border-left: 4px solid #4169E1;",
        "background-color: #f4f6fc;"
      ),
      
      shiny::tags$strong(
        "Scientific interpretation"
      ),
      
      shiny::tags$p(
        "CpG–gene relationships represent statistical associations and do not establish causality."
      ),
      
      shiny::tags$p(
        "ssGSEA scores represent relative enrichment of transcriptomic gene sets. They are not direct measurements of mitochondrial function, metabolic flux or metabolite concentrations."
      ),
      
      shiny::tags$p(
        "This module does not predict clinical immunotherapy response, treatment benefit or survival."
      )
    )
  )
}


# ------------------------------------------------------------
# 3. SERVER
# ------------------------------------------------------------

mitoImmuneSpectrumServer <- function(id,
                                     project_root) {
  
  shiny::moduleServer(
    id,
    function(input, output, session) {
      
      tables_dir <- file.path(
        project_root,
        "results",
        "tables"
      )
      
      models_dir <- file.path(
        project_root,
        "results",
        "models"
      )
      
      robust_file <- file.path(
        tables_dir,
        "TCGA_SKCM_CpG_gene_correlation_robust_pairs.csv"
      )
      
      annotation_file <- file.path(
        tables_dir,
        "TCGA_SKCM_gene_annotation_final.csv"
      )
      
      genesets_file <- file.path(
        models_dir,
        "TCGA_SKCM_MitoImmune_ssGSEA_gene_sets_20261109_v1.rds"
      )

      pathway_scores_file <- file.path(
        tables_dir,
        "TCGA_SKCM_MitoImmune_ssGSEA_pathway_scores_20261109_v1.csv"
      )
      
      required_files <- c(
        robust_file,
        annotation_file,
        genesets_file,
        pathway_scores_file
      )
      
      if (!all(file.exists(required_files))) {
        stop(
          paste(
            "Required Mito-Immune module resources are missing:",
            paste(
              basename(required_files[!file.exists(required_files)]),
              collapse = ", "
            )
          )
        )
      }
      
      robust_pairs <- read.csv(
        robust_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      
      annotation <- read.csv(
        annotation_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      
      gene_sets <- readRDS(
        genesets_file
      )

      pathway_reference <- read.csv(
        pathway_scores_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      component_map <- unique(
        pathway_reference[, c("pathway", "component")]
      )

      if (anyDuplicated(component_map$pathway)) {
        stop("Pathway component reference contains duplicate pathways.")
      }

      if (nrow(component_map) != length(gene_sets)) {
        stop("Pathway component reference does not match the gene sets.")
      }
      
      # ------------------------------------------------------
      # EPIGENETIC PROFILE
      # ------------------------------------------------------
      
      cpg_data <- shiny::reactive({
        
        shiny::req(input$cpg_file)
        
        data <- mito_read_csv(
          input$cpg_file
        )
        
        mito_validate_methylation(
          data,
          require_sample_id = FALSE
        )
      })
      
      output$cpg_status <- shiny::renderText({
        
        shiny::req(input$cpg_file)
        
        tryCatch({
          
          data <- cpg_data()
          
          paste0(
            "Validation: PASSED\n",
            "CpGs: ", nrow(data), "\n",
            "Mean beta: ", round(mean(data$Beta), 4), "\n",
            "Median beta: ", round(median(data$Beta), 4), "\n",
            "Sample ID: ",
            if ("Sample_ID" %in% names(data)) {
              paste(unique(data$Sample_ID), collapse = ", ")
            } else {
              "Not provided"
            }
          )
          
        }, error = function(e) {
          paste("Validation error:", conditionMessage(e))
        })
      })
      
      output$beta_plot <- shiny::renderPlot({
        
        data <- cpg_data()
        
        ggplot2::ggplot(
          data,
          ggplot2::aes(x = Beta)
        ) +
          ggplot2::geom_histogram(
            bins = 30,
            fill = "#4169E1",
            color = "white"
          ) +
          ggplot2::labs(
            x = "Beta value",
            y = "Number of CpGs",
            title = "Methylation beta distribution"
          ) +
          ggplot2::theme_minimal()
      })
      
      output$cpg_association_table <- shiny::renderTable({
        
        data <- cpg_data()
        
        result <- mito_find_cpg_associations(
          data,
          robust_pairs
        )
        
        if (nrow(result) == 0) {
          return(
            data.frame(
              Message = "No uploaded CpGs overlap with the TCGA robust association table."
            )
          )
        }
        
        result[, c(
          "CpG_ID",
          "Gene",
          "Mapping_Type",
          "Pooled_rho",
          "Pooled_FDR",
          "Primary_rho",
          "Metastatic_rho"
        ), drop = FALSE]
        
      }, striped = TRUE, bordered = TRUE)
      
      
      # ------------------------------------------------------
      # TRANSCRIPTOMIC PROFILE
      # ------------------------------------------------------
      
      expression_data <- shiny::reactive({
        
        shiny::req(input$expression_file)
        
        data <- mito_read_csv(
          input$expression_file
        )
        
        mito_map_expression_genes(
          data,
          annotation
        )
      })
      
      expression_scores <- shiny::reactive({
        
        data <- expression_data()
        
        mito_calculate_ssgsea(
          data,
          gene_sets,
          component_map
        )
      })
      
      output$expression_status <- shiny::renderText({
        
        shiny::req(input$expression_file)
        
        tryCatch({
          
          data <- expression_data()
          
          scores <- expression_scores()
          
          paste0(
            "Validation: PASSED\n",
            "Mapped genes: ", nrow(data), "\n",
            "Pathways scored: ", nrow(scores), "\n",
            "Sample ID: ",
            if ("Sample_ID" %in% names(data)) {
              paste(unique(data$Sample_ID), collapse = ", ")
            } else {
              "Not provided"
            }
          )
          
        }, error = function(e) {
          paste("Processing error:", conditionMessage(e))
        })
      })
      
      output$expression_pathway_plot <- shiny::renderPlot({
        
        scores <- expression_scores()
        
        mito_plot_pathways(scores)
      })
      
      output$expression_pathway_table <- shiny::renderTable({
        
        scores <- expression_scores()
        
        scores
        
      }, striped = TRUE, bordered = TRUE)
      
      
      # ------------------------------------------------------
      # INTEGRATED PROFILE
      # ------------------------------------------------------
      
      integrated_data <- shiny::reactive({
        
        shiny::req(
          input$integrated_cpg_file,
          input$integrated_expression_file
        )
        
        cpg <- mito_read_csv(
          input$integrated_cpg_file
        )
        
        expression <- mito_read_csv(
          input$integrated_expression_file
        )
        
        cpg_id <- mito_validate_sample_id(cpg)
        
        expression_id <- mito_validate_sample_id(
          expression
        )
        
        cpg_id_normalized <- mito_normalize_sample_id(
          cpg_id
        )
        
        expression_id_normalized <- mito_normalize_sample_id(
          expression_id
        )
        
        if (!identical(
          cpg_id_normalized,
          expression_id_normalized
        )) {
          stop(
            paste0(
              "Sample ID mismatch. Methylation: ",
              cpg_id,
              "; Expression: ",
              expression_id
            )
          )
        }
        
        cpg <- mito_validate_methylation(
          cpg,
          require_sample_id = TRUE
        )
        
        expression <- mito_map_expression_genes(
          expression,
          annotation
        )
        
        scores <- mito_calculate_ssgsea(
          expression,
          gene_sets,
          component_map
        )
        
        list(
          sample_id = cpg_id_normalized,
          cpg = cpg,
          expression = expression,
          scores = scores
        )
      })
      
      output$integrated_status <- shiny::renderText({
        
        shiny::req(
          input$integrated_cpg_file,
          input$integrated_expression_file
        )
        
        tryCatch({
          
          data <- integrated_data()
          
          paste0(
            "Integrated validation: PASSED\n",
            "Matched Sample_ID: ", data$sample_id, "\n",
            "CpGs: ", nrow(data$cpg), "\n",
            "Mapped genes: ", nrow(data$expression), "\n",
            "Pathways scored: ", nrow(data$scores), "\n",
            "Matching: first 16 characters for TCGA barcodes; exact matching otherwise."
          )
          
        }, error = function(e) {
          paste("Integration error:", conditionMessage(e))
        })
      })
      
      output$integrated_pathway_plot <- shiny::renderPlot({
        
        data <- integrated_data()
        
        mito_plot_pathways(
          data$scores
        )
      })
      
      output$integrated_beta_plot <- shiny::renderPlot({
        
        data <- integrated_data()
        
        ggplot2::ggplot(
          data$cpg,
          ggplot2::aes(x = Beta)
        ) +
          ggplot2::geom_histogram(
            bins = 30,
            fill = "#4169E1",
            color = "white"
          ) +
          ggplot2::labs(
            x = "Beta value",
            y = "Number of CpGs",
            title = paste(
              "Methylation profile:",
              data$sample_id
            )
          ) +
          ggplot2::theme_minimal()
      })
      
      output$integrated_association_table <- shiny::renderTable({
        data <- integrated_data()

        result <- robust_pairs[
          robust_pairs$CpG_ID %in% data$cpg$CpG &
            robust_pairs$Gene %in% data$expression$Gene_Symbol,
          ,
          drop = FALSE
        ]

        if (nrow(result) == 0) {
          return(data.frame(
            Message = paste(
              "No uploaded CpG-gene pairs overlap with the TCGA robust associations and the strictly mapped genes in the uploaded expression profile."
            )
          ))
        }

        result$Uploaded_Beta <- data$cpg$Beta[
          match(result$CpG_ID, data$cpg$CpG)
        ]

        result$Uploaded_expression <- data$expression$Expression[
          match(result$Gene, data$expression$Gene_Symbol)
        ]

        result <- result[
          order(abs(result$Pooled_rho), decreasing = TRUE),
          ,
          drop = FALSE
        ]

        result <- head(result, 30)

        data.frame(
          CpG_ID = result$CpG_ID,
          Uploaded_Beta = result$Uploaded_Beta,
          Gene = result$Gene,
          Uploaded_expression = result$Uploaded_expression,
          Mapping_Type = result$Mapping_Type,
          TCGA_Pooled_rho = result$Pooled_rho,
          TCGA_Pooled_FDR = result$Pooled_FDR,
          TCGA_Primary_rho = result$Primary_rho,
          TCGA_Metastatic_rho = result$Metastatic_rho,
          check.names = FALSE
        )

      }, striped = TRUE, bordered = TRUE)
      
    }
  )
}

