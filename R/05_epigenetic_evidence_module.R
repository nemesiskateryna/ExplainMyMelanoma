
epigeneticEvidenceUI <- function(id) {

  ns <- shiny::NS(id)

  shiny::tags$details(
    class = "epigenetic-evidence-details",

    shiny::tags$summary(
      "Show epigenetic evidence"
    ),

    shiny::div(
      class = "epigenetic-evidence-content",

      # GEO DNA METHYLATION VALIDATION
      shiny::h4("GEO external CpG validation"),

      shiny::p(
        "GSE144487: external assessment of the TCGA-derived",
        " DNA methylation signature."
      ),

      shiny::tableOutput(ns("geo_validation_table")),

      shiny::hr(),

      # TCGA CpG-GENE ASSOCIATIONS
      shiny::h4("CpG–gene expression associations"),

      shiny::p(
        "Robust associations between selected TCGA-SKCM CpG sites and gene expression."
      ),

      shiny::tableOutput(ns("cpg_gene_table")),

      shiny::hr(),

      # GEO IMMUNE SIGNATURES
      shiny::h4("GSE65904: transcriptomic immune signatures"),

      shiny::p(
        "Association of the TIS-17 gene signature with auxiliary immune-related expression signatures."
      ),

      shiny::tableOutput(ns("tis_table")),

      shiny::hr(),

      shiny::div(
        class = "scientific-disclaimer",
        shiny::strong("Scientific limitation"),
        shiny::p(
          "These analyses describe associations and annotation patterns. ",
          "They do not establish that DNA methylation directly causes ",
          "changes in gene expression or immune phenotype."
        )
      )
    )
  )
}
epigeneticEvidenceServer <- function(
  id,
  tables_dir,
  figures_dir
) {

  shiny::moduleServer(id, function(input, output, session) {

    # Safe CSV reader
    read_result <- function(filename) {

      path <- file.path(tables_dir, filename)

      if (!file.exists(path)) {
        warning("Evidence table not found: ", path)
        return(NULL)
      }

      tryCatch(
        read.csv(
          path,
          stringsAsFactors = FALSE,
          check.names = FALSE
        ),
        error = function(e) {
          warning(
            "Could not read ", filename, ": ",
            conditionMessage(e)
          )
          NULL
        }
      )
    }


    # GEO VALIDATION TABLE
    output$geo_validation_table <- shiny::renderTable({

      dat <- read_result(
        "GEO_external_validation_final_summary.csv"
      )

      if (is.null(dat)) return(NULL)

      dat <- dat[
        dat$Section %in% c(
          "GEO cohort",
          "Auxiliary immune signatures",
          "CpG overlap"
        ),
        c("Section", "Metric", "Value"),
        drop = FALSE
      ]

      dat

    }, striped = TRUE, bordered = FALSE, spacing = "s")


    # CpG-GENE ASSOCIATION TABLE
    output$cpg_gene_table <- shiny::renderTable({

      dat <- read_result(
        "TCGA_SKCM_CpG_gene_correlation_robust_pairs.csv"
      )

      if (is.null(dat)) return(NULL)

      required <- c(
        "CpG_ID", "Gene", "Pooled_rho",
        "Pooled_FDR", "Primary_rho", "Metastatic_rho"
      )

      if (!all(required %in% names(dat))) {
        warning("Unexpected CpG-gene table structure.")
        return(NULL)
      }

      selected <- dat[
        (dat$CpG_ID == "cg16890093" & dat$Gene == "TAP1") |
          (dat$CpG_ID == "cg07016276" & dat$Gene == "HLA-F"),
        required,
        drop = FALSE
      ]

      selected

    }, digits = 4, striped = TRUE, bordered = FALSE, spacing = "s")


    # TIS-17 CORRELATIONS
    output$tis_table <- shiny::renderTable({

      dat <- read_result(
        "GEO_TIS17_auxiliary_correlations.csv"
      )

      if (is.null(dat)) return(NULL)

      dat

    }, digits = 4, striped = TRUE, bordered = FALSE, spacing = "s")


    # SAFE IMAGE RENDERING
    render_evidence_image <- function(filename, description) {

      image_path <- file.path(figures_dir, filename)

      if (!file.exists(image_path)) {
        warning("Evidence figure not found: ", image_path)
        return(NULL)
      }

      list(
        src = normalizePath(image_path),
        contentType = "image/png",
        alt = description,
        width = "100%"
      )
    }


    output$tap1_plot <- shiny::renderImage({
      render_evidence_image(
        "TCGA_SKCM_cg16890093_TAP1_correlation_primary_v2.png",
        "Correlation between cg16890093 methylation and TAP1 expression"
      )
    }, deleteFile = FALSE)


    output$hlaf_plot <- shiny::renderImage({
      render_evidence_image(
        "TCGA_SKCM_cg07016276_HLA-F_correlation_primary_v2.png",
        "Correlation between cg07016276 methylation and HLA-F expression"
      )
    }, deleteFile = FALSE)


    output$region_plot <- shiny::renderImage({
      render_evidence_image(
        "GSE144487_CpG_gene_region_distribution.png",
        "Genomic region distribution of GSE144487 CpGs"
      )
    }, deleteFile = FALSE)


    output$island_plot <- shiny::renderImage({
      render_evidence_image(
        "GSE144487_CpG_island_distribution.png",
        "CpG island distribution in GSE144487"
      )
    }, deleteFile = FALSE)

  })
}
