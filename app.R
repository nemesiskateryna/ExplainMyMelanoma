###############################################################################
# ExplainMyMelanoma AI
# Main Shiny application
###############################################################################
# Determine application directory
app_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
# Determine project root
project_root <- app_dir
# Define results directories
tables_dir <- file.path(project_root, "results", "tables")
figures_dir <- file.path(project_root, "results", "figures")

# Validate results directories
if (!dir.exists(tables_dir)) {
  stop("Results tables directory not found: ", tables_dir)
}

if (!dir.exists(figures_dir)) {
  stop("Results figures directory not found: ", figures_dir)
}
# Check packages
required_packages <- c("shiny", "bslib")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}
# Load application modules
source(file.path(app_dir, "R", "01_mlp_inference.R"))
source(file.path(app_dir, "R", "02_csv_input_validation.R"))
source(file.path(app_dir, "R", "03_upload_module.R"))
source(file.path(app_dir, "R", "04_prediction_overview_module.R"))
source(file.path(app_dir, "R", "05_epigenetic_evidence_module.R"))
source(file.path(app_dir, "R", "06_dynamic_epigenetic_visualizations.R"))
source(file.path(app_dir, "R", "07_mitochondrial_immune_spectrum_module.R"))
# Load saved MLP model
model_file <- file.path(
  project_root,
  "results",
  "models",
  "TCGA_SKCM_MLP_final_model.rds"
)
if (!file.exists(model_file)) {
  stop("Saved TCGA-SKCM MLP model not found.")
}
mlp <- readRDS(model_file)
if (length(mlp$selected_cpgs) != 5000) {
  stop("Unexpected number of selected CpG features.")
}

# Load 5-fold out-of-fold MLP permutation importance
importance_file <- file.path(
  project_root,
  "results",
  "models",
  "TCGA_SKCM_MLP_CpG_permutation_importance_5fold.rds"
)

if (!file.exists(importance_file)) {
  stop("MLP permutation importance object not found.")
}

importance_object <- readRDS(importance_file)

if (is.null(importance_object$importance_summary)) {
  stop("Permutation importance summary is missing.")
}

required_importance_columns <- c(
  "CpG",
  "Mean_Delta_LogLoss",
  "Importance_Rank"
)

if (!all(required_importance_columns %in% 
        colnames(importance_object$importance_summary))) {
  stop("Required permutation importance columns are missing.")
}

# Royal palette
royal_blue <- "#4169E1"
gold <- "#D4AF37"
deep_purple <- "#6A0DAD"
navy <- "#191970"
# Application theme
app_theme <- bslib::bs_theme(
  version = 5,
  bg = "#F7F8FC",
  fg = navy,
  primary = royal_blue,
  secondary = deep_purple
)
# USER INTERFACE
ui <- shiny::fluidPage(

  theme = app_theme,

  shiny::tags$head(
    shiny::tags$title("ExplainMyMelanoma AI"),
    shiny::tags$meta(
      name = "viewport",
      content = "width=device-width, initial-scale=1"
    ),
    shiny::tags$style(shiny::HTML("
      body {
        background-color: #F7F8FC;
        font-family: Arial, sans-serif;
      }
      .app-header {
        background: linear-gradient(120deg, #191970, #4169E1);
        color: white;
        padding: 35px 32px;
        border-radius: 16px;
        margin-top: 25px;
        margin-bottom: 25px;
      }
      .app-kicker {
        color: #D4AF37;
        font-size: 12px;
        font-weight: bold;
        letter-spacing: 2px;
        margin-bottom: 12px;
      }
      .app-header h1 {
        font-weight: 800;
        font-size: 34px;
      }
      .app-subtitle {
        color: #E4E9FF;
        font-size: 16px;
      }
      .section-card {
        background: white;
        border: 1px solid #E4E8F2;
        border-radius: 14px;
        padding: 24px;
        margin-bottom: 22px;
        box-shadow: 0 3px 12px rgba(25,25,112,0.04);
      }
      .section-heading {
        color: #191970;
        font-weight: bold;
        margin-bottom: 14px;
      }
      .accent-line {
        width: 55px;
        height: 4px;
        background: #D4AF37;
        border-radius: 5px;
        margin-bottom: 18px;
      }
      .prediction-panel {
        background: #EEF2FF;
        border-left: 5px solid #4169E1;
        border-radius: 10px;
        padding: 22px;
        margin-top: 18px;
      }
      .prediction-panel h2 {
        color: #191970;
        font-weight: 800;
      }
      .scientific-disclaimer {
        background: #FFF9E8;
        border-left: 4px solid #D4AF37;
        border-radius: 8px;
        padding: 16px;
        margin-top: 20px;
        color: #4B4560;
        font-size: 13px;
      }
      .app-footer {
        text-align: center;
        color: #777F95;
        font-size: 12px;
        padding: 25px 10px;
        margin-top: 25px;
        border-top: 1px solid #E4E8F2;
      }
    "))
  ),

  shiny::div(
    class = "container-fluid",

    # HEADER
    shiny::div(
      class = "app-header",
      shiny::div(
        class = "app-kicker",
        "RESEARCH PLATFORM | EPIGENETICS + MACHINE LEARNING"
      ),
      shiny::h1("ExplainMyMelanoma AI"),
      shiny::p(
        class = "app-subtitle",
        "An Explainable Epigenetic Machine Learning Platform for Melanoma Immune Phenotype Prediction"
      )
    ),

    # INTRODUCTION
    shiny::div(
      class = "section-card",
      shiny::h3(class = "section-heading", "About the platform"),
      shiny::div(class = "accent-line"),
      shiny::p(paste(
        "ExplainMyMelanoma AI uses a neural network trained on TCGA-SKCM DNA methylation data to estimate a transcriptionally defined immune phenotype."
      )),
      shiny::p(paste(
        "The model processes 5,000 selected CpG features and applies the original training PCA transformation and MLP model."
      ))
    ),

    # UPLOAD
    shiny::div(
      class = "section-card",
      shiny::h3(class = "section-heading", "Upload methylation profile"),
      shiny::div(class = "accent-line"),
      shiny::p("Upload a CSV file containing two columns: CpG and Beta. Beta values must be between 0 and 1."),
      melanomaUploadUI("upload")
    ),

    # PREDICTION OVERVIEW
    shiny::div(
      class = "section-card",
      shiny::h3(class = "section-heading", "Prediction overview"),
      shiny::div(class = "accent-line"),
      predictionOverviewUI("prediction")
    ),

    # DYNAMIC METHYLATION VISUALIZATIONS
    shiny::div(
      class = "section-card",
      shiny::h3(class = "section-heading", "Dynamic methylation visualizations"),
      shiny::div(class = "accent-line"),
      dynamicEpigeneticUI("dynamic")
    ),
    # EPIGENETIC EVIDENCE
    shiny::div(
      class = "section-card",
      shiny::h3(class = "section-heading", "Epigenetic evidence"),
      shiny::div(class = "accent-line"),
      epigeneticEvidenceUI("epigenetic")
    ),

    # MITOCHONDRIAL–IMMUNE INTEGRATED SPECTRUM
    shiny::div(
      class = "section-card",
      shiny::h3(
        class = "section-heading",
        "Mitochondrial–Immune Integrated Spectrum"
      ),
      shiny::div(class = "accent-line"),
      mitoImmuneSpectrumUI("mito_immune")
    ),

    # SCIENTIFIC INFORMATION
    shiny::div(
      class = "section-card",
      shiny::h3(class = "section-heading", "Scientific interpretation"),
      shiny::div(class = "accent-line"),
      shiny::p(shiny::strong("Model: "), "TCGA-SKCM Multilayer Perceptron (MLP)"),
      shiny::p(shiny::strong("Input: "), "5,000 selected DNA methylation CpG features"),
      shiny::p(shiny::strong("Output: "), "Immune-high or Immune-low phenotype"),
      shiny::p(shiny::strong("Classification threshold: "), "0.5"),
      shiny::p(shiny::strong("Model development: "), "TCGA-SKCM discovery cohort"),
      shiny::p(shiny::strong("GEO DNA methylation validation: "), "GSE144487 — external validation of the CpG signature."),
      shiny::p(shiny::strong("GEO gene-expression analysis: "), "GSE65904 — immune phenotype and TIS-17 analysis."),
      shiny::div(
        class = "scientific-disclaimer",
        shiny::strong("Research-use disclaimer"),
        shiny::p(paste(
          "This platform provides a research-use prediction of a transcriptionally defined immune phenotype. It does not predict clinical response to immunotherapy, patient survival, or treatment benefit."
        )),
        shiny::p(paste(
          "The output is not a clinical diagnosis or a clinical decision-support recommendation."
        ))
      )
    ),

    # FOOTER
    shiny::div(
      class = "app-footer",
      shiny::strong("ExplainMyMelanoma AI"),
      shiny::br(),
      "Research prototype | TCGA-SKCM-derived MLP",
      shiny::br(),
      "Epigenetic machine learning for melanoma research"
    )
  )
)

# SERVER
server <- function(input, output, session) {

  upload_result <- melanomaUploadServer(
    id = "upload",
    mlp = mlp
  )

  predictionOverviewServer(
    id = "prediction",
    prediction = shiny::reactive({
      current_result <- upload_result$result()
      if (is.null(current_result)) {
        return(NULL)
      }
      current_result$prediction
    })
  )

  # Dynamic methylation visualizations
  dynamicEpigeneticServer(
    id = "dynamic",
    profile = shiny::reactive({
      current_result <- upload_result$result()
      if (is.null(current_result)) {
        return(NULL)
      }
      current_result$input_data
    }),
    mlp = mlp,
    importance = importance_object
  )
  # Epigenetic evidence module
  epigeneticEvidenceServer(
    id = "epigenetic",
    tables_dir = tables_dir,
    figures_dir = figures_dir
  )

  # Mitochondrial–Immune Integrated Spectrum
  mitoImmuneSpectrumServer(
    id = "mito_immune",
    project_root = project_root
  )

}

# Shiny application
shiny::shinyApp(
  ui = ui,
  server = server
)

###############################################################################
