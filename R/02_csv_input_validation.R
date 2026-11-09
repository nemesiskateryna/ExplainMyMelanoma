
validate_methylation_csv <- function(file, mlp) {

  # Check uploaded file
  if (length(file) != 1 || !file.exists(file)) {
    stop("CSV file does not exist.")
  }

  # Read CSV
  input_data <- tryCatch(
    read.csv(
      file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    error = function(e) {
      stop("Unable to read CSV: ", conditionMessage(e))
    }
  )

  # Check empty file
  if (nrow(input_data) == 0) {
    stop("CSV file contains no data rows.")
  }

  # Check required columns
  required_columns <- c("CpG", "Beta")

  if (!all(required_columns %in% colnames(input_data))) {
    stop(
      "CSV must contain two required columns: CpG and Beta."
    )
  }

  # Keep required columns
  input_data <- input_data[, required_columns, drop = FALSE]

  # Normalize CpG identifiers
  input_data$CpG <- trimws(as.character(input_data$CpG))

  if (any(is.na(input_data$CpG) | input_data$CpG == "")) {
    stop("CpG identifiers contain missing or empty values.")
  }

  # Check duplicated CpGs
  if (anyDuplicated(input_data$CpG)) {
    stop("CSV contains duplicated CpG identifiers.")
  }

  # Convert Beta values to numeric
  beta_values <- suppressWarnings(
    as.numeric(as.character(input_data$Beta))
  )

  if (any(!is.finite(beta_values))) {
    stop("Beta values contain missing or non-numeric values.")
  }

  if (any(beta_values < 0 | beta_values > 1)) {
    stop("Beta values must be within the interval [0, 1].")
  }

  input_data$Beta <- beta_values

  # Check required model features
  required_cpgs <- mlp$selected_cpgs

  missing_cpgs <- setdiff(
    required_cpgs,
    input_data$CpG
  )

  if (length(missing_cpgs) > 0) {
    stop(
      paste0(
        "Missing ",
        length(missing_cpgs),
        " required CpGs. First missing IDs: ",
        paste(head(missing_cpgs, 5), collapse = ", ")
      )
    )
  }

  extra_cpgs <- setdiff(
    input_data$CpG,
    required_cpgs
  )

  # Run the existing MLP inference function
  prediction <- predict_mlp_profile(
    input_data,
    mlp
  )

  # Return structured result
  list(
    status = "success",
    prediction = prediction,
    input_summary = list(
      n_rows_uploaded = nrow(input_data),
      n_required_cpgs = length(required_cpgs),
      n_extra_cpgs = length(extra_cpgs),
      n_cpgs_used = prediction$n_cpgs_used,
      beta_min = min(input_data$Beta),
      beta_max = max(input_data$Beta)
    ),
    input_data = input_data
  )
}

