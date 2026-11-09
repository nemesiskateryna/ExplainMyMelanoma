###############################################################################
# ExplainMyMelanoma AI — MLP inference
# Uses the previously trained TCGA-SKCM MLP model.
# Does not retrain the model.
###############################################################################

predict_mlp_profile <- function(beta_profile, mlp, threshold = 0.5) {

  if (!is.numeric(threshold) || length(threshold) != 1 ||
      is.na(threshold) || threshold <= 0 || threshold >= 1) {
    stop('Threshold must be a single number between 0 and 1.')
  }

  required_cpgs <- mlp$selected_cpgs

  if (is.data.frame(beta_profile)) {
    if (!all(c('CpG', 'Beta') %in% colnames(beta_profile))) {
      stop("Input data.frame must contain columns 'CpG' and 'Beta'.")
    }
    if (anyDuplicated(beta_profile$CpG)) {
      stop('Input contains duplicated CpG identifiers.')
    }
    if (!is.numeric(beta_profile$Beta)) {
      stop('Beta column must be numeric.')
    }
    values <- beta_profile$Beta
    names(values) <- as.character(beta_profile$CpG)
  } else if (is.numeric(beta_profile) && !is.null(names(beta_profile))) {
    values <- beta_profile
  } else {
    stop('Input must be a named numeric vector or a CpG/Beta data.frame.')
  }

  if (anyDuplicated(names(values))) {
    stop('Input contains duplicated CpG identifiers.')
  }

  missing_cpgs <- setdiff(required_cpgs, names(values))

  if (length(missing_cpgs) > 0) {
    stop(
      paste0(
        'Missing ', length(missing_cpgs),
        ' required CpGs. First missing IDs: ',
        paste(head(missing_cpgs, 10), collapse = ', ')
      )
    )
  }

  values <- values[required_cpgs]

  if (anyNA(values) || any(!is.finite(values))) {
    stop('Beta values contain missing or non-finite values.')
  }

  if (any(values < 0 | values > 1)) {
    stop('Beta values must be within the interval [0, 1].')
  }

  X <- matrix(
    as.numeric(values),
    nrow = 1,
    dimnames = list('Input_sample', required_cpgs)
  )

  X_scaled <- scale(
    X,
    center = mlp$pca$center,
    scale = mlp$pca$scale
  )

  pca_scores <- X_scaled %*% mlp$pca$rotation

  pca_scores_scaled <- scale(
    pca_scores,
    center = mlp$score_means,
    scale = mlp$score_sds
  )

  probabilities <- predict(
    mlp$model,
    newdata = pca_scores_scaled,
    type = 'raw'
  )

  if (!identical(colnames(probabilities), mlp$class_levels)) {
    stop('Prediction class order does not match the model class levels.')
  }

  if (anyNA(probabilities) || any(!is.finite(probabilities))) {
    stop('Model returned invalid probabilities.')
  }

  if (any(abs(rowSums(probabilities) - 1) > 1e-8)) {
    stop('Predicted probabilities do not sum to 1.')
  }

  p_low <- as.numeric(probabilities[1, 'Immune-low'])
  p_high <- as.numeric(probabilities[1, 'Immune-high'])

  predicted_class <- if (p_high >= threshold) {
    'Immune-high'
  } else {
    'Immune-low'
  }

  list(
    predicted_class = predicted_class,
    probability_immune_low = p_low,
    probability_immune_high = p_high,
    threshold = threshold,
    n_cpgs_used = length(required_cpgs),
    model_name = mlp$model_name,
    interpretation = paste(
      'Prediction of transcriptionally defined immune phenotype;',
      'not a prediction of clinical immunotherapy response.'
    )
  )
}

###############################################################################
