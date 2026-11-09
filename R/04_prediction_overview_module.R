###############################################################################
# MODULE 04 — DYNAMIC PREDICTION OVERVIEW
###############################################################################

predictionOverviewUI <- function(id) {

  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::tags$div(
      class = 'prediction-overview-content',
      shiny::uiOutput(ns('overview'))
    )
  )
}


predictionOverviewServer <- function(id, prediction) {

  shiny::moduleServer(id, function(input, output, session) {

    output$overview <- shiny::renderUI({

      pred <- prediction()

      # Initial state before CSV upload
      if (is.null(pred)) {
        return(
          shiny::tags$div(
            class = 'prediction-panel',
            shiny::tags$h4('Prediction overview'),
            shiny::tags$p(
              'Upload a methylation CSV profile to generate a prediction.'
            )
          )
        )
      }

      # Validate prediction fields
      required_fields <- c(
        'predicted_class',
        'probability_immune_low',
        'probability_immune_high',
        'threshold',
        'n_cpgs_used',
        'model_name'
      )

      missing_fields <- setdiff(required_fields, names(pred))

      if (length(missing_fields) > 0) {
        return(
          shiny::tags$div(
            class = 'scientific-disclaimer',
            shiny::tags$strong('Prediction display error'),
            shiny::tags$p(
              paste(
                'Missing prediction fields:',
                paste(missing_fields, collapse = ', ')
              )
            )
          )
        )
      }

      p_low <- as.numeric(pred$probability_immune_low)
      p_high <- as.numeric(pred$probability_immune_high)
      threshold <- as.numeric(pred$threshold)

      if (
        length(p_low) != 1 ||
        length(p_high) != 1 ||
        !is.finite(p_low) ||
        !is.finite(p_high) ||
        p_low < 0 || p_low > 1 ||
        p_high < 0 || p_high > 1 ||
        length(threshold) != 1 ||
        !is.finite(threshold)
      ) {
        return(
          shiny::tags$div(
            class = 'scientific-disclaimer',
            'Prediction contains invalid probability or threshold values.'
          )
        )
      }

      # Format probabilities
      low_percent <- round(p_low * 100, 1)
      high_percent <- round(p_high * 100, 1)

      low_label <- sprintf('%.1f%%', low_percent)
      high_label <- sprintf('%.1f%%', high_percent)

      predicted_class <- as.character(pred$predicted_class)
      model_name <- as.character(pred$model_name)
      n_cpgs <- as.integer(pred$n_cpgs_used)

      class_color <- if (identical(predicted_class, 'Immune-high')) {
        '#4169E1'
      } else {
        '#6A0DAD'
      }

      # Dynamic prediction overview
      shiny::tagList(

        shiny::tags$div(
          class = 'prediction-panel',

          shiny::tags$p(
            style = 'font-size:13px; color:#555F7A; margin-bottom:6px;',
            'PREDICTED IMMUNE PHENOTYPE'
          ),

          shiny::tags$h2(
            predicted_class,
            style = paste0(
              'color:', class_color, ';',
              'font-size:30px; font-weight:800; margin-top:0;'
            )
          ),

          shiny::tags$hr(),

          shiny::tags$h4(
            'Class probabilities',
            style = 'color:#191970; font-weight:700;'
          ),

          # Immune-low probability
          shiny::tags$div(
            style = 'margin-top:16px; margin-bottom:16px;',
            shiny::tags$div(
              style = 'display:flex; justify-content:space-between;',
              shiny::tags$strong('Immune-low'),
              shiny::tags$strong(low_label)
            ),
            shiny::tags$div(
              style = paste0(
                'height:13px; background:#E8EAF2; ',
                'border-radius:8px; overflow:hidden;'
              ),
              shiny::tags$div(
                style = paste0(
                  'height:13px; width:', low_percent, '%;',
                  'background:#6A0DAD; border-radius:8px;'
                )
              )
            )
          ),

          # Immune-high probability
          shiny::tags$div(
            style = 'margin-top:16px; margin-bottom:16px;',
            shiny::tags$div(
              style = 'display:flex; justify-content:space-between;',
              shiny::tags$strong('Immune-high'),
              shiny::tags$strong(high_label)
            ),
            shiny::tags$div(
              style = paste0(
                'height:13px; background:#E8EAF2; ',
                'border-radius:8px; overflow:hidden;'
              ),
              shiny::tags$div(
                style = paste0(
                  'height:13px; width:', high_percent, '%;',
                  'background:#4169E1; border-radius:8px;'
                )
              )
            )
          ),

          shiny::tags$hr(),

          shiny::tags$div(
            style = 'display:flex; flex-wrap:wrap; gap:24px;',
            shiny::tags$div(
              shiny::tags$small('Classification threshold'),
              shiny::tags$p(
                style = 'font-size:18px; font-weight:700; color:#191970;',
                sprintf('%.2f', threshold)
              )
            ),
            shiny::tags$div(
              shiny::tags$small('CpG features used'),
              shiny::tags$p(
                style = 'font-size:18px; font-weight:700; color:#191970;',
                format(n_cpgs, big.mark = ',')
              )
            ),
            shiny::tags$div(
              shiny::tags$small('Model'),
              shiny::tags$p(
                style = 'font-size:18px; font-weight:700; color:#191970;',
                model_name
              )
            )
          )
        ),

        shiny::tags$p(
          class = 'scientific-disclaimer',
          'This result estimates a transcriptionally defined immune phenotype. It is not a prediction of clinical immunotherapy response, survival, or treatment benefit.'
        )
      )
    })

    invisible(NULL)
  })
}

###############################################################################
