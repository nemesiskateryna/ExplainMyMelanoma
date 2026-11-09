melanomaUploadUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::tags$div(
      class = 'upload-panel',
      shiny::tags$h3('Upload methylation profile'),
      shiny::tags$p('Upload a CSV file containing CpG identifiers and Beta values.'),
      shiny::fileInput(
        inputId = ns('profile_file'),
        label = 'Select methylation CSV file',
        accept = c('.csv', 'text/csv'),
        multiple = FALSE
      ),
      shiny::tags$div(
        class = 'csv-format-info',
        shiny::tags$strong('Required format:'),
        shiny::tags$p('Two columns: CpG and Beta.'),
        shiny::tags$p('One CpG identifier per row.'),
        shiny::tags$p('Beta values must be between 0 and 1.')
      )
    ),
    shiny::uiOutput(ns('status')),
    shiny::tableOutput(ns('summary')),
    shiny::tags$p(
      class = 'scientific-disclaimer',
      'Research-use prediction of a transcriptionally defined immune phenotype. This is not a prediction of clinical immunotherapy response and is not a clinical decision-support tool.'
    )
  )
}

melanomaUploadServer <- function(id, mlp) {
  shiny::moduleServer(id, function(input, output, session) {

    result <- shiny::reactiveVal(NULL)
    error_message <- shiny::reactiveVal(NULL)

    shiny::observeEvent(input$profile_file, {
      result(NULL)
      error_message(NULL)

      uploaded_file <- input$profile_file$datapath

      tryCatch({
        validated <- validate_methylation_csv(uploaded_file, mlp)

        result(list(
          prediction = validated$prediction,
          input_summary = validated$input_summary,
          input_data = validated$input_data
        ))
      }, error = function(e) {
        error_message(conditionMessage(e))
      })
    }, ignoreInit = TRUE)

    output$status <- shiny::renderUI({
      if (!is.null(error_message())) {
        return(shiny::tags$div(
          class = 'alert alert-danger',
          shiny::tags$strong('Input validation failed: '),
          error_message()
        ))
      }

      if (!is.null(result())) {
        return(shiny::tags$div(
          class = 'alert alert-success',
          'CSV successfully validated. MLP prediction completed.'
        ))
      }

      shiny::tags$div(
        class = 'alert alert-info',
        'Waiting for a methylation CSV file.'
      )
    })

    output$prediction <- shiny::renderUI({
      current <- result()
      if (is.null(current)) return(NULL)

      prediction <- current$prediction

      shiny::tags$div(
        class = 'prediction-panel',
        shiny::tags$h3('MLP prediction'),
        shiny::tags$h2(prediction$predicted_class),
        shiny::tags$p(sprintf(
          'Immune-high probability: %.2f%%',
          100 * prediction$probability_immune_high
        )),
        shiny::tags$p(sprintf(
          'Immune-low probability: %.2f%%',
          100 * prediction$probability_immune_low
        )),
        shiny::tags$p(paste(
          'Classification threshold:',
          prediction$threshold
        ))
      )
    })

    output$summary <- shiny::renderTable({
      current <- result()
      if (is.null(current)) return(NULL)

      s <- current$input_summary

      data.frame(
        Metric = c(
          'Uploaded rows',
          'Required CpGs',
          'CpGs used by MLP',
          'Additional CpGs',
          'Minimum Beta',
          'Maximum Beta'
        ),
        Value = c(
          s$n_rows_uploaded,
          s$n_required_cpgs,
          s$n_cpgs_used,
          s$n_extra_cpgs,
          round(s$beta_min, 5),
          round(s$beta_max, 5)
        )
      )
    }, rownames = FALSE)

    return(list(
      result = shiny::reactive(result()),
      error = shiny::reactive(error_message())
    ))
  })
}
