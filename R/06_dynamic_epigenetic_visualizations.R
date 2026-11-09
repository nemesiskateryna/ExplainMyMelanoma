# MODULE 06 — DYNAMIC EPIGENETIC VISUALIZATIONS
#
# Three visualizations based on the uploaded methylation profile:
# 1. MLP-selected CpG methylation heatmap
# 2. MLP-selected CpG beta-value distribution
# 3. Uploaded sample projected onto the fixed TCGA training PCA reference
#
# Does not retrain the model or modify training data.
###############################################################################

dynamicEpigeneticUI <- function(id) {

  ns <- shiny::NS(id)

  shiny::tagList(

    shiny::h4("Top 50 CpGs ranked by permutation importance"),

    shiny::p(
      "5-fold out-of-fold MLP importance; colour indicates 
uploaded-sample beta value."
    ),

    shiny::plotOutput(
      ns("heatmap"),
      height = "850px"
    ),

    shiny::hr(),

    shiny::h4("Beta-value distribution"),

    shiny::p(
      "Distribution of beta values across the 5,000 CpG features ",
      "used by the MLP model."
    ),

    shiny::plotOutput(
      ns("distribution"),
      height = "350px"
    ),

    shiny::hr(),

    shiny::h4("TCGA-SKCM reference PCA"),

    shiny::p(
      "The uploaded sample is projected into the original PCA space ",
      "fitted on the 368 TCGA-SKCM training samples. ",
      "Training coordinates remain unchanged."
    ),

    shiny::plotOutput(
      ns("pca"),
      height = "480px"
    ),

    shiny::p(
      class = "scientific-disclaimer",
      "PCA position describes similarity in the model's selected ",
      "methylation feature space. It is not a clinical classification ",
      "or an independent prediction."
    )
  )
}


dynamicEpigeneticServer <- function(id, profile, mlp, importance) {

  shiny::moduleServer(id, function(input, output, session) {

    # Validate and extract the 5,000 MLP features -----------------------------

    selected_profile <- shiny::reactive({

      uploaded <- profile()

      shiny::validate(
        shiny::need(
          !is.null(uploaded),
          "Upload a methylation CSV profile to display the visualizations."
        )
      )

      shiny::validate(
        shiny::need(
          is.data.frame(uploaded),
          "Uploaded profile must be a data frame."
        )
      )

      shiny::validate(
        shiny::need(
          all(c("CpG", "Beta") %in% colnames(uploaded)),
          "Uploaded profile must contain CpG and Beta columns."
        )
      )

      selected_cpgs <- mlp$selected_cpgs

      shiny::validate(
        shiny::need(
          length(selected_cpgs) == 5000,
          "The MLP model does not contain the expected 5,000 CpGs."
        )
      )

      positions <- match(selected_cpgs, uploaded$CpG)

      shiny::validate(
        shiny::need(
          !anyNA(positions),
          "Some required MLP CpGs are missing from the uploaded profile."
        )
      )

      beta <- as.numeric(uploaded$Beta[positions])

      shiny::validate(
        shiny::need(
          all(is.finite(beta)) && all(beta >= 0 & beta <= 1),
          "Invalid beta values detected in the selected CpG profile."
        )
      )

      data.frame(
        feature_index = seq_along(selected_cpgs),
        CpG = selected_cpgs,
        Beta = beta,
        stringsAsFactors = FALSE
      )
    })


    # 1. TOP 50 CpGs BY PERMUTATION IMPORTANCE -------------------------------

    output$heatmap <- shiny::renderPlot({

      dat <- selected_profile()
      shiny::req(dat)
      shiny::req(importance)
      shiny::req(importance$importance_summary)

      importance_table <- importance$importance_summary

      required_columns <- c("CpG", "Importance_Rank")
      shiny::validate(
        shiny::need(
          all(required_columns %in% colnames(importance_table)),
          "Permutation importance data are incomplete."
        )
      )

      dat$CpG <- as.character(dat$CpG)
      importance_table$CpG <- as.character(importance_table$CpG)

      importance_table <- importance_table[
        order(importance_table$Importance_Rank),
        ,
        drop = FALSE
      ]

      ranked_cpgs <- unique(importance_table$CpG)
      ranked_cpgs <- ranked_cpgs[!is.na(ranked_cpgs)]

      top_cpgs <- ranked_cpgs[ranked_cpgs %in% dat$CpG]
      top_cpgs <- head(top_cpgs, 50)

      shiny::validate(
        shiny::need(
          length(top_cpgs) == 50,
          "Fewer than 50 ranked CpGs are available in the uploaded profile."
        )
      )

      plot_data <- dat[dat$CpG %in% top_cpgs, , drop = FALSE]
      plot_data$CpG <- factor(
        plot_data$CpG,
        levels = rev(top_cpgs)
      )

      ggplot2::ggplot(
        plot_data,
        ggplot2::aes(x = 1, y = CpG, fill = Beta)
      ) +
        ggplot2::geom_tile(
          width = 0.8,
          height = 0.9,
          colour = "white",
          linewidth = 0.2
        ) +
        ggplot2::scale_fill_gradientn(
          colours = c("#191970", "#4169E1", "#D4AF37"),
          limits = c(0, 1),
          name = "Beta value"
        ) +
        ggplot2::scale_x_continuous(
          breaks = NULL,
          expand = ggplot2::expansion(mult = c(0.02, 0.02))
        ) +
        ggplot2::scale_y_discrete(
          drop = FALSE
        ) +
        ggplot2::labs(
          title = "Top 50 CpGs ranked by permutation importance",
          subtitle = paste(
            "5-fold out-of-fold MLP importance;",
            "colour indicates uploaded-sample beta value"
          ),
          x = NULL,
          y = "CpG site"
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(
            size = 8,
            colour = "#191970"
          ),
          axis.title.y = ggplot2::element_text(
            colour = "#191970",
            face = "bold"
          ),
          plot.title = ggplot2::element_text(
            colour = "#191970",
            face = "bold"
          ),
          plot.subtitle = ggplot2::element_text(
            colour = "#4169E1"
          ),
          panel.grid = ggplot2::element_blank(),
          legend.position = "right"
        )
    })

    # 2. BETA-VALUE DISTRIBUTION -----------------------------------------------

    output$distribution <- shiny::renderPlot({

      dat <- selected_profile()

      ggplot2::ggplot(
        dat,
        ggplot2::aes(x = Beta)
      ) +

        ggplot2::geom_histogram(
          bins = 40,
          fill = "#4169E1",
          colour = "white",
          linewidth = 0.25
        ) +

        ggplot2::scale_x_continuous(
          breaks = seq(0, 1, by = 0.1)
        ) +

        ggplot2::coord_cartesian(
          xlim = c(0, 1)
        ) +

        ggplot2::labs(
          x = "DNA methylation beta value",
          y = "Number of CpG features",
          title = "Distribution of MLP-selected CpG beta values",
          subtitle = paste(
            "Features:", format(nrow(dat), big.mark = ",")
          )
        ) +

        ggplot2::theme_minimal(base_size = 13) +

        ggplot2::theme(
          plot.title = ggplot2::element_text(
            colour = "#191970",
            face = "bold"
          ),
          panel.grid.minor = ggplot2::element_blank()
        )

    })


    # 3. FIXED TCGA REFERENCE PCA ----------------------------------------------

    output$pca <- shiny::renderPlot({

      dat <- selected_profile()

      pca_ref <- mlp$pca

      shiny::validate(
        shiny::need(
          !is.null(pca_ref$x) &&
            !is.null(pca_ref$rotation) &&
            !is.null(pca_ref$center) &&
            !is.null(pca_ref$scale),
          "Saved PCA reference is incomplete."
        )
      )

      shiny::validate(
        shiny::need(
          nrow(pca_ref$x) == length(mlp$training_labels) &&
            nrow(pca_ref$x) == length(mlp$training_sample_ids),
          "PCA reference and training metadata are inconsistent."
        )
      )

      shiny::validate(
        shiny::need(
          identical(
            rownames(pca_ref$rotation),
            mlp$selected_cpgs
          ),
          "PCA feature order does not match the MLP feature order."
        )
      )

      shiny::validate(
        shiny::need(
          all(c("PC1", "PC2") %in% colnames(pca_ref$x)),
          "The saved PCA reference does not contain PC1 and PC2."
        )
      )

      # Preserve the exact original PCA transformation
      beta <- dat$Beta

      X <- matrix(
        beta,
        nrow = 1,
        dimnames = list(
          "Uploaded sample",
          mlp$selected_cpgs
        )
      )

      X_scaled <- scale(
        X,
        center = pca_ref$center,
        scale = pca_ref$scale
      )

      new_scores <- X_scaled %*% pca_ref$rotation

      shiny::validate(
        shiny::need(
          all(is.finite(new_scores)),
          "The uploaded sample could not be projected into PCA space."
        )
      )

      training_scores <- data.frame(
        PC1 = pca_ref$x[, "PC1"],
        PC2 = pca_ref$x[, "PC2"],
        Immune_phenotype = factor(
          mlp$training_labels,
          levels = c("Immune-low", "Immune-high")
        ),
        Sample = "TCGA training sample"
      )

      new_sample <- data.frame(
        PC1 = as.numeric(new_scores[1, "PC1"]),
        PC2 = as.numeric(new_scores[1, "PC2"])
      )

      ggplot2::ggplot(
        training_scores,
        ggplot2::aes(
          x = PC1,
          y = PC2,
          colour = Immune_phenotype
        )
      ) +

        ggplot2::geom_point(
          alpha = 0.65,
          size = 2.1
        ) +

        ggplot2::geom_point(
          data = new_sample,
          ggplot2::aes(x = PC1, y = PC2),
          inherit.aes = FALSE,
          shape = 23,
          size = 5,
          stroke = 1.2,
          fill = "#D4AF37",
          colour = "#191970"
        ) +

        ggplot2::scale_colour_manual(
          values = c(
            "Immune-low" = "#6A0DAD",
            "Immune-high" = "#4169E1"
          )
        ) +

        ggplot2::labs(
          x = "PC1",
          y = "PC2",
          colour = "Training phenotype",
          title = "Uploaded sample in TCGA-SKCM PCA space",
          subtitle = "Reference: 368 original TCGA-SKCM training samples"
        ) +

        ggplot2::theme_minimal(base_size = 13) +

        ggplot2::theme(
          plot.title = ggplot2::element_text(
            colour = "#191970",
            face = "bold"
          ),
          legend.position = "bottom",
          panel.grid.minor = ggplot2::element_blank()
        ) +

        ggplot2::geom_label(
          data = new_sample,
          ggplot2::aes(x = PC1, y = PC2),
          inherit.aes = FALSE,
          label = "Uploaded sample",
          nudge_y = 0.06 * diff(range(training_scores$PC2)),
          size = 3,
          colour = "#191970",
          fill = "#FFF9E8"
        )

    })

    invisible(NULL)
  })
}

###############################################################################
# END MODULE 06
###############################################################################"
