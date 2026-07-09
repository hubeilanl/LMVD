options(stringsAsFactors = FALSE)

if (!requireNamespace("vegan", quietly = TRUE)) {
  stop("The R package 'vegan' is required for PERMANOVA.")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The R package 'ggplot2' is required for PCoA plotting.")
}

distance_file <- "bray_curtis_dissimilarity_viral_community.tsv"
metadata_file <- "metadata.txt"
coordinate_output_file <- "pcoa_coordinates_viral_community.tsv"
permanova_output_file <- "permanova_viral_community.tsv"
pcoa_figure_output_file <- "pcoa_viral_community.pdf"
permutations <- 999

read_distance_matrix <- function(path) {
  distance_table <- read.table(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE
  )

  sample_ids <- distance_table[[1]]
  distance_matrix <- as.matrix(distance_table[, -1, drop = FALSE])
  storage.mode(distance_matrix) <- "numeric"
  rownames(distance_matrix) <- sample_ids
  colnames(distance_matrix) <- colnames(distance_table)[-1]

  distance_matrix
}

format_p_value <- function(p_value) {
  if (is.na(p_value)) {
    return("NA")
  }
  if (p_value < 0.001) {
    return("<0.001")
  }
  formatC(p_value, format = "f", digits = 3)
}

choose_group_column <- function(metadata) {
  preferred_columns <- c(
    "Manure_type",
    "manure_type",
    "ManureType",
    "manure",
    "Continent",
    "continent",
    "Country",
    "country",
    "Group",
    "group"
  )
  candidate_columns <- unique(c(
    preferred_columns[preferred_columns %in% colnames(metadata)],
    setdiff(colnames(metadata), "SampleID")
  ))

  for (column in candidate_columns) {
    values <- stats::na.omit(metadata[[column]])
    if (length(unique(values)) > 1) {
      return(column)
    }
  }

  NULL
}

read_sample_metadata <- function(path) {
  metadata <- read.table(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    fill = TRUE
  )

  colnames(metadata) <- trimws(gsub("\r", "", colnames(metadata), fixed = TRUE))

  sample_column <- if ("SampleID" %in% colnames(metadata)) {
    "SampleID"
  } else if ("SRA Run" %in% colnames(metadata)) {
    "SRA Run"
  } else {
    stop("The metadata table must contain either 'SampleID' or 'SRA Run'.")
  }

  metadata$SampleID <- trimws(as.character(metadata[[sample_column]]))

  if (!"Manure_type" %in% colnames(metadata) && "Manure" %in% colnames(metadata)) {
    metadata$Manure_type <- metadata[["Manure"]]
  }

  metadata
}

save_pcoa_figure <- function(plot, output_file, width = 4.0, height = 3.5) {
  grDevices::pdf(output_file, width = width, height = height, useDingbats = FALSE)
  print(plot)
  grDevices::dev.off()
}

distance_matrix <- read_distance_matrix(distance_file)

ordination <- cmdscale(as.dist(distance_matrix), eig = TRUE, k = 2)
positive_eigenvalues <- ordination$eig[ordination$eig > 0]
axis_percent <- ordination$eig[seq_len(2)] / sum(positive_eigenvalues) * 100

coordinates <- data.frame(
  SampleID = rownames(ordination$points),
  PCoA1 = ordination$points[, 1],
  PCoA2 = ordination$points[, 2],
  PCoA1_percent = axis_percent[1],
  PCoA2_percent = axis_percent[2],
  stringsAsFactors = FALSE
)

write.table(
  coordinates,
  file = coordinate_output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

figure_data <- coordinates
group_column <- NULL
permanova_label <- NULL

if (file.exists(metadata_file)) {
  metadata <- read_sample_metadata(metadata_file)

  metadata <- metadata[metadata$SampleID %in% rownames(distance_matrix), , drop = FALSE]
  metadata <- metadata[match(rownames(distance_matrix), metadata$SampleID), , drop = FALSE]

  candidate_terms <- c("Manure_type", "Continent", "Country")
  candidate_terms <- candidate_terms[candidate_terms %in% colnames(metadata)]
  candidate_terms <- candidate_terms[
    vapply(metadata[candidate_terms], function(x) length(unique(stats::na.omit(x))) > 1, logical(1))
  ]

  if (length(candidate_terms) > 0) {
    permanova_results <- lapply(candidate_terms, function(term) {
      model_data <- data.frame(Group = metadata[[term]], stringsAsFactors = FALSE)
      keep <- !is.na(model_data$Group)
      distance_subset <- as.dist(distance_matrix[keep, keep, drop = FALSE])
      model_data <- model_data[keep, , drop = FALSE]

      permanova <- vegan::adonis2(
        distance_subset ~ Group,
        data = model_data,
        permutations = permutations
      )

      data.frame(
        Term = term,
        Df = permanova$Df[1],
        SumOfSqs = permanova$SumOfSqs[1],
        R2 = permanova$R2[1],
        F = permanova$F[1],
        P_value = permanova$`Pr(>F)`[1],
        Permutations = permutations,
        stringsAsFactors = FALSE
      )
    })

    permanova_results <- do.call(rbind, permanova_results)

    write.table(
      permanova_results,
      file = permanova_output_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    first_term <- permanova_results[1, , drop = FALSE]
    permanova_label <- paste0(
      first_term$Term,
      ": PERMANOVA R2 = ",
      formatC(first_term$R2, format = "f", digits = 3),
      ", P = ",
      format_p_value(first_term$P_value)
    )
  }

  group_column <- choose_group_column(metadata)
  figure_data <- merge(coordinates, metadata, by = "SampleID", all.x = TRUE, sort = FALSE)
  figure_data <- figure_data[match(coordinates$SampleID, figure_data$SampleID), , drop = FALSE]
}

axis_labels <- c(
  paste0("PCoA1 (", formatC(axis_percent[1], format = "f", digits = 2), "%)"),
  paste0("PCoA2 (", formatC(axis_percent[2], format = "f", digits = 2), "%)")
)

if (!is.null(group_column)) {
  figure_data$GroupForFigure <- factor(figure_data[[group_column]])
} else {
  figure_data$GroupForFigure <- factor("All samples")
}

group_counts <- table(figure_data$GroupForFigure)
ellipse_groups <- names(group_counts[group_counts >= 4])

pcoa_figure <- ggplot2::ggplot(
  figure_data,
  ggplot2::aes(x = PCoA1, y = PCoA2, colour = GroupForFigure)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey70") +
  ggplot2::geom_point(size = 1.8, alpha = 0.85) +
  ggplot2::labs(
    x = axis_labels[1],
    y = axis_labels[2],
    colour = if (!is.null(group_column)) group_column else NULL,
    title = if (!is.null(permanova_label)) permanova_label else NULL
  ) +
  ggplot2::theme_classic(base_size = 7) +
  ggplot2::theme(
    axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "black"),
    legend.title = ggplot2::element_text(size = 6.2),
    legend.text = ggplot2::element_text(size = 5.8),
    plot.title = ggplot2::element_text(size = 7, face = "bold"),
    panel.grid = ggplot2::element_blank()
  )

if (length(ellipse_groups) > 1) {
  ellipse_data <- figure_data[figure_data$GroupForFigure %in% ellipse_groups, , drop = FALSE]
  pcoa_figure <- pcoa_figure +
    ggplot2::stat_ellipse(
      data = ellipse_data,
      ggplot2::aes(fill = GroupForFigure),
      geom = "polygon",
      level = 0.90,
      linetype = "dashed",
      linewidth = 0.25,
      alpha = 0.15,
      show.legend = FALSE
    )
}

save_pcoa_figure(pcoa_figure, pcoa_figure_output_file)
