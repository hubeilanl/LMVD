options(stringsAsFactors = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The R package 'ggplot2' is required for latitudinal-gradient plotting.")
}

metadata_file <- "metadata.txt"
richness_file <- "viral_community_richness_by_sample.tsv"
analysis_input_file <- "viral_latitudinal_richness_input.tsv"
summary_output_file <- "viral_latitudinal_richness_model_summary.tsv"
figure_output_file <- "viral_latitudinal_richness_gradient.pdf"

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

parse_coordinate <- function(coordinate) {
  coordinate <- trimws(as.character(coordinate))
  coordinate <- gsub("\r", "", coordinate, fixed = TRUE)
  coordinate <- gsub("\\s+", "", coordinate)

  if (is.na(coordinate) || coordinate == "") {
    return(c(Latitude = NA_real_, Longitude = NA_real_))
  }

  pattern <- "^([+-]?[0-9]+(?:[.][0-9]+)?)([NnSs])([+-]?[0-9]+(?:[.][0-9]+)?)([EeWw])$"
  match <- regexec(pattern, coordinate, perl = TRUE)
  parts <- regmatches(coordinate, match)[[1]]

  if (length(parts) != 5) {
    return(c(Latitude = NA_real_, Longitude = NA_real_))
  }

  latitude <- as.numeric(parts[2])
  longitude <- as.numeric(parts[4])
  lat_hemi <- toupper(parts[3])
  lon_hemi <- toupper(parts[5])

  if (!is.finite(latitude) || !is.finite(longitude)) {
    return(c(Latitude = NA_real_, Longitude = NA_real_))
  }

  if (latitude >= 0 && lat_hemi == "S") {
    latitude <- -latitude
  }
  if (longitude >= 0 && lon_hemi == "W") {
    longitude <- -longitude
  }

  if (abs(latitude) > 90 || abs(longitude) > 180) {
    return(c(Latitude = NA_real_, Longitude = NA_real_))
  }

  c(Latitude = latitude, Longitude = longitude)
}

model_p_value <- function(model) {
  fstat <- summary(model)$fstatistic
  if (is.null(fstat)) {
    return(NA_real_)
  }
  stats::pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
}

model_aicc <- function(model) {
  n <- stats::nobs(model)
  k <- length(stats::coef(model))
  if (n <= k + 1) {
    return(NA_real_)
  }
  stats::AIC(model) + (2 * k * (k + 1)) / (n - k - 1)
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

metadata <- read_sample_metadata(metadata_file)
coordinate_matrix <- t(vapply(metadata$Coordinate, parse_coordinate, numeric(2)))
metadata$Latitude <- coordinate_matrix[, "Latitude"]
metadata$Longitude <- coordinate_matrix[, "Longitude"]

richness <- read.table(
  richness_file,
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE
)

if (!all(c("SampleID", "ViralRichness") %in% colnames(richness))) {
  stop("The richness table must contain 'SampleID' and 'ViralRichness' columns.")
}

analysis_data <- merge(
  richness,
  metadata,
  by = "SampleID",
  all = FALSE,
  sort = FALSE
)

analysis_data <- analysis_data[
  is.finite(analysis_data$ViralRichness) &
    is.finite(analysis_data$Latitude) &
    is.finite(analysis_data$Longitude),
  ,
  drop = FALSE
]

analysis_data$AbsoluteLatitude <- abs(analysis_data$Latitude)
analysis_data$Hemisphere <- ifelse(analysis_data$Latitude < 0, "Southern", "Northern")
analysis_data <- analysis_data[analysis_data$AbsoluteLatitude <= 60, , drop = FALSE]

if (nrow(analysis_data) < 5) {
  stop("At least five samples with valid coordinates are required.")
}

write.table(
  analysis_data,
  file = analysis_input_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

linear_model <- stats::lm(ViralRichness ~ AbsoluteLatitude, data = analysis_data)
quadratic_model <- stats::lm(
  ViralRichness ~ poly(AbsoluteLatitude, 2, raw = TRUE),
  data = analysis_data
)

summary_table <- data.frame(
  Model = c("First_order_polynomial", "Second_order_polynomial"),
  R2 = c(summary(linear_model)$r.squared, summary(quadratic_model)$r.squared),
  AdjustedR2 = c(summary(linear_model)$adj.r.squared, summary(quadratic_model)$adj.r.squared),
  P_value = c(model_p_value(linear_model), model_p_value(quadratic_model)),
  AICc = c(model_aicc(linear_model), model_aicc(quadratic_model)),
  NumberOfSamples = nrow(analysis_data),
  stringsAsFactors = FALSE
)

write.table(
  summary_table,
  file = summary_output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

prediction_data <- data.frame(
  AbsoluteLatitude = seq(
    min(analysis_data$AbsoluteLatitude),
    max(analysis_data$AbsoluteLatitude),
    length.out = 200
  )
)
prediction_data$FirstOrderFit <- stats::predict(linear_model, newdata = prediction_data)
prediction_data$SecondOrderFit <- stats::predict(quadratic_model, newdata = prediction_data)

subtitle_text <- paste0(
  "First order: R2 = ",
  formatC(summary_table$R2[1], format = "f", digits = 3),
  ", P = ",
  format_p_value(summary_table$P_value[1]),
  "; second order: R2 = ",
  formatC(summary_table$R2[2], format = "f", digits = 3),
  ", P = ",
  format_p_value(summary_table$P_value[2])
)

latitudinal_figure <- ggplot2::ggplot(
  analysis_data,
  ggplot2::aes(x = AbsoluteLatitude, y = ViralRichness)
) +
  ggplot2::geom_point(size = 0.8, alpha = 0.55, colour = "grey35") +
  ggplot2::geom_line(
    data = prediction_data,
    ggplot2::aes(y = FirstOrderFit),
    linewidth = 0.45,
    linetype = "dashed",
    colour = "#2f6f9f"
  ) +
  ggplot2::geom_line(
    data = prediction_data,
    ggplot2::aes(y = SecondOrderFit),
    linewidth = 0.5,
    colour = "black"
  ) +
  ggplot2::labs(
    x = "Absolute latitude",
    y = "Viral richness",
    title = subtitle_text
  ) +
  ggplot2::theme_classic(base_size = 7) +
  ggplot2::theme(
    axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "black"),
    plot.title = ggplot2::element_text(size = 6.5, face = "bold")
  )

grDevices::pdf(figure_output_file, width = 4.0, height = 3.4, useDingbats = FALSE)
print(latitudinal_figure)
grDevices::dev.off()
