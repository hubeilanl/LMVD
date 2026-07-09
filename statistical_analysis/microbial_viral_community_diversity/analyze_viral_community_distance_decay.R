options(stringsAsFactors = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The R package 'ggplot2' is required for distance-decay plotting.")
}

metadata_file <- "metadata.txt"
distance_file <- "bray_curtis_dissimilarity_viral_community.tsv"
pairwise_output_file <- "viral_community_distance_decay_sample_pairs.tsv"
summary_output_file <- "viral_community_distance_decay_summary.tsv"
figure_output_file <- "viral_community_distance_decay.pdf"
figure_max_pairs <- 200000
set.seed(1)

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

haversine_km <- function(latitude_1, longitude_1, latitude_2, longitude_2) {
  radius_km <- 6371.0088
  to_radian <- pi / 180

  lat_1 <- latitude_1 * to_radian
  lat_2 <- latitude_2 * to_radian
  delta_lat <- (latitude_2 - latitude_1) * to_radian
  delta_lon <- (longitude_2 - longitude_1) * to_radian

  a <- sin(delta_lat / 2)^2 +
    cos(lat_1) * cos(lat_2) * sin(delta_lon / 2)^2
  2 * radius_km * asin(pmin(1, sqrt(a)))
}

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

metadata <- read_sample_metadata(metadata_file)
coordinate_matrix <- t(vapply(metadata$Coordinate, parse_coordinate, numeric(2)))
metadata$Latitude <- coordinate_matrix[, "Latitude"]
metadata$Longitude <- coordinate_matrix[, "Longitude"]
metadata <- metadata[is.finite(metadata$Latitude) & is.finite(metadata$Longitude), , drop = FALSE]

distance_matrix <- read_distance_matrix(distance_file)
shared_samples <- intersect(rownames(distance_matrix), metadata$SampleID)
if (length(shared_samples) < 3) {
  stop("At least three samples with valid coordinates are required.")
}

distance_matrix <- distance_matrix[shared_samples, shared_samples, drop = FALSE]
metadata <- metadata[match(shared_samples, metadata$SampleID), , drop = FALSE]

pair_index <- utils::combn(seq_len(nrow(metadata)), 2)
geographic_distance <- haversine_km(
  metadata$Latitude[pair_index[1, ]],
  metadata$Longitude[pair_index[1, ]],
  metadata$Latitude[pair_index[2, ]],
  metadata$Longitude[pair_index[2, ]]
)

pairwise_results <- data.frame(
  SampleID_1 = metadata$SampleID[pair_index[1, ]],
  SampleID_2 = metadata$SampleID[pair_index[2, ]],
  GeographicDistanceKm = geographic_distance,
  Log10GeographicDistanceKmPlusOne = log10(geographic_distance + 1),
  BrayCurtisDissimilarity = distance_matrix[cbind(pair_index[1, ], pair_index[2, ])],
  stringsAsFactors = FALSE
)
pairwise_results$BrayCurtisSimilarity <- 1 - pairwise_results$BrayCurtisDissimilarity

write.table(
  pairwise_results,
  file = pairwise_output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

similarity_test <- suppressWarnings(stats::cor.test(
  pairwise_results$Log10GeographicDistanceKmPlusOne,
  pairwise_results$BrayCurtisSimilarity,
  method = "spearman",
  exact = FALSE
))

dissimilarity_test <- suppressWarnings(stats::cor.test(
  pairwise_results$Log10GeographicDistanceKmPlusOne,
  pairwise_results$BrayCurtisDissimilarity,
  method = "spearman",
  exact = FALSE
))

summary_table <- data.frame(
  Metric = c("BrayCurtisSimilarity", "BrayCurtisDissimilarity"),
  SpearmanRho = c(unname(similarity_test$estimate), unname(dissimilarity_test$estimate)),
  P_value = c(similarity_test$p.value, dissimilarity_test$p.value),
  NumberOfSamples = length(shared_samples),
  NumberOfPairs = nrow(pairwise_results),
  stringsAsFactors = FALSE
)

write.table(
  summary_table,
  file = summary_output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

figure_data <- pairwise_results
if (nrow(figure_data) > figure_max_pairs) {
  figure_data <- figure_data[sample(seq_len(nrow(figure_data)), figure_max_pairs), , drop = FALSE]
}

figure_label <- paste0(
  "Spearman rho = ",
  formatC(unname(similarity_test$estimate), format = "f", digits = 3),
  ", P = ",
  format_p_value(similarity_test$p.value)
)

distance_decay_figure <- ggplot2::ggplot(
  figure_data,
  ggplot2::aes(x = Log10GeographicDistanceKmPlusOne, y = BrayCurtisSimilarity)
) +
  ggplot2::geom_point(size = 0.35, alpha = 0.18, colour = "grey35") +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.5, colour = "black") +
  ggplot2::labs(
    x = "Geographic distance (log10[km + 1])",
    y = "Bray-Curtis similarity",
    title = figure_label
  ) +
  ggplot2::theme_classic(base_size = 7) +
  ggplot2::theme(
    axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "black"),
    plot.title = ggplot2::element_text(size = 7, face = "bold")
  )

grDevices::pdf(figure_output_file, width = 4.0, height = 3.4, useDingbats = FALSE)
print(distance_decay_figure)
grDevices::dev.off()
