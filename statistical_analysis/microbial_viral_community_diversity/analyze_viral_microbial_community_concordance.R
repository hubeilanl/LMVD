options(stringsAsFactors = FALSE)

if (!requireNamespace("vegan", quietly = TRUE)) {
  stop("The R package 'vegan' is required for Procrustes analysis.")
}

viral_distance_file <- "bray_curtis_dissimilarity_viral_community.tsv"
microbial_distance_file <- "bray_curtis_dissimilarity_microbial_community.tsv"
permutations <- 999

read_distance_matrix <- function(distance_file) {
  distance_table <- read.table(
    distance_file,
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

perform_pcoa <- function(distance_matrix) {
  cmdscale(as.dist(distance_matrix), eig = TRUE, k = 2)$points
}

if (!file.exists(viral_distance_file) || !file.exists(microbial_distance_file)) {
  stop("Both viral and microbial Bray-Curtis distance matrices are required.")
}

viral_distance <- read_distance_matrix(viral_distance_file)
microbial_distance <- read_distance_matrix(microbial_distance_file)

shared_samples <- intersect(rownames(viral_distance), rownames(microbial_distance))
if (length(shared_samples) < 3) {
  stop("At least three shared samples are required for Procrustes analysis.")
}

viral_distance <- viral_distance[shared_samples, shared_samples, drop = FALSE]
microbial_distance <- microbial_distance[shared_samples, shared_samples, drop = FALSE]

viral_pcoa <- perform_pcoa(viral_distance)
microbial_pcoa <- perform_pcoa(microbial_distance)

procrustes_fit <- vegan::procrustes(viral_pcoa, microbial_pcoa, symmetric = TRUE)
protest_fit <- vegan::protest(
  viral_pcoa,
  microbial_pcoa,
  permutations = permutations,
  symmetric = TRUE
)

summary_table <- data.frame(
  NumberOfSamples = length(shared_samples),
  ProcrustesCorrelation = unname(protest_fit$t0),
  P_value = protest_fit$signif,
  SumOfSquares = procrustes_fit$ss,
  Permutations = permutations,
  stringsAsFactors = FALSE
)

write.table(
  summary_table,
  file = "viral_microbial_procrustes_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

coordinates <- data.frame(
  SampleID = shared_samples,
  Viral_PCoA1 = viral_pcoa[shared_samples, 1],
  Viral_PCoA2 = viral_pcoa[shared_samples, 2],
  Microbial_PCoA1 = microbial_pcoa[shared_samples, 1],
  Microbial_PCoA2 = microbial_pcoa[shared_samples, 2],
  Procrustes_Microbial_PCoA1 = procrustes_fit$Yrot[shared_samples, 1],
  Procrustes_Microbial_PCoA2 = procrustes_fit$Yrot[shared_samples, 2],
  stringsAsFactors = FALSE
)

write.table(
  coordinates,
  file = "viral_microbial_procrustes_coordinates.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
