options(stringsAsFactors = FALSE)

required_packages <- c("data.table", "vegan")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Required R package(s) are not available: ",
    paste(missing_packages, collapse = ", ")
  )
}

input_file <- "LMVD_vOTU_abundance_table.txt.gz"
output_file <- "bray_curtis_dissimilarity_viral_community.tsv"

command <- if (grepl("\\.gz$", input_file)) {
  paste("gzip -dc", shQuote(input_file))
} else {
  NULL
}

abundance_table <- if (is.null(command)) {
  data.table::fread(
    input_file,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    data.table = FALSE
  )
} else {
  data.table::fread(
    cmd = command,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    data.table = FALSE
  )
}

feature_ids <- abundance_table[[1]]
abundance <- as.matrix(abundance_table[, -1, drop = FALSE])
storage.mode(abundance) <- "numeric"
rownames(abundance) <- make.unique(as.character(feature_ids))
abundance[!is.finite(abundance)] <- 0
abundance[abundance < 0] <- 0

sample_by_feature <- t(abundance)
bray_curtis <- vegan::vegdist(sample_by_feature, method = "bray")
distance_matrix <- as.matrix(bray_curtis)

results <- data.frame(
  SampleID = rownames(distance_matrix),
  distance_matrix,
  check.names = FALSE
)

write.table(
  results,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
