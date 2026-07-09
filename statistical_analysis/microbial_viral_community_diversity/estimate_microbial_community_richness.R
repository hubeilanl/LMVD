options(stringsAsFactors = FALSE)

input_file <- "LMVD_mOTU_abundance_table.txt.gz"
output_file <- "microbial_community_richness_by_sample.tsv"
chunk_lines <- 500

open_abundance_matrix <- function(path) {
  if (grepl("\\.gz$", path)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
}

con <- open_abundance_matrix(input_file)
on.exit(close(con), add = TRUE)

header <- readLines(con, n = 1)
if (length(header) != 1) {
  stop("The microbial abundance table is empty: ", input_file)
}

fields <- strsplit(header, "\t", fixed = TRUE)[[1]]
sample_ids <- fields[-1]
richness <- numeric(length(sample_ids))

repeat {
  lines <- readLines(con, n = chunk_lines)
  if (length(lines) == 0) {
    break
  }

  block <- read.table(
    text = lines,
    sep = "\t",
    header = FALSE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    colClasses = c("character", rep("numeric", length(sample_ids)))
  )

  if (ncol(block) != length(sample_ids) + 1) {
    stop("Unexpected column number in ", input_file)
  }

  abundance <- as.matrix(block[, -1, drop = FALSE])
  richness <- richness + colSums(abundance > 0, na.rm = TRUE)
}

results <- data.frame(
  SampleID = sample_ids,
  MicrobialRichness = as.integer(richness),
  stringsAsFactors = FALSE
)

write.table(
  results,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
