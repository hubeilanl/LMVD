options(stringsAsFactors = FALSE)

input_file <- "LMVD_mOTU_abundance_table.txt.gz"
output_file <- "microbial_community_shannon_diversity_by_sample.tsv"
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
abundance_sum <- numeric(length(sample_ids))
abundance_log_sum <- numeric(length(sample_ids))

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
  abundance[!is.finite(abundance)] <- 0
  abundance[abundance < 0] <- 0

  abundance_sum <- abundance_sum + colSums(abundance, na.rm = TRUE)

  positive <- abundance > 0
  abundance_log <- matrix(0, nrow = nrow(abundance), ncol = ncol(abundance))
  abundance_log[positive] <- abundance[positive] * log(abundance[positive])
  abundance_log_sum <- abundance_log_sum + colSums(abundance_log, na.rm = TRUE)
}

shannon <- numeric(length(sample_ids))
detected <- abundance_sum > 0
shannon[detected] <- log(abundance_sum[detected]) -
  abundance_log_sum[detected] / abundance_sum[detected]

results <- data.frame(
  SampleID = sample_ids,
  MicrobialShannon = shannon,
  MicrobialTotalAbundance = abundance_sum,
  stringsAsFactors = FALSE
)

write.table(
  results,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
