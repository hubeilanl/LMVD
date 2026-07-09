library(ggplot2)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) == 0) {
  wd <- normalizePath(getwd(), mustWork = FALSE)
  root <- if (basename(wd) == "viral_arg_contig_length_distribution") {
    wd
  } else {
    file.path(wd, "viral_arg_contig_length_distribution")
  }
} else {
  script_file <- sub("^--file=", "", script_arg[1])
  script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
  if (!grepl("^/", script_file)) script_file <- file.path(getwd(), script_file)
  root <- dirname(dirname(normalizePath(script_file, mustWork = FALSE)))
}
data_dir <- file.path(root, "data")
results_dir <- file.path(root, "results")
figures_dir <- file.path(root, "figures")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

raw <- read.delim(file.path(data_dir, "viral_arg_contig_lengths_raw.tsv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
required <- c("Viral_ID", "Viral_length")
missing <- setdiff(required, names(raw))
if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

raw$Viral_length <- as.numeric(raw$Viral_length)
viral <- unique(raw[!is.na(raw$Viral_ID) & is.finite(raw$Viral_length),
                    c("Viral_ID", "Viral_length")])
viral$length_kb <- viral$Viral_length / 1000

breaks <- c(5, 10, 20, 30, 40, 50, 100, Inf)
labels <- c("5-10", "10-20", "20-30", "30-40", "40-50", "50-100", ">=100")
viral$length_bin <- cut(viral$length_kb, breaks = breaks, labels = labels,
                        right = FALSE, include.lowest = TRUE)
bins <- as.data.frame(table(factor(viral$length_bin, levels = labels)))
names(bins) <- c("length_bin", "n")
bins$proportion <- bins$n / sum(bins$n)
bins$percent <- bins$proportion * 100
write.csv(bins, file.path(results_dir, "viral_arg_contig_length_bins.csv"),
          row.names = FALSE)

summary_df <- data.frame(
  metric = c("raw_arg_records", "unique_viral_contigs",
             "duplicated_arg_records_removed", "minimum_length_kb",
             "median_length_kb", "mean_length_kb", "maximum_length_kb",
             "proportion_lt10kb", "proportion_lt20kb", "proportion_ge20kb",
             "proportion_ge50kb"),
  value = c(nrow(raw), nrow(viral), nrow(raw) - nrow(viral),
            min(viral$length_kb), median(viral$length_kb),
            mean(viral$length_kb), max(viral$length_kb),
            mean(viral$length_kb < 10), mean(viral$length_kb < 20),
            mean(viral$length_kb >= 20), mean(viral$length_kb >= 50))
)
write.csv(summary_df, file.path(results_dir, "viral_arg_contig_length_summary.csv"),
          row.names = FALSE)
write.csv(viral, file.path(data_dir, "unique_viral_arg_contig_lengths.csv"),
          row.names = FALSE)

bins$label <- paste0(sprintf("%.1f", bins$percent), "%\n", "n=", bins$n)
bins$fill_group <- ifelse(bins$length_bin == "5-10", "5-10 kb", "Other")

p <- ggplot(bins, aes(x = length_bin, y = percent, fill = fill_group)) +
  geom_col(width = 0.76, colour = "white", linewidth = 0.25) +
  geom_text(aes(label = label), vjust = -0.22, size = 2.05,
            lineheight = 0.88) +
  scale_fill_manual(values = c("5-10 kb" = "#E28E2C", Other = "#4E79A7"),
                    guide = "none") +
  scale_y_continuous(limits = c(0, max(bins$percent) * 1.24),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Viral contig length (kb)",
       y = "Proportion of viral contigs (%)") +
  theme_classic(base_size = 7, base_family = "sans")

ggsave(file.path(figures_dir, "viral_arg_contig_length_distribution.pdf"),
       p, width = 89 / 25.4, height = 68 / 25.4,
       units = "in",
       device = function(filename, width, height, ...) {
         grDevices::pdf(file = filename, width = width, height = height,
                        useDingbats = FALSE)
       })

print(summary_df)
