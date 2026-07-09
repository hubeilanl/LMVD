library(ggplot2)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) == 0) {
  wd <- normalizePath(getwd(), mustWork = FALSE)
  root <- if (basename(wd) == "pseudomonadota_arg_profile_sensitivity") {
    wd
  } else {
    file.path(wd, "pseudomonadota_arg_profile_sensitivity")
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

write_tsv <- function(df, filename) {
  write.table(df, file.path(results_dir, filename), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}

prop <- read.delim(file.path(data_dir, "pseudomonadota_arg_type_proportions_by_genome.tsv"),
                   stringsAsFactors = FALSE, check.names = FALSE)
sens <- read.delim(file.path(data_dir, "pseudomonadota_short_contig_sensitivity.tsv"),
                   stringsAsFactors = FALSE, check.names = FALSE)

prop$dataset <- factor(prop$dataset, levels = c("AR-MAG", "Isolate"))
type_order <- c("multidrug", "fluoroquinolone antibiotic", "peptide antibiotic",
                "nucleoside antibiotic", "aminocoumarin antibiotic",
                "aminoglycoside antibiotic", "Other")
prop$ARG_type <- factor(prop$ARG_type, levels = type_order)

summary_rows <- aggregate(ARG_type_percent ~ dataset + ARG_type, prop, function(x) {
  c(mean = mean(x), median = median(x),
    Q1 = as.numeric(quantile(x, 0.25)),
    Q3 = as.numeric(quantile(x, 0.75)))
})
summary_out <- data.frame(
  dataset = summary_rows$dataset,
  ARG_type = summary_rows$ARG_type,
  mean_percent = summary_rows$ARG_type_percent[, "mean"],
  median_percent = summary_rows$ARG_type_percent[, "median"],
  Q1_percent = summary_rows$ARG_type_percent[, "Q1"],
  Q3_percent = summary_rows$ARG_type_percent[, "Q3"]
)
write_tsv(summary_out, "pseudomonadota_arg_type_proportion_summary.tsv")

fig4 <- ggplot(prop, aes(x = ARG_type, y = ARG_type_percent, fill = dataset)) +
  geom_boxplot(width = 0.68, outlier.shape = NA, linewidth = 0.25,
               position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = c("AR-MAG" = "#6F8798", Isolate = "#8A9A5B")) +
  labs(x = "ARG type", y = "ARG type proportion per genome (%)") +
  theme_classic(base_size = 7, base_family = "sans") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.title = element_blank(), legend.position = "top")
ggsave(file.path(figures_dir, "pseudomonadota_arg_type_proportion_comparison.pdf"),
       fig4, width = 150 / 25.4, height = 92 / 25.4,
       units = "in",
       device = function(filename, width, height, ...) {
         grDevices::pdf(file = filename, width = width, height = height,
                        useDingbats = FALSE)
       })

sens$dataset <- factor(sens$dataset,
                       levels = c("Pseudomonadota AR-MAG",
                                  "Livestock isolate Pseudomonadota"))
fig5 <- ggplot(
  sens,
  aes(x = cutoff_kb,
      y = mean_ARG_annotations_per_MAG_original_denominator,
      colour = dataset, shape = dataset)
) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.8, stroke = 0.25) +
  scale_colour_manual(values = c("Pseudomonadota AR-MAG" = "#34495E",
                                 "Livestock isolate Pseudomonadota" = "#7A8A43")) +
  scale_x_continuous(breaks = sort(unique(sens$cutoff_kb)),
                     labels = c("No cutoff", "1", "2", "3", "4", "5")) +
  labs(x = "Contig length cutoff (kb)",
       y = "Mean ARG annotations per MAG/genome") +
  theme_classic(base_size = 7, base_family = "sans") +
  theme(legend.title = element_blank(), legend.position = "top")
ggsave(file.path(figures_dir, "pseudomonadota_short_contig_sensitivity.pdf"),
       fig5, width = 150 / 25.4, height = 92 / 25.4,
       units = "in",
       device = function(filename, width, height, ...) {
         grDevices::pdf(file = filename, width = width, height = height,
                        useDingbats = FALSE)
       })

write_tsv(sens, "pseudomonadota_short_contig_sensitivity.tsv")

print(summary_out)
