library(ggplot2)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) == 0) {
  wd <- normalizePath(getwd(), mustWork = FALSE)
  root <- if (basename(wd) == "votu_prevalence_abundance") {
    wd
  } else {
    file.path(wd, "votu_prevalence_abundance")
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

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

source_file <- file.path(data_dir, "votu_prevalence_abundance_source_data.tsv")
df <- read.delim(source_file, stringsAsFactors = FALSE, check.names = FALSE)

required <- c("vOTU", "type", "prevalence_percent", "mean_abundance_present_samples")
missing <- setdiff(required, names(df))
if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

df$prevalence_percent <- as.numeric(df$prevalence_percent)
df$mean_abundance_present_samples <- as.numeric(df$mean_abundance_present_samples)
df <- df[is.finite(df$prevalence_percent) &
           is.finite(df$mean_abundance_present_samples), , drop = FALSE]

labels <- c(chicken = "Chicken", swine = "Swine", cattle = "Cattle")
colours <- c(chicken = "#4E79A7", swine = "#A05D56", cattle = "#5F8F5B")
df$type <- factor(df$type, levels = names(labels), labels = labels)

stats <- do.call(rbind, lapply(levels(df$type), function(label) {
  sub <- df[df$type == label, , drop = FALSE]
  fit <- lm(mean_abundance_present_samples ~ prevalence_percent, data = sub)
  fs <- summary(fit)
  data.frame(
    type = label,
    n_vOTU_prevalence_gt_1pct = nrow(sub),
    intercept = unname(coef(fit)[1]),
    slope_per_1_percent_prevalence = unname(coef(fit)[2]),
    r_squared = fs$r.squared,
    adjusted_r_squared = fs$adj.r.squared,
    slope_p_value = coef(fs)[2, "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
}))
write.table(stats, file.path(results_dir, "votu_prevalence_abundance_lm_statistics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

theme_set(
  theme_classic(base_size = 7, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(colour = "black"),
      panel.grid = element_blank(),
      legend.position = "none"
    )
)

make_prevalence_abundance_panel <- function(label) {
  sub <- df[df$type == label, , drop = FALSE]
  row <- stats[stats$type == label, , drop = FALSE]
  annotation <- paste0(
    "n = ", format(row$n_vOTU_prevalence_gt_1pct, big.mark = ","), "\n",
    "R2 = ", sprintf("%.3f", row$r_squared), "\n",
    "P = ", format_p(row$slope_p_value)
  )
  ggplot(sub, aes(prevalence_percent, mean_abundance_present_samples)) +
    geom_point(colour = unname(colours[names(labels)[labels == label]]),
               alpha = 0.42, size = 0.7, stroke = 0) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = "black",
                fill = unname(colours[names(labels)[labels == label]]),
                linewidth = 0.55, alpha = 0.18) +
    annotate("text", x = Inf, y = Inf, label = annotation,
             hjust = 1.05, vjust = 1.15, size = 2.15,
             lineheight = 0.95) +
    labs(title = label, x = "Prevalence (%)",
         y = "Mean abundance in present samples")
}

for (label in levels(df$type)) {
  p <- make_prevalence_abundance_panel(label)
  out <- paste0(tolower(label), "_votu_prevalence_abundance.pdf")
  ggsave(file.path(figures_dir, out), p, width = 85 / 25.4,
         height = 70 / 25.4, units = "in",
         device = function(filename, width, height, ...) {
           grDevices::pdf(file = filename, width = width, height = height,
                          useDingbats = FALSE)
         })
}

combined <- ggplot(df, aes(prevalence_percent, mean_abundance_present_samples)) +
  geom_point(aes(colour = type), alpha = 0.42, size = 0.55, stroke = 0) +
  geom_smooth(aes(fill = type), method = "lm", formula = y ~ x, se = TRUE,
              colour = "black", linewidth = 0.45, alpha = 0.16) +
  facet_wrap(~ type, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = unname(colours)) +
  scale_fill_manual(values = unname(colours)) +
  labs(x = "Prevalence (%)", y = "Mean abundance in present samples")
ggsave(file.path(figures_dir, "votu_prevalence_abundance_by_manure_type.pdf"),
       combined, width = 180 / 25.4, height = 62 / 25.4,
       units = "in",
       device = function(filename, width, height, ...) {
         grDevices::pdf(file = filename, width = width, height = height,
                        useDingbats = FALSE)
       })

print(stats)
