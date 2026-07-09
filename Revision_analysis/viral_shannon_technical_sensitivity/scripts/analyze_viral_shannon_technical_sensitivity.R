library(ggplot2)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) == 0) {
  wd <- normalizePath(getwd(), mustWork = FALSE)
  root <- if (basename(wd) == "viral_shannon_technical_sensitivity") {
    wd
  } else {
    file.path(wd, "viral_shannon_technical_sensitivity")
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

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  formatC(p, format = "f", digits = 3)
}

summary_stats <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  c(n = length(x), mean = mean(x), median = median(x), sd = sd(x),
    Q1 = as.numeric(quantile(x, 0.25, names = FALSE)),
    Q3 = as.numeric(quantile(x, 0.75, names = FALSE)),
    minimum = min(x), maximum = max(x))
}

vif_for_predictors <- function(dat, predictors) {
  rows <- lapply(predictors, function(pred) {
    others <- setdiff(predictors, pred)
    fit <- lm(as.formula(paste(pred, "~", paste(others, collapse = " + "))),
              data = dat)
    r2 <- summary(fit)$r.squared
    vif <- ifelse(is.na(r2) || r2 >= 1, Inf, 1 / (1 - r2))
    data.frame(variable = pred, VIF = vif,
               VIF_flag = ifelse(vif > 5, "high_gt_5", "acceptable_le_5"),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

dunn_test_manual <- function(values, groups) {
  ok <- is.finite(values) & !is.na(groups)
  values <- values[ok]
  groups <- droplevels(as.factor(groups[ok]))
  n <- length(values)
  ranks <- rank(values, ties.method = "average")
  tie_sizes <- as.numeric(table(values))
  tie_correction <- 1 - sum(tie_sizes^3 - tie_sizes) / (n^3 - n)
  lev <- levels(groups)
  mean_ranks <- tapply(ranks, groups, mean)
  n_by_group <- table(groups)
  pairs <- combn(lev, 2, simplify = FALSE)
  out <- do.call(rbind, lapply(pairs, function(pair) {
    g1 <- pair[1]
    g2 <- pair[2]
    se <- sqrt((n * (n + 1) / 12) * tie_correction *
                 (1 / as.numeric(n_by_group[g1]) + 1 / as.numeric(n_by_group[g2])))
    z <- (mean_ranks[g1] - mean_ranks[g2]) / se
    p <- 2 * pnorm(-abs(z))
    data.frame(group1 = g1, group2 = g2,
               n_group1 = as.integer(n_by_group[g1]),
               n_group2 = as.integer(n_by_group[g2]),
               z_statistic = as.numeric(z), p_value = as.numeric(p),
               stringsAsFactors = FALSE)
  }))
  out$BH_adjusted_p_value <- p.adjust(out$p_value, method = "BH")
  out
}

raw <- read.delim(file.path(data_dir, "viral_shannon_technical_metrics.tsv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
required <- c("sample", "Shannon index", "nReads", "sum_len", "num_seqs", "Continent")
missing <- setdiff(required, names(raw))
if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

dat <- raw[, required]
for (col in c("Shannon index", "nReads", "sum_len", "num_seqs")) {
  dat[[col]] <- suppressWarnings(as.numeric(dat[[col]]))
}
names(dat)[names(dat) == "Shannon index"] <- "Shannon_index"
dat$Continent <- factor(dat$Continent, levels = unique(raw$Continent))
dat$log_nReads <- log10(dat$nReads + 1)
dat$log_sum_len <- log10(dat$sum_len + 1)
dat$log_num_seqs <- log10(dat$num_seqs + 1)

complete <- dat[complete.cases(dat[, c("Shannon_index", "log_nReads",
                                        "log_sum_len", "log_num_seqs",
                                        "Continent")]), , drop = FALSE]
complete$Continent <- droplevels(complete$Continent)

fit <- lm(Shannon_index ~ log_nReads + log_sum_len + log_num_seqs,
          data = complete)
complete$Shannon_residual <- residuals(fit)

write_tsv(data.frame(sample = complete$sample,
                     Continent = as.character(complete$Continent),
                     Shannon_index = complete$Shannon_index,
                     nReads = complete$nReads,
                     sum_len = complete$sum_len,
                     num_seqs = complete$num_seqs,
                     Shannon_residual = complete$Shannon_residual),
          "viral_shannon_residuals.tsv")

resid_summary <- do.call(rbind, lapply(levels(complete$Continent), function(continent) {
  sub <- complete[complete$Continent == continent, , drop = FALSE]
  s <- summary_stats(sub$Shannon_residual)
  data.frame(Continent = continent, sample_number = as.integer(s["n"]),
             mean = s["mean"], median = s["median"],
             standard_deviation = s["sd"], Q1 = s["Q1"], Q3 = s["Q3"],
             minimum = s["minimum"], maximum = s["maximum"],
             stringsAsFactors = FALSE)
}))
write_tsv(resid_summary, "viral_shannon_residual_summary_by_continent.tsv")

kw <- kruskal.test(Shannon_residual ~ Continent, data = complete)
write_tsv(data.frame(test_name = "Kruskal-Wallis rank sum test",
                     chi_squared_statistic = unname(kw$statistic),
                     degrees_of_freedom = unname(kw$parameter),
                     p_value = kw$p.value),
          "viral_shannon_residual_continent_test.tsv")
if (!is.na(kw$p.value) && kw$p.value < 0.05) {
  write_tsv(dunn_test_manual(complete$Shannon_residual, complete$Continent),
            "viral_shannon_residual_dunn_test.tsv")
}

vif <- vif_for_predictors(complete, c("log_nReads", "log_sum_len", "log_num_seqs"))
write_tsv(vif, "viral_shannon_model_vif.tsv")

fit_summary <- summary(fit)
coef_table <- as.data.frame(coef(fit_summary))
coef_table$term <- rownames(coef_table)
coef_table <- coef_table[, c("term", "Estimate", "Std. Error", "t value", "Pr(>|t|)")]
names(coef_table) <- c("term", "estimate", "standard_error", "t_value", "p_value")
write_tsv(coef_table, "viral_shannon_model_coefficients.tsv")

raw_median <- aggregate(Shannon_index ~ Continent, complete, median)
resid_median <- aggregate(Shannon_residual ~ Continent, complete, median)
median_compare <- merge(raw_median, resid_median, by = "Continent")
names(median_compare) <- c("Continent", "raw_Shannon_median", "residual_Shannon_median")
raw_order <- median_compare[order(-median_compare$raw_Shannon_median), "Continent"]
resid_order <- median_compare[order(-median_compare$residual_Shannon_median), "Continent"]
median_compare$raw_rank_descending <- match(median_compare$Continent, raw_order)
median_compare$residual_rank_descending <- match(median_compare$Continent, resid_order)
write_tsv(median_compare, "viral_shannon_raw_vs_residual_median_order.tsv")

sink(file.path(results_dir, "viral_shannon_model_summary.txt"))
cat("Technical residual model\n")
cat("Formula: Shannon index ~ log_nReads + log_sum_len + log_num_seqs\n\n")
cat("Number of samples used:", nobs(fit), "\n\n")
cat("Coefficient table:\n")
print(coef_table, row.names = FALSE)
cat("\nR-squared:", fit_summary$r.squared, "\n")
cat("Adjusted R-squared:", fit_summary$adj.r.squared, "\n")
cat("Residual standard error:", fit_summary$sigma, "\n")
cat("\nVIF:\n")
print(vif, row.names = FALSE)
sink()

interpretation <- c(
  "Continental differences remained statistically detectable after accounting for read count and assembly-output metrics.",
  paste0("Raw Shannon median order, high to low: ", paste(raw_order, collapse = " > "), "."),
  paste0("Residual Shannon median order, high to low: ", paste(resid_order, collapse = " > "), "."),
  paste0("Kruskal-Wallis P value for residual continent differences: ", format_p(kw$p.value), "."),
  "This sensitivity analysis does not perform full catalogue-wide read recruitment."
)
writeLines(interpretation, file.path(results_dir, "viral_shannon_interpretation.txt"),
           useBytes = TRUE)

theme_set(theme_classic(base_size = 7, base_family = "sans"))
p <- ggplot(complete, aes(x = Continent, y = Shannon_residual)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35,
             colour = "grey45") +
  geom_boxplot(width = 0.56, fill = "#D9E2E8", colour = "#34495E",
               outlier.shape = NA, linewidth = 0.35) +
  geom_jitter(width = 0.16, height = 0, size = 0.35, alpha = 0.35,
              colour = "#34495E") +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Kruskal-Wallis P = ", format_p(kw$p.value)),
           hjust = 1.04, vjust = 1.35, size = 2.25) +
  labs(x = "Continent", y = "Residual viral Shannon diversity") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(figures_dir, "viral_shannon_residual_by_continent.pdf"),
       p, width = 155 / 25.4, height = 95 / 25.4,
       units = "in",
       device = function(filename, width, height, ...) {
         grDevices::pdf(file = filename, width = width, height = height,
                        useDingbats = FALSE)
       })

print(interpretation)
