# Random-forest analysis of environmental predictors associated with viral alpha diversity

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) == 0) {
  analysis_dir <- normalizePath(getwd(), mustWork = FALSE)
  if (basename(analysis_dir) != "random_forest_environmental_predictors") {
    analysis_dir <- file.path(
      analysis_dir,
      "Revision_analysis",
      "statistical_analysis",
      "random_forest_environmental_predictors"
    )
  }
} else {
  script_file <- sub("^--file=", "", script_arg[1])
  script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
  if (!grepl("^/", script_file)) script_file <- file.path(getwd(), script_file)
  analysis_dir <- dirname(normalizePath(script_file, mustWork = FALSE))
}

input_file <- file.path(analysis_dir, "shannon_environmental_predictors.tsv")
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

data <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)

required_columns <- c(
  "Shannon",
  "Antimicrobial_usage_in_livestock",
  "Scale_of_livestock_breeding",
  "Agricultural_production",
  "Economic_level",
  "Livestock_feces_production",
  "Pesticides_usage",
  "Livestock_meat_production",
  "Fertilization_of_livestock_manure",
  "Human_population"
)

missing_columns <- setdiff(required_columns, names(data))
if (length(missing_columns) > 0) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

data <- data[, required_columns]
for (column in required_columns) {
  data[[column]] <- suppressWarnings(as.numeric(data[[column]]))
}

complete_rows <- stats::complete.cases(data)
if (!all(complete_rows)) {
  data <- data[complete_rows, , drop = FALSE]
}

if (nrow(data) < 10) {
  stop("Too few complete observations for random-forest analysis.")
}

if (!requireNamespace("randomForest", quietly = TRUE)) {
  stop(
    "The R package 'randomForest' is required for the random-forest model."
  )
}

ntree <- 500
permutation_replicates <- 1000
random_seed <- 123

set.seed(random_seed)
rf_model <- randomForest::randomForest(
  Shannon ~ .,
  data = data,
  importance = TRUE,
  ntree = ntree,
  na.action = stats::na.omit
)

importance_table <- as.data.frame(randomForest::importance(rf_model), check.names = FALSE)
importance_table$Predictor <- rownames(importance_table)
rownames(importance_table) <- NULL

if ("%IncMSE" %in% names(importance_table)) {
  importance_table <- importance_table[
    order(importance_table[["%IncMSE"]], decreasing = TRUE),
    ,
    drop = FALSE
  ]
}

utils::write.table(
  importance_table,
  file.path(analysis_dir, "random_forest_variable_importance.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

model_summary <- data.frame(
  Metric = c(
    "n_samples",
    "n_predictors",
    "ntree",
    "random_seed",
    "mean_squared_residual",
    "variance_explained_percent"
  ),
  Value = c(
    nrow(data),
    length(required_columns) - 1,
    ntree,
    random_seed,
    tail(rf_model$mse, 1),
    tail(rf_model$rsq, 1) * 100
  )
)

utils::write.table(
  model_summary,
  file.path(analysis_dir, "random_forest_model_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

set.seed(111)
train_index <- sample(seq_len(nrow(data)), floor(nrow(data) * 0.7))
training_data <- data[train_index, , drop = FALSE]
testing_data <- data[-train_index, , drop = FALSE]

set.seed(random_seed)
rf_train <- randomForest::randomForest(
  Shannon ~ .,
  data = training_data,
  importance = TRUE,
  ntree = ntree,
  na.action = stats::na.omit
)

prediction_training <- stats::predict(rf_train, newdata = training_data)
prediction_testing <- stats::predict(rf_train, newdata = testing_data)

performance <- data.frame(
  Dataset = c("training", "testing"),
  N = c(nrow(training_data), nrow(testing_data)),
  RMSE = c(
    sqrt(mean((training_data$Shannon - prediction_training)^2)),
    sqrt(mean((testing_data$Shannon - prediction_testing)^2))
  ),
  R_squared = c(
    stats::cor(training_data$Shannon, prediction_training)^2,
    stats::cor(testing_data$Shannon, prediction_testing)^2
  )
)

utils::write.table(
  performance,
  file.path(analysis_dir, "random_forest_train_test_performance.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (requireNamespace("rfPermute", quietly = TRUE)) {
  set.seed(random_seed)
  rf_permutation <- rfPermute::rfPermute(
    Shannon ~ .,
    data = data,
    importance = TRUE,
    ntree = ntree,
    nrep = permutation_replicates,
    num.cores = 1
  )

  permutation_importance <- as.data.frame(
    randomForest::importance(rf_permutation, scale = TRUE),
    check.names = FALSE
  )
  permutation_importance$Predictor <- rownames(permutation_importance)
  rownames(permutation_importance) <- NULL

  p_values <- rf_permutation$pval
  if (length(dim(p_values)) == 3) {
    p_values <- p_values[, , dim(p_values)[3], drop = FALSE]
    p_values <- as.data.frame(p_values[, , 1], check.names = FALSE)
  } else {
    p_values <- as.data.frame(p_values, check.names = FALSE)
  }
  p_values$Predictor <- rownames(p_values)
  rownames(p_values) <- NULL

  permutation_results <- merge(
    permutation_importance,
    p_values,
    by = "Predictor",
    all = TRUE,
    suffixes = c("_importance", "_p_value")
  )

  utils::write.table(
    permutation_results,
    file.path(analysis_dir, "random_forest_permutation_importance.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
} else {
  message(
    "Package 'rfPermute' is not installed; permutation-based significance ",
    "testing was skipped."
  )
}

print(importance_table)
print(model_summary)
print(performance)
