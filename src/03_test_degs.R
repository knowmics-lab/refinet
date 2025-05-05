library(keras)
library(tensorflow)
library(progress)
library(dplyr)
library(openxlsx)
library(edgeR)
library(limma)
source("commons.R")

model_path <- "./models_positional_encoding_with_classifier_single"

# # Set up output capture
# output_file <- file.path(model_path, "degs_test_output_closest_encoded_control_samples.txt")
# if (file.exists(output_file)) {
#   file.remove(output_file)
# }
# sink(file = output_file, append = TRUE, split = TRUE) # split=TRUE shows output in console and file

# Load Factor Levels
load(file.path(model_path, "factor_levels.rda"))

class_map <- readr::read_csv(file.path(model_path, "biopsy_site_categories_with_id.csv"),
  col_names = TRUE
)
class_map_with_site <- class_map
class_map$site <- NULL
class_map <- unique(class_map)
colnames(class_map) <- c(
  "predicted_system_level_1",
  "predicted_system_level_2",
  "predicted_organ_level_3",
  "class_id"
)

# Load Models

full_model <- load_model_hdf5(file.path(model_path, "final_full_model.h5"))
encoder_model <- load_model_hdf5(file.path(model_path, "final_encoder.h5"))

train_data_env <- new.env()
load("new_dataset.rda", envir = train_data_env)
load("new_dataset_count.rda", envir = train_data_env)
train_data_env$metadata <- train_data_env$phenoDF[train_data_env$phenoDF$sample.type %in% c("normal", "adjacent"), ] #
train_data_env$tpm_matrix <- train_data_env$tpm_matrix[, train_data_env$metadata$sample.id]
train_data_env$count_matrix <- train_data_env$count_matrix[, train_data_env$metadata$sample.id]
input_data <- t(as.array(train_data_env$tpm_matrix))
train_data_env$data_encoded <- predict(encoder_model, input_data)
rownames(train_data_env$data_encoded) <- rownames(input_data)
train_data_env$metadata <- train_data_env$metadata %>% inner_join(class_map_with_site, by = c("biopsy.site" = "site"))
rm(input_data)

# Load TPM data
ipf_data <- readRDS("datasets/ipf.rds")
other_data <- readRDS("datasets/other_datasets_tpm.rds")
tpm_data <- cbind(ipf_data, other_data)

# Load RAW counts data
ipf_count <- readRDS("datasets/ipf_count.rds")
other_count <- readRDS("datasets/other_datasets_count.rds")
count_data <- cbind(ipf_count, other_count)

# Load Metedata
ipf_metadata <- readRDS("ipf_metadata.rds")
ipf_metadata <- ipf_metadata[, c("sample_id", "type", "dataset")]
other_datasets_metadata <- readRDS("other_datasets_metadata.rds")
other_datasets_metadata <- other_datasets_metadata[, c("sample_id", "type", "dataset")]
metadata <- rbind(ipf_metadata, other_datasets_metadata)
metadata <- na.omit(metadata)

# Cleanup
rm(
  ipf_data,
  other_data,
  ipf_count,
  other_count,
  ipf_metadata,
  other_datasets_metadata
)

# Keep TPM and count data only for samples with metadata
tpm_data <- tpm_data[, metadata$sample_id]
count_data <- count_data[, metadata$sample_id]

# Classify and encode all case samples in one batch to speed up the process
input_data <- t(as.array(tpm_data[, metadata$type == "case"]))
sample_predictions <- predict(full_model, input_data)[[2]]
sample_predicted_classes <- get_class_names(sample_predictions, sample_factors)
names(sample_predicted_classes) <- rownames(input_data)
sample_encoded <- predict(encoder_model, input_data)
rownames(sample_encoded) <- rownames(input_data)

compute_samples_tissue_class <- function(predicted_classes,
                                         all_classes,
                                         samples_data,
                                         train_data_env,
                                         method = c("random", "octad", "refinet"),
                                         num_marker_genes = 1000) {
  method <- match.arg(method)
  if (method == "random") {
    # Randomly sample one of the classes
    sample(all_classes, 1)
  } else if (method == "octad") {
    # Compute the marker genes and select the closest sample
    # according to the spearman correlation between the samples
    # and the train data
    train_count_matrix <- train_data_env$count_matrix
    train_median <- apply(train_count_matrix, 1, median)
    train_expressed_gene <- names(train_median)[train_median > 1]
    train_count_matrix <- train_count_matrix[train_expressed_gene, ]
    rank_train <- apply(train_count_matrix, 2, rank)
    rank_train_sd <- apply(rank_train, 1, sd)
    train_marker_genes <- names(sort(rank_train_sd, decreasing = TRUE))[seq_len(num_marker_genes)]
    train_marker_genes <- intersect(rownames(samples_data), train_marker_genes)
    correlation_matrix <- cor(
      samples_data[train_marker_genes, ],
      train_count_matrix[train_marker_genes, ],
      method = "spearman"
    )
    correlation_matrix[is.na(correlation_matrix)] <- 0
    cell_line_median_cor <- apply(correlation_matrix, 2, median) %>%
      sort(decreasing = TRUE)
    best_sample <- names(cell_line_median_cor)[1]
    train_data_env$metadata$class_id[train_data_env$metadata$sample.id == best_sample]
  } else if (method == "refinet") {
    # Use the result of the classifier for each sample to get the majority class
    names(sort(table(predicted_classes), decreasing = TRUE))[1]
  } else {
    stop("Invalid method")
  }
}


compute_control_samples <- function(samples_data,
                                    tissue_class,
                                    train_data_env,
                                    method = c("random", "octad", "refinet", 
                                               "encoded"),
                                    n_samples = nrow(samples_data),
                                    n_var_genes = 500,
                                    var_method = c("IQR", "MAD"),
                                    mode = c("closest", "furthest")) {
  method <- match.arg(method)
  var_method <- match.arg(var_method)
  mode <- match.arg(mode)
  if (length(mode) > 1) {
    stop("Only one mode can be used at a time")
  }
  mode <- mode == "closest"
  if (is.null(tissue_class) && method == "encoded_all") {
    stop("You must provide a tissue class when using the encoded_all method")
  }
  # Get the train data for the tissue class
  sample_ids <- train_data_env$metadata$sample.id[train_data_env$metadata$class_id == tissue_class]
  count_train_data <- t(train_data_env$count_matrix[, sample_ids])
  encoded_train_data <- train_data_env$data_encoded[sample_ids, ]
  control_samples <- NULL
  if (method == "random") {
    # Randomly sample n_samples from the train sample names
    control_samples <- sample(sample_ids, min(n_samples, length(sample_ids)))
  } else if (method == "refinet") {
    real_n_samples <- nrow(samples_data)
    if (real_n_samples > length(sample_ids)) {
      return (sample_ids) # If there are not enough samples in the train data, return all the samples
    }
    control_samples <- character(0)
    in_train_data   <- encoded_train_data
    in_samples_data <- samples_data
    while (length(control_samples) < real_n_samples) {
      # Compute the spearman correlation between all pairs of rows in in_train_data and in_samples_data
      cor_matrix <- cor(t(in_train_data), t(in_samples_data), method = "spearman")
      # Get the row and column with the highest (mode = "closest") or lowest (mode = "furthest") correlation
      max_cor_idx <- if (mode) which.max(cor_matrix) else which.min(cor_matrix)
      max_cor_row <- (max_cor_idx - 1) %% nrow(cor_matrix) + 1
      max_cor_col <- (max_cor_idx - 1) %/% nrow(cor_matrix) + 1
      control_samples <- c(control_samples, rownames(cor_matrix)[max_cor_row])
      # Remove the row and column from the matrices
      in_train_data <- in_train_data[-max_cor_row, , drop = FALSE]
      in_samples_data <- in_samples_data[-max_cor_col, , drop = FALSE]
    }
  } else {
    if (method == "octad") {
      # Compute the control samples by using the most variable genes
      # Most variable genes are computed on the train data using the IQR or MAD
      common_genes <- intersect(colnames(samples_data), colnames(count_train_data))
      if (length(common_genes) == 0) {
        stop("No common genes found between the samples and the train data")
      }
      if (length(common_genes) < n_var_genes) {
        stop("Not enough common genes found between the samples and the train data")
      }
      count_train_data <- count_train_data[, common_genes]
      samples_data <- samples_data[, common_genes]
      if (var_method == "IQR") {
        iqr_gene <- apply(count_train_data, 2, stats::IQR) # get the IQR per gene
        varying_genes <- order(iqr_gene, decreasing = TRUE)[seq_len(min(n_var_genes, length(iqr_gene)))]
      } else if (var_method == "MAD") {
        mad_gene <- apply(count_train_data, 2, stats::mad) # get the MAD per gene
        varying_genes <- order(mad_gene, decreasing = TRUE)[seq_len(min(n_var_genes, length(mad_gene)))]
      }

      in_train_data <- count_train_data[, varying_genes]
      in_samples_data <- samples_data[, varying_genes]
    } else if (method == "encoded") {
      # Compute the control samples using the encoded train data
      if (ncol(encoded_train_data) != ncol(samples_data)) {
        stop("The encoded train data and the samples data have different number of columns")
      }
      in_train_data <- encoded_train_data
      in_samples_data <- samples_data
    } else {
      stop("Invalid method")
    }
    # Compute the spearman correlation between all pairs of rows in in_train_data and in_samples_data
    cor_matrix <- cor(t(in_train_data), t(in_samples_data), method = "spearman")
    # Compute the median correlations for each train sample
    median_correlations <- apply(cor_matrix, 1, median)
    # Get the n_samples (or less) with the highest median correlations
    real_n_samples <- min(n_samples, length(median_correlations))

    sample_idx <- order(median_correlations, decreasing = mode)[seq_len(real_n_samples)]
    control_samples <- names(median_correlations)[sample_idx]
  }


  return(control_samples)
}

get_control_samples <- function(samples_data_count,
                                samples_data_encoded,
                                train_data_env,
                                predicted_classes,
                                all_classes,
                                method = c("random", "octad", "refinet", "encoded"),
                                n_samples = nrow(samples_data_count),
                                n_marker_genes = 1000,
                                n_var_genes = 500,
                                mode = c("closest", "furthest")) {
  method <- match.arg(method)
  # var_method <- match.arg(var_method)
  mode <- match.arg(mode)

  control_samples <- NULL

  while (is.null(control_samples) ||
    (length(control_samples) == 0 && method == "random")) {
    # compute the dataset tissue
    dataset_tissue_type <- compute_samples_tissue_class(
      predicted_classes = predicted_classes,
      all_classes = all_classes,
      samples_data = samples_data_count,
      train_data_env = train_data_env,
      method = if (method == "encoded") "octad" else method,
      num_marker_genes = n_marker_genes
    )

    # Get the control samples
    control_samples <- compute_control_samples(
      samples_data = if (method == "octad") {
        t(samples_data_count)
      } else {
        samples_data_encoded
      },
      tissue_class = dataset_tissue_type,
      train_data_env = train_data_env,
      method = method,
      mode = mode,
      n_samples = n_samples,
      n_var_genes = n_var_genes # , var_method = var_method
    )
  }

  train_data_env$count_matrix[, control_samples]
}

compute_degs <- function(count_matrix, sample_types) {
  design <- model.matrix(~ 0 + sample_types)
  colnames(design) <- c("case", "control")
  dge <- DGEList(count_matrix)
  keep <- filterByExpr(dge, design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  contrasts <- makeContrasts(case - control, levels = design)
  fit2 <- suppressWarnings(eBayes(contrasts.fit(fit, contrasts)))
  topTable(fit2, coef = "case - control", number = Inf)
}

count <- function(x) {
  (length(which(x)))
}

vector_frobenius_norm <- function(u, v) (sqrt(sum((u - v)^2)))

conditional_cor <- function(v1, v2, cond = NULL, ...) {
  if (is.null(cond)) {
    n <- intersect(names(v1), names(v2))
    cor(v1[n], v2[n], ...)
  } else {
    n <- names(v1[cond])
    cor(v1[cond], v2[n], ...)
  }
}

factorize <- function(x) {
  tmp <- character(length(x))
  tmp[x < 0] <- "low"
  tmp[x > 0] <- "high"
  tmp[x == 0] <- NA
  factor(tmp, levels = c("high", "low"))
}

evaluate_as_classifier <- function(p, o) {
  library(caret)
  cm <- caret::confusionMatrix(factorize(p), factorize(o))
  unname(cm$byClass[c("Pos Pred Value", "Neg Pred Value", "Balanced Accuracy")])
}

my_score <- function(real_signature, predicted_signature, is_upregulated = TRUE) {
  num_genes <- length(predicted_signature)
  ks_score <- 0
  predicted_signature <- rank(predicted_signature)
  real_rank <- predicted_signature[as.vector(real_signature)]
  real_position <- sort(real_rank)
  num_tags <- length(real_position)
  if (num_tags > 1) {
    a <- 0
    b <- 0
    a <- max(vapply(seq_len(num_tags), function(j) {
      j / num_tags - real_position[j] / num_genes
    }, FUN.VALUE = numeric(1)))
    b <- max(vapply(seq_len(num_tags), function(j) {
      real_position[j] / num_genes - (j - 1) / num_tags
    }, FUN.VALUE = numeric(1)))
    if (a > b) {
      ks_score <- a
    } else {
      ks_score <- -b
    }
  } else {
    ks_score <- 0
  }
  return(ks_score)
}

compute_degs_results <- function(degs_original, degs_predicted) {
  common_genes <- intersect(rownames(degs_original), rownames(degs_predicted))
  degs_original <- degs_original[common_genes, ]
  degs_predicted <- degs_predicted[common_genes, ]

  degs_original_fdr_005 <- degs_original[degs_original$adj.P.Val < 0.05, ]
  degs_predicted_fdr_005 <- degs_predicted[degs_predicted$adj.P.Val < 0.05, ]

  degs_original_fdr_001 <- degs_original[degs_original$adj.P.Val < 0.01, ]
  degs_predicted_fdr_001 <- degs_predicted[degs_predicted$adj.P.Val < 0.01, ]

  # Get the DEGs that are in both datasets
  degs_both_datasets_005 <- intersect(rownames(degs_original_fdr_005), rownames(degs_predicted_fdr_005))
  degs_both_datasets_001 <- intersect(rownames(degs_original_fdr_001), rownames(degs_predicted_fdr_001))

  common_original_fdr_005 <- setNames(degs_original_fdr_005[degs_both_datasets_005, "logFC"], degs_both_datasets_005)
  common_predicted_fdr_005 <- setNames(degs_predicted_fdr_005[degs_both_datasets_005, "logFC"], degs_both_datasets_005)

  common_original_fdr_001 <- setNames(degs_original_fdr_001[degs_both_datasets_001, "logFC"], degs_both_datasets_001)
  common_predicted_fdr_001 <- setNames(degs_predicted_fdr_001[degs_both_datasets_001, "logFC"], degs_both_datasets_001)

  common_percentage_fdr_005 <- length(degs_both_datasets_005) / nrow(degs_original_fdr_005)
  common_percentage_fdr_001 <- length(degs_both_datasets_001) / nrow(degs_original_fdr_001)

  frobenius_norm_005 <- vector_frobenius_norm(common_original_fdr_005, common_predicted_fdr_005)
  frobenius_norm_001 <- vector_frobenius_norm(common_original_fdr_001, common_predicted_fdr_001)

  cor_001 <- conditional_cor(
    common_original_fdr_001,
    common_predicted_fdr_001
  )
  cor_005 <- conditional_cor(
    common_original_fdr_005,
    common_predicted_fdr_005
  )
  eval_001 <- evaluate_as_classifier(
    common_predicted_fdr_001,
    common_original_fdr_001
  )
  eval_005 <- evaluate_as_classifier(
    common_predicted_fdr_005,
    common_original_fdr_005
  )

  return(c(
    common_percentage_fdr_005 = common_percentage_fdr_005,
    common_percentage_fdr_001 = common_percentage_fdr_001,
    frobenius_norm_005 = frobenius_norm_005,
    frobenius_norm_001 = frobenius_norm_001,
    cor_001 = cor_001,
    cor_005 = cor_005,
    ppv_001 = eval_001[1],
    npv_001 = eval_001[2],
    ba_001 = eval_001[3],
    ppv_005 = eval_005[1],
    npv_005 = eval_005[2],
    ba_005 = eval_005[3]
  ))
}


random_repetitions <- 10

parameters_combinations <- data.frame(
  method = c("octad", "refinet"),
  mode = c("closest", "closest"),
  stringsAsFactors = FALSE
)

datasets <- unique(metadata$dataset)

random_results <- data.frame(
  dataset = character(),
  method = character(),
  common_percentage_fdr_005 = numeric(),
  frobenius_norm_005 = numeric(),
  cor_005 = numeric(),
  ppv_005 = numeric(),
  npv_005 = numeric(),
  ba_005 = numeric(),
  stringsAsFactors = FALSE
)

results <- data.frame(
  dataset = character(),
  method = character(),
  common_percentage_fdr_005 = numeric(),
  frobenius_norm_005 = numeric(),
  cor_005 = numeric(),
  ppv_005 = numeric(),
  npv_005 = numeric(),
  ba_005 = numeric(),
  stringsAsFactors = FALSE
)

for (dataset in datasets) {
  cat("Testing dataset:", dataset, "\n")
  dataset_samples <- metadata$sample_id[metadata$dataset == dataset]
  dataset_sample_types <- metadata$type[metadata$dataset == dataset]

  if (length(which(dataset_sample_types == "control")) == 0) {
    cat("No control samples found for dataset:", dataset, "\n")
    next()
  }
  if (length(which(dataset_sample_types == "case")) == 0) {
    cat("No case samples found for dataset:", dataset, "\n")
    next()
  }

  dataset_case_samples <- intersect(dataset_samples[dataset_sample_types == "case"], names(sample_predicted_classes))
  dataset_predicted_classes <- sample_predicted_classes[dataset_case_samples]
  dataset_encoded <- sample_encoded[dataset_case_samples, ]

  # get the TPM data for this dataset
  dataset_tpm <- tpm_data[, dataset_samples]
  dataset_count <- count_data[, dataset_samples]

  degs_original <- compute_degs(dataset_count, dataset_sample_types)

  for (repetition in seq_len(random_repetitions)) {
    cat("  Random repetition:", repetition, "\r")
    predicted_control_count <- get_control_samples(
      samples_data_count = dataset_count[, dataset_case_samples],
      samples_data_encoded = dataset_encoded,
      train_data_env = train_data_env,
      predicted_classes = dataset_predicted_classes,
      n_samples = count(dataset_sample_types == "control"),
      all_classes = levels(sample_factors),
      method = "random",
      n_marker_genes = 1000,
      n_var_genes = 500,
      mode = "closest"
    )
    # Create a new dataset with the original case samples and the predicted control samples
    new_dataset <- cbind(
      dataset_count[, dataset_sample_types == "case"],
      predicted_control_count
    )
    new_dataset_sample_types <- c(
      dataset_sample_types[dataset_sample_types == "case"],
      rep("control", ncol(predicted_control_count))
    )
    degs_new_dataset <- compute_degs(new_dataset, new_dataset_sample_types)

    comparison <- compute_degs_results(degs_original, degs_new_dataset)
    random_results <- rbind(random_results, data.frame(
      dataset = dataset,
      method = "random",
      common_percentage_fdr_005 = comparison["common_percentage_fdr_005"],
      frobenius_norm_005 = comparison["frobenius_norm_005"],
      cor_005 = comparison["cor_005"],
      ppv_005 = comparison["ppv_005"],
      npv_005 = comparison["npv_005"],
      ba_005 = comparison["ba_005"],
      stringsAsFactors = FALSE
    ))
  }
  for (i in seq_len(nrow(parameters_combinations))) {
    cat("  Parameters combination:", i, " of ", nrow(parameters_combinations), "\r")
    method <- parameters_combinations$method[i]
    mode <- parameters_combinations$mode[i]
    predicted_control_count <- get_control_samples(
      samples_data_count = dataset_count[, dataset_case_samples],
      samples_data_encoded = dataset_encoded,
      train_data_env = train_data_env,
      predicted_classes = dataset_predicted_classes,
      n_samples = count(dataset_sample_types == "control"),
      all_classes = levels(sample_factors),
      method = method,
      n_marker_genes = 1000,
      n_var_genes = 500,
      mode = mode
    )
    # Create a new dataset with the original case samples and the predicted control samples
    new_dataset <- cbind(
      dataset_count[, dataset_sample_types == "case"],
      predicted_control_count
    )
    new_dataset_sample_types <- c(
      dataset_sample_types[dataset_sample_types == "case"],
      rep("control", ncol(predicted_control_count))
    )
    degs_new_dataset <- compute_degs(new_dataset, new_dataset_sample_types)

    comparison <- compute_degs_results(degs_original, degs_new_dataset)
    results <- rbind(results, data.frame(
      dataset = dataset,
      method = method,
      common_percentage_fdr_005 = comparison["common_percentage_fdr_005"],
      frobenius_norm_005 = comparison["frobenius_norm_005"],
      cor_005 = comparison["cor_005"],
      ppv_005 = comparison["ppv_005"],
      npv_005 = comparison["npv_005"],
      ba_005 = comparison["ba_005"],
      stringsAsFactors = FALSE
    ))
  }
  cat("\n")

}

save(random_results, results, file = "test_degs_results.rda")
final_results <- rbind(results, random_results)

library(openxlsx)
write.xlsx(random_results, file = "test_degs_results_random.xlsx")
write.xlsx(results, file = "test_degs_results_methods.xlsx")
write.xlsx(final_results, file = "test_degs_results_all.xlsx")
