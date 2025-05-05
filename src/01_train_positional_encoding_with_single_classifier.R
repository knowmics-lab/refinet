# Setup console size
width <- 160
if (Sys.getenv("COLUMNS") != "") {
  width <- as.integer(Sys.getenv("COLUMNS"))
}
options(width = width)
options(setWidthOnResize = TRUE)

# Set up output capture
output_dir <- "models_positional_encoding_with_classifier_single"
output_file <- file.path(output_dir, "training_output.txt")
dir.create(output_dir, showWarnings = FALSE)
if (file.exists(output_file)) {
  file.remove(output_file)
}
sink(file = output_file, append = TRUE, split = TRUE) # split=TRUE shows output in console and file

gather_site_map <- function(site_map_file, site_map_with_id_file,
                            category_columns = c("system_level_1", "system_level_2", "organ_level_3"),
                            id_column = "class_id") {
  site_map <- readr::read_csv(site_map_file, col_names = TRUE)
  site_map$tmp_id <- do.call(paste, c(site_map[category_columns], sep = "_"))
  tmp_ids <- unique(site_map$tmp_id)
  ids <- paste0("class_", seq_along(tmp_ids))
  site_map[[id_column]] <- ids[match(site_map$tmp_id, tmp_ids)]
  site_map$tmp_id <- NULL
  readr::write_csv(site_map, site_map_with_id_file)
  return(site_map)
}

# Load the training data
load("datasets/new_dataset.rda")

# Print data dimensions for debugging
cat("Data dimensions:\n")
cat(paste("TPM matrix dimensions:", paste(dim(tpm_matrix), collapse = " x ")), "\n")
cat(paste("PhenoDF dimensions:", paste(dim(phenoDF), collapse = " x ")), "\n")

n_samples <- ncol(tpm_matrix)
n_genes <- nrow(tpm_matrix)

sample_ids <- phenoDF$sample.id
biopsy_site <- trimws(setNames(phenoDF$biopsy.site, phenoDF$sample.id))
sample_type <- setNames(phenoDF$sample.type, phenoDF$sample.id)
site_map <- gather_site_map(
  site_map_file = "datasets/biopsy_site_categories.csv",
  site_map_with_id_file = file.path(output_dir, "biopsy_site_categories_with_id.csv"),
  category_columns = c("system_level_1", "system_level_2", "organ_level_3"),
  id_column = "class_id"
)
class_map <- setNames(site_map$class_id, site_map$site)
sample_class <- setNames(class_map[biopsy_site], names(biopsy_site))
sample_system_level_1 <- setNames(site_map$system_level_1[match(biopsy_site, site_map$site)], names(biopsy_site))
sample_system_level_2 <- setNames(site_map$system_level_2[match(biopsy_site, site_map$site)], names(biopsy_site))
sample_organ_level_3 <- setNames(site_map$organ_level_3[match(biopsy_site, site_map$site)], names(biopsy_site))

cat("Number of classes:", length(unique(sample_class)), "\n")
cat("Number of system level 1 classes:", length(unique(sample_system_level_1)), "\n")
cat("Number of system level 2 classes:", length(unique(sample_system_level_2)), "\n")
cat("Number of organ level 3 classes:", length(unique(sample_organ_level_3)), "\n")

# Ensure tpm_matrix has the correct column order
tpm_matrix <- tpm_matrix[, sample_ids]

# Part 1: Data preparation (same as train.R)
set.seed(42)

create_strata_set <- function(strata, biopsy_site, sample_type,
                              exclude_sample_types = NULL,
                              exclude_biopsy_sites = NULL) {
  unique_strata <- unique(strata)
  if (is.null(exclude_sample_types) && is.null(exclude_biopsy_sites)) {
    return(list(
      partitioned_strata = unique_strata,
      excluded_strata = character(0)
    ))
  }
  if (is.null(exclude_sample_types)) {
    exclude_sample_types <- unique(sample_type)
  }
  if (is.null(exclude_biopsy_sites)) {
    exclude_biopsy_sites <- unique(biopsy_site)
  }
  excluded_table <- expand.grid(exclude_biopsy_sites, exclude_sample_types)
  excluded_strata <- apply(excluded_table, 1, paste, collapse = "_")
  unique_strata <- unique_strata[!unique_strata %in% excluded_strata]
  return(list(
    partitioned_strata = unique_strata,
    excluded_strata = excluded_strata
  ))
}


# Function to perform stratified partitioning of samples (same as train.R)
create_stratified_partition <- function(sample_ids,
                                        biopsy_site,
                                        sample_type,
                                        train_fraction = 0.9,
                                        exclude_sample_types = NULL,
                                        exclude_biopsy_sites = NULL) {
  # Ensure all inputs have the same length
  if (length(sample_ids) != length(biopsy_site) || length(sample_ids) != length(sample_type)) {
    stop("All input vectors must have the same length")
  }
  
  # Create stratification groups based on biopsy site and sample type combinations
  strata <- paste(biopsy_site, sample_type, sep = "_")
  
  # Initialize vectors to store sample IDs
  train_test_samples <- c()
  holdout_samples <- c()
  
  all_strata_sets <- create_strata_set(strata, biopsy_site, sample_type, exclude_sample_types, exclude_biopsy_sites)
  partitioned_strata <- all_strata_sets$partitioned_strata
  excluded_strata <- all_strata_sets$excluded_strata
  
  # Stratified sampling for each combination
  for (stratum in partitioned_strata) {
    # Get indices for current stratum
    stratum_indices <- which(strata == stratum)
    stratum_samples <- sample_ids[stratum_indices]
    n_samples_stratum <- length(stratum_samples)
    
    # If stratum has less than 5 samples, add all to training set
    if (n_samples_stratum < 5) {
      train_test_samples <- c(train_test_samples, stratum_samples)
      next
    }
    
    # For larger strata, randomly sample train_fraction for train/test set
    n_train_test <- ceiling(train_fraction * n_samples_stratum)
    
    # Randomly select samples for train/test set
    selected_samples <- sample(stratum_samples, size = n_train_test)
    
    # Add to respective vectors
    train_test_samples <- c(train_test_samples, selected_samples)
    holdout_samples <- c(
      holdout_samples,
      setdiff(stratum_samples, selected_samples)
    )
  }
  if (length(excluded_strata) > 0) {
    # Add excluded strata to the training set
    stratum_indices <- which(strata %in% excluded_strata)
    stratum_samples <- sample_ids[stratum_indices]
    train_test_samples <- unique(c(train_test_samples, stratum_samples))
  }
  
  # Sort sample IDs for consistency
  train_test_samples <- sort(train_test_samples)
  holdout_samples <- sort(holdout_samples)
  
  # Return list with both sets of sample IDs
  return(list(train_test = train_test_samples, holdout = holdout_samples))
}

# Create the holdout set using only primary, metastatic, and recurrent sample types
partition <- create_stratified_partition(
  sample_ids, biopsy_site, sample_type,
  exclude_sample_types = NULL # c("normal", "adjacent")
)
train_test_samples <- partition$train_test
holdout_samples <- partition$holdout

# Create train/test and holdout sets
train_test_set <- tpm_matrix[, train_test_samples]
holdout_set <- tpm_matrix[, holdout_samples]

# Create a train/test partition for model training
biopsy_site_train <- biopsy_site[train_test_samples]
sample_type_train <- sample_type[train_test_samples]
train_partition <- create_stratified_partition(
  train_test_samples, biopsy_site_train, sample_type_train, 0.9
)
train_samples <- train_partition$train_test
test_samples <- train_partition$holdout

# Create train/test and holdout sets for cross-validation
train_set <- tpm_matrix[, train_samples]
test_set <- tpm_matrix[, test_samples]

# Print final data dimensions for debugging
cat("Final data dimensions:\n")
cat(paste("Train set:", paste(dim(train_set), collapse = " x ")), "\n")
cat(paste("Test set:", paste(dim(test_set), collapse = " x ")), "\n")
cat(paste("Holdout set:", paste(dim(holdout_set), collapse = " x ")), "\n")

train_sample_class <- sample_class[train_samples]
test_sample_class <- sample_class[test_samples]
holdout_sample_class <- sample_class[holdout_samples]

# Save the sample IDs of the train, test, and holdout sets
save(train_samples, test_samples, holdout_samples, file = file.path(output_dir, "sample_ids_positional_encoding.rda"))

# Convert data to arrays and transpose
cat("Converting and transposing data...\n")
train_data <- t(as.array(train_set))
test_data <- t(as.array(test_set))
holdout_data <- t(as.array(holdout_set))

# Print transposed data dimensions for debugging
cat("Transposed data dimensions:\n")
cat(paste("Train data:", paste(dim(train_data), collapse = " x ")), "\n")
cat(paste("Test data:", paste(dim(test_data), collapse = " x ")), "\n")
cat(paste("Holdout data:", paste(dim(holdout_data), collapse = " x ")), "\n")

# Part 2: Positional Encoding Implementation

# Load required libraries
library(keras)
library(tensorflow)

# Define positional encoding layer
positional_encoding_layer <- function(seq_len, d_model) {
  function(inputs) {
    # Cast dimensions to integers
    seq_len <- as.integer(seq_len)
    d_model <- as.integer(d_model)
    
    # Create position indices
    even_indices <- tf$range(0L, d_model, 2L)
    odd_indices <- tf$range(1L, d_model, 2L)
    
    # Initialize encoding matrix with zeros
    encoding_matrix <- tf$zeros(shape = c(seq_len, d_model))
    
    # Calculate angles for even positions
    angles_even <- tf$cast(
      tf$tile(tf$range(0L, seq_len), tf$reshape(tf$size(even_indices), shape = c(1L))),
      tf$float32
    ) / tf$pow(
      10000.0,
      tf$cast(tf$tile(even_indices, tf$reshape(seq_len, shape = c(1L))), tf$float32) / tf$cast(d_model, tf$float32)
    )
    
    # Calculate angles for odd positions
    angles_odd <- tf$cast(
      tf$tile(tf$range(0L, seq_len), tf$reshape(tf$size(odd_indices), shape = c(1L))),
      tf$float32
    ) / tf$pow(
      10000.0,
      tf$cast(tf$tile(odd_indices, tf$reshape(seq_len, shape = c(1L))), tf$float32) / tf$cast(d_model, tf$float32)
    )
    
    # Fill even positions with sine values
    indices_even <- tf$stack(list(
      tf$tile(tf$range(0L, seq_len), tf$reshape(tf$size(even_indices), shape = c(1L))),
      tf$tile(even_indices, tf$reshape(seq_len, shape = c(1L)))
    ), axis = 1L)
    encoding_matrix <- tf$tensor_scatter_nd_update(
      encoding_matrix,
      indices_even,
      tf$sin(angles_even)
    )
    
    # Fill odd positions with cosine values
    indices_odd <- tf$stack(list(
      tf$tile(tf$range(0L, seq_len), tf$reshape(tf$size(odd_indices), shape = c(1L))),
      tf$tile(odd_indices, tf$reshape(seq_len, shape = c(1L)))
    ), axis = 1L)
    encoding_matrix <- tf$tensor_scatter_nd_update(
      encoding_matrix,
      indices_odd,
      tf$cos(angles_odd)
    )
    
    # Add batch dimension and broadcast
    encoding_matrix <- tf$expand_dims(encoding_matrix, axis = 0L)
    batch_size <- tf$shape(inputs)[1L]
    encoding_matrix <- tf$tile(encoding_matrix, tf$stack(list(batch_size, 1L, 1L)))
    
    # Add positional encoding to inputs
    inputs + encoding_matrix
  }
}

# Define improved transformer encoder block with positional encoding
transformer_encoder_block <- function(inputs, d_model = 512, num_heads = 8, ff_dim = 512, dropout_rate = 0.1) {
  # Layer normalization before attention
  x <- layer_layer_normalization(epsilon = 1e-6)(inputs)
  
  # Add positional encoding
  pos_encoding <- positional_encoding_layer(1L, d_model)(x)
  x <- layer_add()(list(x, pos_encoding))
  
  # Multi-head attention
  attention_output <- layer_multi_head_attention(
    num_heads = num_heads,
    key_dim = d_model %/% num_heads
  )(x, x)
  
  # Add & Norm with residual connection
  x <- layer_add()(list(attention_output, x))
  x <- layer_layer_normalization(epsilon = 1e-6)(x)
  x <- layer_dropout(rate = dropout_rate)(x)
  
  # Feed-forward network
  ffn_output <- layer_dense(units = ff_dim)(x)
  ffn_output <- layer_batch_normalization()(ffn_output)
  ffn_output <- layer_activation_leaky_relu(alpha = 0.2)(ffn_output)
  ffn_output <- layer_dropout(rate = dropout_rate)(ffn_output)
  
  ffn_output <- layer_dense(units = d_model)(ffn_output)
  ffn_output <- layer_batch_normalization()(ffn_output)
  
  # Add & Norm with residual connection
  x <- layer_add()(list(ffn_output, x))
  x <- layer_layer_normalization(epsilon = 1e-6)(x)
  x <- layer_dropout(rate = dropout_rate)(x)
  
  return(x)
}

# Build improved encoder with balanced capacity and positional encoding
build_encoder <- function(input_shape, latent_dim = 64) {
  inputs <- layer_input(shape = input_shape)
  
  # Initial dense layers
  x <- layer_dense(units = 512)(inputs)
  x <- layer_batch_normalization()(x)
  x <- layer_activation_leaky_relu(alpha = 0.2)(x)
  x <- layer_dropout(rate = 0.1)(x)
  
  # Reshape for transformer blocks
  x <- layer_reshape(target_shape = c(1, 512))(x)
  
  # Multiple transformer encoder blocks
  for (i in 1:3) {
    x <- transformer_encoder_block(x, d_model = 512)
  }
  
  # Flatten and final dense layers
  x <- layer_flatten()(x)
  
  x <- layer_dense(units = 256)(x)
  x <- layer_batch_normalization()(x)
  x <- layer_activation_leaky_relu(alpha = 0.2)(x)
  x <- layer_dropout(rate = 0.1)(x)
  
  x <- layer_dense(units = 128)(x)
  x <- layer_batch_normalization()(x)
  x <- layer_activation_leaky_relu(alpha = 0.2)(x)
  x <- layer_dropout(rate = 0.1)(x)
  
  outputs <- layer_dense(units = latent_dim)(x)
  
  return(keras_model(inputs, outputs))
}

# Build improved decoder with balanced capacity and positional encoding
build_decoder <- function(latent_dim = 64, output_shape) {
  inputs <- layer_input(shape = latent_dim)
  
  # Initial dense layers
  x <- layer_dense(units = 128)(inputs)
  x <- layer_batch_normalization()(x)
  x <- layer_activation_leaky_relu(alpha = 0.2)(x)
  x <- layer_dropout(rate = 0.1)(x)
  
  x <- layer_dense(units = 256)(x)
  x <- layer_batch_normalization()(x)
  x <- layer_activation_leaky_relu(alpha = 0.2)(x)
  x <- layer_dropout(rate = 0.1)(x)
  
  # Reshape for transformer blocks
  x <- layer_dense(units = 512)(x)
  x <- layer_reshape(target_shape = c(1, 512))(x)
  
  # Multiple transformer encoder blocks (using encoder blocks for simplicity)
  for (i in 1:3) {
    x <- transformer_encoder_block(x, d_model = 512)
  }
  
  # Flatten and final dense layers
  x <- layer_flatten()(x)
  
  x <- layer_dense(units = 1024)(x)
  x <- layer_batch_normalization()(x)
  x <- layer_activation_leaky_relu(alpha = 0.2)(x)
  x <- layer_dropout(rate = 0.1)(x)
  
  outputs <- layer_dense(units = output_shape)(x)
  
  return(keras_model(inputs, outputs))
}

# Build classification heads with exact number of classes
build_classifier_head <- function(encoded_input, num_classes, name) {
  x <- layer_dense(units = 64)(encoded_input)
  x <- layer_batch_normalization()(x)
  x <- layer_activation_leaky_relu(alpha = 0.2)(x)
  x <- layer_dropout(rate = 0.1)(x)
  
  x <- layer_dense(units = 32)(x)
  x <- layer_batch_normalization()(x)
  x <- layer_activation_leaky_relu(alpha = 0.2)(x)
  x <- layer_dropout(rate = 0.1)(x)
  
  outputs <- layer_dense(units = num_classes, activation = "softmax", name = name)(x)
  return(outputs)
}

# Build the improved autoencoder with positional encoding and classification heads
build_autoencoder_with_classifiers <- function(input_shape, latent_dim = 64, num_classes = 132) {
  # Build encoder and decoder
  encoder <- build_encoder(input_shape, latent_dim)
  decoder <- build_decoder(latent_dim, input_shape)
  
  # Build autoencoder
  inputs <- layer_input(shape = input_shape)
  encoded <- encoder(inputs)
  decoded <- decoder(encoded)
  
  # Build classification head with exact number of classes
  sample_classifier <- build_classifier_head(encoded, num_classes, "sample_classifier")
  
  # Create the full model
  full_model <- keras_model(
    inputs = inputs,
    outputs = list(
      decoded,
      sample_classifier
    )
  )
  
  return(list(
    full_model = full_model,
    encoder = encoder,
    decoder = decoder
  ))
}

# Get number of classes for each classification task from the complete dataset
num_classes <- length(unique(sample_class))

# Create and compile the improved model with positional encoding and classifiers
models <- build_autoencoder_with_classifiers(input_shape = n_genes, latent_dim = 64, num_classes = num_classes)
full_model <- models$full_model
encoder <- models$encoder
decoder <- models$decoder

# Create improved learning rate schedule with better warmup and decay
lr_schedule <- function(epoch) {
  initial_lr <- 0.0005 # Reduced initial learning rate
  warmup_epochs <- 20 # Longer warmup
  decay_epochs <- 380
  min_lr <- 1e-6
  
  if (epoch < warmup_epochs) {
    # Linear warmup from initial_lr/20 to initial_lr
    return(initial_lr * 0.05 + (initial_lr * 0.95 * epoch / warmup_epochs))
  } else {
    # Cosine decay with better minimum learning rate handling
    progress <- (epoch - warmup_epochs) / decay_epochs
    cosine_decay <- 0.5 * (1 + cos(pi * progress))
    return(min_lr + (initial_lr - min_lr) * cosine_decay)
  }
}

encode_labels <- function(labels, factor) {
  factor_levels <- levels(factor)
  factor_labels <- factor(labels, levels = factor_levels)
  return(to_categorical(as.integer(factor_labels) - 1, num_classes = length(factor_levels)))
}

get_labels <- function(map) (unique(map))

# Shuffle the labels to avoid any bias in the training data
sample_classes <- get_labels(sample_class)

# Create factor levels from all samples before splitting
sample_factors <- factor(sample_class, levels = sample_classes)

# Save the factor levels for later use
save(sample_factors, file = file.path(output_dir, "factor_levels.rda"))

# Print number of levels for verification
cat("Number of classes:", nlevels(sample_factors), "\n")

# Convert labels using consistent factor levels
train_labels_sample <- encode_labels(train_sample_class, sample_factors)

test_labels_sample <- encode_labels(test_sample_class, sample_factors)

holdout_labels_sample <- encode_labels(holdout_sample_class, sample_factors)

# Compile the model with gradient clipping and improved optimizer settings
full_model %>% compile(
  optimizer = optimizer_adam(
    learning_rate = lr_schedule(0),
    clipnorm = 1.0,
    beta_1 = 0.9,
    beta_2 = 0.999,
    epsilon = 1e-7
  ),
  loss = list(
    "mse", # Reconstruction loss
    "categorical_crossentropy" # Sample classifier
  ),
  loss_weights = list(
    1.0, # Reconstruction loss weight
    0.5 # Sample classifier weight
  ),
  metrics = list(
    model_1 = list("mae"),
    sample_classifier = list("categorical_accuracy") # Sample classifier
  )
)

# Print model summary
print(full_model)

# Define improved callbacks
callbacks <- list(
  callback_early_stopping(
    monitor = "val_loss",
    patience = 50,
    restore_best_weights = TRUE,
    min_delta = 1e-5
  ),
  callback_model_checkpoint(
    filepath = file.path(output_dir, "best_model.h5"),
    monitor = "val_loss",
    save_best_only = TRUE
  ),
  callback_tensorboard(
    log_dir = file.path(output_dir, "logs"),
    histogram_freq = 1,
    write_graph = TRUE
  ),
  callback_reduce_lr_on_plateau(
    monitor = "val_loss",
    factor = 0.5,
    patience = 15,
    min_lr = 1e-6,
    min_delta = 1e-5
  )
)

# Train the improved model with better batch size
history <- full_model %>% fit(
  x = train_data,
  y = list(
    train_data, # Reconstruction target
    train_labels_sample # Sample classification target
  ),
  epochs = 400,
  batch_size = 64, # Increased batch size for better stability
  validation_data = list(
    test_data, # Input data
    list(
      test_data, # Reconstruction target
      test_labels_sample # Sample classification target
    )
  ),
  callbacks = callbacks,
  verbose = 1,
  shuffle = TRUE
)

# Create training plots
library(ggplot2)
library(tidyr)
library(dplyr)

# Convert history to a data frame
history_df <- as.data.frame(history$metrics)
history_df$epoch <- seq_len(nrow(history_df))

# Reshape data for plotting
plot_data <- history_df %>%
  pivot_longer(
    cols = -epoch,
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    set = ifelse(grepl("val_", metric), "Validation", "Training"),
    metric = gsub("val_", "", metric)
  )

# Create the plot with white background
p <- ggplot(plot_data, aes(x = epoch, y = value, color = set)) +
  geom_line(linewidth = 1) +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Training Statistics (Positional Encoding Model with Classifiers)",
    x = "Epoch",
    y = "Value",
    color = "Dataset"
  ) +
  theme_minimal()
  # theme(
  #   plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
  #   axis.title = element_text(size = 12),
  #   axis.text = element_text(size = 10),
  #   legend.title = element_text(size = 12),
  #   legend.text = element_text(size = 10),
  #   strip.text = element_text(size = 12, face = "bold"),
  #   panel.background = element_rect(fill = "white"),
  #   plot.background = element_rect(fill = "white")
  # )

# Save the plot
ggsave(file.path(output_dir, "training_curves.pdf"), p, width = 10, height = 6, bg = "white")
ggsave(
  file.path(output_dir, "training_curves.png"), p,
  width = 10, height = 6, dpi = 300, bg = "white"
)

# Print final metrics
cat("\nFinal Training Metrics:\n")
cat(paste("Loss:", tail(history_df$loss, 1)), "\n")
cat(paste("MAE:", tail(history_df$model_1_mae, 1)), "\n")
cat("\nFinal Validation Metrics:\n")
cat(paste("Loss:", tail(history_df$val_loss, 1)), "\n")
cat(paste("MAE:", tail(history_df$val_model_1_mae, 1)), "\n")

# Save the final models
save_model_hdf5(encoder, file.path(output_dir, "final_encoder.h5"))
save_model_hdf5(decoder, file.path(output_dir, "final_decoder.h5"))
save_model_hdf5(full_model, file.path(output_dir, "final_full_model.h5"))

# Save individual classifiers
# Create separate models for each classifier
sample_classifier_model <- keras_model(
  inputs = full_model$input,
  outputs = full_model$get_layer("sample_classifier")$output
)

# Save individual classifiers
save_model_hdf5(sample_classifier_model, file.path(output_dir, "final_sample_classifier.h5"))

# Evaluate on holdout set
holdout_evaluation <- full_model %>% evaluate(
  holdout_data,
  list(
    holdout_data,
    holdout_labels_sample
  ),
  verbose = 0
)
cat("\nHoldout set evaluation:\n")
print(holdout_evaluation)

# Calculate and print detailed performance metrics for each classifier
library(caret)

# Function to calculate and print metrics for a classifier
print_classifier_metrics <- function(predictions, true_labels, classifier_name, classes_factor) {
  classes_names <- levels(classes_factor)
  # Convert predictions to class labels
  pred_classes <- classes_names[apply(predictions, 1, which.max)]
  true_classes <- classes_names[apply(true_labels, 1, which.max)]
  
  # Create confusion matrix
  cm <- confusionMatrix(
    factor(pred_classes, levels = classes_names),
    factor(true_classes, levels = classes_names)
  )
  
  # Print overall metrics
  cat("\nPerformance metrics for", classifier_name, "classifier:\n")
  cat("Overall Accuracy:", cm$overall["Accuracy"], "\n")
  
  # Print per-class metrics
  cat("\nPer-class metrics:\n")
  for (i in seq_len(ncol(true_labels))) {
    class_name <- classes_names[i]
    precision <- cm$byClass[i, "Precision"]
    recall <- cm$byClass[i, "Recall"]
    f1_score <- cm$byClass[i, "F1"]
    if (is.na(precision) || precision == 0) {
      if (is.na(recall) && is.na(f1_score)) {
        next()
      }
    }
    cat("\nClass:", class_name, "\n")
    cat("Precision:", precision, "\n")
    cat("Recall:", recall, "\n")
    cat("F1-score:", f1_score, "\n")
  }
  
  # Print confusion matrix
  cat("\nConfusion Matrix:\n")
  print(cm$table)
}

# Get predictions for holdout set
holdout_predictions <- full_model %>% predict(holdout_data)

# Calculate metrics for each classifier
print_classifier_metrics(
  holdout_predictions[[2]], # Sample classifier predictions
  holdout_labels_sample,
  "Sample",
  sample_factors
)

# Save training history
save(history, file = file.path(output_dir, "training_history.rda"))

# Load required packages for visualization
library(Rtsne)
library(ggplot2)
library(dplyr)
library(tidyr)

# Function to create t-SNE plot
create_tsne_plot <- function(encoded_data, metadata, color_by, title) {
  # Remove duplicates from encoded data and keep track of unique indices
  unique_indices <- !duplicated(encoded_data)
  unique_encoded_data <- encoded_data[unique_indices, ]
  unique_metadata <- metadata[unique_indices]
  
  # Perform t-SNE
  tsne_result <- Rtsne(unique_encoded_data, perplexity = 30, theta = 0.5, dims = 2)
  
  # Create data frame for plotting
  plot_data <- data.frame(
    tSNE1 = tsne_result$Y[, 1],
    tSNE2 = tsne_result$Y[, 2],
    color = unique_metadata
  )
  
  # Create the plot with improved legend layout and white background
  p <- ggplot(plot_data, aes(x = tSNE1, y = tSNE2, color = color)) +
    geom_point(alpha = 0.6, size = 1) +
    theme_minimal() +
    labs(
      title = title,
      color = color_by
    ) +
    theme_minimal() +
    guides(color = guide_legend(
      ncol = 1,
      override.aes = list(size = 2),
      keywidth = 0.5,
      keyheight = 0.5
    ))
  
  # theme(
  #   plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
  #   axis.title = element_text(size = 12),
  #   axis.text = element_text(size = 10),
  #   legend.title = element_text(size = 12),
  #   legend.text = element_text(size = 8),
  #   legend.position = "right",
  #   legend.box.just = "left",
  #   legend.margin = margin(0, 0, 0, 0),
  #   legend.box.margin = margin(0, 0, 0, 0),
  #   legend.key.size = unit(0.5, "lines"),
  #   legend.spacing.y = unit(0.2, "lines"),
  #   plot.margin = margin(10, 10, 10, 10),
  #   panel.background = element_rect(fill = "white"),
  #   plot.background = element_rect(fill = "white")
  # ) 
  
  return(p)
}

# Get encoded representations for all datasets
train_encoded <- predict(encoder, train_data)
test_encoded <- predict(encoder, test_data)
holdout_encoded <- predict(encoder, holdout_data)

# Combine all encoded data
all_encoded <- rbind(train_encoded, test_encoded, holdout_encoded)

# save the encoded data
save(all_encoded, file = file.path(output_dir, "encoded_data.rda"))

# Create metadata vectors for coloring
all_samples <- c(train_samples, test_samples, holdout_samples)
all_biopsy_sites <- biopsy_site[all_samples]
all_sample_types <- sample_type[all_samples]
all_system_level_1 <- sample_system_level_1[all_samples]
all_system_level_2 <- sample_system_level_2[all_samples]
all_organ_level_3 <- sample_organ_level_3[all_samples]

# Create t-SNE plots with adjusted dimensions
tsne_biopsy <- create_tsne_plot(
  all_encoded,
  all_biopsy_sites,
  "Biopsy Site",
  "t-SNE of Encoded Data (Positional Encoding Model, Colored by Biopsy Site)"
)

tsne_sample_type <- create_tsne_plot(
  all_encoded,
  all_sample_types,
  "Sample Type",
  "t-SNE of Encoded Data (Positional Encoding Model, Colored by Sample Type)"
)

tsne_system_level_1 <- create_tsne_plot(
  all_encoded,
  all_system_level_1,
  "System Level 1",
  "t-SNE of Encoded Data (Positional Encoding Model, Colored by System Level 1)"
)

tsne_system_level_2 <- create_tsne_plot(
  all_encoded,
  all_system_level_2,
  "System Level 2",
  "t-SNE of Encoded Data (Positional Encoding Model, Colored by System Level 2)"
)

tsne_organ_level_3 <- create_tsne_plot(
  all_encoded,
  all_organ_level_3,
  "Organ Level 3",
  "t-SNE of Encoded Data (Positional Encoding Model, Colored by Organ Level 3)"
)

# Save the plots with adjusted dimensions and white background
ggsave(
  file.path(output_dir, "tsne_biopsy_site.png"), tsne_biopsy,
  width = 12, height = 8, dpi = 300, bg = "white"
)
ggsave(
  file.path(output_dir, "tsne_sample_type.png"), tsne_sample_type,
  width = 10, height = 8, dpi = 300, bg = "white"
)
ggsave(
  file.path(output_dir, "tsne_system_level_1.png"), tsne_system_level_1,
  width = 12, height = 8, dpi = 300, bg = "white"
)
ggsave(
  file.path(output_dir, "tsne_system_level_2.png"), tsne_system_level_2,
  width = 12, height = 8, dpi = 300, bg = "white"
)
ggsave(
  file.path(output_dir, "tsne_organ_level_3.png"), tsne_organ_level_3,
  width = 12, height = 8, dpi = 300, bg = "white"
)

save.image(file.path(output_dir, "all_data.rda"))

# Close the output capture at the end of the script
sink()
cat("Console output has been saved to:", output_file, "\n")
