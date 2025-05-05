library(keras)
library(tensorflow)
library(progress)
library(dplyr)
library(openxlsx)

# Path to the model
# sample_ids_path <- "./sample_ids_positional_encoding.rda"
model_path <- "./models_positional_encoding_with_classifier"

# Set up output capture
output_file <- file.path(model_path, "independent_test_output.txt")
if (file.exists(output_file)) {
  file.remove(output_file)
}
sink(file = output_file, append = TRUE, split = TRUE) # split=TRUE shows output in console and file
warn <- options(warn = -1)
on.exit(options(warn))

load(file.path(model_path, "factor_levels.rda"))

get_class_names <- function(predictions, factor_levels) {
  classes_names <- levels(factor_levels)
  classes_names[apply(predictions, 1, which.max)]
}

classifier_level_1 <- load_model_hdf5(file.path(model_path, "final_classifier_level_1.h5"))
classifier_level_2 <- load_model_hdf5(file.path(model_path, "final_classifier_level_2.h5"))
classifier_level_3 <- load_model_hdf5(file.path(model_path, "final_classifier_level_3.h5"))

# Load data in ipf.rds
ipf_data <- readRDS("datasets/ipf.rds")

# Convert data to arrays and transpose
print("Converting and transposing data...")
ipf_data <- t(as.array(ipf_data))

ipf_data_level_1_predictions <- predict(classifier_level_1, ipf_data)
ipf_data_level_2_predictions <- predict(classifier_level_2, ipf_data)
ipf_data_level_3_predictions <- predict(classifier_level_3, ipf_data)

ipf_data_level_1_predicted_classes <- get_class_names(ipf_data_level_1_predictions, level_1_factors)
ipf_data_level_2_predicted_classes <- get_class_names(ipf_data_level_2_predictions, level_2_factors)
ipf_data_level_3_predicted_classes <- get_class_names(ipf_data_level_3_predictions, level_3_factors)

# Create a data frame with the predicted labels
ipf_predicted_df <- data.frame(
  sample_id = rownames(ipf_data),
  level_1_prediction = ipf_data_level_1_predicted_classes,
  level_2_prediction = ipf_data_level_2_predicted_classes,
  level_3_prediction = ipf_data_level_3_predicted_classes
)

# Save the predicted labels to a file
write.table(ipf_predicted_df, "ipf_predicted_labels_from_classifier.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)

ipf_predicted_df$sample_id <- gsub("_GeneCount", "", ipf_predicted_df$sample_id)
meta <- readRDS("datasets/ipf_metadata.rds")

ipf_predicted_df <- left_join(ipf_predicted_df, meta, by = "sample_id")

library(openxlsx)
write.xlsx(ipf_predicted_df, file = "ipf_predicted_labels.xlsx")

level_1_accuracy <- length(which(ipf_predicted_df$level_1_prediction == "Respiratory System")) / nrow(ipf_predicted_df)
cat("\n\nLevel 1 accuracy: ", level_1_accuracy, "\n")
cat("Level 1 classification:\n")
print(
  round((table(ipf_predicted_df$level_1_prediction) / nrow(ipf_predicted_df)) * 100, 2)
)

level_2_accuracy <- length(which(ipf_predicted_df$level_2_prediction == "Lung")) / nrow(ipf_predicted_df)
cat("\n\nLevel 2 accuracy: ", level_2_accuracy, "\n")
cat("Level 2 classification:\n")
print(
  round((table(ipf_predicted_df$level_2_prediction) / nrow(ipf_predicted_df)) * 100, 2)
)

level_3_accuracy <- length(which(ipf_predicted_df$level_3_prediction == "Lung")) / nrow(ipf_predicted_df)
cat("\n\nLevel 3 accuracy: ", level_3_accuracy, "\n")
cat("Level 3 classification:\n")
print(
  round((table(ipf_predicted_df$level_3_prediction) / nrow(ipf_predicted_df)) * 100, 2)
)

###############################################################################
# Load the second independent test set

other_datasets_metadata <- readRDS("datasets/other_datasets_metadata.rds")
other_datasets_data     <- readRDS("datasets/other_datasets_tpm.rds")

# Convert data to arrays and transpose
print("Converting and transposing data...")
other_data <- t(as.array(other_datasets_data))

other_data_level_1_predictions <- predict(classifier_level_1, other_data)
other_data_level_2_predictions <- predict(classifier_level_2, other_data)
other_data_level_3_predictions <- predict(classifier_level_3, other_data)

other_data_level_1_predicted_classes <- get_class_names(other_data_level_1_predictions, level_1_factors)
other_data_level_2_predicted_classes <- get_class_names(other_data_level_2_predictions, level_2_factors)
other_data_level_3_predicted_classes <- get_class_names(other_data_level_3_predictions, level_3_factors)

# Create a data frame with the predicted labels
other_predicted_df <- data.frame(
  sample_id = rownames(other_data),
  level_1_prediction = other_data_level_1_predicted_classes,
  level_2_prediction = other_data_level_2_predicted_classes,
  level_3_prediction = other_data_level_3_predicted_classes
)

# Save the predicted labels to a file
write.table(other_predicted_df, "other_data_predicted_labels_from_classifier.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE
)


other_predicted_df <- left_join(other_predicted_df, other_datasets_metadata, by = "sample_id")

purity <- readRDS("datasets/other_datasets_purity.rds")

other_predicted_df <- left_join(other_predicted_df, purity, by = "sample_id")

write.xlsx(other_predicted_df, file = "other_predicted_labels.xlsx")


print_classifier_metrics <- function(predictions, true_labels, classifier_name) {
  library(caret)
  all_levels <- unique(c(predictions, true_labels))
  # Create confusion matrix
  cm <- confusionMatrix(
    factor(predictions, levels = all_levels),
    factor(true_labels, levels = all_levels),
  )
  # Print overall metrics
  cat("\nPerformance metrics for", classifier_name, "classifier:\n")
  cat("Overall Accuracy:", cm$overall["Accuracy"], "\n")
  # Print per-class metrics
  cat("\nPer-class metrics:\n")
  for (i in seq_len(length(all_levels))) {
    class_name <- all_levels[i]
    precision <- cm$byClass[i, "Precision"]
    recall <- cm$byClass[i, "Recall"]
    f1_score <- cm$byClass[i, "F1"]
    if (is.na(precision) || precision == 0) {
      if (is.na(recall) && is.na(f1_score)) {
        next()
      }
    }
    cat("\nClass:", class_name, "\n")
    cat("Precision:", cm$byClass[i, "Precision"], "\n")
    cat("Recall:", cm$byClass[i, "Recall"], "\n")
    cat("F1-score:", cm$byClass[i, "F1"], "\n")
  }
  
  # Print confusion matrix
  cat("\nConfusion Matrix:\n")
  print(cm$table)
}

print_classifier_metrics(
  predictions = other_predicted_df$level_1_prediction,
  true_labels = other_predicted_df$system_level_1,
  classifier_name = "Level 1"
)
print_classifier_metrics(
  predictions = other_predicted_df$level_2_prediction,
  true_labels = other_predicted_df$system_level_2,
  classifier_name = "Level 2"
)
print_classifier_metrics(
  predictions = other_predicted_df$level_3_prediction,
  true_labels = other_predicted_df$organ_level_3,
  classifier_name = "Level 3"
)

# Close the output capture at the end of the script
sink()
