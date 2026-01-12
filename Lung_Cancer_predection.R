# -------------------------------------------#
# Section 0: Load Required Libraries---
# -------------------------------------------#
# If you wanna save time run the code from section 5, line 141 by reading the final dataset attached to this script

library(dplyr) # Data wrangling (filtering, merging, summarizing)
library(tidyr) # Data wrangling 
library(RCurl) # Downloading and unzipping GEO SOFT & matrix files
library(R.utils) # Downloading and unzipping GEO SOFT & matrix files
library(reshape2) # Long-to-wide format conversion for plotting
library(ggplot2) # Visualization (PCA, boxplot, ROC, calibration)
library(limma) # Differential expression analysis (moderated t-test)
library(glmnet) # LASSO regression for feature selection
library(caret) # Model training, cross-validation, performance metrics
library(pROC) # ROC curve and AUC evaluation
library(pheatmap)
library(RColorBrewer)
library(splitstackshape) # For stratified sampling
library(ggrepel) # Labeling volcano plot points
library(plotly) # Interactive 3D PCA plotting
library(e1071) # SVM model backend for caret
library(randomForest) # Random Forest classifier
library(cowplot) # Visualization (3D-PCA)
# -------------------------------------------#
# Section 1: Build miRNA Mapping from SOFT File----
# -------------------------------------------#
# Download and decompress the SOFT file (if not already done)
soft_file_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137140/soft/GSE137140_family.soft.gz"
#download.file(soft_file_url, destfile = "GSE137140_family.soft.gz", mode = "wb")
gunzip("GSE137140_family.soft.gz", remove = FALSE)

# Read the entire SOFT file as text
soft_data <- readLines("GSE137140_family.soft")

# Locate the platform table (miRNA mapping table)
table_start <- grep("!platform_table_begin", soft_data) + 1

# Read the mapping table from the SOFT file, skipping metadata lines
miRNA_data <- read.delim("GSE137140_family.soft", header = FALSE, sep = "\t", skip = table_start)

# Assign column names: first column = Probe ID (MIMAT), third column = real miRNA names
colnames(miRNA_data) <- c("ID", "miRNA", "miRNA_ID_LIST")
head(miRNA_data)

# If a probe maps to multiple miRNAs, keep only the first miRNA name:
miRNA_data$miRNA_ID_LIST <- sapply(strsplit(miRNA_data$miRNA_ID_LIST, ", "), `[`, 1)

# Create mapping table: keys are probe IDs and values are real miRNA names
miRNA_map <- setNames(miRNA_data$miRNA_ID_LIST, miRNA_data$ID)
head(miRNA_map)

# Save miRNA list
write.table(miRNA_data, file = "miRNA_data.txt", sep = "\t", row.names = FALSE, quote = FALSE)
# -------------------------------------------#
# Section 2: Extract Sample Metadata ----
# -------------------------------------------#
## 1-Extract Sample IDs ----
# Sample IDs are stored in lines starting with "!Series_sample_id"
sample_id_lines <- soft_data[grep("^!Series_sample_id", soft_data)]
head(sample_id_lines)

# Remove the "!Series_sample_id = " prefix
sample_ids <- gsub("^!Series_sample_id = ", "", sample_id_lines)

# Trim any extra spaces
sample_ids <- trimws(sample_ids)

# Verify extracted sample IDs
head(sample_ids)
length(sample_ids) # 3924

## 2-Extract Disease State Information ----
# Download the series matrix file (if not already done)
matrix_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137140/matrix/GSE137140_series_matrix.txt.gz"
download.file(matrix_url, destfile = "GSE137140_series_matrix.txt.gz", mode = "wb")

# Read the series matrix file as text for metadata extraction
matrix_lines <- readLines("GSE137140_series_matrix.txt.gz")
# Disease state info is stored in "!Sample_characteristics_ch1" in the series matrix
disease_state_line <- matrix_lines[grep("^!Sample_characteristics_ch1", matrix_lines)]
disease_states <- unlist(strsplit(disease_state_line, "\t"))[-1]
disease_states <- gsub("\"", "", disease_states)
disease_states <- gsub("disease state: ", "", disease_states)
disease_states <- trimws(disease_states)
# Filter out irrelevant metadata (e.g., "tissue: serum") by keeping only entries that contain key words
filtered_disease_states <- disease_states[grep("Non-cancer control|Lung cancer", disease_states)]
head(filtered_disease_states)

# Adjust length if necessary to match sample_ids:
if (length(filtered_disease_states) > length(sample_ids)) {
  filtered_disease_states <- filtered_disease_states[1:length(sample_ids)]
} else if (length(filtered_disease_states) < length(sample_ids)) {
  warning("Fewer disease state entries than sample IDs; please verify the metadata extraction.")
}

# Create metadata dataframe
metadata_df <- data.frame(Sample_ID = sample_ids, Cancer_Status = filtered_disease_states, stringsAsFactors = FALSE)

# Simplify labels 
metadata_df$Cancer_Status <- gsub("Lung cancer, pre-operation", "Lung_cancer", metadata_df$Cancer_Status)
metadata_df$Cancer_Status <- gsub("Non-cancer control", "Non_cancer", metadata_df$Cancer_Status)
str(metadata_df)
head(metadata_df)

# -------------------------------------------#
# Section 3: Read and Process Expression Data from Series Matrix File----
# -------------------------------------------#
# Read the series matrix file (expression data) with comment lines removed
expr_data <- read.delim("GSE137140_series_matrix.txt.gz", header = TRUE, sep = "\t", comment.char = "!")
# Remove the first column (ID_REF) as it contains probe IDs (MIMAT IDs)
probe_ids <- expr_data[, 1]
rownames(expr_data) <- probe_ids
expr_data <- expr_data[, -1]

# Replace row names (MIMAT IDs) with real miRNA names using mapping table
rownames(expr_data) <- miRNA_map[rownames(expr_data)]

# For any NA mappings, retain the original probe IDs
na_ids <- is.na(rownames(expr_data))
rownames(expr_data)[na_ids] <- probe_ids[na_ids]

# Transpose expression data so that rows are samples and columns are miRNAs
expr_data_t <- as.data.frame(t(expr_data))
expr_data_t$Sample_ID <- rownames(expr_data_t)


# -------------------------------------------#
# Section 4: Merge Expression Data with Metadata-----
# -------------------------------------------#
#Merge Expression Data with Metadata
final_dataset <- merge(expr_data_t, metadata_df, by = "Sample_ID")
# Remove unwanted samples (e.g., "Lung cancer, post-operation")
final_dataset <- subset(final_dataset, Cancer_Status != "Lung cancer, post-operation")
print(table(final_dataset$Cancer_Status))
#Lung_cancer  Non_cancer 
#1566        2178 
#head(final_dataset)

# -------------------------------------------#
# Section 5: Save and Re-read the Final Dataset----
# -------------------------------------------#
write.table(final_dataset, file = "GSE137140_final_dataset.txt", sep = "\t", row.names = FALSE, quote = FALSE)
final_dataset <- read.delim("GSE137140_final_dataset.txt", sep = "\t", header = TRUE)

####################################################################################################################################-
# -------------------------------------------#
# Section 6: Data Preparation ----
# -------------------------------------------#
# Load and preprocess data (from your previous steps)
final_dataset$Cancer_Status <- as.factor(final_dataset$Cancer_Status)

# 1. Identify feature columns safely
feature_cols <- setdiff(
  colnames(final_dataset),
  c("Sample_ID", "Cancer_Status")
)

# Check variance structure
feature_variances <- apply(final_dataset[, feature_cols], 2, var)
cat("\nVariance summary:\n")
print(summary(feature_variances))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.02489  3.27322  6.01594  5.38907  7.30307 20.54887
# Some miRNAs have very low variance (e.g., min = 0.024), meaning they barely change across samples — potentially uninformative.
#Others vary more widely — which is expected and useful for classification.

# Visualize variance distribution
ggplot(data.frame(variance = feature_variances), aes(x = variance)) +
  geom_histogram(bins = 50) +
  labs(title = "Distribution of Feature Variances")

# -------------------------------------------#
# Section 7: Filter Low-Expressed miRNAs----
# -------------------------------------------#
# Step 0: Ensure Cancer_Status is properly set
final_dataset$Cancer_Status <- as.character(final_dataset$Cancer_Status)
final_dataset$Cancer_Status[final_dataset$Cancer_Status == "Lung cancer"] <- "Lung_cancer"
final_dataset$Cancer_Status[final_dataset$Cancer_Status == "Non-cancer"] <- "Non_cancer"
final_dataset$Cancer_Status <- factor(final_dataset$Cancer_Status, levels = c("Non_cancer", "Lung_cancer"))


# Step 1: Extract only expression values (excluding Sample_ID and Cancer_Status)
expr_values <- final_dataset[, !(colnames(final_dataset) %in% c("Sample_ID", "Cancer_Status"))]

# Also extract sample group labels
group1 <- final_dataset$Cancer_Status == "Lung_cancer"
group2 <- final_dataset$Cancer_Status == "Non_cancer"

# Convert expression data to matrix
expr_matrix <- as.matrix(expr_values)

# Step 2: Define threshold parameters
expression_threshold <- 1.0   # Expression must exceed this value
min_fraction <- 0.5           # At least 50% of samples in a group must exceed the threshold

# Step 3: Keep miRNAs expressed in at least one group
keep <- apply(expr_matrix, 2, function(x) {
  mean(x[group1] > expression_threshold) >= min_fraction ||
    mean(x[group2] > expression_threshold) >= min_fraction
})

# Step 4: Filter expression matrix
expr_values_filtered <- expr_matrix[, keep]

# Convert back to data frame
expr_values_filtered <- as.data.frame(expr_values_filtered)

# Check result
cat("Original number of miRNAs:", ncol(expr_matrix), "\n")  
# Original number of miRNAs: 2565 
cat("Number of miRNAs after filtering:", ncol(expr_values_filtered), "\n")
# Number of miRNAs after filtering: 1504 


# -------------------------------------------#
# Section 8: Quantile Normalization----
#--------------------------------------------#
# Quantile normalization
final_dataset_nor <- normalizeBetweenArrays(as.matrix(expr_values_filtered), method = "quantile")

# Compute the median expression value for each sample after quantile normalization
sample_medians_normalized <- apply(final_dataset_nor, 2, median)
print(summary(sample_medians_normalized)) # Should be very similar across samples
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#3.909   3.913   3.913   3.913   3.913   3.916 

#They are now nearly identical across samples — a sign that quantile normalization worked well.
#This greatly reduces technical variation between samples, which is crucial for fair comparisons and modeling.
# We still want to remove miRNAs with low variance (via nearZeroVar) after normalization to get rid of uninformative features.


# Add back Sample_ID and Cancer_Status
final_dataset_nor <- cbind(
  final_dataset[, c("Sample_ID", "Cancer_Status")],
  expr_values_filtered
)

#----------------------------------------#
# Section 9: PCA plot on original dataset-----
#----------------------------------------#
# Step 1: Remove ID and class label to get expression matrix
expr_numeric1 <- final_dataset[, !(colnames(final_dataset) %in% c("Sample_ID", "Cancer_Status"))]

# Step 2: Perform PCA
pca_result_original <- prcomp(expr_numeric1, scale. = TRUE)

## 3D PCA Plot on the original data------
# Step 3: Create a data frame for the first 3 principal components
pca_df3d_original <- data.frame(
  PC1 = pca_result_original$x[, 1],
  PC2 = pca_result_original$x[, 2],
  PC3 = pca_result_original$x[, 3],
  Cancer_Status = final_dataset$Cancer_Status
)

# Step 4: Plot in 3D using plotly
plot_ly(
  data = pca_df3d_original,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~Cancer_Status,
  colors = c("Non_cancer" = "skyblue", "Lung_cancer" = "tomato"),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3, opacity = 0.7)
) %>%
  layout(
    title = "3D PCA Plot: Original Dataset",
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    )
  )

#-----------------------------------#
# Section 10: Perform PCA on normalized data-----------
#----------------------------------#
# Step 1: Remove ID and class label to get expression matrix
expr_numeric_nor <- final_dataset_nor[, !(colnames(final_dataset_nor) %in% c("Sample_ID", "Cancer_Status"))]

# Step 2: Apply PCA
pca_result_normalized <- prcomp(expr_numeric_nor, scale. = TRUE)

##3D PCA plot Filtering + Quantile-Normalized miRNA-----

# Step 3: Create PCA dataframe with top 3 PCs
pca_df3d_nor <- data.frame(
  PC1 = pca_result_normalized$x[, 1],
  PC2 = pca_result_normalized$x[, 2],
  PC3 = pca_result_normalized$x[, 3],
  Cancer_Status = final_dataset_nor$Cancer_Status
)

# Step 4: Interactive 3D PCA plot using plotly
plot_ly(
  data = pca_df3d_nor,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~Cancer_Status,
  colors = c("Non_cancer" = "skyblue", "Lung_cancer" = "tomato"),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3, opacity = 0.7)
) %>%
  layout(
    title = "3D PCA of miRNA Expression (Filtered + Quantile Normalized)",
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    )
  )

# -------------------------------------------#
# Section 11: Proper Train-Test Split ----
# -------------------------------------------#
# First, ensure your data is properly formatted
#final_dataset_nor <- as.data.frame(final_dataset_nor)
#final_dataset_nor$Cancer_Status <- as.factor(final_dataset_nor$Cancer_Status)
# Verify the structure
#str(final_dataset_nor)

# split data into training and test sets using caret
set.seed(123)
train_index <- createDataPartition(final_dataset_nor$Cancer_Status, p = 0.8, list = FALSE)

trainData <- final_dataset_nor[train_index, ]
testData  <- final_dataset_nor[-train_index, ]

# Verify the split
cat("Training set size:", nrow(trainData), "\n")
# Training set size: 2996
cat("Test set size:", nrow(testData), "\n")
# Test set size: 748
cat("\nClass distribution in training set:\n")
print(prop.table(table(trainData$Cancer_Status)))
# Lung cancer  Non-cancer
# 0.4182243  (~1095 samples)   0.5817757  (~15256 samples)
cat("\nClass distribution in test set:\n")
print(prop.table(table(testData$Cancer_Status)))
# Lung cancer  Non-cancer
# 0.4184492  (~470 samples)  0.5815508 (~653 samples)

# -------------------------------------------#
# Section 12: Feature Selection Using  LASSO regression with 10-fold CV (ONLY on training data) ----
# -------------------------------------------#
#Feature selection was performed using LASSO regression with 10-fold stability selection. LASSO was applied to 80% random subsamples of the training data using the glmnet R package. miRNAs selected in over 50% of the iterations were retained as stable features for downstream classification.
train_feature_selection <- function(train_data) {
  # Prepare predictor matrix and response vector
  X_train <- as.matrix(train_data[, !(colnames(train_data)) %in% c("Sample_ID", "Cancer_Status")])
  y_train <- train_data$Cancer_Status

  # Stability selection with repeated LASSO
  stable_features <- replicate(10, {
    subsample_idx <- sample(nrow(X_train), size = 0.8 * nrow(X_train))
    cv_fit <- cv.glmnet(
      X_train[subsample_idx, ],
      y_train[subsample_idx],
      family = "binomial",
      alpha = 1
    )
    coefs <- coef(cv_fit, s = "lambda.min")
    rownames(coefs)[coefs[, 1] != 0][-1] # Exclude intercept
  }) %>%
    unlist() %>%
    table()

  # Select features chosen in >50% of iterations
  names(stable_features[stable_features > 5])
}

selected_miRNAs <- train_feature_selection(trainData)
cat("Stably selected features:", length(selected_miRNAs), "\n")
# Stably selected features: 41

# -------------------------------------------#
# Section 13: Model Training and Comparison for SVM, GLM, RF ----
# -------------------------------------------#
# Ensure Sample_ID is removed before training
trainData_clean <- trainData[, !(colnames(trainData) %in% "Sample_ID")]

# Train control
ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)

# Train Models
# SVM model
svm_model <- train(
  Cancer_Status ~ ., 
  data = trainData_clean, 
  method = "svmRadial",
  preProcess = c("center", "scale"),
  metric = "ROC",
  trControl = ctrl, 
  tuneLength = 5
)

# Generalized Logistic regression model              
glm_model <- train(Cancer_Status ~ .,
                   data = trainData_clean, 
                   method = "glm",
                   family = "binomial", 
                   metric = "ROC",
                   trControl = ctrl)

# Random forest model
rf_model <- train(Cancer_Status ~ ., 
                  data = trainData_clean,
                  method = "rf",
                  metric = "ROC",
                  trControl = ctrl,
                  ntree = 100)

# -------------------------------------------#
# Section 14: Model Evaluation ----
# -------------------------------------------#
# ROC Plot Function
plot_roc_panel <- function(probs, labels, model_name) {
  roc_obj <- roc(labels, probs, levels = rev(levels(labels)))
  auc_val <- round(auc(roc_obj), 3)
  coords_best <- coords(roc_obj, "best", ret = c("sensitivity", "specificity"))
  sens <- round(coords_best["sensitivity"] * 100, 1)
  spec <- round(coords_best["specificity"] * 100, 1)
  
  ggroc(roc_obj, color = "red", size = 1.2) +
    geom_abline(linetype = "dotted", color = "darkgreen") +
    theme_minimal() +
    labs(title = model_name, x = "1 - Specificity", y = "Sensitivity") +
    annotate("text", x = 0.65, y = 0.2,
             label = paste0("AUC: ", auc_val,
                            "\nSensitivity: ", sens,
                            "\nSpecificity: ", spec),
             size = 4.5, hjust = 0)
}

# Generate Predictions
probs_svm <- predict(svm_model, newdata = testData, type = "prob")
probs_glm <- predict(glm_model, newdata = testData, type = "prob")
probs_rf  <- predict(rf_model,  newdata = testData, type = "prob")

# Generate ROC Panels
p1 <- plot_roc_panel(probs_svm$Lung_cancer, testData$Cancer_Status, "SVM Model")
#p1 <- plot_roc_panel(probs_svm$Lung_cancer, testData$Cancer_Status, "SVM (Radial)")
p2 <- plot_roc_panel(probs_glm$Lung_cancer, testData$Cancer_Status, "G-Logistic Regression")
p3 <- plot_roc_panel(probs_rf$Lung_cancer, testData$Cancer_Status, "Random Forest")

# Combine and display plots
plot_grid(p1, p2, p3, nrow = 1)

# -------------------------------------------#
# Section 15: Random Forest Visualization Plots ----
# -------------------------------------------#
## 1. Feature Importance Plot----
#The used code performs feature selection, which helps identify the most important miRNAs based on how strongly they are associated with the outcome (Cancer_Status) based on RF model.
var_imp <- varImp(rf_model)$importance
var_imp$miRNA <- rownames(var_imp)
var_imp <- var_imp[order(-var_imp$Overall), ]

ggplot(var_imp[1:30, ], aes(x = reorder(miRNA, Overall), y = Overall)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 30 Important miRNAs",
    x = "miRNA",
    y = "Importance"
  ) +
 theme_minimal()

# save the top 30 miRNA
top_miRNAs <- rownames(var_imp)[1:30]

## 2- PCA Plot with RF Model Predictions-----
pca_data <- prcomp(final_dataset_nor[, top_miRNAs], scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca_data$x[, 1],
  PC2 = pca_data$x[, 2],
  Status = final_dataset_nor$Cancer_Status,
  Prediction = predict(rf_model, final_dataset_nor)
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Status)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  labs(title = "PCA Plot with RF Model Predictions") +
  theme_minimal()

##3- Calibration Plot-----
calib_data <- data.frame(
  Predicted = predict(rf_model, newdata = testData, type = "prob")$Lung_cancer,
  Actual = as.numeric(testData$Cancer_Status == "Lung_cancer")
)

ggplot(calib_data, aes(x = Predicted, y = Actual)) +
  geom_point(alpha = 0.3, color = "blue") +
  geom_smooth(method = "loess", se = TRUE, color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(
    title = "Calibration Plot: Random Forest",
    x = "Predicted Probability (Lung cancer)",
    y = "Observed Frequency"
  ) +
  theme_minimal()

##4- Boxplot of Top 30 miRNAs (Random Forest)------

boxplot_data <- final_dataset_nor[, c("Sample_ID", "Cancer_Status", top_miRNAs)]
boxplot_long <- melt(boxplot_data, id.vars = c("Sample_ID", "Cancer_Status"),
                     variable.name = "miRNA", value.name = "Expression")
boxplot_long$Cancer_Status <- gsub(" ", "_", boxplot_long$Cancer_Status)

ggplot(boxplot_long, aes(x = miRNA, y = Expression, fill = Cancer_Status)) +
  geom_boxplot(outlier.size = 0.7, outlier.alpha = 0.3, position = position_dodge(width = 0.8)) +
  theme_minimal() +
  scale_fill_manual(values = c("Lung_cancer" = "#E69F00", "Non_cancer" = "lightblue")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(
    title = "Expression Levels of Top 30 miRNAs (Random Forest)",
    x = "miRNA",
    y = "Normalized Expression (log2)",
    fill = "Group"
  )
# save the top 30 miRNA
top_miRNAs <- rownames(var_imp)[1:30]

## 5 - Volcano Plot Highlighting RF Top miRNAs ----

# Step 1: Ensure proper factor levels (Non_cancer as reference)
final_dataset_nor$Cancer_Status <- relevel(final_dataset_nor$Cancer_Status, ref = "Non_cancer")

# Step 2: Save the top 30 miRNAs from Random Forest importance
top_miRNAs <- rownames(var_imp)[1:30]

# Step 3: Prepare expression matrix and design
expr_matrix <- as.matrix(final_dataset_nor[, !(colnames(final_dataset_nor) %in% c("Sample_ID", "Cancer_Status"))])
design <- model.matrix(~ Cancer_Status, data = final_dataset_nor)

# Step 4: Differential expression analysis
fit <- lmFit(t(expr_matrix), design)
fit <- eBayes(fit)
results <- topTable(fit, coef = "Cancer_StatusLung_cancer", number = Inf, adjust.method = "fdr")

# Step 5: Annotate significance
results$Significant <- "Not Significant"
results$Significant[results$adj.P.Val < 0.05 & results$logFC > 1] <- "Upregulated"
results$Significant[results$adj.P.Val < 0.05 & results$logFC < -1] <- "Downregulated"

# Step 6: Highlight Top RF features
results$Top_Feature <- ifelse(rownames(results) %in% top_miRNAs, "Top miRNA", "Other")

# Step 7: Label only significant top miRNAs
top_labeled <- results[rownames(results) %in% top_miRNAs & results$Significant != "Not Significant", ]
top_labeled <- top_labeled[order(-abs(top_labeled$logFC)), ][1:min(15, nrow(top_labeled)), ]
results$Label <- ifelse(rownames(results) %in% rownames(top_labeled), rownames(results), "")

# Step 8: Volcano plot
ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = Significant, shape = Top_Feature)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_text_repel(aes(label = Label), size = 3, box.padding = 0.4, max.overlaps = Inf) +
  scale_color_manual(values = c(
    "Upregulated" = "#E41A1C",
    "Downregulated" = "#377EB8",
    "Not Significant" = "gray70"
  )) +
  scale_shape_manual(values = c("Top miRNA" = 17, "Other" = 16)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Volcano Plot (Random Forest Top miRNAs)",
    x = "Log2 Fold Change (Lung Cancer vs. Non-Cancer)",
    y = "-Log10 Adjusted P-Value",
    color = "Regulation",
    shape = "Feature"
  ) +
  theme_minimal()

#-------------------------------------------#
# Section 16: GLR Visualization Plots ----
# -------------------------------------------#
## 1. Feature Importance Plot----
#The used code performs feature selection, which helps identify the most important miRNAs based on how strongly they are associated with the outcome (Cancer_Status) based on RF model.
var_imp_glm <- varImp(glm_model)$importance
var_imp_glm$miRNA <- rownames(var_imp)
var_imp_glm <- var_imp[order(-var_imp$Overall), ]

ggplot(var_imp_glm[1:30, ], aes(x = reorder(miRNA, Overall), y = Overall)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 30 Important miRNAs - GLR",
    x = "miRNA",
    y = "Importance"
  ) +
  theme_minimal()

# save the top 30 miRNA
top_miRNAs_glm <- rownames(var_imp)[1:30]

## 2- PCA Plot with GLR Model Predictions-----
pca_data_glm <- prcomp(final_dataset_nor[, top_miRNAs_glm], scale. = TRUE)
pca_df_glm <- data.frame(
  PC1 = pca_data$x[, 1],
  PC2 = pca_data$x[, 2],
  Status = final_dataset_nor$Cancer_Status,
  Prediction = predict(glm_model, final_dataset_nor)
)

ggplot(pca_df_glm, aes(x = PC1, y = PC2, color = Status)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  labs(title = "PCA Plot with GLR Model Predictions") +
  theme_minimal()

##3- Calibration Plot GLR-----
calib_data_glm <- data.frame(
  Predicted = predict(glm_model, newdata = testData, type = "prob")$Lung_cancer,
  Actual = as.numeric(testData$Cancer_Status == "Lung_cancer")
)

ggplot(calib_data_glm, aes(x = Predicted, y = Actual)) +
  geom_point(alpha = 0.3, color = "blue") +
  geom_smooth(method = "loess", se = TRUE, color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(
    title = "Calibration Plot: GLM",
    x = "Predicted Probability (Lung cancer)",
    y = "Observed Frequency"
  ) +
  theme_minimal()

##4- Boxplot of Top 30 miRNAs (GLR)------

boxplot_data_glm <- final_dataset_nor[, c("Sample_ID", "Cancer_Status", top_miRNAs_glm)]
boxplot_long_glm <- melt(boxplot_data_glm, id.vars = c("Sample_ID", "Cancer_Status"),
                     variable.name = "miRNA", value.name = "Expression")
boxplot_long_glm$Cancer_Status <- gsub(" ", "_", boxplot_long_glm$Cancer_Status)

ggplot(boxplot_long_glm, aes(x = miRNA, y = Expression, fill = Cancer_Status)) +
  geom_boxplot(outlier.size = 0.7, outlier.alpha = 0.3, position = position_dodge(width = 0.8)) +
  theme_minimal() +
  scale_fill_manual(values = c("Lung_cancer" = "#E69F00", "Non_cancer" = "lightblue")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(
    title = "Expression Levels of Top 30 miRNAs (GLM)",
    x = "miRNA",
    y = "Normalized Expression (log2)",
    fill = "Group"
  )

## 5 - Volcano Plot Highlighting GLR Top miRNAs ----

# Step 1: Ensure proper factor levels (Non_cancer as reference)
final_dataset_nor$Cancer_Status <- relevel(final_dataset_nor$Cancer_Status, ref = "Non_cancer")

# Step 2: Save the top 30 miRNAs from Random Forest importance
top_miRNAs_glm <- rownames(var_imp_glm)[1:30]

# Step 3: Prepare expression matrix and design
expr_matrix <- as.matrix(final_dataset_nor[, !(colnames(final_dataset_nor) %in% c("Sample_ID", "Cancer_Status"))])
design <- model.matrix(~ Cancer_Status, data = final_dataset_nor)

# Step 4: Differential expression analysis
fit <- lmFit(t(expr_matrix), design)
fit <- eBayes(fit)
results <- topTable(fit, coef = "Cancer_StatusLung_cancer", number = Inf, adjust.method = "fdr")

# Step 5: Annotate significance
results$Significant <- "Not Significant"
results$Significant[results$adj.P.Val < 0.05 & results$logFC > 1] <- "Upregulated"
results$Significant[results$adj.P.Val < 0.05 & results$logFC < -1] <- "Downregulated"

# Step 6: Highlight Top RF features
results$Top_Feature <- ifelse(rownames(results) %in% top_miRNAs_glm, "Top miRNA", "Other")

# Step 7: Label only significant top miRNAs
top_labeled <- results[rownames(results) %in% top_miRNAs_glm & results$Significant != "Not Significant", ]
top_labeled <- top_labeled[order(-abs(top_labeled$logFC)), ][1:min(15, nrow(top_labeled)), ]
results$Label <- ifelse(rownames(results) %in% rownames(top_labeled), rownames(results), "")

# Step 8: Volcano plot
ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = Significant, shape = Top_Feature)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_text_repel(aes(label = Label), size = 3, box.padding = 0.4, max.overlaps = Inf) +
  scale_color_manual(values = c(
    "Upregulated" = "#E41A1C",
    "Downregulated" = "#377EB8",
    "Not Significant" = "gray70"
  )) +
  scale_shape_manual(values = c("Top miRNA" = 17, "Other" = 16)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Volcano Plot (GLR model Top miRNAs)",
    x = "Log2 Fold Change (Lung Cancer vs. Non-Cancer)",
    y = "-Log10 Adjusted P-Value",
    color = "Regulation",
    shape = "Feature"
  ) +
  theme_minimal()


















