#==============================================================================
# PARPi7 Prediction using lncRNA Expression
#==============================================================================
# Author: [Your Name]
# Date: [Current Date]
# Description: This script analyzes the predictive power of 29 lncRNAs for 
#              PARPi7 sensitivity in ovarian cancer data. It includes data 
#              preprocessing, model training using various machine learning
#              algorithms, and model evaluation.
#
# Data source: TCGA DDR Data Resources from Cell paper
#==============================================================================

#==============================================================================
# 1. SETUP AND DEPENDENCIES
#==============================================================================

# Clear workspace
rm(list = ls())

# Load required libraries
required_packages <- c(
  "tidyverse",      # For data manipulation and visualization
  "caret",          # For machine learning workflows
  "randomForest",   # For random forest models
  "rpart",          # For decision trees
  "rpart.plot",     # For visualizing decision trees
  "multiROC",       # For ROC curve analysis
  "C50",            # For C5.0 decision tree models
  "kernlab",        # For SVM 
  "mlbench",        # For benchmark ML algorithms
  "caretEnsemble",  # For ensemble models
  "ggplot2",        # For visualization
  "doSNOW",         # For parallel processing
  "RColorBrewer"    # For color palettes
)

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, dependencies = TRUE)

# Load all packages
lapply(required_packages, library, character.only = TRUE)

# Load libraries
invisible(lapply(required_packages, library, character.only = TRUE))

# Set up parallel processing
cl <- makeCluster(min(10, parallel::detectCores()-1), type = "SOCK")
registerDoSNOW(cl)

# Set seed for reproducibility
set.seed(54321)

#==============================================================================
# 2. DATA LOADING AND PREPROCESSING
#==============================================================================

# Set working directory (consider using here::here() instead for project organization)
# setwd("~/Desktop/Machine Learning/HDR-prediction/TCGA_DDR_Data_Resources_cell paper")

# Load preprocessed data with 29 lncRNAs
# Update the file path to where your data is stored
#lncRNAs29 <- read.csv("~/Desktop/Machine Learning/HDR-prediction/Data/29lncRNAs_with_all_features.csv")

# Load necessary library
library(dplyr)

# Read the dataset
lncRNA_allData_ov <- read.csv("C:/Users/kd6/Google Drive/Work_MA/5 HDR_paper/lncRNA_allData_ov.csv")

# Define column names explicitly as a vector
lncRNA_columns <- c(
  "PARPi7_bin",  
  "ENSG00000255651.2", "ENSG00000234264.1", "ENSG00000272172.1", "ENSG00000260369.2",
  "ENSG00000233029.3", "ENSG00000261183.1", "ENSG00000261824.2", "ENSG00000267575.2",
  "ENSG00000272635.1", "ENSG00000267498.1", "ENSG00000255471.1", "ENSG00000229874.2",
  "ENSG00000248925.1", "ENSG00000257671.1", "ENSG00000231187.2", "ENSG00000270137.1",
  "ENSG00000260954.1", "ENSG00000264885.1", "ENSG00000266999.1",
  "ENSG00000261798.1", "ENSG00000250271.1", "ENSG00000254031.1", "ENSG00000260322.1",
  "ENSG00000244040.1", "ENSG00000242078.1", "ENSG00000266126.1", "ENSG00000258413.1",
  "ENSG00000242540.2", "ENSG00000231966.1"
)

# Convert to dataframe and filter
plot2 <- lncRNA_allData_ov %>%
  filter(type == "OV") %>%
  select(all_of(lncRNA_columns)) %>%
  drop_na()

# Check the output
head(plot2)

# Rename target variable to something more interpretable
plot2 <- plot2 %>% rename(CE = PARPi7_bin)
plot2$CE <- as.factor(plot2$CE)  # Convert target to factor

levels(plot2$CE) <- make.names(levels(plot2$CE))
# Check the new factor levels
print(levels(plot2$CE))


# Display summary information about the dataset
cat("Dataset dimensions:", dim(plot2), "\n")
cat("Target variable distribution:\n")
print(table(plot2$CE))

# Save the preprocessed data
write.csv(plot2, "lncRNA_processed_data.csv", row.names = FALSE)

#==============================================================================
# 3. DATA SPLITTING
#==============================================================================

# Split data into training (70%) and testing (30%) sets, maintaining class distributions
indexes <- createDataPartition(as.factor(plot2$CE),
                              times = 1,
                              p = 0.7,
                              list = FALSE)
CCNE1train <- plot2[indexes,]
CCNE1test <- plot2[-indexes,]

# Verify class proportions are maintained across splits
cat("Full dataset proportions:\n")
print(prop.table(table(plot2$CE)))
cat("Training set proportions:\n")
print(prop.table(table(CCNE1train$CE)))
cat("Test set proportions:\n")
print(prop.table(table(CCNE1test$CE)))

#==============================================================================
# 4. EXPLORATORY DATA ANALYSIS
#==============================================================================

# Create feature plots to visualize relationships between predictors and outcome
transparentTheme(trans = .9)

# Box plots - shows distribution of each predictor by outcome class
featurePlot(x = plot2[, -1],  # all columns except target 
            y = plot2$CE, 
            plot = "box", 
            scales = list(y = list(relation="free"),
                         x = list(rot = 90)))

# Density plots - shows distribution shape of each predictor by outcome class
featurePlot(x = plot2[, -1], 
            y = plot2$CE, 
            plot = "density", 
            scales = list(x = list(relation="free"), 
                         y = list(relation="free")), 
            adjust = 1.5)

#==============================================================================
# 5. MODEL TRAINING SETUP
#==============================================================================

# Training control parameters for 10-fold cross-validation
train_control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 3,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  search = "grid"
)

# For XGBoost models, we'll use a grid search for hyperparameter tuning
xgb_grid <- expand.grid(
  eta = c(0.05, 0.075, 0.1),
  nrounds = c(50, 75, 100),
  max_depth = 6:8,
  min_child_weight = c(2.0, 2.25, 2.5),
  colsample_bytree = c(0.3, 0.4, 0.5),
  gamma = 0,
  subsample = 1
)

#==============================================================================
# 6. MODEL TRAINING
#==============================================================================

# Function to train and return model
train_model <- function(data, method_name, tune_grid = NULL, train_control) {
  cat("Training", method_name, "model...\n")
  
  if (is.null(tune_grid)) {
    model <- train(
      CE ~ .,
      data = data,
      method = method_name,
      trControl = train_control,
      metric = "ROC"
    )
  } else {
    model <- train(
      CE ~ .,
      data = data,
      method = method_name,
      tuneGrid = tune_grid,
      trControl = train_control,
      metric = "ROC"
    )
  }
  
  return(model)
}

# Train multiple models
model_list <- list()

# Random Forest
model_list$rf <- train_model(CCNE1train, "rf", NULL, train_control)

# XGBoost Tree
model_list$xgbTree <- train_model(CCNE1train, "xgbTree", xgb_grid, train_control)

# XGBoost DART
model_list$xgbDART <- train_model(CCNE1train, "xgbDART", NULL, train_control)

# C5.0
model_list$C5.0 <- train_model(CCNE1train, "C5.0", NULL, train_control)

# Support Vector Machine with Radial Kernel
model_list$svmRadial <- train_model(CCNE1train, "svmRadial", NULL, train_control)

# K-Nearest Neighbors
model_list$kknn <- train_model(CCNE1train, "kknn", NULL, train_control)

#==============================================================================
# 7. MODEL COMPARISON
#==============================================================================

# Compare models using resamples
results <- resamples(model_list)

# Summarize and visualize model comparison
summary(results)

# Boxplots of model performance
bwplot(results, main="Model Performance Comparison")

# Dot plots of model performance
dotplot(results, main="Model Performance Comparison")

# Save model comparison results
write.csv(results$values, "model_comparison_results.csv")

#==============================================================================
# 8. BEST MODEL EVALUATION
#==============================================================================

# Identify the best model based on ROC
model_metrics <- summary(results)$statistics$ROC
best_model_name <- names(model_list)[which.max(model_metrics[,"Mean"])]
best_model <- model_list[[best_model_name]]

cat("Best model:", best_model_name, "\n")
print(best_model)

# Make predictions on test set
predictions <- predict(best_model, CCNE1test)
pred_probs <- predict(best_model, CCNE1test, type = "prob")

# Calculate confusion matrix
conf_matrix <- confusionMatrix(predictions, as.factor(CCNE1test$CE))
print(conf_matrix)

# Save confusion matrix
write.csv(conf_matrix$table, "confusion_matrix.csv")
write.csv(conf_matrix$byClass, "classification_metrics.csv")

#==============================================================================
# 9. FEATURE IMPORTANCE
#==============================================================================

# Extract variable importance from the best model
var_imp <- varImp(best_model, scale = TRUE)

# Plot variable importance
plot(var_imp, main = "Variable Importance")

# Save variable importance
write.csv(var_imp$importance, "variable_importance.csv")

#==============================================================================
# 10. ROC CURVE ANALYSIS
#==============================================================================
# Check what's in pred_probs
dim(pred_probs)
head(pred_probs)

# Create a clearer dataframe
true_class <- as.factor(CCNE1test$CE)
true_labels <- model.matrix(~ 0 + true_class) # Creates one-hot encoded matrix
colnames(true_labels) <- paste0(levels(true_class), "_true")

# Make sure pred_probs has matching column names
# If pred_probs columns don't have proper names, fix them
if(is.null(colnames(pred_probs))) {
  colnames(pred_probs) <- paste0(levels(true_class), "_pred")
} else {
  # Ensure they have _pred suffix
  colnames(pred_probs) <- paste0(colnames(pred_probs), "_pred")
}

# Combine them properly
final_df <- cbind(data.frame(true_labels), data.frame(pred_probs))

# Now try the multi_roc function
# Check the structure first
str(final_df)
head(final_df)


# For binary classification
if(length(levels(as.factor(CCNE1test$CE))) == 2) {
  roc_obj <- roc(CCNE1test$CE, pred_probs[,1])
  plot(roc_obj)
  auc(roc_obj)
}

#roc_res <- multi_roc(final_df, force_diag=TRUE)

# For multi_roc, you'll need the appropriate package
# If using pROC package:
library(pROC)

# Generate ROC curves - approach depends on which package you're using
# If using pROC for binary classification:
# roc_obj <- roc(CCNE1test$CE, pred_probs[,1])

# For multiclass ROC (check if you have the multiROC package)
# install.packages("multiROC")  # if needed
library(multiROC)  # or equivalent package
roc_res <- multi_roc(final_df, force_diag=TRUE)
plot_roc_data <- plot_roc_data(roc_res)

# Save ROC data
write.csv(plot_roc_data, "roc_curve_data.csv")

# Plot ROC curves
library(ggplot2)
ggplot(plot_roc_data, aes(x = 1-Specificity, y=Sensitivity)) +
  geom_path(aes(color = Group), size=1.5) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               colour='grey', linetype = 'dotdash') +
  theme_bw() + 
  labs(title = "ROC Curves for lncRNA-based PARPi7 Prediction",
       x = "False Positive Rate (1-Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.justification=c(1, 0), legend.position=c(.95, .05),
        legend.title=element_blank(), 
        legend.background = element_rect(fill=NULL, size=0.5, 
                                         linetype="solid", colour ="black"))

#==============================================================================
# 11. ENSEMBLE MODEL (OPTIONAL)
#==============================================================================
library(caret)
library(caretEnsemble)

# First, create the proper trainControl for all individual models
model_control <- trainControl(
  method = "cv",
  number = 10,
  savePredictions = "all",  # This is critical for later ensemble creation
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

# Recreate all your individual models with this trainControl
# For example (adjust these to match your actual models):
model1 <- train(CE ~ ., data = CCNE1train, 
                method = "rf", 
                trControl = model_control,
                metric = "ROC")

model2 <- train(CE ~ ., data = CCNE1train, 
                method = "glmnet", 
                trControl = model_control,
                metric = "ROC")

model3 <- train(CE ~ ., data = CCNE1train, 
                method = "xgbTree", 
                trControl = model_control,
                metric = "ROC")

model_list <- list(model1, model2, model3)

# Set up control for ensemble
ensemble_control <- trainControl(
  method = "cv",
  number = 10,
  savePredictions = "final",  
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)


# Train ensemble model
ensemble_model <- caretEnsemble(
  model_list,  # Directly pass list of models
  metric = "ROC",
  trControl = ensemble_control
)

# Summarize ensemble
summary(ensemble_model)

)

# Train the ensemble model
library(caretEnsemble)
ensemble_model <- caretEnsemble(
  model_list,
  metric = "ROC",
  trControl = ensemble_control
)



# Get the column name with the maximum probability for each row (class label)
ensemble_preds_class <- ifelse(ensemble_preds$X0 > ensemble_preds$X1, "X0", "X1")

# Convert to factor
ensemble_preds_class <- factor(ensemble_preds_class, levels = c("X0", "X1"))

# Ensure both 'ensemble_preds_class' and 'CCNE1test$CE' have the same levels
levels_true <- levels(CCNE1test$CE)

# Relevel ensemble predictions to match the true label levels
ensemble_preds_class <- factor(ensemble_preds_class, levels = levels_true)

# Now compute the confusion matrix
ensemble_conf_matrix <- confusionMatrix(ensemble_preds_class, CCNE1test$CE)

# Print confusion matrix
print(ensemble_conf_matrix)








#==============================================================================
# 12. CLEANUP
#==============================================================================

# Stop parallel cluster
stopCluster(cl)

# Save important objects
save(model_list, best_model, ensemble_model, file = "lncRNA_models.RData")

cat("Analysis complete!\n")
