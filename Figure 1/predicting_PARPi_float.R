# Load necessary libraries
library(tidyverse)
library(tidymodels)
library(ranger)
library(kernlab)
library(vip)
library(workflows)
library(tune)
tidymodels_prefer()

# Load data (adapt path as needed)
DR <- read.csv("C:/Users/kd6/Google Drive/Work_MA/5 HDR_paper/lncRNA_allData_ov.csv")

# Filter breast cancer samples and select relevant columns
hrd_data <- DR %>% 
  filter(type == "OV") %>% 
  select( PARPi7, 
          ENSG00000255651.2,
          ENSG00000234264.1,
          ENSG00000272172.1,
          ENSG00000260369.2,
          ENSG00000233029.3,
          ENSG00000261183.1,
          ENSG00000261824.2,
          ENSG00000267575.2,
          ENSG00000272635.1,
          ENSG00000267498.1,
          ENSG00000255471.1,
          ENSG00000229874.2,
          ENSG00000248925.1,
          ENSG00000257671.1,
          ENSG00000231187.2,
          ENSG00000270137.1,
          ENSG00000260954.1,
          ENSG00000264885.1,
          ENSG00000266999.1,
          ENSG00000261798.1,
          ENSG00000250271.1,
          ENSG00000254031.1,
          ENSG00000260322.1,
          ENSG00000244040.1,
          ENSG00000242078.1,
          ENSG00000266126.1,
          ENSG00000258413.1,
          ENSG00000242540.2,
          ENSG00000231966.1) %>% drop_na()

# Check dimensions and data structure
dim(hrd_data)
head(hrd_data)
summary(hrd_data$PARPi7)

# Implement a three-way split: train (60%), validation (20%), test (20%)
set.seed(123)
# Initial split to get the test set
initial_split_obj <- initial_split(hrd_data, prop = 0.8, strata = PARPi7)
training_validation <- training(initial_split_obj)
testing <- testing(initial_split_obj)

# Split the training data into training and validation
set.seed(234)
train_val_split <- initial_split(training_validation, prop = 0.75, strata = PARPi7)
training <- training(train_val_split)
validation <- testing(train_val_split)

# Check dimensions of all sets
dim(training)
dim(validation)
dim(testing)

# Create a proper preprocessing recipe
hrd_recipe <- recipe(PARPi7 ~ ., data = training) %>%
  step_zv(all_predictors()) %>%
  step_normalize(all_predictors()) %>%
  step_impute_knn(all_predictors())


# Create cross-validation folds for proper tuning
set.seed(345)
cv_folds <- vfold_cv(training, v = 10, strata = PARPi7)

# Define model specifications with tuning parameters
# 1. Linear Regression with regularization
lm_spec <- linear_reg(penalty = tune(), mixture = tune()) %>%
  set_engine("glmnet") %>%
  set_mode("regression")

# 2. Random Forest
rf_spec <- rand_forest(
  mtry = tune(),
  trees = tune(),
  min_n = tune()
) %>%
  set_engine("ranger", importance = "permutation") %>%
  set_mode("regression")

# 3. SVM
svm_spec <- svm_rbf(
  cost = tune(),
  rbf_sigma = tune()
) %>%
  set_engine("kernlab") %>%
  set_mode("regression")

# Create workflows
lm_workflow <- workflow() %>%
  add_recipe(hrd_recipe) %>%
  add_model(lm_spec)

rf_workflow <- workflow() %>%
  add_recipe(hrd_recipe) %>%
  add_model(rf_spec)

svm_workflow <- workflow() %>%
  add_recipe(hrd_recipe) %>%
  add_model(svm_spec)

# Define tuning grids
lm_grid <- grid_regular(
  penalty(range = c(-3, 0), trans = log10_trans()),
  mixture(range = c(0, 1)),
  levels = c(10, 5)
)

rf_grid <- grid_regular(
  mtry(range = c(5, 25)),
  trees(range = c(500, 1500)),
  min_n(range = c(2, 10)),
  levels = c(5, 3, 3)
)

svm_grid <- grid_regular(
  cost(range = c(-1, 3), trans = log10_trans()),
  rbf_sigma(range = c(-3, 0), trans = log10_trans()),
  levels = c(10, 10)
)

# Set up parallel processing if available
# library(doParallel)
# all_cores <- parallel::detectCores() - 1
# cl <- makeCluster(all_cores)
# registerDoParallel(cl)
library(doParallel)
registerDoParallel(cores = 4) 
# Tuning control
ctrl <- control_grid(
  save_pred = TRUE,
  save_workflow = TRUE,
  verbose = TRUE,
  parallel_over = "resamples"  # Specify what to parallelize over
)

# Tune models
set.seed(456)
lm_tuned <- tune_grid(
  lm_workflow,
  resamples = cv_folds,
  grid = lm_grid,
  metrics = metric_set(rmse, rsq, mae),
  control = ctrl
)

set.seed(567)
rf_tuned <- tune_grid(
  rf_workflow,
  resamples = cv_folds,
  grid = rf_grid,
  metrics = metric_set(rmse, rsq, mae),
  control = ctrl
)

set.seed(678)
svm_tuned <- tune_grid(
  svm_workflow,
  resamples = cv_folds,
  grid = svm_grid,
  metrics = metric_set(rmse, rsq, mae),
  control = ctrl
)

# Compare tuning results
autoplot(lm_tuned)
autoplot(rf_tuned)
autoplot(svm_tuned)

# Select best models
best_lm <- select_best(lm_tuned, metric = "rmse")
best_rf <- select_best(rf_tuned, metric = "rmse")
best_svm <- select_best(svm_tuned, metric = "rmse")

# Finalize workflows with best parameters
final_lm <- finalize_workflow(lm_workflow, best_lm)
final_rf <- finalize_workflow(rf_workflow, best_rf)
final_svm <- finalize_workflow(svm_workflow, best_svm)

# Fit final models
lm_fit <- fit(final_lm, training)
rf_fit <- fit(final_rf, training)
svm_fit <- fit(final_svm, training)

# Evaluate on validation set
evaluate_model <- function(model_fit, new_data, model_name) {
  predictions <- predict(model_fit, new_data) %>%
    bind_cols(new_data %>% select(PARPi7))
  
  metrics <- metric_set(rmse, rsq, mae)
  results <- metrics(predictions, truth = PARPi7, estimate = .pred)
  
  # Create visualization
  p <- ggplot(predictions, aes(x = PARPi7, y = .pred)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(
      title = paste("Model:", model_name),
      x = "Actual PARPi7",
      y = "Predicted PARPi7",
      subtitle = paste(
        "RMSE =", round(results %>% filter(.metric == "rmse") %>% pull(.estimate), 3),
        "| R² =", round(results %>% filter(.metric == "rsq") %>% pull(.estimate), 3)
      )
    ) +
    theme_minimal()
  
  list(metrics = results, plot = p, predictions = predictions)
}

# Validation set evaluation
lm_val <- evaluate_model(lm_fit, validation, "Linear Regression")
rf_val <- evaluate_model(rf_fit, validation, "Random Forest")
svm_val <- evaluate_model(svm_fit, validation, "SVM")

# Compare validation results
bind_rows(
  lm_val$metrics %>% mutate(model = "Linear Regression"),
  rf_val$metrics %>% mutate(model = "Random Forest"),
  svm_val$metrics %>% mutate(model = "SVM")
) %>%
  pivot_wider(names_from = .metric, values_from = .estimate) %>%
  arrange(rmse)

# Extract best model based on validation performance
# (Assuming RF is best based on your original code)
best_model <- rf_fit

# Evaluate on test set (final evaluation)
final_results <- evaluate_model(best_model, testing, "Best Model (Random Forest)")
final_results$metrics

# Visualize final results
final_results$plot

# Feature importance for random forest
if (inherits(best_model, "workflow")) {
  rf_importance <- extract_fit_parsnip(best_model) %>%
    vip(num_features = 20)
  print(rf_importance)
}

# Save model and results
saveRDS(best_model, "best_hrd_prediction_model.rds")
write_csv(final_results$predictions, "hrd_test_predictions.csv")

# Model explanation (partial dependence plots for important features)
if (inherits(best_model, "workflow")) {
  library(DALEXtra)
  explainer <- explain_tidymodels(
    best_model,
    data = testing %>% select(-PARPi7),
    y = testing$PARPi7,
    label = "Best Model"
  )
  
  # Get top 5 important features
  top_features <- rf_importance$data %>%
    arrange(desc(Importance)) %>%
    head(5) %>%
    pull(Variable)
  
  # Create partial dependence plots
  pdp_plots <- lapply(top_features, function(feature) {
    pdp <- model_profile(explainer, variables = feature)
    plot(pdp, geom = "profiles") +
      ggtitle(paste("Partial Dependence Plot for", feature))
  })
  
  print(pdp_plots)
}

# Combine the metrics from all models into a single data frame
comparison_results <- bind_rows(
  lm_val$metrics %>% mutate(model = "Linear Regression"),
  rf_val$metrics %>% mutate(model = "Random Forest"),
  svm_val$metrics %>% mutate(model = "SVM")
) %>%
  pivot_wider(names_from = .metric, values_from = .estimate) %>%
  arrange(rmse)

# Plot comparison for RMSE and R²
library(ggplot2)

# Combine RMSE and R² into one long format for comparison
comparison_long <- comparison_results %>%
  select(model, rmse, rsq) %>%
  pivot_longer(cols = c(rmse, rsq), names_to = "metric", values_to = "value")

# Create a bar plot comparing RMSE and R² for each model
ggplot(comparison_long, aes(x = model, y = value, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Model Comparison: RMSE and R²",
    y = "Metric Value",
    x = "Model",
    fill = "Metric"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("rmse" = "lightblue", "rsq" = "lightgreen"))


# Extract cross-validation metrics for each model
lm_cv_results <- lm_tuned %>%
  collect_metrics() %>%
  filter(.metric %in% c("rmse", "rsq")) %>%
  mutate(model = "Linear Regression")

rf_cv_results <- rf_tuned %>%
  collect_metrics() %>%
  filter(.metric %in% c("rmse", "rsq")) %>%
  mutate(model = "Random Forest")

svm_cv_results <- svm_tuned %>%
  collect_metrics() %>%
  filter(.metric %in% c("rmse", "rsq")) %>%
  mutate(model = "SVM")

# Bind results into one dataframe
cv_comparison <- bind_rows(lm_cv_results, rf_cv_results, svm_cv_results)

# Ensure that the cross-validation fold information is present in the data
cv_comparison <- cv_comparison %>%
  mutate(fold = rep(1:10, length.out = nrow(cv_comparison)))  # Add fold information if missing

# Adjust the cv_comparison data for plotting
cv_comparison_plot_data <- cv_comparison %>%
  filter(.metric %in% c("rmse", "rsq")) %>%
  select(fold, .metric, mean, model)  # Extract necessary columns

# Check the structure again
str(cv_comparison_plot_data)

# Plot the cross-validation results with the iteration (fold) on the x-axis
ggplot(cv_comparison_plot_data, aes(x = fold, y = mean, color = model, group = interaction(model, fold))) +
  geom_line() +
  geom_point() +
  facet_wrap(~ .metric, scales = "free_y") +  # Separate RMSE and R² in different facets
  labs(
    title = "10x Cross-Validation Results for Models",
    x = "Fold",
    y = "Metric Value",
    color = "Model"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Violin plot with individual values (overlay points on the plot)
ggplot(cv_comparison_plot_data, aes(x = model, y = mean, fill = model)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  # Create the violin plot
  geom_jitter(aes(color = model), width = 0.1, alpha = 0.6) +  # Overlay individual data points
  facet_wrap(~ .metric, scales = "free_y") +  # Separate RMSE and R² in different facets
  labs(
    title = "Distribution of Cross-Validation Results by Model",
    x = "Model",
    y = "Metric Value",
    fill = "Model",
    color = "Model"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Now, compare validation set performance for each model
comparison_results <- bind_rows(
  lm_val$metrics %>% mutate(model = "Linear Regression"),
  rf_val$metrics %>% mutate(model = "Random Forest"),
  svm_val$metrics %>% mutate(model = "SVM")
) %>%
  pivot_wider(names_from = .metric, values_from = .estimate) %>%
  arrange(rmse)

# Combine Cross-validation and Validation Data
comparison_long <- bind_rows(
  cv_comparison %>% mutate(data_type = "Cross-validation"),
  comparison_results %>% select(model, rmse, rsq) %>% pivot_longer(cols = c(rmse, rsq), names_to = "metric", values_to = "value") %>% mutate(data_type = "Validation")
)

# Plot combined comparison of cross-validation and validation results
ggplot(comparison_long, aes(x = model, y = value, fill = interaction(metric, data_type))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Model Comparison: Cross-validation vs. Validation",
    y = "Metric Value",
    x = "Model",
    fill = "Metric & Data Type"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("rmse.Cross-validation" = "lightblue", "rsq.Cross-validation" = "lightgreen", "rmse.Validation" = "orange", "rsq.Validation" = "yellow"))


# Calculate Pearson correlation for each model
lm_pearson <- cor(lm_val$predictions$PARPi7, lm_val$predictions$.pred)
rf_pearson <- cor(rf_val$predictions$PARPi7, rf_val$predictions$.pred)
svm_pearson <- cor(svm_val$predictions$PARPi7, svm_val$predictions$.pred)

# Create a data frame to bind predictions and model labels
all_predictions <- bind_rows(
  lm_val$predictions %>% mutate(model = "Linear Regression"),
  rf_val$predictions %>% mutate(model = "Random Forest"),
  svm_val$predictions %>% mutate(model = "SVM")
)

# Plot actual vs predicted values with smoothing lines
ggplot(all_predictions, aes(x = PARPi7, y = .pred, color = model)) +
  geom_point(alpha = 0.5) +  # Scatter plot for predictions
  geom_smooth(method = "lm", se = FALSE, aes(group = model), linetype = "solid") +  # Smoothing line for each model
  facet_wrap(~ model) +  # Separate by model
  labs(
    title = "Comparison of Models: Actual vs Predicted (Validation Data)",
    subtitle = paste(
      "Linear Regression Pearson: ", round(lm_pearson, 3), "\n",
      "Random Forest Pearson: ", round(rf_pearson, 3), "\n",
      "SVM Pearson: ", round(svm_pearson, 3)
    ),
    x = "Actual PARPi7",
    y = "Predicted PARPi7",
    color = "Model"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

