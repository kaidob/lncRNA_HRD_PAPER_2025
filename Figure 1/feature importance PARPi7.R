# Clear workspace
rm(list = ls())

# Load necessary libraries
library(tidyverse)
library(tidymodels)
library(ranger)
library(vip)
library(finetune)
library(tune)
library(shapviz)

tidymodels_prefer()

# Load data (adapt path as needed)
DR <- read.csv("lncRNA_allData_ov.csv")

# Filter ovarian cancer samples and select relevant columns
hrd_data <- DR %>% 
  filter(type == "OV") %>% 
  select( PARPi7,  # Change target to HRD_Score
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

# Split data into training (80%) and testing (20%)
set.seed(123)
data_split <- initial_split(hrd_data, prop = 0.8)
train_data <- training(data_split)
test_data  <- testing(data_split)

# Define Random Forest model
rf_model <- rand_forest(mode = "regression", trees = 500) %>%
  set_engine("ranger", importance = "permutation")

# Create a workflow
rf_workflow <- workflow() %>%
  add_model(rf_model) %>%
  add_formula(PARPi7 ~ .)

# Fit the model
rf_fit <- fit(rf_workflow, data = train_data)

# Extract feature importance
rf_importance <- vip::vi(rf_fit)

# Print feature importance
print(rf_importance)

# Plot feature importance
rf_importance %>%
  ggplot(aes(x = Importance, y = fct_reorder(Variable, Importance))) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  labs(title = "Feature Importance (Random Forest)", x = "Importance", y = "Features")

# ---- SHAP Analysis ----
# Fit a SHAP model
# Install and load necessary packages
#install.packages("fastshap")
library(fastshap)

# Extract the fitted model from the workflow
rf_extracted <- extract_fit_parsnip(rf_fit)$fit

# Define a correct prediction wrapper for fastshap
pred_wrapper <- function(object, newdata) {
  predict(object, data = newdata)$predictions  # Change "new_data" to "data"
}

# Compute SHAP values using fastshap
set.seed(123)
shap_values <- fastshap::explain(
  rf_extracted, 
  X = as.data.frame(test_data[-1]), 
  pred_wrapper = pred_wrapper,  # Use corrected wrapper
  nsim = 50
)

# Convert SHAP values to a tidy format
shap_df <- as_tibble(shap_values) %>% 
  mutate(Sample = row_number()) %>% 
  pivot_longer(-Sample, names_to = "Feature", values_to = "SHAP_Value")

# Plot SHAP summary
ggplot(shap_df, aes(x = SHAP_Value, y = fct_reorder(Feature, SHAP_Value, .fun = median))) +
  geom_boxplot(fill = "steelblue") +
  theme_minimal() +
  labs(title = "SHAP Feature Importance", x = "SHAP Value", y = "Feature")

# ---- Statistical Validation ----










# Extract the fitted ranger model from the tidymodels workflow
rf_extracted <- extract_fit_parsnip(rf_fit)$fit

# Load vip package
library(vip)


pred_wrapper <- function(object, newdata) {
  predict(object, data = newdata)$predictions  # Correct format for ranger
}

# Compute permutation-based feature importance
set.seed(123)
perm_test <- vip::vi_permute(
  rf_extracted, 
  train = as.data.frame(train_data),  # Provide full training data
  target = "PARPi7",  # Specify the target column
  metric = "rmse", 
  nsim = 100, 
  pred_wrapper = pred_wrapper  # Pass the defined prediction wrapper
)


# Calculate mean and standard deviation of the importance values
perm_test <- perm_test %>%
  mutate(
    Importance_StDev = apply(perm_test[, -1], 1, sd)  # Compute StDev for each feature
  )

# Print results
print(perm_test)

# Plot permutation importance
ggplot(perm_test, aes(x = Importance, y = fct_reorder(Variable, Importance))) +
  geom_col(fill = "red") +
  theme_minimal() +
  labs(title = "Permutation Test Feature Importance", x = "Importance", y = "Feature")

ggplot(perm_test, aes(x = Importance, y = fct_reorder(Variable, Importance))) +
  geom_col(aes(fill = Importance), width = 0.7) +
  geom_errorbarh(aes(xmin = Importance - Importance_StDev, xmax = Importance + Importance_StDev), height = 0.3) +
  theme_minimal() +
  labs(
    title = "Permutation Test Feature Importance with StDev", 
    x = "Importance (with StDev)", 
    y = "Feature"
  ) +
  theme(legend.position = "none")



# Permutation test to check feature importance significance
set.seed(123)
perm_test <- vip::vi_permute(rf_fit, train_data[-1], train_data$HRD_Score, metric = "rmse", nsim = 100)

# Print permutation test results
print(perm_test)

# Plot permutation importance
perm_test %>%
  ggplot(aes(x = Importance, y = fct_reorder(Variable, Importance))) +
  geom_col(fill = "red") +
  theme_minimal() +
  labs(title = "Permutation Test Feature Importance", x = "Importance", y = "Features")

