# Load necessary libraries
library(tidyverse)
library(caret)
library(randomForest)
library(reprtree)

# Set working directory and load data

PANCAN_HDR <- read_csv("DDR_lncRNA_Cell_paper")

# Check dimensions of the dataset
dim(PANCAN_HDR)

# Display the first few rows of selected columns
head(PANCAN_HDR[, 1:49])

# Filter data for type "OV" and select relevant columns
OV3 <- PANCAN_HDR %>%
  filter(type == "OV") %>%
  select(CCNE1_SNP6_array, 46:10087) %>%
  drop_na()

# Check dimensions after filtering
dim(OV3)

# Fit a random forest model
fit_rf <- randomForest(CCNE1_SNP6_array ~ ., data = OV3, ntree = 5)

# Calculate variable importance using mean decrease in Gini
RF_importance <- importance(fit_rf, type = 2)

# Save and read the importance scores
write.csv(RF_importance, "RF_CCNE1_SNP6_array_OV")
RF <- read_csv("RF_CCNE1_SNP6_array_OV")

# Sort importance scores by MeanDecreaseGini
RF <- RF %>% arrange(desc(MeanDecreaseGini))

# Save the sorted importance scores
write.csv(RF, "RF_CCNE1_SNP6_array_OV")

# Display the top 20 important features
head(RF, 20)

# Plot the first tree from the random forest model
plot(getTree(fit_rf, 1, labelVar = TRUE))

# Plot variable importance
varImpPlot(fit_rf)

# Display summary and structure of the importance dataframe
summary(RF)
str(RF)
