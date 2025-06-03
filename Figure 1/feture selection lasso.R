# Load necessary libraries
library(glmnet)
library(dplyr)
library(psych)
library(tidyverse)

# Set seed for reproducibility
set.seed(123)

# Set working directory and load data

PANCAN_HDR <- read_csv("DDR_lncRNA_Cell_paper_wo_highCor")

# Prepare data
y <- PANCAN_HDR %>%
  drop_na() %>%
  select(CCNE1_SNP6_array) %>%
  scale(center = TRUE, scale = FALSE) %>%
  as.matrix()

X <- PANCAN_HDR %>%
  drop_na() %>%
  select(46:10087) %>%
  as.matrix()

# Check dimensions
dim(y)
dim(X)

# Perform 10-fold cross-validation to select lambda for Ridge Regression
lambdas_to_try <- 10^seq(-3, 5, length.out = 100)
ridge_cv <- cv.glmnet(X, y, alpha = 0, lambda = lambdas_to_try,
                      standardize = TRUE, nfolds = 10, family = "gaussian")

# Plot cross-validation results
plot(ridge_cv)

# Best cross-validated lambda
lambda_cv <- ridge_cv$lambda.min

# Fit final Ridge Regression model
model_cv <- glmnet(X, y, alpha = 0, lambda = lambda_cv, standardize = TRUE, family = "gaussian")

# Extract coefficients
ridge_coefs <- as.matrix(coef(model_cv))
write.csv(ridge_coefs, "RR_CCNE1_SNP6_array")

# Read and sort coefficients
RF <- read_csv("RR_CCNE1_SNP6_array") %>%
  arrange(desc(s0))

# Save sorted coefficients
write.csv(RF, "RR_CCNE1_SNP6_array")

# Display top coefficients
head(RF)

# Lasso Regression
model_lasso <- glmnet(X, y, alpha = 1, lambda = lambda_cv, standardize = TRUE)

# Extract coefficients
lasso_coefs <- as.matrix(coef(model_lasso))
write.csv(lasso_coefs, "LASSO_HRD_Score")

# Read and sort coefficients
LS <- read_csv("LASSO_HRD_Score") %>%
  arrange(s0)

# Save sorted coefficients
write.csv(LS, "LASSO_HRD_Score")

# Display top coefficients
head(LS)

# Predict and calculate residuals for Ridge Regression
y_hat_cv <- predict(model_cv, X)
ssr_cv <- t(y - y_hat_cv) %*% (y - y_hat_cv)
rsq_ridge_cv <- cor(y, y_hat_cv)^2

# Use information criteria to select lambda
X_scaled <- scale(X)
aic <- bic <- c()

for (lambda in seq_along(lambdas_to_try)) {
  model <- glmnet(X, y, alpha = 0, lambda = lambdas_to_try[lambda], standardize = TRUE)
  betas <- as.vector(as.matrix(coef(model))[-1, ])
  resid <- y - (X_scaled %*% betas)
  ld <- lambdas_to_try[lambda] * diag(ncol(X_scaled))
  H <- X_scaled %*% solve(t(X_scaled) %*% X_scaled + ld) %*% t(X_scaled)
  df <- tr(H)
  aic[lambda] <- nrow(X_scaled) * log(t(resid) %*% resid) + 2 * df
  bic[lambda] <- nrow(X_scaled) * log(t(resid) %*% resid) + 2 * df * log(nrow(X_scaled))
}

# Plot information criteria against lambda values
plot(log(lambdas_to_try), aic, col = "orange", type = "l",
     ylim = c(190, 260), ylab = "Information Criterion")
lines(log(lambdas_to_try), bic, col = "skyblue3")
legend("bottomright", lwd = 1, col = c("orange", "skyblue3"), legend = c("AIC", "BIC"))

# Lasso with alpha = 1
cv1 <- cv.glmnet(X, y, family = "gaussian", nfold = 10, type.measure = "deviance", parallel = TRUE, alpha = 1)
md1 <- glmnet(X, y, family = "gaussian", lambda = lambda_cv, alpha = 1)
coef(md1)

# Ridge with alpha = 0
cv2 <- cv.glmnet(X, y, family = "gaussian", nfold = 10, type.measure = "deviance", parallel = TRUE, alpha = 0)
md2 <- glmnet(X, y, family = "gaussian", lambda = lambda_cv, alpha = 0)
coef(md2)

# Elastic Net with 0 < alpha < 1
a <- seq(0.1, 0.9, 0.05)
search <- foreach(i = a, .combine = rbind) %dopar% {
  cv <- cv.glmnet(X, y, family = "gaussian", nfold = 10, type.measure = "deviance", parallel = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
}

cv3 <- search[search$cvm == min(search$cvm), ]
md3 <- glmnet(X, y, family = "gaussian", lambda = cv3$lambda.1se, alpha = cv3$alpha)
coef(md3)
