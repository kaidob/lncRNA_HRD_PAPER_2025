# Benötigte Pakete laden
library(survival)
library(survminer)
library(dplyr)

# Angenommen, df2 enthält die Daten und ENSG00000272172.1 ist eine numerische Variable
# Wir unterteilen die Expression des Gens in 3 Gruppen (hoch, mittel, niedrig)
df2$expression_group <- cut(df2$ENSG00000272172.1,
                            breaks = quantile(df2$ENSG00000272172.1, probs = 0:3 / 3, na.rm = TRUE),
                            labels = c("Low", "Medium", "High"),
                            include.lowest = TRUE)

# Überlebensobjekt erstellen
surv_obj <- Surv(time = df2$OS_time, event = df2$OS)

# Kaplan-Meier Analyse durchführen
km_fit <- survfit(surv_obj ~ expression_group, data = df2)

# Kaplan-Meier Kurve plotten
ggsurvplot(km_fit, data = df2,
           pval = TRUE,               # p-Wert anzeigen
           risk.table = TRUE,         # Risikotabelle anzeigen
           pval.method = TRUE,        # p-Wert-Methode anzeigen
           legend.title = "Gene Expression",
           legend.labs = c("Low", "Medium", "High"),  # Beschriftungen für die Gruppen
           xlab = "Time (Days)",       # X-Achse beschriften
           ylab = "Survival Probability",  # Y-Achse beschriften
           title = "Kaplan-Meier Curve for ENSG00000272172.1 Expression",
           palette = c("red", "blue", "green"))  # Farben für die Gruppen


# Load necessary libraries
library(survival)
library(survminer)

# Assuming df2 has the data, and ENSG00000272172.1 is the gene expression variable

# Create the expression groups (Low, Medium, High) based on quantiles
df2$expression_group <- cut(df2$ENSG00000272172.1,
                            breaks = quantile(df2$ENSG00000272172.1, probs = 0:3 / 3, na.rm = TRUE),
                            labels = c("Low", "Medium", "High"),
                            include.lowest = TRUE)

# Create a survival object
surv_obj <- Surv(time = df2$OS_time, event = df2$OS)

# Fit the Cox proportional hazards model
cox_model <- coxph(surv_obj ~ expression_group, data = df2)

# Summarize the model to get the Cox regression statistics
summary(cox_model)

# Check the proportional hazards assumption
cox_ph_assum <- cox.zph(cox_model)

# Print the proportional hazards test results
cox_ph_assum

#Code for Cox Proportional Hazards Model with Confounding Effects:
# Creating the survival object
surv_obj <- Surv(time = df2$OS_time, event = df2$OS)

# Fitting the Cox model with confounders (age, purity)
cox_model <- coxph(surv_obj ~ expression_group + age + purity, data = df2)

# Summarizing the model
summary(cox_model)


# Load necessary libraries
library(survival)
library(survminer)
library(ggplot2)

# Ensure numeric variables
df2$OS_time <- as.numeric(df2$OS_time)
df2$OS <- as.numeric(df2$OS)
df2$age <- as.numeric(df2$age)
df2$purity <- as.numeric(df2$purity)

# Convert 'expression_group ' to factor (categorical variable)
df2$expression_group  <- as.factor(df2$expression_group )

# Fit a Cox proportional hazards model adjusting for age and purity
cox_model <- coxph(Surv(OS_time, OS) ~ expression_group  + age + purity, data = df2)

# Perform pairwise comparisons for all groups in "expression_group "
pairwise_results <- pairwise_survdiff(Surv(OS_time, OS) ~ expression_group , data = df2, p.adjust.method = "BH")

# Print pairwise comparison results
print(pairwise_results)

# Print model summary
summary(cox_model)

# Plot adjusted survival curves
ggadjustedcurves(cox_model, data = df2, variable = "expression_group", method = "conditional") +
  ggtitle("Adjusted Survival Curves by ENSG00000272172.1 Expression") +
  theme_minimal()

cox.zph(cox_model)


############################################################
##############           PFI #############
###############################################


# Überlebensobjekt erstellen
surv_obj <- Surv(time = df2$PFI_time, event = df2$PFI)

# Kaplan-Meier Analyse durchführen
km_fit <- survfit(surv_obj ~ expression_group, data = df2)

# Kaplan-Meier Kurve plotten
ggsurvplot(km_fit, data = df2,
           pval = TRUE,               # p-Wert anzeigen
           risk.table = TRUE,         # Risikotabelle anzeigen
           pval.method = TRUE,        # p-Wert-Methode anzeigen
           legend.title = "Gene Expression",
           legend.labs = c("Low", "Medium", "High"),  # Beschriftungen für die Gruppen
           xlab = "Time (Days)",       # X-Achse beschriften
           ylab = "Survival Probability",  # Y-Achse beschriften
           title = "Kaplan-Meier Curve for ENSG00000272172.1 Expression",
           palette = c("red", "blue", "green"))  # Farben für die Gruppen


# Create a survival object
surv_obj <- Surv(time = df2$PFI_time, event = df2$PFI)

# Fit the Cox proportional hazards model
cox_model <- coxph(surv_obj ~ expression_group, data = df2)

# Summarize the model to get the Cox regression statistics
summary(cox_model)

# Check the proportional hazards assumption
cox_ph_assum <- cox.zph(cox_model)

# Print the proportional hazards test results
cox_ph_assum

#Code for Cox Proportional Hazards Model with Confounding Effects:
# Creating the survival object
surv_obj <- Surv(time = df2$PFI_time, event = df2$PFI)

# Fitting the Cox model with confounders (age, purity)
cox_model <- coxph(surv_obj ~ expression_group + age + purity, data = df2)

# Summarizing the model
summary(cox_model)


# Load necessary libraries
library(survival)
library(survminer)
library(ggplot2)

# Ensure numeric variables
df2$PFI_time <- as.numeric(df2$PFI_time)
df2$PFI <- as.numeric(df2$PFI)
df2$age <- as.numeric(df2$age)
df2$purity <- as.numeric(df2$purity)

# Convert 'expression_group ' to factor (categorical variable)
df2$expression_group  <- as.factor(df2$expression_group )

# Fit a Cox proportional hazards model adjusting for age and purity
cox_model <- coxph(Surv(PFI_time, PFI) ~ expression_group  + age + purity, data = df2)

# Perform pairwise comparisons for all groups in "expression_group "
pairwise_results <- pairwise_survdiff(Surv(PFI_time, PFI) ~ expression_group , data = df2, p.adjust.method = "BH")

# Print pairwise comparison results
print(pairwise_results)

# Print model summary
summary(cox_model)

# Plot adjusted survival curves
ggadjustedcurves(cox_model, data = df2, variable = "expression_group", method = "conditional") +
  ggtitle("Adjusted PFI Curves by ENSG00000272172.1 Expression") +
  theme_minimal()

cox.zph(cox_model)







library(ggplot2)
library(ggpubr)

# Create age group variable based on age ranges
df2$age_group <- cut(df2$age,
                     breaks = c(20, 40, 60, Inf),
                     labels = c("20-40", "41-60", "61+"),
                     right = FALSE)

# Violin plot with individual points and p-value from ANOVA
ggplot(df2, aes(x = age_group, y = ENSG00000272172.1)) +
  geom_violin(trim = FALSE, fill = "lightblue", alpha = 0.5) + # Violin plot
  geom_jitter(width = 0.1, color = "darkred", alpha = 0.6) + # Individual points
  stat_compare_means(method = "anova", label = "p.signif") + # Add p-value (ANOVA)
  labs(title = "Violin Plot: ENSG00000272172.1 Expression by Age Group",
       x = "Age Group",
       y = "ENSG00000272172.1 Expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        title = element_text(size = 14))




library(ggplot2)

# Scatter plot with linear regression line
ggplot(df2, aes(x = age, y = ENSG00000272172.1)) +
  geom_point(color = "darkred", alpha = 0.6) +  # Scatter plot of individual points
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  # Linear regression line
  labs(title = "Correlation Between Age and ENSG00000272172.1 Expression",
       x = "Age",
       y = "ENSG00000272172.1 Expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        title = element_text(size = 14)) +
  stat_cor(aes(label = ..r.label..), method = "pearson", label.x = 40, label.y = 5)  # Add Pearson correlation coefficient







