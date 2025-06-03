library(tidyverse)


 
 df2 <- suvival_cluster_age_purity
 
 
 # Load necessary libraries
 library(survival)
 library(survminer)
 library(ggplot2)
 library(dplyr)
 
 df2$six <- as.factor(df2$six)
 # Ensure OS_time and OS are numeric
 df2$OS_time <- as.numeric(df2$OS_time)
 df2$OS <- as.numeric(df2$OS)
 # Define the survival object
 surv_obj <- Surv(time = df2$OS_time, event = df2$OS)
 
 # Fit the Kaplan-Meier model
 fit <- survfit(surv_obj ~ six, data = df2)
 
 ggsurvplot(
   fit,
   data = df2,
   pval = TRUE,                   # Show p-value for significance
   conf.int = FALSE,               # Show confidence intervals
   risk.table = TRUE,             # Show number at risk table
   palette = "Dark2",             # Color palette for different groups
   ggtheme = theme_minimal(),      # Apply a clean theme
   legend.title = "Group (six)",   # Title for the legend
   legend.labs = levels(df2$six)   # Labels for the legend
 )
 
 
 df2$six <- as.factor(df2$six)
 
 # Ensure PFI_time and PFI are numeric
 df2$PFI_time <- as.numeric(df2$PFI_time)
 df2$PFI <- as.numeric(df2$PFI)
 
 # Define the survival object for PFI
 surv_obj_pfi <- Surv(time = df2$PFI_time, event = df2$PFI)
 
 # Fit the Kaplan-Meier model for PFI
 fit_pfi <- survfit(surv_obj_pfi ~ six, data = df2)
 
 # Plot Kaplan-Meier curves for PFI
 ggsurvplot(
    fit_pfi,
    data = df2,
    pval = TRUE,                    # Show p-value for significance
    conf.int = FALSE,               # Show confidence intervals
    risk.table = TRUE,             # Show number at risk table
    palette = "Dark2",             # Color palette for different groups
    ggtheme = theme_minimal(),     # Apply a clean theme
    legend.title = "Group (six)",  # Title for the legend
    legend.labs = levels(df2$six)  # Labels for the legend
 )
 
 ########### R Code for Cox Model with Confounders:
 
 # Load necessary libraries
 library(survival)
 library(survminer)
 library(ggplot2)
 
 # Ensure numeric variables
 df2$OS_time <- as.numeric(df2$OS_time)
 df2$OS <- as.numeric(df2$OS)
 df2$age <- as.numeric(df2$age)
 df2$purity <- as.numeric(df2$purity)
 
 # Convert 'six' to factor (categorical variable)
 df2$six <- as.factor(df2$six)
 
 # Fit a Cox proportional hazards model adjusting for age and purity
 cox_model <- coxph(Surv(OS_time, OS) ~ six + age + purity, data = df2)
 
 # Perform pairwise comparisons for all groups in "six"
 pairwise_results <- pairwise_survdiff(Surv(OS_time, OS) ~ six, data = df2, p.adjust.method = "BH")
 
 # Print pairwise comparison results
 print(pairwise_results)
 
 # Print model summary
 summary(cox_model)
 
 # Plot adjusted survival curves
 ggadjustedcurves(cox_model, data = df2, variable = "six", method = "conditional") +
   ggtitle("Adjusted Survival Curves by Six Groups") +
   theme_minimal()
 
 cox.zph(cox_model)
 
 
 ############### PFI
 
 # Ensure numeric variables
 df2$PFI_time <- as.numeric(df2$PFI_time)
 df2$PFI <- as.numeric(df2$PFI)
 df2$age <- as.numeric(df2$age)
 df2$purity <- as.numeric(df2$purity)
 
 # Convert 'six' to factor (categorical variable)
 df2$six <- as.factor(df2$six)
 
 # Fit a Cox proportional hazards model adjusting for age and purity
 cox_model <- coxph(Surv(OS_time, OS) ~ six + age + purity, data = df2)
 
 # Perform pairwise comparisons for all groups in "six"
 pairwise_results <- pairwise_survdiff(Surv(OS_time, OS) ~ six, data = df2, p.adjust.method = "BH")
 
 # Print pairwise comparison results
 print(pairwise_results)
 
 # Print model summary
 summary(cox_model)
 
 # Plot adjusted survival curves
 ggadjustedcurves(cox_model, data = df2, variable = "six", method = "conditional") +
    ggtitle("Adjusted Survival Curves by Six Groups") +
    theme_minimal()
 
 cox.zph(cox_model)
 
 
 
 library(multcomp)
 
 # Fit Cox model
 cox_model_OS <- coxph(Surv(OS_time, OS) ~ six + age + purity, data = df2)
 
 # Perform multiple comparisons for "six" categories
 summary(glht(cox_model, linfct = mcp(six = "Tukey")))
 
 
 # Fit Cox model
 cox_model_PFI <- coxph(Surv(PFI_time, PFI) ~ six + age + purity, data = df2)
 
 # Perform multiple comparisons for "six" categories
 summary(glht(cox_model, linfct = mcp(six = "Tukey")))
 
 
 # Load required package
 library(forestmodel)
 
 # Generate forest plot of Cox model results
 forest_plot_OS <- forest_model(cox_model_OS)
 forest_plot_PFI <- forest_model(cox_model_PFI)
 # Print the plot
 
 print(forest_plot_OS )
 print(forest_plot_PFI )
 
 
 # Load required libraries
 library(dplyr)
 
 # Ensure that PFI_time is numeric
 df2$PFI_time <- as.numeric(df2$PFI_time)
 
 # Summarizing key variables excluding gender
 summary_table <- df2 %>%
    dplyr::summarise(
       N = n(),
       Age_Mean = mean(age, na.rm = TRUE),
       Age_Median = median(age, na.rm = TRUE),
       Age_Range = paste(min(age, na.rm = TRUE), "to", max(age, na.rm = TRUE)),
       Tumor_Purity_Mean = mean(purity, na.rm = TRUE),
       Tumor_Purity_Median = median(purity, na.rm = TRUE),
       Tumor_Purity_Range = paste(min(purity, na.rm = TRUE), "to", max(purity, na.rm = TRUE)),
       OS_Median = median(OS_time, na.rm = TRUE),
       OS_Range = paste(min(OS_time, na.rm = TRUE), "to", max(OS_time, na.rm = TRUE)),
       PFS_Median = median(PFI_time, na.rm = TRUE),
       PFS_Range = paste(min(PFI_time, na.rm = TRUE), "to", max(PFI_time, na.rm = TRUE)),
       Death_Count = sum(OS == 1),
       Alive_Count = sum(OS == 0),
       Progression_Count = sum(PFI == 1),
       No_Progression_Count = sum(PFI == 0)
    )
 
 # Print the summary table
 print(summary_table)
 # Load necessary library
 library(knitr)
 
 # Create a data frame with the summarized information
 summary_data <- data.frame(
   N = 164,
   Age_Mean = 58.81707,
   Age_Median = 57,
   Age_Range = "34 to 87",
   Tumor_Purity_Mean = 0.8005625,
   Tumor_Purity_Median = 0.81,
   Tumor_Purity_Range = "0.41 to 1",
   OS_Median = 1004,
   OS_Range = "11 to 5481",
   PFS_Median = 426,
   PFS_Range = "11 to 5481",
   Death_Count = 92,
   Alive_Count = 72,
   Progression_Count = 119,
   No_Progression_Count = 45
 )
 
 # Print the table in a nice format using knitr
 kable(summary_data, caption = "Patient Summary Table", format = "markdown", digits = 2)
 
