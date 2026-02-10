############################################################################################################# 2.
##### 2. Metabolites - Imputation (Table S4; Fig.S6&7)
setwd("~/Desktop")
rm(list=ls())
library(tidyverse)
library(dplyr)
library(ggplot2)
library(reshape2)
library(imputeLCMD) # 2.2 - QRILC imputation
library(glmnet) # 2.3.1 - LASSO feature selection
library(randomForest) # 2.3.2 - Boruta feature selection
library(Boruta) # 2.3.2 - Boruta feature selection
metabo <- read.csv("metabo_raw.csv")
############################################################################################################# 2.1
######################### 2.1 Investigating missingness ########################
metabo <- metabo[,-1]
missing_counts <- colSums(is.na(metabo))
missing_counts[missing_counts > 0] # values that were missing

## Count complete vs. incomplete cases (# of Participants containing missing values)
tibble(Status = c("Incomplete", "Complete"),
       Count = as.numeric(table(complete.cases(metabo))),
       Percent = round(100 * table(complete.cases(metabo)) / nrow(metabo), 2)) # 12.33% incomplete

## Number and proportion of missing values per variable
cbind("# NA" = sort(colSums(is.na(metabo))),
      "% NA" = round(sort(colMeans(is.na(metabo))) * 100, 2)) # Fig. S6
### COMMENTS:  Citrate_metabolit (1), His_metabolite (1), 
### bOHbutyrate_metabolite (29), and Acetoacetate_metabolite (57).

## Overall missingness
cat("\nOverall missingness:", round(100 * sum(is.na(metabo)) / prod(dim(metabo)), 2), "%\n\n") # 0.05%

############################################################################################################# 2.2
################### 2.2 Imputation, Diagnostics - Fig.S6 #######################
### Make sure data are NATURAL log transformed before plotting density plots and box plots
## Find columns that had NAs in the original data
ln_data <- log(metabo) # log transformation before imputation

set.seed(1234)
imp_data <- impute.QRILC(ln_data, tune.sigma = 1)[[1]] # QRILC imputation

na_mask <- is.na(ln_data) # masking the data

imp_long <- imp_data %>%
  as.data.frame() %>%
  mutate(row = row_number()) %>%
  pivot_longer(-row, names_to = "Variable", values_to = "Value")

mask_long <- na_mask %>%
  as.data.frame() %>%
  mutate(row = row_number()) %>%
  pivot_longer(-row, names_to = "Variable", values_to = "Data")

## Combine the imputed data and the mask
plot_data <- left_join(imp_long, mask_long, by = c("row", "Variable"))

plot_data$Data <- ifelse(plot_data$Data, "Imputed", "Observed")

## Select the variables that were missing
plot_imp <- plot_data %>%
  filter(Variable %in% c("His_metabolite", "Citrate_metabolite", 
                         "bOHbutyrate_metabolite", "Acetoacetate_metabolite")) %>%
  mutate(Variable = dplyr::recode(as.character(Variable),
                                  "His_metabolite" = "Histidine",
                                  "bOHbutyrate_metabolite" = "3-hydroxybutyrate",
                                  "Acetoacetate_metabolite" = "Acetoacetate",
                                  "Citrate_metabolite" = "Citrate"))

ggplot(plot_imp, aes(x = Variable, y = Value, color = Data)) +
  geom_jitter(width = 0.2, size = 2) +
  labs(title = "Stripplot of imputed metabolite data using QRILC imputation",
       y = "Log-transformed value") +
  scale_color_manual(values = c("Observed" = "darkblue", "Imputed" = "orange")) +
  theme_minimal() # Fig.S6

############################################################################################################# 2.3.1
################### 2.3.1 LASSO feature selection - Fig.S7a ####################
bmi <- read.csv("all data_preclean_ID_recoded.csv") 
x <- as.matrix(imp_data) # variables, already natural log transformed 
y <- bmi$BMI # continuous outcome
 
set.seed(1234)
lasso_cv <- cv.glmnet(x,y, alpha = 1, family = "gaussian")  # the optimal lambda based on CV error
plot(lasso_cv) # MSE across different lambda; Fig.S7a_B
lasso_coefs <- coef(lasso_cv, s = "lambda.min")  # leftmost line (minimum MSE during CV)
# lasso_coefs <- coef(lasso_cv, s = "lambda.1se") # rightmost line (the lambda within one SD of the minimum MSE)
### NOTES: if use rightmost line, intersect values will be less (lasso + boruta)

lambda_opt <- lasso_cv$lambda.min # optimal lambda value
lasso_results <- glmnet(x, y, alpha = 1, family = "gaussian") 
plot(lasso_results, xvar = 'lambda', label = TRUE, ) # lambdas in relations to the coefficients
abline(v = log(lambda_opt), col = "red", lty = 2) # Fig.S7a.A

## Filter only non-zero coefficients (selected features)
lasso_selected <- data.frame(
  Feature = rownames(lasso_coefs)[lasso_coefs[, 1] != 0][-1],
  Coefficient = lasso_coefs[lasso_coefs[, 1] != 0, 1][-1]) # 36 variables
rownames(lasso_selected) <- NULL # reset row index

############################################################################################################# 2.3.2
################## 2.3.2 Boruta feature selection - Fig.S7b ####################
boruta <- exp(imp_data) # convert it back to original data for Boruta 
boruta$BMI <- bmi$BMI

### NOTES: log is not needed for boruta algorithm. Use original data is more accurate
set.seed(1234)
boruta_results <- Boruta(BMI ~ ., data = boruta, doTrace = 1)

## Finalize feature selection
final_boruta <- TentativeRoughFix(boruta_results)

## Convert Boruta results to long format
df_long <- as.data.frame(final_boruta$ImpHistory) %>%
  pivot_longer(cols = everything(), names_to = "Feature", values_to = "Importance") %>%
  left_join(attStats(final_boruta) %>% rownames_to_column(var = "Feature") %>% select(Feature, decision), by = "Feature") %>%
  mutate(decision = if_else(is.na(decision), "Reference", as.character(decision))) %>%
  filter(is.finite(Importance))

## Plot Boruta feature importance; 
boruta_colors <- c("Confirmed" = "darkgreen", "Rejected" = "orange", "Tentative" = "gray", "Reference" = "darkblue")

ggplot(df_long, aes(x = reorder(Feature, Importance, FUN = median), y = Importance, fill = decision)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 1, outlier.alpha = 0.6) +
  scale_fill_manual(values = boruta_colors) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, size = 3)) +
  labs(title = "Boruta Feature Importance", x = "Features", y = "Importance Score", fill = "Decision") # Fig.S7b

############################################################################################################# 2.3.3
###################### 2.3.3 Intersect feature - Table S4 ######################
## Get confirmed important features
boruta_selected <- getSelectedAttributes(final_boruta, withTentative = FALSE)

## Extract meanImp scores
meanImp_scores <- attStats(final_boruta) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Feature") %>%
  select(Feature, meanImp)

## Intersect features 
final_features <- intersect(lasso_selected$Feature, boruta_selected)
final_features 
final_features_df <- data.frame(Feature = final_features) %>%
  left_join(lasso_selected, by = "Feature") %>%
  left_join(meanImp_scores %>% select(Feature, meanImp), by = "Feature") %>%
  arrange(desc(meanImp)) # Table S4

met_clean <- imp_data[,final_features] # still log
met_clean <- exp(met_clean)
met_clean$ID <- metabo$ID

#write.csv(met_clean,"~/Desktop/met_clean.csv", row.names=FALSE)