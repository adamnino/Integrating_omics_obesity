############################################################################################################# 10.
##### 10. Reduced ML modelling without no covariates for WC (Table S15)
setwd("~/Desktop")
rm(list=ls()) 
library(dplyr)
library(tidyverse)
library(foreach) # 10.1 - Reduced ML models
library(doParallel) # 10.1 - Reduced ML models
library(doRNG) # 10.1 - Reduced ML models
df <- read.csv("final data_17SNP.csv")
df$S_LDL_P_metabolite <- df$S_LDL_P_metabolite*1000 # changing units
df$InsulinpmolL <- log(df$InsulinpmolL) # only variable needs to be log-transformed based on 6_demo_PCA&VIF
############################################################################################################# 10.1
################# 10.1 Reduced ML models (no covariates) for WC #################
df <- df %>% 
  mutate(across(c(2,3,6,8,9,11,13,15,30,32,33), as.factor)) %>% 
  mutate(across(c(4,5,7,16:29,31,34:55), as.numeric)) %>% 
  mutate(across(c(1,10,12,14), as.character))

## Reduced ML models for WC - Table S15
model_list <- list(
  model1 = df[, c(16,37,38,44,45)], # model 1 = 4-Metabo.
  model1.1 = df[, c(16,37,38,41,44,45)], # model 1.1 = 4-Metabo. + Albumin
  model2 = df[, c(16,17,27)], # model 2 = 2-Clinc.
  model3 = df[, c(16,17,27,37,38)],  # model 3 = 2-Clinc. + AAA
  model3.1 = df[, c(16,17,27,37,38,41)],  # model 3.1 = 2-Clinc. + AAA + Albumin
  model3.2 = df[, c(16,17,27,44,45)],  # model 3.2 = 2-Clinc. + Lipoproteins
  model3.3 = df[, c(16,17,27,41,44,45)],  # model 3.3 = 2-Clinc. + Lipoproteins + Albumin
  model3.4 = df[, c(16,17,27,37,38,44,45)], # model 3.4 = 2-Clinc. + 4-Metabo.
  model3.5 = df[, c(16,17,27,37,38,41,44,45)] # model 3.5 = 2-Clinc. + 4-Metabo. + Albumin
)

set.seed(1053) # repeated 10 times * 5 outer folds * 3 inner folds
repeats_outer <- 10
cores <- parallel::detectCores() - 1 # use all cores except one (7 cores)
registerDoParallel(cores)
registerDoRNG(seed = 1234)

## Reproducible bootstrap CI function; Additionally, this function handles edge cases 
## (e.g., missing or constant values in single-predictor models). If bootstrapping was 
## not possible, the mean was reported and the confidence interval set equal to it. If 
##there are no NAs, this function produces the same outcomes as the one in 8_demo_ML model.
boot_ci <- function(x, seed = 1000) {
  x <- x[!is.na(x)]
  if (length(x) < 2 || length(unique(x)) == 1) return(rep(unique(x), 3))
  
  set.seed(seed)
  b <- try(boot(x, function(d, i) mean(d[i]), R = 1000), silent = TRUE)
  ci <- try(boot.ci(b, type = "perc")$percent[4:5], silent = TRUE)
  
  c(mean = mean(x), lower = ifelse(inherits(ci, "try-error"), NA, ci[1]),
    upper = ifelse(inherits(ci, "try-error"), NA, ci[2]))
} 

## Run models with reproducibility (Seed & RNG); NOTE: uncomment "noise data" for single-models
results_list <- foreach(model_name = names(model_list), 
                        .packages = c("caret", "boot"), 
                        .options.RNG = 1234) %dorng% {
                          
  data <- model_list[[model_name]]
  importance_list <- list() 
  hyperparams_list <- list()
  metric_storage <- list()
  calibration_list <- list()
                          
  for (r in 1:repeats_outer) {
    outer_folds <- createFolds(data$WC, k = 5)
                            
    for (f in seq_along(outer_folds)) {
      test_idx <- outer_folds[[f]]
      train_data <- data[-test_idx, ]
      test_data <- data[test_idx, ]
                              
      tuned_model <- train(
        WC ~ ., data = train_data,
        method = "glmnet", # regularized linear regression (GLMNET)
        preProcess = c("center", "scale"), # GLMNET models
        tuneLength = 10, # 10 tuning combinations
        trControl = trainControl(method = "cv", number = 3, allowParallel = TRUE, search = "random"))
                              
      ## Predictions
      pred_train <- predict(tuned_model, newdata = train_data)
      pred_test  <- predict(tuned_model, newdata = test_data)
      actual_train <- train_data$WC
      actual_test  <- test_data$WC
                              
      ## Evaluation metrics
      metrics <- data.frame(
        Set = c("Train", "Test"),
        R2 = c(R2(pred_train, actual_train), R2(pred_test, actual_test)),
        RMSE = c(RMSE(pred_train, actual_train), RMSE(pred_test, actual_test)),
        MAE = c(MAE(pred_train, actual_train), MAE(pred_test, actual_test)),
        Repeat = r,
        Fold = f)
                              
      metric_storage[[length(metric_storage) + 1]] <- metrics
                              
      ## Variable importance
      importance_vals <- as.matrix(coef(tuned_model$finalModel, s = tuned_model$bestTune$lambda)) # GLMNET model
      importance_vals <- importance_vals[rownames(importance_vals) != "(Intercept)", , drop = FALSE] # GLMNET model
                              
      importance_df <- data.frame(
        Variable = rownames(importance_vals),
        Importance = as.numeric(importance_vals), # GLMNET model
        Repeat = r,
        Fold = f)
      importance_list[[length(importance_list) + 1]] <- importance_df
                              
      ## Hyperparameter tracking
      best_params <- tuned_model$bestTune
      best_params$Repeat <- r
      best_params$Fold <- f
      hyperparams_list[[length(hyperparams_list) + 1]] <- best_params
                              
      ## Calibration data
      calibration_data <- data.frame(
        Actual = actual_test,
        Predicted = pred_test,
        Repeat = r,
        Fold = f)
      calibration_list[[length(calibration_list) + 1]] <- calibration_data
    }
  }
                          
  ## Combine all metrics
  all_metrics <- do.call(rbind, metric_storage)
                          
  ## Calculate bootstrapped CI for each metric per dataset (train & test)
  ci_results <- all_metrics %>%
    group_by(Set) %>%
    summarise(
      R2_mean = mean(R2),
      R2_lower = boot_ci(R2)[2],
      R2_upper = boot_ci(R2)[3],
                              
      RMSE_mean = mean(RMSE),
      RMSE_lower = boot_ci(RMSE)[2],
      RMSE_upper = boot_ci(RMSE)[3],
                              
      MAE_mean = mean(MAE),
      MAE_lower = boot_ci(MAE)[2],
      MAE_upper = boot_ci(MAE)[3])
                          
  list(
    Summary = ci_results,
    All_Metrics = all_metrics,
    Importance = do.call(rbind, importance_list),
    Hyperparams = do.call(rbind, hyperparams_list),
    Calibration = do.call(rbind, calibration_list))
}

stopImplicitCluster() # shut down parallel backend

names(results_list) <- names(model_list) # assign model names
 
############################################################################################################# 10.2.1
############## 10.2.1 Data extractions - Evaluation metrics 95% CI ##############
########################## Step 1 - Extract all metrics ########################
all_metrics <- bind_rows(lapply(names(results_list), function(model) {
  metrs <- results_list[[model]]$All_Metrics
  metrs$Model <- model
  metrs}))

########## Step 2 - Replace NAs with group mean (by Model and Dataset) ##########
all_metrics <- all_metrics %>%
  group_by(Model, Set) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
  ungroup()

########## Step 3 - Calculate mean values for each variable-model pair #########
mean_metrics <- all_metrics %>%
  group_by(Model, Set) %>%
  summarise(R2 = mean(R2, na.rm = TRUE), 
            RMSE = mean(RMSE, na.rm = TRUE),
            MAE = mean(MAE, na.rm = TRUE),
            .groups = "drop")

############ Step 4 - Extract summary of evaluation metrics & 95% CI ###########
all_avg <- bind_rows(
  lapply(names(results_list), function(model_name) {
    results_list[[model_name]]$Summary %>%
      pivot_longer(cols = -Set, 
                   names_to = c("Metric", ".value"),
                   names_pattern = "(.*)_(mean|lower|upper)") %>%
      rename(Mean = mean, CI_Lower = lower, CI_Upper = upper) %>%
      mutate(Model = model_name)}))

#######################  Step 5 - Formatting the results #######################
format_num <- function(x) {
  ifelse(is.na(x), NA, format(round(x, 3), nsmall = 3))} # to ensure 3 decimal places in the results
df_wide <- all_avg %>%
  mutate(
    across(c(Mean, CI_Lower, CI_Upper), as.numeric),
    metric_set = paste0(Metric, "_", tolower(Set))) %>%
  pivot_wider(
    id_cols = Model,
    names_from = metric_set,
    values_from = c(Mean, CI_Lower, CI_Upper))

## Calculate differences (Train - Test) for each metric (Optional)
df_wide <- df_wide %>%
  mutate(
    ## R2 difference
    diff_R2 = Mean_R2_train - Mean_R2_test,
    diff_R2_CI_lower = CI_Lower_R2_train - CI_Upper_R2_test, 
    diff_R2_CI_upper = CI_Upper_R2_train - CI_Lower_R2_test,
    
    ## RMSE difference
    diff_RMSE = Mean_RMSE_train - Mean_RMSE_test,
    diff_RMSE_CI_lower = CI_Lower_RMSE_train - CI_Upper_RMSE_test,
    diff_RMSE_CI_upper = CI_Upper_RMSE_train - CI_Lower_RMSE_test,
    
    ## MAE difference
    diff_MAE = Mean_MAE_train - Mean_MAE_test,
    diff_MAE_CI_lower = CI_Lower_MAE_train - CI_Upper_MAE_test,
    diff_MAE_CI_upper = CI_Upper_MAE_train - CI_Lower_MAE_test)

## Final formatted result
df_final <- df_wide %>%
  mutate(
    across(where(is.numeric), format_num), # 3 decimal places
    
    R2_train = paste(Mean_R2_train, " [", CI_Lower_R2_train, ", ", CI_Upper_R2_train, "]", sep = ""),
    RMSE_train = paste(Mean_RMSE_train, " [", CI_Lower_RMSE_train, ", ", CI_Upper_RMSE_train, "]", sep = ""),
    MAE_train = paste(Mean_MAE_train, " [", CI_Lower_MAE_train, ", ", CI_Upper_MAE_train, "]", sep = ""),
    
    R2_test = paste(Mean_R2_test, " [", CI_Lower_R2_test, ", ", CI_Upper_R2_test, "]", sep = ""),
    RMSE_test = paste(Mean_RMSE_test, " [", CI_Lower_RMSE_test, ", ", CI_Upper_RMSE_test, "]", sep = ""),
    MAE_test = paste(Mean_MAE_test, " [", CI_Lower_MAE_test, ", ", CI_Upper_MAE_test, "]", sep = ""),
    
    diff_R2 = paste(diff_R2, " [", diff_R2_CI_lower, ", ", diff_R2_CI_upper, "]", sep = ""),
    diff_RMSE = paste(diff_RMSE, " [", diff_RMSE_CI_lower, ", ", diff_RMSE_CI_upper, "]", sep = ""),
    diff_MAE = paste(diff_MAE, " [", diff_MAE_CI_lower, ", ", diff_MAE_CI_upper, "]", sep = "")) %>%
  select(
    Model,
    R2_train, RMSE_train, MAE_train,
    R2_test, RMSE_test, MAE_test,
    diff_R2, diff_RMSE, diff_MAE)
view(df_final) # Table S15
