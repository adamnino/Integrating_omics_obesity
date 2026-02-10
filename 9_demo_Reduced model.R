############################################################################################################# 9.
##### 9. Reduced ML modelling with/without no covariates (Table 2; Fig.3; Table S13,14,16,17; Fig.S13&14)
setwd("~/Desktop") 
rm(list=ls())
library(dplyr)
library(tidyverse)
library(foreach) # 9.1 - Reduced ML models
library(doParallel) # 9.1 - Reduced ML models
library(doRNG) # 9.1 - Reduced ML models
library(patchwork) # 9.3.1 - Figures
library(ggpattern) # 9.3.1 - Figures
library(ggpointdensity) # 9.3.2 - Figures (Density&calibration)
library(extrafont) # 9.3.2 - Figures (Font)
library(ggcorrplot) # 9.3.4 - Correlation
df <- read.csv("final data_17SNP.csv") 
df$S_LDL_P_metabolite <- df$S_LDL_P_metabolite*1000 # changing units
df$InsulinpmolL <- log(df$InsulinpmolL) # only variable needs to be log-transformed based on 6_demo_PCA&VIF
############################################################################################################# 9.1
##################### 9.1 Reduced ML models (no covariates) ####################
df <- df %>% 
  mutate(across(c(2,3,6,8,9,11,13,15,30,32,33), as.factor)) %>% 
  mutate(across(c(4,5,7,16:29,31,34:55), as.numeric)) %>% 
  mutate(across(c(1,10,12,14), as.character))

#### NOTE: for simplicity and clarity, the code below is designed for Table 2; Fig.3; 
#### Table S13&14; and Fig.S13. To generate Table S16&17 and Fig.S14, simply including 
#### covariates (age, sex, districts) in the models. Additionally, single-predictor models 
#### were run independently to avoid potential numerical instabilities in parallelized training 
#### of GLMNET. A small noise variable was included to prevent algorithmic errors in models with 
#### only one predictor.

## Reduced ML models - Fig.3, Fig.S13&14
model_list <- list(
  model1 = df[, c(5,37,38,44,45)], # model 1 = 4-Metabo.
  model1.1 = df[, c(5,37,38,41,44,45)], # model 1.1 = 4-Metabo. + Albumin
  model2 = df[, c(5,17,27)], # model 2 = 2-Clinc.
  model3 = df[, c(5,17,27,37,38)],  # model 3 = 2-Clinc. + AAA
  model3.1 = df[, c(5,17,27,37,38,41)],  # model 3.1 = 2-Clinc. + AAA + Albumin
  model3.2 = df[, c(5,17,27,44,45)],  # model 3.2 = 2-Clinc. + Lipoproteins
  model3.3 = df[, c(5,17,27,41,44,45)],  # model 3.3 = 2-Clinc. + Lipoproteins + Albumin
  model3.4 = df[, c(5,17,27,37,38,44,45)], # model 3.4 = 2-Clinc. + 4-Metabo.
  model3.5 = df[, c(5,17,27,37,38,41,44,45)], # model 3.5 = 2-Clinc. + 4-Metabo. + Albumin

  # model1a = df[, c(5,41)], # model 1a = Albumin
  # model1b = df[, c(5,44)],  # model 1b = S_LDL_P
  # model1c = df[, c(5,37)], # model 1c = Phe
  # model1d = df[, c(5,38)], # model 1d = Tyr
  # model1e = df[, c(5,45)],  # model 1f = L_HDL_CE
  # model1f = df[, c(5,42)],  # model 1e = GlycA
  # model2a = df[, c(5,17)], # model 2a = SBP
  # model2b = df[, c(5,27)], # model 2b = Insulin

  model3a = df[, c(5,17,27,41)],  # model 3a = 2-Clinc. + Albumin
  model3b = df[, c(5,17,27,42)],  # model 3b = 2-Clinc. + GlycA
  model3c = df[, c(5,17,27,41,42)],  # model 3c = 2-Clinc. + GlycA + Albumin
  model3d = df[, c(5,17,27,37,38,42)],  # model 3d = 2-Clinc. + AAA + GlycA
  model3e = df[, c(5,17,27,42,44,45)],  # model 3e = 2-Clinc. + Lipoprotein + GlycA
  model3f = df[, c(5,17,27,37,38,42,44,45)], # model 3f = 2-Clinc. + 4-Metabo. + GlycA
  model3g = df[, c(5,17,27,37,38,41,42,44,45)], # model 3g = 2-Clinc. + 4-Metabo. + GlycA + Albumin
  model3h = df[, c(5,17,27,37,38,41,42)],  # model 3h = 2-Clinc. + AAA + GlycA + Albumin
  model3i = df[, c(5,17,27,41,42,44,45)],  # model 3i = 2-Clinc. + Lipoprotein + GlycA + Albumin

  model4a = df[, c(5,37,38,42,44,45)], # model 4a = 4-Metabo. + GlycA
  model4b = df[, c(5,37,38,41,42,44,45)], # model 4b = 4-Metabo. + GlycA + Albumin

  model5a = df[, c(5,16,17,27)],  # model 5a = 2-Clinc. + WC
  model5b = df[, c(5,16,37,38,41,44,45)],  # model 5b = 4-Metabo. + Albumin + WC
  model5c = df[, c(5,16,17,27,37,38,41,44,45)]  # model 5c = 2-Clinc. + 4-Metabo. + Albumin + WC
)

set.seed(1053) # repeated 10 times * 5 outer folds * 3 inner folds
repeats_outer <- 10
cores <- parallel::detectCores() - 1 # use all cores except one (7 cores)
registerDoParallel(cores)
registerDoRNG(seed = 1234)

## Reproducible bootstrap CI function; Additionally, this function handles edge cases 
## (e.g., missing or constant values in single-predictor models 1a–1e, 2a–2b). If 
## bootstrapping was not possible, the mean was reported and the confidence interval 
## set equal to it. If there are no NAs, this function produces the same outcomes 
## as the one in 8_demo_ML model.
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
    outer_folds <- createFolds(data$BMI, k = 5)

    for (f in seq_along(outer_folds)) {
      test_idx <- outer_folds[[f]]
      train_data <- data[-test_idx, ]
      test_data <- data[test_idx, ]
      
      # train_data$Noise <- rnorm(nrow(train_data)) # small noise for single-models; w/o covariates
      # test_data$Noise  <- rnorm(nrow(test_data)) # small noise for single-models; w/o covariates only
      
      tuned_model <- train(
        BMI ~ ., data = train_data,
        method = "glmnet", # regularized linear regression (GLMNET)
        preProcess = c("center", "scale"), # GLMNET models
        tuneLength = 10, # 10 tuning combinations
        trControl = trainControl(method = "cv", number = 3, allowParallel = TRUE, search = "random"))

      ## Predictions
      pred_train <- predict(tuned_model, newdata = train_data)
      pred_test  <- predict(tuned_model, newdata = test_data)
      actual_train <- train_data$BMI
      actual_test  <- test_data$BMI

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

############################################################################################################# 9.2.1
############## 9.2.1 Data extractions - Evaluation metrics 95% CI ##############

#### NOTE: For reduced GLMNET model, 6 NAs were observed for model1a.This was 
#### expected as R2 were extremely small in those models. Run Steps 1-6.

########################## Step 1 - Extract all metrics ########################
all_metrics <- bind_rows(lapply(names(results_list), function(model) {
  metrs <- results_list[[model]]$All_Metrics
  metrs$Model <- model
  metrs}))
view(all_metrics) # check which model contain NAs
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
view(mean_metrics) # copy averaged values to Step 5 after running Step 4
############ Step 4 - Extract summary of evaluation metrics & 95% CI ###########
all_avg <- bind_rows(
  lapply(names(results_list), function(model_name) {
    results_list[[model_name]]$Summary %>%
      pivot_longer(cols = -Set, 
                   names_to = c("Metric", ".value"),
                   names_pattern = "(.*)_(mean|lower|upper)") %>%
      rename(Mean = mean, CI_Lower = lower, CI_Upper = upper) %>%
      mutate(Model = model_name)}))
view(all_avg) # for GLMNET model, check which row and column the NA is in
###### Step 5 - Replace NAs in Step 4 with the means calculated in Step 3 ######
# all_avg[1, 3] <- 0.004147951 # GLMNET model 1a, test (W/O covariates only)
# all_avg[4, 3] <- 0.009949075 # GLMNET model 1a, train (W/O covariates only)

# all_avg[55, 3] <- 0.004147951 # GLMNET model 1a, test (WITH covariates only)
# all_avg[58, 3] <- 0.009949075 # GLMNET model 1a, train (WITH covariates only)
#######################  Step 6 - Formatting the results #######################
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
view(df_final) # Table S13 or Table S16

############################################################################################################# 9.2.2
################### 9.2.2 Data extractions - hyperparameters ###################
## Check if R2 (test set) are consistent across different hyperparameter combinations
all_hyperparams <- bind_rows(lapply(names(results_list), function(model) {
  hpams <- results_list[[model]]$Hyperparams
  hpams$Model <- model
  hpams}))

labels_reduced_model <- c(
  model1 = "4-Metabo.",
  model1.1 = "4-Metabo. + Albumin",
  model2 = "2-Clinc.",
  model3 = "2-Clinc. + AAA",
  model3.1 = "2-Clinc. + AAA + Albumin",
  model3.2 = "2-Clinc. + Lipoproteins",
  model3.3 = "2-Clinc. + Lipoproteins + Albumin",
  model3.4 = "2-Clinc. + 4-Metabo.",
  model3.5 = "2-Clinc. + 4-Metabo. + Albumin",
  # model1a = "Albumin", # fig.s13a
  # model1b = "S_LDL_P", # fig.s13a
  # model1c = "Phenylalanine", # fig.s13a
  # model1d = "Tyrosine", # fig.s13a
  # model1e = "L_HDL_CE", # fig.s13a
  # model1f = "GlycA", # fig.s13a
  # model2a = "SBP", # fig.s13a
  # model2b = "Insulin", # fig.s13a
  model3a = "2-Clinc. + Albumin",
  model3b = "2-Clinc. + GlycA",
  model3c = "2-Clinc. + GlycA + Albumin",
  model3d = "2-Clinc. + AAA + GlycA",
  model3e = "2-Clinc. + Lipoproteins + GlycA",
  model3f = "2-Clinc. + 4-Metabo. + GlycA",
  model3g = "2-Clinc. + 4-Metabo. + GlycA + Albumin",
  model3h = "2-Clinc. + AAA + GlycA + Albumin",
  model3i =  "2-Clinc. + Lipoprotein + GlycA + Albumin",
  model4a = "4-Metabo. + GlycA",
  model4b = "4-Metabo. + GlycA + Albumin",
  model5a = "2-Clinc. + WC",
  model5b = "4-Metabo. + Albumin + WC",
  model5c = "2-Clinc. + 4-Metabo. + Albumin + WC"
) # Fig.3; fig.s13a

all_hyper_metric <- all_hyperparams %>%
  mutate(
    R2 = all_metrics$R2[all_metrics$Set == "Test"],
    Model = factor(Model, levels = names(labels_reduced_model)))

ggplot(all_hyper_metric, aes(x = lambda, y = alpha, color = R2)) +
  geom_point() +
  facet_wrap(~ Model, labeller = as_labeller(labels_reduced_model)) + 
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_minimal() +
  labs(
    # title = "Hyperparameters Combinations For Each Model Colored By R2",
    x = "Lambda (log10 scale)", y = "Alpha") +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(size = 6)) # Fig.S13a or Fig.S14a

############################################################################################################# 9.2.3
############### 9.2.3 Data extractions - importance/coefficients ###############
all_importance <- bind_rows(lapply(names(results_list), function(model) {
  iprte <- results_list[[model]]$Importance
  iprte$Model <- model
  iprte}))

## Mean importance for each variable-model pair
avg_importance <- all_importance %>%
  group_by(Model, Variable) %>%
  summarise(AvgImportance = mean(Importance, na.rm = TRUE), .groups = "drop")

############################################################################################################# 9.2.4
## 9.2.4 Data extractions - p values b/t models on test sets - Table 2,S14&16 ##
df_test <- all_metrics %>%
  filter(Set == "Test") %>% 
  select(R2, Repeat, Fold, Model) %>% # change R2/RMSE/MAE
  pivot_wider(names_from = Model, values_from = R2) %>% # change R2/RMSE/MAE
  mutate(id = paste0("Rep", Repeat, "_Fold", Fold))

long_df_T1_3 <- df_test %>% 
  pivot_longer(cols = c("model3a","model3b","model3c","model3d","model3e","model3f",
                        "model3g","model3h","model3i","model4a","model4b","model3",
                        "model3.1","model3.2","model3.3","model3.4","model3.5","model1",
                        "model1.1","model2"),
               names_to = "Model", values_to = "R2") # change R2/RMSE/MAE

## Pairwise Wilcoxon signed-rank test
p_value <- pairwise.wilcox.test(long_df_T1_3$R2,  # change R2/RMSE/MAE
                     long_df_T1_3$Model, paired = TRUE, p.adjust.method = "none")  # Table 2,S14&16
view(p_value$p.value)

############################################################################################################# 9.3.1
######################## 9.3.1 Bar plots Fig.3 - GLM ###########################
labels_fig.3 <- c(
  model1 = "4-Metabo.",
  model1.1 = "4-Metabo. + Albumin",
  model2 = "2-Clinc.",
  model3 = "AAA",
  model3.1 = "AAA + Albumin",
  model3.2 = "Lipoproteins",
  model3.3 = "Lipoproteins + Albumin",
  model3.4 = "4-Metabo.",
  model3.5 = "4-Metabo. + Albumin"
) # For Fig.3

model_fig.3 <- data.frame(
  Model_Name = names(labels_fig.3),
  Display_Name = unname(labels_fig.3),
  Group = c(
    rep("Single model", 3), # model1, model1.1, model2
    rep("2-Clinc. +", 6)) # model3, model3.1....model3.5
) # For Fig.3

all_avg_fig.3 <- all_avg %>%
  filter(Model %in% c("model1", "model1.1", "model2","model3","model3.1","model3.2",
                      "model3.3","model3.4","model3.5"))

plot_data <- all_avg_fig.3 %>%
  filter(Set == "Test") %>%
  mutate(Model = factor(Model, levels = names(labels_fig.3))) %>%
  left_join(model_fig.3, by = c("Model" = "Model_Name")) %>%
  mutate(Group = factor(Group, levels = unique(model_fig.3$Group)))

plot_metric <- function(metric_name, y_label_expr) {
  plot_data %>%
    filter(Metric == metric_name) %>%
    ggplot(aes(x = Display_Name, y = Mean)) +
    geom_bar(stat = "identity", width = 0.6, fill = "white", color = "black") +
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, color = "black") +
    #geom_text(aes(label = sprintf("%.3f", Mean), y = CI_Upper + 0.01), size = 3.5, vjust = -0.5) + # decimal dot
    geom_text(aes(label = gsub("\\.", "\u00B7", sprintf("%.3f", Mean)), y = CI_Upper + 0.01),size = 3.5, vjust = -0.5) + # middle dot
    scale_y_continuous(labels = function(x) gsub("\\.", "\u00B7", sprintf("%.3f", x))) +  # middle dot
    labs(x = NULL, y = y_label_expr) +
    facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
    theme_bw(base_size = 15.5, base_family = "Times") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "white", color = "black"),
          strip.text = element_text(face = "bold"),
          panel.spacing.x = unit(0, "pt"))} # Fig.3

## Combine plots
p1 <- plot_metric("R2", expression(R^2))
p1 # Fig.3
#p2 <- plot_metric("RMSE", "RMSE") # optional figure
#p3 <- plot_metric("MAE", "MAE") # optional figure

############################################################################################################# 9.3.2
############## 9.3.2 Calibration plots - GLM - Fig.S13b&c;S14b&c ###############
## Density & calibration plots
all_calibration_data <- lapply(names(results_list), function(model) {
  df <- results_list[[model]]$Calibration
  df$Model <- model
  df}) %>% bind_rows()

calibration_data_fig.S13 <- all_calibration_data %>%
  filter(Model %in% c("model1", "model1.1", "model2","model3","model3.1","model3.2",
                      "model3.3","model3.4","model3.5")) %>%
  mutate(Model = forcats::fct_relevel(Model, "model1", "model1.1", "model2",
                                      "model3","model3.1","model3.2","model3.3","model3.4","model3.5")) # Fig.S13
labels_fig.S13 <- c(
  model1 = "4-Metabo.",
  model1.1 = "4-Metabo. + Albumin",
  model2 = "2-Clinc.",
  model3 = "2-Clinc. + AAA",
  model3.1 = "2-Clinc. + AAA + Albumin",
  model3.2 = "2-Clinc. + Lipoproteins",
  model3.3 = "2-Clinc. + Lipoproteins + Albumin",
  model3.4 = "2-Clinc. + 4-Metabo.",
  model3.5 = "2-Clinc. + 4-Metabo. + Albumin"
) # Fig.S13b&c;S14b&c

metrics <- df_final%>% # df_final is in 9.3.1 - Step 6; for Fig.S13b&c;S14b&c
  select(Model, R2_test, RMSE_test, MAE_test)%>%
  mutate(
    label = sprintf(
      "R² = %s\nRMSE = %s\nMAE = %s",
      R2_test,
      RMSE_test,
      MAE_test)) 

metrics <- metrics %>% # For Fig.S13b&14b
  filter(Model %in% c("model1", "model1.1","model2","model3",
                      "model3.1","model3.2","model3.3","model3.4","model3.5")) # reduced model

metrics$Model <- as.factor(metrics$Model)

pos_df <- calibration_data_fig.S13 %>%
  group_by(Model) %>%
  summarize(
    x_pos = min(Predicted, na.rm = TRUE), 
    y_pos = max(Actual, na.rm = TRUE))

metrics <- metrics %>%
  left_join(pos_df, by = "Model")

ggplot(calibration_data_fig.S13, aes(x = Predicted, y = Actual)) + 
  geom_pointdensity(size = 0.1, alpha = 0.5) +
  scale_color_viridis_c(option = "viridis", name = "Point Density") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.5) +
  facet_wrap(~ Model, scales = "free", labeller = as_labeller(labels_fig.S13)) + 
  geom_text(
    data = metrics, aes(x = x_pos, y = y_pos, label = label), inherit.aes = FALSE,
    hjust = 0, vjust = 1, size = 2, family = "Times New Roman") +
  labs(
    # title = "Calibration Plots: Predicted vs. Observed BMI",
    x = "Predicted BMI", y = "Observed BMI") +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(face = "bold"),
    strip.text = element_text(size = 8),
    legend.title = element_text(),
    legend.text = element_text()) # Fig.S13b&S14b

## Residual plots for all models  # Fig.S13c&14c
all_residual_data <- calibration_data_fig.S13 %>% 
  mutate(Residual = Actual - Predicted) %>%
  select(Model, Residual, Predicted, Actual)

ggplot(all_residual_data, aes(x = Predicted, y = Residual)) +
  geom_point(alpha = 0.5, shape = 21, color = "gray30", fill = "white", size = 0.01) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  facet_wrap(~ Model, scales = "free", labeller = as_labeller(labels_fig.S13)) + 
  labs(
    #title = "Residual Plots: Residuals vs. Predicted BMI",
    x = "Predicted BMI",y = "Residual (Actual - Predicted)") +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(family = "Times New Roman"),
    axis.title = element_text(family = "Times New Roman"),
    axis.text = element_text(family = "Times New Roman"),
    strip.text = element_text(family = "Times New Roman", size = 8)) # Fig.S13c&14c