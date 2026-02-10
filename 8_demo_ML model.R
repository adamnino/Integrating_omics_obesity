############################################################################################################# 8.
##### 8. ML modelling (Fig.2; Table S10-12; Fig.S11&12)
setwd("~/Desktop")
rm(list=ls()) 
library(dplyr)
library(tidyverse)
library(foreach) # 8.1 - ML models
library(doParallel) # 8.1 - ML models
library(doRNG) # 8.1 - ML models
library(patchwork) # 8.3.1 - Figures
library(ggpattern) # 8.3.1 - Figures
library(ggpointdensity) # 8.3.2 - Figures (Density&calibration)
library(extrafont) # 8.3.2 - Figures (Font)
library(ggcorrplot) # 8.3.4 - Correlation
df <- read.csv("final data_17SNP.csv")
df$S_LDL_P_metabolite <- df$S_LDL_P_metabolite*1000 # changing units

vars_to_log <- c(
  "TriacylglycerolmmolL", "hsCRPmgL", "AdiponectinμgL", "InsulinpmolL",
  "Acetate_metabolite", "XXL_VLDL_C_metabolite","XL_HDL_CE_pct_metabolite",
  "Animalprotein_gram_perday", "Animalfat_gram_perday", "Vegetablefat_gram_perday",
  "Totalcarb_gram_perday", "Fiber_gram_perday", "Calcium_gram_perday")

for (var in vars_to_log) {
  if (all(var %in% names(df)) && is.numeric(df[[var]])) {df[[var]] <- log(df[[var]])} 
  else {warning(paste("Variable", var, "is missing or not numeric"))}}

df$PA_all <- log(df$PA_all+1) # log (x+1) is used for PA because PA contains true value "0"; log(0+1)=0
############################################################################################################# 8.1
################################ 8.1 ML models #################################
df <- df %>% 
  mutate(across(c(2,3,6,8,9,11,13,15,30,32,33), as.factor)) %>% 
  mutate(across(c(4,5,7,16:29,31,34:56), as.numeric)) %>% 
  mutate(across(c(1,10,12,14), as.character))

## Define all models - Fig.2; Table S10; Fig.S11
model_list <- list(
  model0 = df[, c(2:5)],# model 0 = base
  model1 = df[, c(2:5,16)], # model 1 = WC
  model2 = df[, c(2:5,7,16:33,56)], # model 2 = clinical
  model3 = df[, c(2:5,7,17:33,56)],  # model 3 = clinical (no WC)
  model4 = df[, c(2:5,34:47)],  # model 4 = metabolites
  model5 = df[, c(2:5,48:54)],  # model 5 = diet
  model6 = df[, c(2:5,55)],  # model 6 = SNPs
  
  model7 = df[, c(2:5,7,17:47,56)], # model 7 = clinical + metabolites (no WC)
  model8 = df[, c(2:5,7,17:33,48:54,56)],  # model 8 = clinical + diet (no WC)
  model9 = df[, c(2:5,7,17:33,55,56)],  # model 9 = clinical + SNP (no WC)
  model10 = df[, c(2:5,7,17:56)],  # model 10 = clinical + metabolites + SNP + diet (no WC)
  
  model11 = df[, c(2:5,7,16:47,56)], # model 11 = clinical + metabolites
  model12 = df[, c(2:5,7,16:33,48:54,56)],  # model 12 = clinical + diet
  model13 = df[, c(2:5,7,16:33,55,56)],  # model 13 = clinical + SNP
  model14 = df[, c(2:5,7,16:56)]  # model 14 = clinical + metabolites + SNP + diet 
)

set.seed(1053) # repeated 10 times * 5 outer folds * 3 inner folds
repeats_outer <- 10
cores <- parallel::detectCores() - 1 # use all cores except one (7 cores)
registerDoParallel(cores)
registerDoRNG(seed = 1234)

## Reproducible bootstrap CI function
boot_ci <- function(metric_values, seed = 1000) {
  set.seed(seed)
  boot_out <- boot(data = metric_values, statistic = function(data, i) mean(data[i]), R = 1000)
  ci <- boot.ci(boot_out, type = "perc")
  return(c(mean = mean(boot_out$t0), lower = ci$percent[4], upper = ci$percent[5]))}

## Run models with reproducibility (Seed & RNG)
results_list <- foreach(model_name = names(model_list), 
                        .packages = c("caret", "boot", "NeuralNetTools"), 
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
      
      tuned_model <- train(
        BMI ~ ., data = train_data,
        method = "glmnet", # regularized linear regression (GLMNET)
        
        # method = "xgbTree", # extreme gradient boosting (XGBoost)

        # method = "ranger", # random forest (RF)
        # importance = 'permutation', # RF model
        
        # method = "nnet", # neural network (NN)
        # linout = TRUE, # NN model
        # trace = FALSE, # NN model
        # maxit = 500, # NN model
        preProcess = c("center", "scale"), # NN & GLMNET models
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
      # importance_vals <- varImp(tuned_model, scale = FALSE)$importance # RF & XGBoost models

      # importance_vals <- olden(tuned_model$finalModel, bar_plot = FALSE) # NN model using olden method

      importance_vals <- as.matrix(coef(tuned_model$finalModel, s = tuned_model$bestTune$lambda)) # GLMNET model
      importance_vals <- importance_vals[rownames(importance_vals) != "(Intercept)", , drop = FALSE] # GLMNET model
      
      importance_df <- data.frame(
        Variable = rownames(importance_vals),
        # Importance = importance_vals$Overall, # RF & XGBoost models
        # Importance = importance_vals$importance, # NN model
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

#### NOTES: GLMNET(4 minutes); RF(65 minutes); XGBoost(29 minutes); NN(80 minutes).

############################################################################################################# 8.2.1
######### 8.2.1 Data extractions - Evaluation metrics 95% CI - Table S10 ########

#### NOTE: For GLMNET model, 14 NAs (7NAs for train, 7NAs for test) were observed 
### (8 NAs in model0, 2 NAs in model5, and 4 NAs in model6) in different folds and
### repeats. This was expected as R2 were extremely small in those models. 

### For all other models, run Steps 4&6 only. For GLMNET model, run Steps 1-6.
########################## Step 1 - Extract all metrics ########################
all_metrics <- bind_rows(lapply(names(results_list), function(model) {
  metrs <- results_list[[model]]$All_Metrics
  metrs$Model <- model
  metrs}))
view(all_metrics) # Check which model contain NAs
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
# all_avg[1, 3] <- 0.02717022 # GLMNET  model 0, test
# all_avg[4, 3] <- 0.05331103 # GLMNET  model 0, train
# 
# all_avg[31, 3] <- 0.03101666 # GLMNET model 5, test
# all_avg[34, 3] <- 0.07228506 # GLMNET model 5, train
# 
# all_avg[37, 3] <- 0.02907036 # GLMNET model 6, test
# all_avg[40, 3] <- 0.06628059 # GLMNET model 6, train

######################  Step 6 - Formatting the results #######################
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
view(df_final) # Table S10

############################################################################################################# 8.2.2
################### 8.2.2 Data extractions - hyperparameters ###################
## Check if R2 (test set) are consistent across different hyperparameter combinations
all_hyperparams <- bind_rows(lapply(names(results_list), function(model) {
  hpams <- results_list[[model]]$Hyperparams
  hpams$Model <- model
  hpams}))

labels_all_model <- c(
  model0 = "Baseline*",
  model1 = "WC",
  model2 = "Clinical (with WC)",
  model3 = "Clinical (no WC)",
  model4 = "Metabolites",
  model5 = "Diet",
  model6 = "GRS",
  model7 = "Clinical (no WC) + Metabolites",
  model8 = "Clinical (no WC) + Diet",
  model9 = "Clinical (no WC) + GRS",
  model10 = "Clinical (no WC) + Metabolites + Diet + GRS",
  model11 = "Clinical (with WC) + Metabolites",
  model12 = "Clinical (with WC) + Diet",
  model13 = "Clinical (with WC) + GRS",
  model14 = "Clinical (with WC) + Metabolites + Diet + GRS"
) # Fig. 1, Fig.S11

all_hyper_metric <- all_hyperparams %>%
  mutate(
    R2 = all_metrics$R2[all_metrics$Set == "Test"],
    Model = factor(Model, levels = names(labels_all_model))) 

ggplot(all_hyper_metric, aes(x = lambda, y = alpha, color = R2)) +
  geom_point() +
  facet_wrap(~ Model, labeller = as_labeller(labels_all_model)) + 
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_minimal() +
  labs(
    #title = "Hyperparameters Combinations For Each Model Colored By R2",
    x = "Lambda (log10 scale)", y = "Alpha") +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(size = 6)) # Fig.S11a

############################################################################################################# 8.2.3
############### 8.2.3 Data extractions - importance/coefficients ###############
all_importance <- bind_rows(lapply(names(results_list), function(model) {
  iprte <- results_list[[model]]$Importance
  iprte$Model <- model
  iprte}))

## Mean importance for each variable-model pair
avg_importance <- all_importance %>%
  group_by(Model, Variable) %>%
  summarise(AvgImportance = mean(Importance, na.rm = TRUE), .groups = "drop")

############################################################################################################# 8.2.4
#### 8.2.4 Data extractions - p values b/t models on test sets - Table S10 #####
df_test <- all_metrics %>%
  filter(Set == "Test") %>% 
  select(R2, Repeat, Fold, Model) %>% # change R2/RMSE/MAE
  pivot_wider(names_from = Model, values_from = R2) %>% # change R2/RMSE/MAE
  mutate(id = paste0("Rep", Repeat, "_Fold", Fold))

long_df_T1_3 <- df_test %>% 
  pivot_longer(cols = c("model0","model1","model2","model3","model4","model5",
                        "model6","model7","model8","model9","model10","model11",
                        "model12","model13","model14"),
               names_to = "Model", values_to = "R2") # change R2/RMSE/MAE

## Pairwise Wilcoxon signed-rank test
p_value <- pairwise.wilcox.test(long_df_T1_3$R2,  # change R2/RMSE/MAE
                                long_df_T1_3$Model, paired = TRUE, p.adjust.method = "none")  # Table S10
view(p_value$p.value)

############################################################################################################# 8.3.1
######################## 8.3.1 Bar plots Fig. 2 - GLM ##########################
labels_fig.2 <- c(
  model0 = "Baseline*",
  model1 = "WC",
  model2 = "Clinical (with WC)",
  model3 = "Clinical (no WC)",
  model4 = "Metabolites",
  model5 = "Diet",
  model6 = "GRS",
  model7 = "Metabolites",
  model8 = "Diet",
  model9 = "GRS",
  model10 = "Metabolites + Diet + GRS",
  model11 = "Metabolites",
  model12 = "Diet",
  model13 = "GRS",
  model14 = "Metabolites + Diet + GRS") # For Fig.2
 
model_fig.2 <- data.frame(
  Model_Name = names(labels_fig.2),
  Display_Name = unname(labels_fig.2),
  Group = c(
    rep("Base model", 1), # model0
    rep("Single model", 6), # model1, model2, model3, model4, model5, model6
    rep("Clinical (no WC) +", 4), # model7, model8, model9, model10
    rep("Clinical (with WC) +", 4)) # model11, model12, model13, model14
) # For Fig.2

plot_data <- all_avg %>%
  filter(Set == "Test") %>%
  mutate(Model = factor(Model, levels = names(labels_fig.2))) %>%
  left_join(model_fig.2, by = c("Model" = "Model_Name")) %>%
  mutate(Group = factor(Group, levels = unique(model_fig.2$Group)))

## Plotting function - Fig. 2
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
          panel.spacing.x = unit(0, "pt"))} # Fig. 2

## Combine plots
p1 <- plot_metric("R2", expression(R^2))
p1 # Fig. 2
#p2 <- plot_metric("RMSE", "RMSE") # optional figure
#p3 <- plot_metric("MAE", "MAE") # optional figure

############################################################################################################# 8.3.2
################# 8.3.2 Calibration plots - GLM - Fig.S11b&c ###################
## Density & calibration plots
all_calibration_data <- lapply(names(results_list), function(model) {
  df <- results_list[[model]]$Calibration
  df$Model <- model
  df}) %>% 
  bind_rows()

calibration_data_fig.S11 <- all_calibration_data %>%
  mutate(Model = forcats::fct_relevel(Model, 'model0', 'model1', 'model2','model3'
                                      ,'model4','model5','model6','model7','model8',
                                      'model9','model10','model11','model12','model13',
                                      'model14')) # Fig.S11b&c

metrics <- df_final%>% # df_final is in 8.3.1 - Step 6; for Fig.S11b
  select(Model, R2_test, RMSE_test, MAE_test)%>%
  mutate(
    label = sprintf(
      "R² = %s\nRMSE = %s\nMAE = %s",
      R2_test,
      RMSE_test,
      MAE_test)) 

metrics$Model <- as.factor(metrics$Model)

pos_df <- calibration_data_fig.S11 %>%
  group_by(Model) %>%
  summarize(
    x_pos = min(Predicted, na.rm = TRUE), 
    y_pos = max(Actual, na.rm = TRUE))

metrics <- metrics %>%
  left_join(pos_df, by = "Model")

ggplot(calibration_data_fig.S11, aes(x = Predicted, y = Actual)) +
  geom_pointdensity(size = 0.1, alpha = 0.5) +
  scale_color_viridis_c(option = "viridis", name = "Point Density") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.5) +
  facet_wrap(~ Model, scales = "free", labeller = as_labeller(labels_all_model)) +
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
    legend.text = element_text()) # Fig.S11b

## Residual plots for all models  # Fig.S11c
all_residual_data <- calibration_data_fig.S11 %>% 
  mutate(Residual = Actual - Predicted) %>%
  select(Model, Residual, Predicted, Actual)

ggplot(all_residual_data, aes(x = Predicted, y = Residual)) +
  geom_point(alpha = 0.5, shape = 21, color = "gray30", fill = "white", size = 0.01) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  facet_wrap(~ Model, scales = "free", labeller = as_labeller(labels_all_model)) + 
  labs(
    # title = "Residual Plots: Residuals vs. Predicted BMI",
    x = "Predicted BMI", y = "Residual (Actual - Predicted)") +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(family = "Times New Roman"),
    axis.title = element_text(family = "Times New Roman"),
    axis.text = element_text(family = "Times New Roman"),
    strip.text = element_text(family = "Times New Roman", size = 8)) # Fig.S11c

############################################################################################################# 8.3.3
##### 8.3.3 Figures for the best model - importance - Fig.S12a; Table S12 #######
## Plot all variable importance (mean ± SD) for each model
plot_var_importance <- function(model) {
  df_summary <- all_importance %>%
    filter(Model == model) %>%
    group_by(Variable) %>%
    summarize(
      AvgImportance = mean(Importance),
      SDImportance = sd(Importance),
      .groups = "drop") %>%
    arrange(desc(AvgImportance))
  
  model_title <- labels_all_model[[model]] %||% model
  
  ggplot(df_summary, aes(x = reorder(Variable, AvgImportance), y = AvgImportance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_errorbar(
      aes(ymin = AvgImportance - SDImportance, ymax = AvgImportance + SDImportance),
      width = 0.3, color = "black") +
    coord_flip() +
    labs(
      title = paste("Variable Coefficient -", model_title),
      x = "Variable", y = "Mean Coefficient ± SD") +
    theme_minimal()}

invisible(lapply(unique(all_importance$Model), function(model) {
  print(plot_var_importance(model))})) # Fig.S12a

## Find standardized beta coefficients
find_var_coef <- function(variable_names, model_name = NULL) {
  if (is.null(model_name)) {
    # check across all models for multiple variables
    result <- all_importance %>%
      filter(Variable %in% variable_names) %>%
      group_by(Model, Variable) %>%
      summarize(
        Mean_Coefficient = mean(Importance),
        SD = sd(Importance),
        .groups = "drop"
      ) %>%
      arrange(Variable, desc(abs(Mean_Coefficient)))
  } else {
    # check in specific model for multiple variables
    result <- all_importance %>%
      filter(Model == model_name, Variable %in% variable_names) %>%
      group_by(Variable) %>%
      summarize(
        Mean_Coefficient = mean(Importance),
        SD = sd(Importance),
        N = n(),
        .groups = "drop"
      ) %>%
      arrange(desc(abs(Mean_Coefficient)))
  }
  
  return(result)
}

find_var_coef(c("InsulinpmolL", "SBP"), "model3") # Table S12
find_var_coef(c("Tyr_metabolite", "Phe_metabolite","L_HDL_CE_metabolite",
                "S_LDL_P_metabolite","GlycA_metabolite","Albumin_metabolite"), "model4")  # Table S12
############################################################################################################# 8.3.4
################# 8.3.4 Correlation plot - Fig.S12b; Table S12 ##################
cor_matrix <- df %>% 
  select(BMI,InsulinpmolL,SBP,Albumin_metabolite,Tyr_metabolite,L_HDL_CE_metabolite,
         Phe_metabolite,S_LDL_P_metabolite,GlycA_metabolite)

cor_matrix <- cor(cor_matrix, use = "complete.obs", method = "pearson")

ggcorrplot(cor_matrix, lab = TRUE, type = "upper", colors = c("darkblue", "white", "orange")) 

variables <- c("InsulinpmolL", "SBP", "Albumin_metabolite", "Tyr_metabolite", 
               "L_HDL_CE_metabolite", "Phe_metabolite", "S_LDL_P_metabolite", "GlycA_metabolite")

results <- data.frame(
  Variable = variables,
  Correlation = numeric(length(variables)),
  P_value = numeric(length(variables)),
  N = numeric(length(variables)))

## Calculate correlations with BMI
for (i in seq_along(variables)) {
  test_result <- cor.test(df$BMI, df[[variables[i]]], 
                          method = "pearson", 
                          use = "complete.obs")
  
  results$Correlation[i] <- test_result$estimate
  results$P_value[i] <- test_result$p.value
  results$N[i] <- sum(complete.cases(df$BMI, df[[variables[i]]]))}
view(results) # Table S12