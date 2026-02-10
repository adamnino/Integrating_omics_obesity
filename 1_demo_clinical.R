############################################################################################################# 1.
##### 1. Clinical - Imputation (Table S2&3; Fig.S3&4)
setwd("~/Desktop")
rm(list=ls()) 
library(tidyverse)
library(dplyr)
library(ggplot2) 
library(reshape2)
library(mice) # 1.3 - MICE imputation
library(broom) # 1.5 - pool the estimates and extract results
cli <- read.csv("clinical_raw.csv")
cli <- cli %>% 
  mutate(across(c(2,3,5,8,18,20,22:28,30,32,34,36,38,40,41,57,62:71), as.factor)) %>% 
  mutate(across(c(4,42:56,59,), as.numeric)) %>% 
  mutate(across(c(1,6,7,9:17,19,21,29,31,33,35,37,39,58,60,61,), as.character))

############################################################################################################# 1.1
################## 1.1 Investigating missingness - Table S2 ####################
missing_counts <- colSums(is.na(cli))
missing_counts[missing_counts > 0] # values that were missing

## Count complete vs. incomplete cases (# of Participants containing missing values)
tibble(Status = c("Incomplete", "Complete"),
       Count = as.numeric(table(complete.cases(cli))),
       Percent = round(100 * table(complete.cases(cli)) / nrow(cli), 2)) # 16.13% incomplete. Thus, use 17 imputed data sets during MICE imputation.

## Number and proportion of missing values per variable
cbind("# NA" = sort(colSums(is.na(cli))),
      "% NA" = round(sort(colMeans(is.na(cli))) * 100, 2)) # Table S2

## Overall missingness
cat("\nOverall missingness:", round(100 * sum(is.na(cli)) / prod(dim(cli)), 2), "%\n\n") # 1.08%
############################################################################################################# 1.2
######################### 1.2 Distribution of the data #########################
plot_fig <- function(start_col, end_col, nr, nc, mgp_vals, mar_vals) {
  par(mfrow = c(nr, nc), mgp = mgp_vals, mar = mar_vals)
  
  for (i in start_col:end_col) {
    if (is.numeric(cli[, i])) {
      hist(cli[, i], nclass = 30, xlab = "",
           main = paste0(names(cli[i]), " (", 
                         round(mean(is.na(cli[, i])) * 100, 2), "% NA)"),
           cex.main = 0.7, cex.axis = 0.6, cex.lab = 0.8)
    } else {
      barplot(table(cli[, i]), ylab = "Frequency",
              main = paste0(names(cli[i]), " (", 
                            round(mean(is.na(cli[, i])) * 100, 2), "% NA)"),
              cex.main = 0.7, cex.axis = 0.6, cex.lab = 0.8)
    }
  }
} # function to plot distribution

nc1 <- nc2 <- 6 # columns for first & second figures
nr1 <- nr2 <- 5 # rows for first & second figures
split_idx <- ceiling(ncol(cli) / 2) # split index

plot_fig(2, split_idx, 5, 6, c(1.5, 0.5, 0), c(2, 3, 3, 0.5)) # first 27 variables
plot_fig(split_idx + 1, ncol(cli), nr2, nc2, c(1, 0.5, 0), c(2, 3, 3, 0.5)) # another 28 variables

############################################################################################################# 1.3
### NOTES: There are no additional observations that can be included in the 
### imputation stage to improve model stability. Most excluded observations lack 
### clinical data, with only a few (< 5) containing clinical data. Thus, 
### including them is unnecessary for model stability during the imputation stage.
############################# 1.3 Predictor matrix #############################
## Make the predictor matrix. This matrix indicates which variables will be used 
## as predictors for missing values
imp0 <- mice(cli, maxit = 0, 
             defaultMethod = c("norm", "logreg", "polyreg", "polr")) # setup-run

meth <- imp0$meth
cli.mat <- imp0$predictorMatrix # adjust predictor matrix

### Logreg doesn't  work with variables with only 1 NA (the first 5 variables);
### and adiponectin has slightly skewed distribution. Smoking.statues also has 
### skewed distribution. Additionally, most variables are categorical. Thus, 
### use rf. Except "Polycystic_Ovary" because it gives false prediction.
vars_rf <- c("Diabetes_other_endocrine","Renal_disease","Cancer",
             "Thrombosis","AdiponectinμgL","GlucosemmolL","hsCRPmgL",
             "Heart_stroke_vascular","Asthma","High_cholesterol",
             "Hypertension","Thyroid_disorder","Food_allergy","WC","HbA1c",
             "Smoking_YN","Diabetes_mother","Hypertension_mother",
             "Obesity_mother","Heartattack_mother","Stroke_mother","Diabetes_father",
             "Heartattack_father","Obesity_father","Hypertension_father","Stroke__father")
meth[vars_rf] <- "rf" # use random forest for these variables 
meth

############################################################################################################# 1.4
############### 1.4 Imputation, Convergence, Diagnostics - Fig.S3 ##############
## MICE imputation
imp <- mice(cli, method = meth, predictorMatrix = cli.mat, 
            m = 17, maxit = 30,
            seed = 1730)

## Convergence plots
imputed_vars <- names(meth[meth != ""]) # get only variables that were imputed

plot_imp <- function(data, y_label) {
  data <- melt(data)
  data <- data[data$Var1 %in% imputed_vars, ]  # filter imputed variables
  
  ggplot(data, aes(x = Var2, y = value, color = Var3)) +
    geom_line() +
    facet_wrap("Var1", scales = 'free') +
    theme(legend.position = 'none') +
    xlab("iteration") +
    ylab(y_label)
} # function to plot imputed data

plot_imp(imp$chainMean, "imputed value") # plot mean of imputed values; Fig.S3a.A
plot_imp(imp$chainVar, "SD of imputed value") # plot SD of imputed values; Fig.S3a.B

## Diagnostics for continuous & categorical variables
propplot <- function(x, formula, facet = "wrap", ...) {
  library(ggplot2)
  
  cd <- data.frame(mice::complete(x, "long", include = TRUE))
  cd$.imp <- factor(cd$.imp)
  
  r <- as.data.frame(is.na(x$data))
  
  impcat <- x$meth != "" & sapply(x$data, is.factor)
  vnames <- names(impcat)[impcat]
  
  if (missing(formula)) {
    formula <- as.formula(paste(paste(vnames, collapse = "+",
                                      sep = ""), "~1", sep = ""))
  }
  
  tmsx <- terms(formula[-3], data = x$data)
  xnames <- attr(tmsx, "term.labels")
  xnames <- xnames[xnames %in% vnames]
  
  if (paste(formula[3]) != "1") {
    wvars <- gsub("[[:space:]]*\\|[[:print:]]*", "", paste(formula)[3])
    # wvars <- all.vars(as.formula(paste("~", wvars)))
    wvars <- attr(terms(as.formula(paste("~", wvars))), "term.labels")
    if (grepl("\\|", formula[3])) {
      svars <- gsub("[[:print:]]*\\|[[:space:]]*", "", paste(formula)[3])
      svars <- all.vars(as.formula(paste("~", svars)))
    } else {
      svars <- ".imp"
    }
  } else {
    wvars <- NULL
    svars <- ".imp"
  }
  
  for (i in seq_along(xnames)) {
    xvar <- xnames[i]
    select <- cd$.imp != 0 & !r[, xvar]
    cd[select, xvar] <- NA
  }
  
  
  for (i in which(!wvars %in% names(cd))) {
    cd[, wvars[i]] <- with(cd, eval(parse(text = wvars[i])))
  }
  
  meltDF <- reshape2::melt(cd[, c(wvars, svars, xnames)], id.vars = c(wvars, svars))
  meltDF <- meltDF[!is.na(meltDF$value), ]
  
  
  wvars <- if (!is.null(wvars)) paste0("`", wvars, "`")
  
  a <- plyr::ddply(meltDF, c(wvars, svars, "variable", "value"), plyr::summarize,
                   count = length(value))
  b <- plyr::ddply(meltDF, c(wvars, svars, "variable"), plyr::summarize,
                   tot = length(value))
  mdf <- merge(a,b)
  mdf$prop <- mdf$count / mdf$tot
  
  plotDF <- merge(unique(meltDF), mdf)
  plotDF$value <- factor(plotDF$value,
                         levels = unique(unlist(lapply(x$data[, xnames], levels))),
                         ordered = T)
  
  p <- ggplot(plotDF, aes(x = value, fill = get(svars), y = prop)) +
    geom_bar(position = "dodge", stat = "identity") +
    theme(legend.position = "bottom", ...) +
    ylab("proportion") +
    scale_fill_manual(name = "",
                      values = c("darkblue",
                                 colorRampPalette(
                                   RColorBrewer::brewer.pal(9, "Oranges"))(x$m + 3)[1:x$m + 3])) +
    guides(fill = guide_legend(nrow = 1))
  
  if (facet == "wrap")
    if (length(xnames) > 1) {
      print(p + facet_wrap(c("variable", wvars), scales = "free"))
    } else {
      if (is.null(wvars)) {
        print(p)
      } else {
        print(p + facet_wrap(wvars, scales = "free"))
      }
    }
  
  if (facet == "grid")
    if (!is.null(wvars)) {
      print(p + facet_grid(paste(paste(wvars, collapse = "+"), "~ variable"),
                           scales = "free"))
    }
} # function similar to density plot

propplot(imp) # for categorical variables; Fig.S3b.A

densityplot(imp,~WC+HbA1c+GlucosemmolL,
            col = c("darkblue", "orange")) # for continuous variables; Fig.S3b.B

stripplot(imp, AdiponectinμgL,
          pch = c(1, 19), 
          col = c("darkblue", "orange"), 
          cex = 0.8) # Fig.S3b.C

stripplot(imp, hsCRPmgL,
          pch = c(1, 19), 
          col = c("darkblue", "orange"), 
          cex = 0.8) # Fig.S3b.D

############################################################################################################# 1.5
########################### 1.5 Analysis - Table S3 ############################
## Linear regression of imputed datasets
models <- with(imp, lm(BMI~Sitename + Sex + Age + Food_supplements + Medication + 
                       Illness_current + Disease_Other + Heart_stroke_vascular + 
                       Arthritis + Asthma + Renal_disease + Cancer + Malaria + Thrombosis +
                       Polycystic_Ovary + High_cholesterol + Hypertension + Thyroid_disorder +
                       Diabetes_other_endocrine + Food_allergy + Menstrual_cycle_regular +
                       WC + SBP + DBP + Pulse + Heamoglobin_gdl + GlucosemmolL +
                       CholesterolmmolL + TriacylglycerolmmolL + HDLCmmolL + LDLCmmolL +
                       hsCRPmgL + HbA1c + AdiponectinμgL + InsulinpmolL + Smoking_YN +
                       PA_all + Diabetes_mother + Hypertension_mother +
                       Obesity_mother + Heartattack_mother + Stroke_mother + Diabetes_father +
                       Hypertension_father + Obesity_father + Heartattack_father +
                       Stroke__father))

summary(models$analyses[[1]])

pooled_df <- tidy(pool(models)) # pool the estimates and extract results

## Calculating 95% CI based on Rubin's rules 
pooled_df$lower <- pooled_df$estimate - qt(0.975, pooled_df$df) * pooled_df$std.error
pooled_df$upper <- pooled_df$estimate + qt(0.975, pooled_df$df) * pooled_df$std.error
#write.csv(pooled_df,"~/Desktop/mice_MICE1.csv", row.names=FALSE) # Table S3 (MICE1)

############################################################################################################# 1.6.1
##################### 1.6.1 Sensitivity analysis - Table S3 ####################
### Sensitivity analysis 1 - complete case by removing observations
## Fit the linear model for complete cases only (Complete case by default)
lm_cc <- lm(BMI~Sitename + Sex + Age + Food_supplements + Medication + 
              Illness_current + Disease_Other + Heart_stroke_vascular + 
              Arthritis + Asthma + Renal_disease + Cancer + Malaria + Thrombosis +
              Polycystic_Ovary + High_cholesterol + Hypertension + Thyroid_disorder +
              Diabetes_other_endocrine + Food_allergy + Menstrual_cycle_regular +
              WC + SBP + DBP + Pulse + Heamoglobin_gdl + GlucosemmolL +
              CholesterolmmolL + TriacylglycerolmmolL + HDLCmmolL + LDLCmmolL +
              hsCRPmgL + HbA1c + AdiponectinμgL + InsulinpmolL + Smoking_YN +
              PA_all + Diabetes_mother + Hypertension_mother +
              Obesity_mother + Heartattack_mother + Stroke_mother + Diabetes_father +
              Hypertension_father + Obesity_father + Heartattack_father +
              Stroke__father, data = cli)

CI <- as.data.frame(confint(lm_cc)) # calculating 95% CI
CI <- CI[-c(27,34),] # remove PCOS and menstrual as they were for female only

cc_results <- data.frame(summary(lm_cc)$coefficients, CI) # remove 106 obs out of 657 containing NAs
#write.csv(cc_results,"~/Desktop/mice_CC1.csv", row.names=FALSE) # Table S3 (Complete case 1)

############################################################################################################# 1.6.2
##################### 1.6.2 Sensitivity analysis - Table S3 ####################
### Sensitivity analysis 2 - complete case by removing medical history of father
## Remove medical history of father as a predictor since it contains a lot of NAs
cli_nf <- cli[,c(1:60,62:66)] # nf = no father
100 * nic(cli_nf) / nrow(cli_nf) # incomplete case is now 7.61%

## Repeat 1.3 to 1.6
imp0_nf <- mice(cli_nf, maxit = 0) # set-up run
imp0_nf <- mice(cli_nf, maxit = 0, 
                defaultMethod = c("norm", "logreg", "polyreg", "polr"))

cli.mat_nf <- imp0_nf$cliictorMatrix # adjust predictor matrix
meth_nf <- imp0_nf$meth
meth_nf 

## Remove variables names for medical hisotry of father
vars_rf_nf <- vars_rf[-((length(vars_rf) - 4):length(vars_rf))] 

meth_nf[vars_rf_nf] <- "rf" # set the rest with 'rf' method
meth_nf # check meth_nf, ensure it's all 'rf' except PCOS

## Run imputation; m = 17 is used to be consistent with imputation model 1
imp_nf <- mice(cli_nf, method = meth_nf, cliictorMatrix = cli.mat_nf, 
               m = 17, maxit = 30,
               seed = 1730) 

## Linear regression of imputed datasets
models_nf <- with(imp_nf, lm(BMI~Sitename + Sex + Age + Food_supplements + Medication + 
                             Illness_current + Disease_Other + Heart_stroke_vascular + 
                             Arthritis + Asthma + Renal_disease + Cancer + Malaria + Thrombosis + 
                             Polycystic_Ovary + High_cholesterol + Hypertension + Thyroid_disorder + 
                             Diabetes_other_endocrine + Food_allergy + Menstrual_cycle_regular + 
                             WC + SBP + DBP + Pulse + Heamoglobin_gdl + GlucosemmolL + 
                             CholesterolmmolL + TriacylglycerolmmolL + HDLCmmolL + LDLCmmolL + 
                             hsCRPmgL + HbA1c + AdiponectinμgL + InsulinpmolL + Smoking_YN + 
                             PA_all + Diabetes_mother + Hypertension_mother + 
                             Obesity_mother + Heartattack_mother + Stroke_mother))

summary(models_nf$analyses[[1]])

pooled_df_nf <- tidy(pool(models_nf)) # pool the estimates and extract results

## Calculating 95 CI for imputed datasets based on Rubin's rules
pooled_df_nf$lower <- pooled_df_nf$estimate - qt(0.975, pooled_df_nf$df) * pooled_df_nf$std.error
pooled_df_nf$upper <- pooled_df_nf$estimate + qt(0.975, pooled_df_nf$df) * pooled_df_nf$std.error
#write.csv(pooled_df_nf,"~/Desktop/mice_MICE2.csv", row.names=FALSE) # Table S3 (MICE 2)

lm_cc_nf <- lm(BMI~Sitename + Sex + Age + Food_supplements + Medication + 
                 Illness_current + Disease_Other + Heart_stroke_vascular + 
                 Arthritis + Asthma + Renal_disease + Cancer + Malaria + Thrombosis +
                 Polycystic_Ovary + High_cholesterol + Hypertension + Thyroid_disorder +
                 Diabetes_other_endocrine + Food_allergy + Menstrual_cycle_regular +
                 WC + SBP + DBP + Pulse + Heamoglobin_gdl + GlucosemmolL +
                 CholesterolmmolL + TriacylglycerolmmolL + HDLCmmolL + LDLCmmolL +
                 hsCRPmgL + HbA1c + AdiponectinμgL + InsulinpmolL + Smoking_YN +
                 PA_all + Diabetes_mother + Hypertension_mother +
                 Obesity_mother + Heartattack_mother + Stroke_mother, data = cli_nf) # complete case analysis 2
summary(lm_cc_nf)

CI_nf <- as.data.frame(confint(lm_cc_nf)) # calculating 95% CI for complete case
CI_nf <- CI_nf[-c(27,34),] # remove PCOS and menstrual as they were for female only
cc_results_nf <- data.frame(summary(lm_cc_nf)$coefficients, CI_nf) 
#write.csv(cc_results_nf,"~/Desktop/mice_CC2.csv", row.names=FALSE) # Table S3 (CC 2)

############################################################################################################# 1.6.3
################ 1.6.3 Forest plots for estimates and CI - Fig.S4 ##############
## Coefficients and 95% CI of MICE1, CC1, MICE2, and CC2 are combined in excel to 
## create "forest plot.csv"
df_plot <- read.csv("Forest plot.csv")

order <- rev(names(df_plot)[!names(df_plot)%in% c("term","Method")])
df_long <- df_plot %>%
  pivot_longer(cols = -c(term, Method), names_to = "Variable", values_to = "Value") %>%
  pivot_wider(names_from = term, values_from = Value) %>%
  na.omit() %>%
  mutate(Variable = factor(Variable, levels = order))

## Define general method levels
met.levls <- c("CC 1", "CC 2", "MICE 1", "MICE 2")

## Define different color mappings for different comparisons
color_maps <- list(
  "Overall comparison" = c("Com 1" = "#ff7f0e", "CC 2" = "#0096FF", "MICE 1" = "black", "MICE 2" = "#33a02c"),
  "MICE imputation 1 vs. Complete case 1" = c("MICE 1" = "black", "CC 1" = "#ff7f0e"),
  "MICE imputation 1 vs. Complete case 2" = c("MICE 1" = "black", "CC 2" = "#0096FF"),
  "MICE imputation 1 vs. MICE imputation 2" = c("MICE 1" = "black", "MICE 2" = "#33a02c")
)

shapes <- c("CC 1" = 15, "CC 2" = 16, "MICE 1" = 17, "MICE 2" = 18) 

## Function to create forest plots with different color mappings
forest_plot <- function(data, title, colors) {
  ggplot(data, aes(x = estimate, y = Variable, color = Method, shape = Method)) +
    geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, position = position_dodge(width = 0.6)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    scale_color_manual(values = colors) + 
    scale_shape_manual(values = shapes) +
    theme_minimal() +
    labs(title = title, x = "Estimate (95% CI)", y = "Variables") +
    theme(legend.position = "bottom")
}

## Groups for filtering
groups <- list(
  "Overall comparison" = met.levls,
  "MICE imputation 1 vs. Complete case 1" = c("MICE 1", "CC 1"),
  "MICE imputation 1 vs. Complete case 2" = c("MICE 1", "CC 2"),
  "MICE imputation 1 vs. MICE imputation 2" = c("MICE 1", "MICE 2")
)

plots <- lapply(names(groups), function(title) {
  forest_plot(df_long %>% filter(Method %in% groups[[title]]), title, color_maps[[title]])
})
plots # Fig.S4

############################################################################################################# 1.7
##################### 1.7 Extract averaged imputed datasets ####################
### Imputed values have similar estimates and 95CI; thus, choose MICE model 1
imp_df <- lapply(1:17, function(i) complete(imp, i))

## Averaging numeric imputed variables
num_vars <- names(Filter(is.numeric, imp_df[[1]][imputed_vars]))  # select only numeric imputed vars
avg_num <- Reduce("+", lapply(imp_df, function(df) df[, num_vars, drop=FALSE])) / length(imp_df)

## Categorical imputed variables
cat_vars <- setdiff(imputed_vars, num_vars)

## Convert factor to numeric for imputed categorical variables
cli[, cat_vars] <- lapply(cli[, cat_vars], function(x) as.numeric(as.character(x)))

## Compute mode for categorical variables, ignoring rows with "2"
cat_imp <- as.data.frame(lapply(cat_vars, function(var) {
  imp_values <- apply(sapply(imp_df, `[[`, var), 2, as.numeric)
  imp_values[imp_values == 2] <- NA  
  as.integer(apply(imp_values, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA)
    as.integer(names(sort(table(x), decreasing=TRUE)[1]))  
  }))
}))

names(cat_imp) <- cat_vars

cat_imp$Polycystic_Ovary[is.na(cat_imp$Polycystic_Ovary)] <- 2

## Combine numeric and categorical imputed data, convert categorical back to factor
final_imp_data <- cbind(ID = cli$ID, avg_num, lapply(cat_imp, function(x) factor(as.character(x))))

## Add order column and order
cli_imp <- cli %>%
  mutate(order = 1:n()) %>%
  arrange(order) %>%
  select(-order) # remove the order column after ordering

final_imp_data <- final_imp_data %>%
  mutate(order = 1:n()) %>%
  arrange(order) %>%
  select(-order)

## Create extracted_df with the ordered data 
extracted_df <- cli_imp
extracted_df[, names(final_imp_data)] <- final_imp_data

#write.csv(extracted_df, "~/Desktop/clinical_imp.csv", row.names = FALSE)
### NOTES: scores for medical history, OMH, OE, MH_parents need to be calculated
### in Excel before combining other models. Named as "clinical_clean.csv".