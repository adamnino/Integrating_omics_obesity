############################################################################################################# 7.
##### 7. Classic linear model (Table S9) & Study demographic characteristics (Table 1&S8)
setwd("~/Desktop")
rm(list=ls()) 
library(dplyr)
library(tidyverse)
library(broom) # 7.3 - classic linear model (tidy coefficient, SE, p-value)
df <- read.csv("final data_17SNP.csv")
df$S_LDL_P_metabolite <- df$S_LDL_P_metabolite*1000 # mmol/L --> umol/L
df <- df %>% 
  mutate(across(c(2,3,6,8,9,11,13,15,30,32,33), as.factor)) %>% 
  mutate(across(c(4,5,7,16:29,31,34:56), as.numeric)) %>% 
  mutate(across(c(1,10,12,14), as.character))
############################################################################################################# 7.1
############# 7.1 Basic demographic characteristics - Table 1&S8 ###############
basic <- df[,-c(1,8:15)]

## Function to summarize numeric and factor variables
summarize_group <- function(data, group_name) {
  ## Numeric group summary
  num_summary <- data %>%
    select(where(is.numeric)) %>%
    summarise(across(everything(),
                     list(mean = ~mean(.x, na.rm = TRUE),
                          sd   = ~sd(.x, na.rm = TRUE)),
                     .names = "{.col}_{.fn}")) %>%
    pivot_longer(everything(),
                 names_to = c("Variable", ".value"),
                 names_pattern = "(.*)_(mean|sd)") %>%
    mutate(Value = sprintf("%.3f ± %.3f", mean, sd), # change 2f to 3f for metabolites values
           Group = group_name) %>%
    select(Variable, Value, Group, mean, sd)
  
  ## Factor group summary (loop over factors to compute percentages)
  fac_vars <- data %>% select(where(is.factor))
  fac_summary <- lapply(names(fac_vars), function(var) {
    tab <- prop.table(table(fac_vars[[var]])) * 100
    data.frame(Variable = var,
               Level = names(tab),
               Percentage = as.numeric(tab),
               Value = sprintf("%.1f%%", as.numeric(tab)),
               Group = group_name)}) %>% 
    bind_rows()
  bind_rows(num_summary, fac_summary)}

## Format p-values 
format_p <- function(p) {
  if (is.na(p)) return(NA)
  if (p < 0.001) "<0.001"
  else if (p < 0.01) formatC(p, format = "e", digits = 1)
  else round(p, 3)}

## Summaries by sex & overall
male_summary   <- summarize_group(filter(basic, Sex == 0), "Male")
female_summary <- summarize_group(filter(basic, Sex == 1), "Female")
overall_summary <- summarize_group(basic, "Overall")

summary_df <- bind_rows(male_summary, female_summary, overall_summary)

## P-values (Mann-Whitney test)
p_value_df <- basic %>%
  select(where(is.numeric)) %>%
  names() %>%
  setNames(., .) %>%
  sapply(function(var) {
    male_vals   <- basic %>% filter(Sex == 0) %>% pull(var)
    female_vals <- basic %>% filter(Sex == 1) %>% pull(var)
    tryCatch(wilcox.test(male_vals, female_vals)$p.value, error = function(e) NA)}) %>%
  (\(x) data.frame(Variable = names(x),
                   P_Value = sapply(x, format_p),
                   P_Value_Numeric = x))()

## Merge with summary
summary_wide <- summary_df %>%
  pivot_wider(id_cols = c(Variable, Level),
              names_from = Group,
              values_from = Value,
              values_fill = list(Value = "N/A")) %>%
  left_join(p_value_df, by = "Variable") %>%
  mutate(P_Value = ifelse(duplicated(Variable), NA, P_Value)) # Table 1&S8

############################################################################################################# 7.2
########################### 7.2 Log transformation #############################
plot_fig <- function(start_col, end_col, nr, nc, mgp_vals, mar_vals) {
  par(mfrow = c(nr, nc), mgp = mgp_vals, mar = mar_vals)
  
  for (i in start_col:end_col) {
    if (is.numeric(df[, i])) {
      hist(df[, i], nclass = 30, xlab = "",
           main = paste0(names(df[i]), " (", 
                         round(mean(is.na(df[, i])) * 100, 2), "% NA)"),
           cex.main = 0.7, cex.axis = 0.6, cex.lab = 0.8)} 
    else {
      barplot(table(df[, i]), ylab = "Frequency",
              main = paste0(names(df[i]), " (", 
                            round(mean(is.na(df[, i])) * 100, 2), "% NA)"),
              cex.main = 0.7, cex.axis = 0.6, cex.lab = 0.8)}}} 
nc1 <- nc2 <- 5 # columns for first & second figures
nr1 <- nr2 <- 6 # rows for first & second figures
split_idx <- ceiling(ncol(df) / 2) # split index

plot_fig(2, split_idx, 5, 6, c(1.5, 0.5, 0), c(1, 1, 2, 0.5)) # first 27 variables
plot_fig(split_idx + 1, ncol(df), nr2, nc2, c(1, 0.5, 0), c(1, 1, 2, 0.5)) # another 28 variables

vars_to_log <- c(
  "TriacylglycerolmmolL", "hsCRPmgL", "AdiponectinμgL", "InsulinpmolL",
  "Acetate_metabolite", "XXL_VLDL_C_metabolite","XL_HDL_CE_pct_metabolite",
  "Animalprotein_gram_perday", "Animalfat_gram_perday", "Vegetablefat_gram_perday",
  "Totalcarb_gram_perday", "Fiber_gram_perday", "Calcium_gram_perday")

for (var in vars_to_log) {
  if (all(var %in% names(df)) && is.numeric(df[[var]])) {df[[var]] <- log(df[[var]])} 
  else {warning(paste("Variable", var, "is missing or not numeric"))}}

df$PA_all <- log(df$PA_all+1) # log (x+1) is used for PA because PA contains true value "0"; log(0+1)=0

############################################################################################################# 7.3
################### 7.3 Classic linear model - Table S9a-c #####################
### Notes: Replacing "XXL_VLDL_C_metabolite","S_LDL_P_metabolite" with
### "XXL_VLDL_CE_metabolite","S_LDL_PL_metabolite" to assess their impact on BMI 
### variation. Remember to log transfer "XXL_VLDL_CE_metabolite". Report the results
### in Table S9d.
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

format_coef_se <- function(estimate, std_error) {
  format_value <- function(x) {
    ifelse(abs(x) < 0.1,
           format(round(x, 3), nsmall = 3),
           ifelse(abs(x) < 1,
                  format(round(x, 2), nsmall = 2),
                  format(round(x, 1), nsmall = 1)))}
  est_fmt <- format_value(estimate)
  se_fmt <- format_value(std_error)
  paste0(est_fmt, " (", se_fmt, ")")}

format_pval <- function(pval) {
  format(pval, scientific = TRUE, digits = 3)
}

models <- imap(model_list, ~ {
  data <- .x
  lm(BMI ~ ., data = data)}) # LM model

model_summaries <- imap(models, ~ {
  tidy(.x) %>%
    transmute(
      term,
      !!.y := format_coef_se(estimate, std.error),
      !!paste0(.y, "_p") := format_pval(p.value))})

#summary(models$model14)
 
final_summary <- reduce(model_summaries, full_join, by = "term") # results; Table S9a-c

model_stats <- imap_dfr(models, ~ {
  s <- summary(.x)
  f <- s$fstatistic
  
  tibble(
    model = .y,
    adj_r2 = round(s$adj.r.squared, 3),
    model_p = signif(pf(f[1], f[2], f[3], lower.tail = FALSE), 3))}) # R2 and p values

############################################################################################################# 7.4
################### 7.4 Sensitivity analysis - Table S9d #######################
df_sa <- read.csv("model_clean.csv")
ext_ID <- c("228", "629", "408")
df_sa <- df_sa %>% filter(!ID %in% ext_ID)
df_sa$Sex <- as.factor(df_sa$Sex)
df_sa$Sitename <- as.factor(df_sa$Sitename)
df_sa$LDL_to_HDL <- df$LDL_to_HDL # to check whether combining LDL and HDL improves R2
 
model_list_sa <- list(
  model1 = df_sa[, c(2:5,24)], # model 1 = HDL
  model2 = df_sa[, c(2:5,25)], # model 2 = LDL
  model1.5 = df_sa[, c(2:5,63)], # model 1.5 = LDL_to_HDL
  model3 = df_sa[, c(2:5,46)], # model 3 = XXL_VLDL_CE
  model4 = df_sa[, c(2:5,48)], # model 4 = S_LDL-PL
  model5 = df_sa[, c(2:5,52)], # model 5 = Total energy
  model6 = df_sa[, c(2:5,53)], # model 6 = Total protein
  model7 = df_sa[, c(2:5,56)]  # model 7 = Total fat
)

models_sa <- imap(model_list_sa, ~ {
  data <- .x
  outcome <- names(data)[4]
  predictors <- names(data)[-4]
  lm(as.formula(paste(outcome, "~", paste(predictors, collapse = " + "))), data = data)})

model_summaries_sa <- imap(models_sa, ~ {
  tidy(.x) %>% 
    transmute(
      term,
      !!.y := format_coef_se(estimate, std.error),
      !!paste0(.y, "_p") := format_pval(p.value))})

final_summary_sa <- reduce(model_summaries_sa, full_join, by = "term") # Tablse S9d

# summary(models_sa$model7)

model_stats_sa <- imap_dfr(models_sa, ~ {
  s <- summary(.x)
  f <- s$fstatistic
  
  tibble(model = .y, adj_r2 = round(s$adj.r.squared, 3),
         model_p = signif(pf(f[1], f[2], f[3], lower.tail = FALSE), 3))}) # Tablse S9d
