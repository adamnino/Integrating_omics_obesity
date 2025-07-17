############################################################################################################# 7.
##### 7. Classic linear model (Table S10) & Study demographic characteristics (Table S9)
setwd("~/Desktop")
rm(list=ls()) #clear all variables in workplace
library(dplyr)
library(tidyverse)
library(broom) # 7.3 - Classic linear model (tidy coefficient, SE, p-value)
df <- read.csv("final data_17SNP.csv")
df$S_LDL_P_metabolite <- df$S_LDL_P_metabolite*1000 # mmol/L --> umol/L
df <- df %>% 
  mutate(across(c(2,3,6,8,9,11,13,15,30,32,33), as.factor)) %>% 
  mutate(across(c(4,5,7,16:29,31,34:57), as.numeric)) %>% 
  mutate(across(c(1,10,12,14), as.character))

############################################################################################################# 7.1
############## 7.1 Basic demographic characteristics - Table S9 ################
basic <- df[,-c(1,8:15)]

factor_summary <- function(x) {
  prop <- prop.table(table(x)) * 100
  data.frame(Level = names(prop), Percentage = as.numeric(prop))
}

## Split numeric and factor vars
num_vars <- basic %>% select(where(is.numeric))
fac_vars <- basic %>% select(where(is.factor))

numeric_summary <- basic %>%
  select(where(is.numeric)) %>%
  summarise(across(
    everything(),
    list(mean = ~round(mean(.x, na.rm = TRUE), 3),
         sd   = ~round(sd(.x, na.rm = TRUE), 3))
  )) %>%
  pivot_longer(
    everything(),
    names_to = c("Variable", "Statistic"),
    names_pattern = "^(.*)_(mean|sd)"
  ) %>%
  pivot_wider(names_from = Statistic, values_from = value) # Numeric summary

fac_summary <- lapply(names(fac_vars), function(var) {
  temp <- factor_summary(fac_vars[[var]])
  temp$Variable <- var
  temp
}) %>%
  bind_rows() %>%
  select(Variable, Level, Percentage) # Factor summary

## Combine both; Table S9
summary_df <- bind_rows(numeric_summary, fac_summary) 

############################################################################################################# 7.2
########################### 7.2 Log transformation #############################
plot_fig <- function(start_col, end_col, nr, nc, mgp_vals, mar_vals) {
  par(mfrow = c(nr, nc), mgp = mgp_vals, mar = mar_vals)
  
  for (i in start_col:end_col) {
    if (is.numeric(df[, i])) {
      hist(df[, i], nclass = 30, xlab = "",
           main = paste0(names(df[i]), " (", 
                         round(mean(is.na(df[, i])) * 100, 2), "% NA)"),
           cex.main = 0.7, cex.axis = 0.6, cex.lab = 0.8)
    } else {
      barplot(table(df[, i]), ylab = "Frequency",
              main = paste0(names(df[i]), " (", 
                            round(mean(is.na(df[, i])) * 100, 2), "% NA)"),
              cex.main = 0.7, cex.axis = 0.6, cex.lab = 0.8)
    }
  }
} 
nc1 <- nc2 <- 5 # Columns for first & second figures
nr1 <- nr2 <- 6 # Rows for first & second figures
split_idx <- ceiling(ncol(df) / 2) # Split index

plot_fig(2, split_idx, 5, 6, c(1.5, 0.5, 0), c(1, 1, 2, 0.5)) # First 27 variables
plot_fig(split_idx + 1, ncol(df), nr2, nc2, c(1, 0.5, 0), c(1, 1, 2, 0.5)) # Another 28 variables

vars_to_log <- c(
  "TriacylglycerolmmolL", "hsCRPmgL", "AdiponectinμgL", "InsulinpmolL",
  "Acetate_metabolite", "XXL_VLDL_C_metabolite","XL_HDL_CE_pct_metabolite",
  "Animalprotein_gram_perday", "Animalfat_gram_perday", "Vegetablefat_gram_perday",
  "Totalcarb_gram_perday", "Calcium_gram_perday")

for (var in vars_to_log) {
  if (all(var %in% names(df)) && is.numeric(df[[var]])) {
    df[[var]] <- log(df[[var]])
  } else {
    warning(paste("Variable", var, "is missing or not numeric"))
  }
}

############################################################################################################# 7.3
#################### 7.3 Classic linear model - Table S10 ######################
model_list <- list(
  model0 = df[, c(2:5)],# model0 = Base
  model1 = df[, c(2:5,16)], # model 1 = WC
  model2 = df[, c(2:5,7,16:33,57)], # model 2 = clinical
  model3 = df[, c(2:5,7,17:33,57)],  # model 3 = clinical (no WC)
  model4 = df[, c(2:5,34:47)],  # model 4 = Metabolites
  model5 = df[, c(2:5,48:55)],  # model 5 = Diet
  model6 = df[, c(2:5,56)],  # model 6 = SNPs
  
  model7 = df[, c(2:5,7,17:47,57)], # model 7 = clinical + metabolites (no WC)
  model8 = df[, c(2:5,7,17:33,48:55,57)],  # model 8 = clinical + diet (no WC)
  model9 = df[, c(2:5,7,17:33,56,57)],  # model 9 = clinical + SNP (no WC)
  model10 = df[, c(2:5,7,17:57)],  # model 10 = clinical + metabolites + SNP + diet (no WC)
  
  model11 = df[, c(2:5,7,16:47,57)], # model 11 = clinical + metabolites
  model12 = df[, c(2:5,7,16:33,48:55,57)],  # model 12 = clinical + diet
  model13 = df[, c(2:5,7,16:33,56,57)],  # model 13 = clinical + SNP
  model14 = df[, c(2:5,7,16:57)]  # model 14 = clinical + metabolites + SNP + diet 
)

format_coef_se <- function(estimate, std_error) {
  format_value <- function(x) {
    ifelse(abs(x) < 0.1,
           format(round(x, 3), nsmall = 3),
           ifelse(abs(x) < 1,
                  format(round(x, 2), nsmall = 2),
                  format(round(x, 1), nsmall = 1)))
  }
  est_fmt <- format_value(estimate)
  se_fmt <- format_value(std_error)
  paste0(est_fmt, " (", se_fmt, ")")
}

format_pval <- function(pval) {
  ifelse(
    pval < 0.0001,
    "< 0.0001",
    ifelse(
      pval < 0.05,
      paste0(format(round(pval, 6), nsmall = 6), "*"),
      format(round(pval, 3), nsmall = 3)
    )
  )
}

models <- imap(model_list, ~ {
  data <- .x
  lm(BMI ~ ., data = data)  
}) # LM model

model_summaries <- imap(models, ~ {
  tidy(.x) %>%
    transmute(
      term,
      !!.y := format_coef_se(estimate, std.error),
      !!paste0(.y, "_p") := format_pval(p.value)
    )
})

#summary(models$model14)

final_summary <- reduce(model_summaries, full_join, by = "term") # Results; Table S10

model_stats <- imap_dfr(models, ~ {
  s <- summary(.x)
  f <- s$fstatistic
  
  tibble(
    model = .y,
    adj_r2 = round(s$adj.r.squared, 3),
    model_p = signif(pf(f[1], f[2], f[3], lower.tail = FALSE), 3)
  )
}) # R2 and p values
