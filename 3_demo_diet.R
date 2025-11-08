############################################################################################################# 3.
##### 3. Diet - Imputation
setwd("~/Desktop")
rm(list=ls())
library(tidyverse)
library(dplyr)
library(ggplot2)
library(VIM) # 3.2 - RF imputation 
diet <- read.csv("diet_raw.csv")
############################################################################################################# 3.1
####################### 3.1 Investigating missingness ##########################
diet %>%
  summarise(across(everything(), ~ sum(is.na(.))))
Ncc_diet <- tibble(
  Status = c("Incomplete", "Complete"),
  Count = table(complete.cases(diet)),
  Percent = round(100 * table(complete.cases(diet)) / nrow(diet), 2))
Ncc_diet # 0.15% incomplete

############################################################################################################# 3.2
########################## 3.2 Imputation, Diagnostics #########################
set.seed(5)
imp <- kNN(diet, k = 5)

imp_clean <- imp[, 1:ncol(diet)] # imputed value is 6 (ID464)

imp_clean %>%
  summarise(across(everything(), ~ sum(is.na(.))))

##  Diagnostic plot
imp_flags <- imp[, grep("_imp$", names(imp))] # Identify imputed values

imp_clean$Data <- imp_flags[, "Eat_salty_imp"] # Mark imputed value

imp_clean$Data <- ifelse(imp_clean$Data, "Imputed", "Observed")

ggplot(imp_clean, aes(x = "", y = Eat_salty, color = Data)) +
  geom_jitter(width = 0.2, size = 3) +
  labs(title = "Stripplot of imputed diet data using KNN imputation",
       y = "Eat salty", x = "") +
  scale_color_manual(values = c("Observed" = "darkblue", "Imputed" = "orange")) +
  theme_minimal()

#write.csv(imp_clean, "~/Desktop/diet_clean.csv", row.names = FALSE)