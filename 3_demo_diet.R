############################################################################################################# 3.
##### 3. Diet - Imputation (Fig.S8) 
setwd("~/Desktop")
rm(list=ls())
library(tidyverse)
library(dplyr)
library(ggplot2)
library(VIM) # 3.2 - RF imputation 
diet <- read.csv("diet_raw.csv")
############################################################################################################# 3.1
####################### 3.1 Investigating missingness ##########################
missing_counts <- colSums(is.na(diet))
missing_counts[missing_counts > 0] # values that were missing

## Count complete vs. incomplete cases (# of Participants containing missing values)
tibble(Status = c("Incomplete", "Complete"),
       Count = as.numeric(table(complete.cases(diet))),
       Percent = round(100 * table(complete.cases(diet)) / nrow(diet), 2)) # 1.22% incomplete

## Number and proportion of missing values per variable
cbind("# NA" = sort(colSums(is.na(diet))),
      "% NA" = round(sort(colMeans(is.na(diet))) * 100, 2)) # Fig. S8

## Overall missingness
cat("\nOverall missingness:", round(100 * sum(is.na(diet)) / prod(dim(diet)), 2), "%\n\n") # 0.11%

############################################################################################################# 3.2
###################### 3.2 Imputation, Diagnostics - Fig.S8 ####################
set.seed(5)
imp <- kNN(diet, k = 5)

imp_clean <- imp[, 1:ncol(diet)]

imp_clean %>%
  summarise(across(everything(), ~ sum(is.na(.))))
 
##  Diagnostic plot
imp_flags <- imp[, grep("_imp$", names(imp))] # identify imputed values

imp_clean$Data <- imp_flags[, "Fiber_gram_perday_imp"] # mark imputed value

imp_clean$Data <- ifelse(imp_clean$Data, "Imputed", "Observed")

ggplot(imp_clean, aes(x = "", y = Fiber_gram_perday, color = Data)) +
  geom_jitter(width = 0.2, size = 3) +
  labs(title = "Stripplot of imputed diet data using KNN imputation",
       y = "Fiber", x = "") +
  scale_color_manual(values = c("Observed" = "darkblue", "Imputed" = "orange")) +
  theme_minimal() # Fig.S8

imp_clean <- imp_clean[,-12]
#write.csv(imp_clean, "~/Desktop/diet_clean.csv", row.names = FALSE)