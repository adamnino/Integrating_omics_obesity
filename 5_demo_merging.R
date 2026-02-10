############################################################################################################# 5.
##### 5. Merging all data
setwd("~/Desktop")
rm(list=ls())
library(tidyverse)
library(dplyr)
cli <- read.csv("clinical_clean.csv") 
met <- read.csv("met_clean.csv")
diet <- read.csv("diet_clean.csv") 
snp <- read.csv("snp_clean.csv") # 17 SNPs
snp <- snp[,c(18,19)] # 1 GRS  
   
############################################################################################################# 5.1
############################ 5.1 Merging all data ##############################
## Merged data containing clinical, diet, metabolites, and SNPs
dfs <- list(cli, met, diet, snp)
merged_data <- reduce(dfs, ~ left_join(.x, .y, by = "ID"))

merged_data <- merged_data %>% select(1:5, BMI_group, everything())

merged_data  %>% summarise(across(everything(), ~ sum(is.na(.)))) %>% sum()

#write.csv(merged_data,"~/Desktop/model_clean.csv", row.names=FALSE)