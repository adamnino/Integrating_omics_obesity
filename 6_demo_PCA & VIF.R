############################################################################################################# 6.
##### 6. PCA outliers (Fig. S9) & VIF check (Table S7)
setwd("~/Desktop") 
rm(list=ls()) 
library(tidyverse)
library(dplyr)
library(ggplot2)
library(car) # 6.1 - VIF
library(moments) # 6.3 - data normality
df <- read.csv("model_clean.csv") 
############################################################################################################# 6.1
############# 6.1 Investigating multicollenarity (VIF) - Table S7 ##############
### Step 1, check VIF before combining/deleting variables (51 variables, n = 657) 
### Step 2, check VIF after combining/deleting variables (45 variables, n = 657)
### Step 3, check VIF after removing outliers based on PCA (45 variables, n = 654)
df <- df %>% 
  mutate(across(c(2,3,6,8,9,11,13,15,30,32,34,35), as.factor)) %>% 
  mutate(across(c(4,5,7,16:31,33,36:62), as.numeric)) %>%  
  mutate(across(c(1,10,12,14), as.character)) # steps 1-3

# ext_ID <- c("228", "629", "408") # step 3 only
# df <- df %>% filter(!ID %in% ext_ID) # remove ID after PCA (Step 3 only)

df_new <- df # steps 1-3
df_new$LDL_to_HDL <- df_new$LDLCmmolL / df_new$HDLCmmolL # steps 2&3
df_new <- df_new %>%
  select(-c("LDLCmmolL","HDLCmmolL","XXL_VLDL_CE_metabolite","S_LDL_PL_metabolite",
            "Totalprotein_gram_perday","TotalEnergy_kcal_gram_perday","Totalfat_gram_perday"))  # steps 2$3

lm_test <- lm(BMI ~ ., df_new %>% 
                select(-c("ID","BMI_group", "Category","ME_notes","Sub_notes","OMH_notes",
                          "OE_notes", "ME_code","OMH_code","OE_code")))  # steps 1-3
summary(lm_test) # steps 1-3

vif <- vif(lm_test) # steps 1-3; Table S7
view(vif) # steps 1-3; after obtaining VIF following the step, rerun the code for step 2 from line 10

## After running step 2, run the following code
df <- df_new # step 3; Convert it back to df, prepare it for PCA

## For classic LM & ML model; after obtaining VIF following step 3, save df as final data_17SNP.csv
# write.csv(df,"~/Desktop/final data_17SNP.csv", row.names=FALSE)  

########################################################################################## 6.2
################ 6.2 PCA to visualize outliers - Fig. S9 #######################
### Only continuous outcomes from clinical, metabolite, and genetics. Data such as 
### FFQ, PA, and medical history are self-reported, which could be biased. 

## Contain measurable (reliable) values only; age is not included because they are half 20, half 21.
pca <- df[,-c(1:16,28:33,48:56)] %>%  # remove "WC" for PCA; using df from Step2&3 only
  log() # Log10 for PCA

pca_selected <- prcomp(pca, center = TRUE, scale = TRUE) # PCA
exp_var_selected <- summary(pca_selected)$importance[2, ] * 100 # explained variance

## Prepare results with metadata
pca_selected_results <- as.data.frame(pca_selected$x) %>%
  mutate(
    category = df$Category,
    bmi = df$BMI,
    ID = df$ID,
    ID_new = seq_len(nrow(.)),
    bmi_group = df$BMI_group,
    ME = ifelse(df$ME_notes == "No ME", "", as.character(df$ME_notes)),
    OE = df$OE_notes,
    OMH = df$OMH_notes,
    Sub = df$Sub_notes)

## Centroid calculation for each cluster; Threshold for outliers (beyond 3 SD in this study)!!
centroid <- colMeans(pca_selected_results[, c("PC1", "PC2")])
pca_selected_results <- pca_selected_results %>%
  mutate(
    distance = sqrt((PC1 - centroid["PC1"])^2 + (PC2 - centroid["PC2"])^2),
    is_outlier = distance > mean(distance) + 3 * sd(distance))

## Plot with labels for outliers only (check category, BMI group, and ME)
ggplot(pca_selected_results, aes(PC1, PC2, color = as.factor(bmi_group))) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95, linetype = "dashed", size = 1) +
  theme_minimal() +
  labs(
    title = "PCA - with outliers", # change with outliers/without outliers
    x = sprintf("PC1 (%.1f%%)", exp_var_selected[1]),
    y = sprintf("PC2 (%.1f%%)", exp_var_selected[2]),
    color = "BMI") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(
    aes(label = ifelse(is_outlier, paste(ID_new, category, OE, OMH, ME, sep = "-"), "")),
    vjust = -1, hjust = 0.6, size = 1.9) # Fig.S9

### Exclude the following participants; rerun the code from line 1
### 228: over, female, irregular, 2, D, Digestive system - Biliary - fatty liver
### 629: obese, female,            1, D, Digestive system - Biliary - fatty liver
### 408: over, female, irregular, 3, I, Endocrine      - thyroid_4 -  hyperthyroid (current), medication


############################################################################################################# 6.3
########################### 6.3 Data normality #################################
## Skewness and kurtosis; excess kurtosis = kurtosis−3
skew_vals <- sapply(df, function(x) if (is.numeric(x)) skewness(x, na.rm = TRUE) else NA)
kurt_vals <- sapply(df, function(x) if (is.numeric(x)) kurtosis(x, na.rm = TRUE) - 3 else NA)
normality_summary <- data.frame(
  Variable = names(df),
  Skewness = skew_vals,
  Abs_Skewness = abs(skew_vals),
  Kurtosis_Excess = kurt_vals,
  Abs_Kurtosis_Excess = abs(kurt_vals))

non_normal_vars <- normality_summary[
  (normality_summary$Abs_Skewness >= 2 | normality_summary$Abs_Kurtosis_Excess >= 4) &
    !is.na(normality_summary$Abs_Skewness),] # ignore MHF (treated as discrete data).
