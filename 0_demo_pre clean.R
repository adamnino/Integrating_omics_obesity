############################################################################################################# 0.
##### 0. Pre cleaning - clinical, metabolites, diet, and genetic data
setwd("~/Desktop")
rm(list=ls()) 
library(tidyverse)
library(dplyr) 
df <- read.csv("all data_preclean_ID_recoded.csv")
############################################################################################################# 0.1
####################### 0.1 Clinical data pre cleaning #########################
cli <- df[,1:71]

site_map <- setNames(1:12, c(
  "Site 01", "Site 02", "Site 03", "Site 04", "Site 05", "Site 06", "Site 07",
  "Site 08", "Site 09", "Site 10","Site 11", "Site 12"))

cli$Sitename <- site_map[cli$Sitename]

cli %>% summarise(across(everything(), ~ sum(is.na(.))))

cli <- cli %>% mutate(foodsuppinfo_count = replace_na(foodsuppinfo_count, 0),
                      Medication = replace_na(Medication, 0))

### NOTES: Duration for Polycystic_Ovary_when, Highcholesterol_duration, 
### Hypertension_duration, Throid_disorder_duration, and Diabetes_duration
### needs to be replaced with 0. NAs for PCOS need to be imputed based on sex.

cli_sup <- cli[1:637,] # ppl with complete medical history

cli_sup <- cli_sup %>% mutate(across(c(Highcholesterol_duration, Hypertension_duration,
                              Thyroid_disorder_duration, Diabetes_duration),
                              ~ replace_na(., 0)))

#cli_sup %>% summarise(across(everything(), ~ sum(is.na(.))))

sub_NA <- cli[638:657, ] %>% # ppl with missing medical history
  mutate(Highcholesterol_duration = High_cholesterol,
         Hypertension_duration = Hypertension,
         Thyroid_disorder_duration = Thyroid_disorder,
         Diabetes_duration = Diabetes_other_endocrine,
         Thrombosis_details = as.character(Thrombosis_details))

cli <- bind_rows(cli_sup, sub_NA)

# sub_NA %>% summarise(across(everything(), ~ sum(is.na(.)))) # 7 NAs for PCOS (12-5=7; 5 male)

cli <- cli %>% mutate(across(3, ~ case_when(. == 1 ~ 0, . == 2 ~ 1, TRUE ~ .))) # male (0), female (1)

cli %>% summarise(across(everything(), ~ sum(is.na(.)))) # final check

## Coding PCOS as 0 for female (no PCOS), 1 for female (yes PCOS), 2 for male  (n/a)
cli$Polycystic_Ovary[cli$Sex == 0] <- 2
cli$Polycystic_Ovary_when[cli$Sex == 0] <- 2

## Coding Menstrual_cycle_regular as 0 for female (regular menstrual), 1 for female (irregular menstrual), 2 for male  (n/a)
cli$Menstrual_cycle_regular[cli$Sex == 0] <- 2

#write.csv(cli,"~/Desktop/clinical_raw.csv", row.names=FALSE)

############################################################################################################# 0.2
####################### 0.2 Metabolite data pre cleaning #######################
### COMMENTS: Fixed in excel (percentages for 5 metabolites were not calculated): 
### XL_VLDL_PL_pct_metabolite (4), XL_VLDL_C_pct_metabolite (4), 
### XL_VLDL_CE_pct_metabolite (4), XL_VLDL_FC_pct_metabolite (4), 
### XL_VLDL_TG_pct_metabolite (4). 
metabo <- df[,c(1,42,179:449)]
metabo %>% summarise(across(everything(), ~ sum(is.na(.))))

metabo <- metabo %>% select(-Glycerol_metabolite) # bad quality according to Nightingale report, 201 NAs

metabo <- metabo %>% mutate(across(where(is.character), ~ na_if(., "NA")))

metabo <- metabo[,c(1,24:272)]

metabo[metabo == 0] <- NA # replace 0s with NA

metabo %>% summarise(across(everything(), ~ any(is.na(.)))) %>%
  select(where(~ .)) %>%   # keep only columns where TRUE (contains NA)
  ncol()  # sum of missing metabolites (4)

metabo %>% summarise(across(everything(), ~ sum(is.na(.)))) %>% sum() # sum of all NAs (88)

#write.csv(metabo,"~/Desktop/metabo_raw.csv", row.names=FALSE)

############################################################################################################# 0.3
####################### 0.3 Diet data pre cleaning #########################
diet <- df[,c(1,72:81)]
diet$Fiber_gram_perday <- replace(diet$Fiber_gram_perday, diet$Fiber_gram_perday < 0, NA)

diet %>% summarise(across(everything(), ~ sum(is.na(.)))) # 8 NAs for fiber

#write.csv(diet,"~/Desktop/diet_raw.csv", row.names=FALSE)

############################################################################################################# 0.4
####################### 0.4 Genetic data pre cleaning ##########################
snp <- df[,c(1,84:178)]
snp %>%
  summarise(across(everything(), ~ sum(. %in% c("?", "Bad"), na.rm = TRUE)))

snp <- snp %>%
  mutate(across(where(is.character), ~ na_if(., "?")))

snp <- snp %>%
  mutate(across(where(is.character), ~ na_if(., "Bad")))

snp %>%
  summarise(across(everything(), ~ sum(is.na(.))))

## Count complete vs. incomplete cases (# of Participants containing missing values)
tibble(Status = c("Incomplete", "Complete"),
       Count = as.numeric(table(complete.cases(snp))),
       Percent = round(100 * table(complete.cases(snp)) / nrow(snp), 2)) # 17.05% incomplete. (112 out of 657)

column <- (names(snp)[sapply(snp, function(x) !any(is.na(x)))])

snp_no_na <- snp %>% select(all_of(column))

snp_no_na %>% summarise(across(everything(), ~ sum(is.na(.))))

snp_no_na <- as.data.frame(lapply(snp_no_na, function(col) gsub(":", "", col))) %>% 
  select(-c("rs405697","rs12255372")) # rs405697 (Alzheimer); rs12255372(T2D)

write.csv(snp_no_na,"~/Desktop/snp_raw.csv", row.names=FALSE)
