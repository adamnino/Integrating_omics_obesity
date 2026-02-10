############################################################################################################# 4.
##### 4. Genetic - GRS calculation (Table S6)
setwd("~/Desktop")
rm(list=ls()) 
library(dplyr)
library(tidyr)  
library(purrr)
library(HardyWeinberg) # 4.2 - HWE
bmi <- read.csv("all data_preclean_ID_recoded.csv")
snp <- read.csv("snp_raw.csv") 
############################################################################################################# 4.1
############# 4.1 Genotype count & Coefficient p-value - Table S6 ##############
snp <- snp [,-1]
snps_num <- c("rs543874_2429", "rs7138803", "rs1558902", "rs7903146", "rs2112347", 
          "rs10182181", "rs13021737", "rs16851483", "rs17094222", "rs17405819", 
          "rs7141420", "rs10733682", "rs11583200", "rs3101336", "rs3736485", 
          "rs205262", "rs12940622")
risk_alleles <- c("G", "A", "A", "C", "T", "G", "G", "T", "C", "T", "T", "A", "C", "C", "A", "G", "G") # Locke et al.
 
snp <- snp %>% select(any_of(snps_num)) # only selecting 17 SNPs
 
## Genotype counts
genotype_counts <- map_dfr(names(snp), ~ {
  genotype_table <- table(snp[[.x]])
  total <- sum(genotype_table)
  
  data.frame(
    SNP = .x,
    Genotype = names(genotype_table),
    Count = as.integer(genotype_table),
    Percent = round(100 * as.integer(genotype_table) / total),
    stringsAsFactors = FALSE
  )
}) %>%
  mutate(Count_Percent = paste0(Count, " (", Percent, "%)"))


## Coefficient p-value (FDR)
coef_p_snps_num <- function(df, outcome, risk_alleles) {
  stopifnot(ncol(df) == length(risk_alleles))
  
  map2_dfr(df, risk_alleles, function(genotype_col, ra) {
    allele_recode <- vapply(strsplit(genotype_col, ""), function(alleles) sum(alleles == ra), numeric(1))
    model <- lm(outcome ~ allele_recode)
    coef_summary <- summary(model)$coefficients["allele_recode", ]
    
    tibble(
      beta = coef_summary["Estimate"],
      p_value = coef_summary["Pr(>|t|)"]
    )
  }) %>%
    mutate(
      SNP = colnames(df),
      risk_allele = risk_alleles,
      FDR = p.adjust(p_value, method = "BH")
    ) %>%
    select(SNP, risk_allele, beta, p_value, FDR)
}

results <- coef_p_snps_num(df = snp, outcome = bmi$BMI, risk_alleles = risk_alleles) # Table S6

############################################################################################################# 4.2
#################### 4.2 HWE & MAF calculations - Table S6 #####################
## Allele frequencies & p-values for HWE
hwe <- genotype_counts

run_hwe_all <- function(hwe, snp_id) {
  snp_data <- hwe %>%
    filter(SNP == snp_id) %>%
    mutate(Genotype_std = sapply(strsplit(as.character(Genotype), ""), function(x) paste0(sort(x), collapse = ""))) %>%
    group_by(Genotype_std) %>%
    summarise(Count = sum(Count), .groups = "drop")
  
  allele_counts <- snp_data %>%
    rowwise() %>%
    mutate(alleles = list(strsplit(Genotype_std, "")[[1]])) %>%
    unnest(cols = c(alleles)) %>%
    group_by(alleles) %>%
    summarise(count = sum(Count), .groups = "drop")
  
  if (nrow(allele_counts) != 2) {
    return(data.frame(SNP = snp_id, p_value = NA, ref_freq = NA, alt_freq = NA))
  }
  
  alleles_ordered <- allele_counts %>%
    arrange(desc(count)) %>%
    pull(alleles)
  
  ref <- alleles_ordered[1]
  alt <- alleles_ordered[2]
  
  AA <- paste0(ref, ref)
  AB <- paste0(sort(c(ref, alt)), collapse = "")
  BB <- paste0(alt, alt)
  
  x <- c(
    AA = snp_data$Count[snp_data$Genotype_std == AA],
    AB = snp_data$Count[snp_data$Genotype_std == AB],
    BB = snp_data$Count[snp_data$Genotype_std == BB]
  )
  x[is.na(x)] <- 0
  
  total_alleles <- sum(x) * 2
  ref_count <- 2 * x["AA"] + x["AB"]
  alt_count <- 2 * x["BB"] + x["AB"]
  freq_ref <- round(ref_count / total_alleles, 4)
  freq_alt <- round(alt_count / total_alleles, 4)
  
  result <- HWChisq(x, cc = 0, verbose = FALSE)
  
  return(data.frame(
    SNP = snp_id,
    p_value = result$pval,
    ref_freq = paste0(ref, " (", freq_ref, ")"),
    alt_freq = paste0(alt, " (", freq_alt, ")")
  ))
} 

snp_list <- unique(hwe$SNP) # 17 SNPs

hwe_results <- do.call(rbind, lapply(snp_list, function(snp) run_hwe_all(hwe, snp))) # Table S6

############################################################################################################# 4.3
############################ 4.3 GRS calculation ###############################
## Count risk alleles in each genotype
names(risk_alleles) <- snps_num
count_risk <- function(genotype, risk) stringr::str_count(genotype, risk)
 
## Add GRS calculation
grs_df <- snp %>%
  rowwise() %>%
  mutate(GRS = sum(
    c_across(all_of(snps_num)) %>%
      purrr::map2_dbl(risk_alleles[snps_num], ~ count_risk(.x, .y)),
    na.rm = TRUE
  )) %>%
  ungroup()

grs_df$ID <- bmi$ID

#write.csv(grs_df, "~/Desktop/snp_clean.csv", row.names = FALSE)