
# Load required packages
library(dplyr)
library(readr)
library(nutriR)
library(ggplot2)
library(tidyr)
library(purrr)
library(countrycode)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
options(scipen=999)
rm(list = ls())

# Load the pre-processed data from the previous script
prepped_data <- readRDS("./output/gdd_prepped_data_for_caloric_adjustment.rds")

# Load the aggregated population data
population_data <- readRDS("./output/wpp2024_population_aggregated_2018.rds")


# All intakes are reported adjusted to 700 kcal/day for ages <1 year, 
# 1000 kcal/day for ages 1–<2 years, 1300 kcal/day for ages 2–5 years, 
# 1700 kcal/day for ages 6–10 years, 2000 kcal/day for ages 11–74 years, 
# and 1700 kcal/day for ages ≥75 years.

#source: Miller V, Reedy J, Cudhea F, Zhang J, Shi P, Erndt-Marino J, Coates J, 
#Micha R, Webb P, Mozaffarian D; Global Dietary Database. 
# Global, regional, and national consumption of animal-source foods between 
#1990 and 2018: findings from the Global Dietary Database. 
# Lancet Planet Health. 2022 Mar;6(3):e243-e256. 
# doi: 10.1016/S2542-5196(21)00352-1. 

# GDD's caloric standardization values
gdd_kcal_standards_old <- tibble::tribble(
  ~age_group, ~standard_kcal,
  "0-0.99",      700,
  "1-4",        1300, # Note: You need to decide how to handle the 1-1.99 and 2-4 split.
  # A simple weighted average is fine: (1*1000 + 3*1300)/4 = 1225. Let's use that.
  "5-9",        1700,
  "10-14",      2000,
  "15-19",      2000,
  "20-24",      2000,
  "25-29",      2000,
  "30-34",      2000,
  "35-39",      2000,
  "40-44",      2000,
  "45-49",      2000,
  "50-54",      2000,
  "55-59",      2000,
  "60-64",      2000,
  "65-69",      2000,
  "70-74",      2000,
  "75-79",      1700,
  "80-84",      1700,
  "85-89",      1700,
  "90-94",      1700,
  "95-99",      1700
)



# GDD's caloric standardization values
gdd_kcal_standards <- tibble::tribble(
  ~age_group, ~standard_kcal,
  "0-0.99",     700,
  "1-4",        1300, # Note: You need to decide how to handle the 1-1.99 and 2-4 split.
  # A simple weighted average is fine: (1*1000 + 3*1300)/4 = 1225. Let's use that.
  "5-9",        1700,
  "10-14",      2000,
  "15-19",      2000,
  "20-24",      2000,
  "25-29",      2000,
  "30-34",      2000,
  "35-39",      2000,
  "40-44",      2000,
  "45-49",      2000,
  "50-54",      2000,
  "55-59",      2000,
  "60-64",      2000,
  "65-69",      2000,
  "70-74",      2000,
  "75-79",      2000,
  "80-84",      2000,
  "85-89",      2000,
  "90-94",      2000,
  "95-99",      2000
)
#question: is it really not sex dependent? seems strange

# Correcting the 1-4 group based on the components:
# GDD: 1–<2 years = 1000 kcal; 2–5 years = 1300 kcal.
# Your 1-4 group is a 1-year band and a 3-year band.
# Weighted avg: (1 * 1000 + 3 * 1300) / 4 = 1225 kcal.
#gdd_kcal_standards$standard_kcal[gdd_kcal_standards$age_group == "1-4"] <- 1225


# Join with standardization values and calculate protein's share of calories
prepped_data <- prepped_data %>%
  left_join(gdd_kcal_standards, by = "age_group") %>%
  mutate(
    # Protein provides 4 kcal/gram
    protein_kcal_share_mean  = (gdd_mean * 4) / standard_kcal,
    protein_kcal_share_lower = (gdd_lower * 4) / standard_kcal,
    protein_kcal_share_upper = (gdd_upper * 4) / standard_kcal
  )


# ==========================================================================
# SCRIPT 2: "OPTIMAL CALORIES" SCENARIO
# ==========================================================================

# --- Step 1: Load and Prepare Marco's Country-Specific EER Data ---
# We apply the same rigorous cleaning logic used for the weight data in Script 1.

cat("\n--- Loading and preparing Marco's EER data... ---\n")

library(readxl)
marco_eer_raw <- read_excel("./data/data_extract_040725.xlsx", sheet = 1)

# First, pivot the data from long to wide, creating separate columns for mean, low, high
marco_eer_wide <- marco_eer_raw %>%
  pivot_wider(
    names_from = Stats,
    values_from = Value,
    names_prefix = "eer_"
  )

# Now, apply the full cleaning and harmonization process
marco_eer_processed <- marco_eer_wide %>%
  # Filter for relevant Year (data is only for 2018)
  filter(Year == 2018) %>%
  # Filter out aggregate regions (e.g., "WLD", "HIC") by keeping only 3-char country codes
  filter(nchar(Region) == 3) %>%
  # Filter for single-sex data ("MLE", "FML"), excluding "BTH"
  filter(Sex %in% c("MLE", "FML")) %>%
  # Filter out aggregate age groups
  filter(!Age %in% c("all-a", "20+")) %>%
  # Rename Region to iso3
  rename(iso3 = Region) %>%
  # Map Sex codes to the standard "Males"/"Females"
  mutate(
    sex = case_when(
      Sex == "MLE"   ~ "Males",
      Sex == "FML" ~ "Females"
    )
  ) %>%
  # Harmonize Age Group labels
  mutate(
    age_group = case_when(
      Age == "<1" ~ "0-0.99",
      Age == "95+" ~ "95-99",
      TRUE ~ Age
    )
  ) %>%
  # Select and rename the final, clean columns
  select(
    iso3,
    sex,
    age_group,
    eer_kcal_marco_mean = eer_mean,
    eer_kcal_marco_low  = eer_low,
    eer_kcal_marco_high = eer_high
  ) %>%
  distinct() # Ensure no accidental duplicates

cat("--- Marco's EER data prepared. ---\n")


# --- Step 2: Create the Master Data Frame for THIS ANALYSIS ---
cat("\n--- Assembling data for the Optimal Calories analysis... ---\n")

analysis_data <- prepped_data %>%
  left_join(population_data, by = c("iso3", "sex", "age_group")) %>%
  left_join(marco_eer_processed, by = c("iso3", "sex", "age_group"))

# --- Step 3 and beyond (the rest of the analysis) follows from here... ---

# For now, let's verify the join and data structure.
cat("--- Verification of the data assembly ---\n")
glimpse(analysis_data)
coverage_pct <- round(100 * sum(!is.na(analysis_data$eer_kcal_marco_mean)) / nrow(analysis_data), 1)
cat("Coverage of Marco's EER data:", coverage_pct, "%\n")
if (coverage_pct < 95) {
    cat("⚠️ WARNING: Coverage from Marco's EER data is low. Check column names and join keys.\n")
} else {
    cat("✅ EER data join successful.\n")
}


# ==========================================================================
# STEP 2.1: IMPUTE MISSING EER DATA TO ACHIEVE 100% COVERAGE
# ==========================================================================
cat("\n--- Imputing missing EER data based on diagnosis... ---\n")

# --- 2.1.1: Impute South Sudan (SSD) using Sudan (SDN) data ---

# First, create the lookup table from our existing analysis_data containing SDN's EER.
# This is the most reliable source.
sudan_eer_lookup <- analysis_data %>%
  filter(iso3 == "SDN") %>%
  select(
    sex, age_group, 
    impute_mean = eer_kcal_marco_mean,
    impute_low  = eer_kcal_marco_low,
    impute_high = eer_kcal_marco_high
  ) %>%
  distinct()

# Now, perform the imputation.
# We will also add a tracking column to document our changes.
analysis_data <- analysis_data %>%
  # Join the Sudan lookup data. This will add the 'impute' columns to every row,
  # but they will only be used for SSD.
  left_join(sudan_eer_lookup, by = c("sex", "age_group")) %>%
  
  # Use mutate to fill in the gaps for SSD
  mutate(
    # Create the tracking column first, initially NA for all rows
    eer_imputed_source = NA_character_,
    
    # Flag if this specific imputation happens
    imputed_from_sdn = (iso3 == "SSD" & is.na(eer_kcal_marco_mean)),
    
    # Impute the values
    eer_kcal_marco_mean = if_else(imputed_from_sdn, impute_mean, eer_kcal_marco_mean),
    eer_kcal_marco_low  = if_else(imputed_from_sdn, impute_low,  eer_kcal_marco_low),
    eer_kcal_marco_high = if_else(imputed_from_sdn, impute_high, eer_kcal_marco_high),
    
    # Update the tracking column based on the flag
    eer_imputed_source = if_else(imputed_from_sdn, "Imputed from SDN", eer_imputed_source)
  ) %>%
  
  # Remove the temporary helper columns
  select(-starts_with("impute"), -imputed_from_sdn)

cat("--- SSD imputation complete. ---\n")


# --- 2.1.2: Impute KIR & MHL for age_group "95-99" using their "90-94" data ---

# Create the lookup table using data from the 90-94 age group for these specific countries.
# We use the already-updated `analysis_data` as our source.
impute_90_94_lookup <- analysis_data %>%
  filter(iso3 %in% c("KIR", "MHL") & age_group == "90-94") %>%
  select(
    iso3, sex,
    lookup_mean = eer_kcal_marco_mean,
    lookup_low  = eer_kcal_marco_low,
    lookup_high = eer_kcal_marco_high
  ) %>%
  # Make sure we don't use NA values as a source
  filter(!is.na(lookup_mean))

# Join this lookup data and fill in the final gaps
analysis_data <- analysis_data %>%
  left_join(impute_90_94_lookup, by = c("iso3", "sex")) %>%
  mutate(
    # Flag the specific rows to be imputed in this step
    imputed_from_90_94 = (iso3 %in% c("KIR", "MHL") & age_group == "95-99" & is.na(eer_kcal_marco_mean)),
    
    # Impute the values
    eer_kcal_marco_mean = if_else(imputed_from_90_94, lookup_mean, eer_kcal_marco_mean),
    eer_kcal_marco_low  = if_else(imputed_from_90_94, lookup_low,  eer_kcal_marco_low),
    eer_kcal_marco_high = if_else(imputed_from_90_94, lookup_high, eer_kcal_marco_high),
    
    # Use case_when to safely update the tracking column without overwriting the SDN flag
    eer_imputed_source = case_when(
      imputed_from_90_94 ~ "Imputed from 90-94",
      TRUE ~ eer_imputed_source # Keep existing value (NA or "Imputed from SDN") otherwise
    )
  ) %>%
  # Clean up
  select(-starts_with("lookup"), -imputed_from_90_94)

cat("--- 95-99 age group imputation complete. ---\n")


# ==========================================================================
# FINAL VERIFICATION
# ==========================================================================
cat("\n--- Final Verification of EER data coverage ---\n")

# Check the final coverage percentage
final_coverage_pct <- round(100 * sum(!is.na(analysis_data$eer_kcal_marco_mean)) / nrow(analysis_data), 1)
cat("Final EER Coverage:", final_coverage_pct, "%\n")

# Count remaining NAs (should be 0)
remaining_nas <- sum(is.na(analysis_data$eer_kcal_marco_mean))
cat("Number of remaining NAs for EER mean:", remaining_nas, "\n")
if (remaining_nas == 0) {
  cat("✅ SUCCESS: All missing EER values have been imputed.\n")
}

# Display the rows that were imputed to confirm the logic worked
cat("\n--- Review of imputed rows: ---\n")
analysis_data %>%
  filter(!is.na(eer_imputed_source)) %>%
  select(iso3, sex, age_group, eer_kcal_marco_mean, eer_imputed_source) %>%
  print(n=50)

# ==========================================================================
#"OPTIMAL CALORIES" SCENARIO - CALCULATION & RESULTS
# ==========================================================================

# --- Step 3: Calculate De-standardized Protein Intake for the Optimal Scenario ---
cat("\n--- Calculating de-standardized protein intake for the optimal scenario... ---\n")

optimal_calories_data <- analysis_data %>%
  # The core assumption of this scenario: Caloric Intake = MEAN EER
  mutate(
    kcal_consumed_optimal = eer_kcal_marco_mean, 
    # ----------------------
    
    # De-standardize protein intake using this optimal caloric value
    protein_grams_optimal = (kcal_consumed_optimal * protein_kcal_share_mean) / 4
  )

cat("--- De-standardization complete. 'optimal_calories_data' created. ---\n")


# --- Step 4: Calculate Protein Inadequacy ---
# This uses the same robust functions from your previous work.
cat("\n--- Calculating prevalence of protein inadequacy... ---\n")

# Re-define the function just in case the session was cleared
calculate_inadequacy <- function(mean_intake, cv_intake, distribution_type, requirement) {
  if (is.na(requirement) || requirement <= 0 || is.na(mean_intake) || is.na(cv_intake) || mean_intake <= 0 || cv_intake <= 0) return(NA_real_)
  if (distribution_type == "gamma") {
    shape_k <- 1 / (cv_intake^2); if (is.infinite(shape_k)) return(NA_real_)
    scale_theta <- mean_intake * cv_intake^2
    pgamma(requirement, shape = shape_k, scale = scale_theta)
  } else if (distribution_type == "log-normal") {
    meanlog <- log(mean_intake) - 0.5 * log(1 + cv_intake^2); sdlog <- sqrt(log(1 + cv_intake^2))
    if(is.na(sdlog) || sdlog <= 0) return(NA_real_)
    plnorm(requirement, meanlog = meanlog, sdlog = sdlog)
  } else { NA_real_ }
}


#####REMOVE 0-1 age group from the data
#optimal_calories_data <- optimal_calories_data %>% dplyr::filter(age_group != "0-0.99")

# --- Step 4: Calculate Protein Inadequacy (EAR and OPTIMAL) ---
cat("\n--- Calculating prevalence of protein inadequacy (EAR & OPTIMAL)... ---\n")

optimal_calories_results <- optimal_calories_data %>%
  rowwise() %>%
  mutate(
    # (A) standard EAR-based inadequacy
    prevalence_inadequate_EAR = calculate_inadequacy(
      mean_intake = protein_grams_optimal,
      cv_intake = cv,
      distribution_type = best_dist,
      requirement = ear_mean_g_day
    ),
    # (B) "optimal" inadequacy per Stu thresholds
    prevalence_inadequate_OPT = calculate_inadequacy(
      mean_intake = protein_grams_optimal,
      cv_intake = cv,
      distribution_type = best_dist,
      requirement = opt_mean_g_day   # <-- new variable from script 1
    )
  ) %>%
  ungroup()

cat("--- Inadequacy calculation complete (EAR & OPTIMAL). ---\n")


# --- Step 5: Summarize and Display the Final Result ---
cat("\n--- Final Summary for 'Optimal Calories' Scenario ---\n")

global_summary_optimal <- optimal_calories_results %>%
  summarise(
    global_inad_EAR = weighted.mean(prevalence_inadequate_EAR, w = population, na.rm = TRUE),
    global_inad_OPT = weighted.mean(prevalence_inadequate_OPT, w = population, na.rm = TRUE)
  )

summary_by_sex_optimal <- optimal_calories_results %>%
  group_by(sex) %>%
  summarise(
    inad_EAR = weighted.mean(prevalence_inadequate_EAR, w = population, na.rm = TRUE),
    inad_OPT = weighted.mean(prevalence_inadequate_OPT, w = population, na.rm = TRUE)
  ) %>%
  mutate(
    EAR_pct = scales::percent(inad_EAR, accuracy = 0.1),
    OPT_pct = scales::percent(inad_OPT, accuracy = 0.1)
  )

cat("\n=======================================================================\n")
cat("          FINAL RESULT: OPTIMAL CALORIES SCENARIO\n")
cat("=======================================================================\n")
cat("Global prevalence of protein inadequacy:\n")
cat(">>>  EAR-based: ", scales::percent(global_summary_optimal$global_inad_EAR, accuracy = 0.1), "\n")
cat(">>>  OPTIMAL-based (Stu): ", scales::percent(global_summary_optimal$global_inad_OPT, accuracy = 0.1), "\n")
cat("=======================================================================\n\n")

cat("--- Breakdown by Sex ---\n")
print(summary_by_sex_optimal)





# ==========================================================================
# FINAL STEP OF SCRIPT 2: SAVE THE COMPLETE RESULTS
# ==========================================================================

# The 'optimal_calories_results' object contains EVERYTHING from this script:
#  - The original 'prepped_data' (GDD intakes, body weights, EARs, etc.)
#  - The population data
#  - Marco's EER data (cleaned and imputed)
#  - The calculated inadequacy for the "Optimal Calories" scenario

# This is the single, complete object we need for future work.

cat("\n--- Saving the complete 'optimal_calories_results' object... ---\n")

saveRDS(optimal_calories_results, file = "./output/script2_final_results.rds")

cat("--- Data saved successfully to ./output/script2_final_results.rds ---\n")

