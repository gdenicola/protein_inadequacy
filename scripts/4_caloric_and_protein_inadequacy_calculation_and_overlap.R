##load needed elements from script 1,2, 3
#calculate protein share of caloric intake by de-standardizing GDD data
#calculate absolute protein levels by multiplying protein share for actual calories
#calculate absolute protein inadequacy caloric inadequacy and check overlap


# Load required packages
library(dplyr)
library(readr)
library(nutriR)
library(ggplot2)
library(tidyr)
library(purrr)
library(countrycode)
library(scales) # For the percent() function


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
options(scipen=999)
rm(list = ls())

# ==========================================================================
# SCRIPT 4: THE GRAND FINALE - TRUE INADEQUACY ESTIMATES
#
# PURPOSE:
# 1. Loads the final "Protein World" and "Calorie World" datasets.
# 2. Joins them into a master analytical dataframe.
# 3. Calculates the "true" protein intake by de-standardizing it against
#    the final, modeled caloric intake.
# 4. Calculates the final prevalence of both PROTEIN and CALORIE inadequacy.
# 5. Reports the final, global headline numbers.
# ==========================================================================


# --- Step 1: Load the Two Essential Datasets ---
# We load the final, definitive outputs from our previous work.

cat("--- Loading the spoils of victory... ---\n")

# Load the "Protein World" from Script 2's final output.
# It contains everything needed for protein: shares, distribution shapes, and requirements (EARs).
protein_data_full <- readRDS("./output/script2_final_results.rds")

# Load the "Calorie World" from Script 3's final output.
# It contains our best estimate of actual caloric intake and its distribution shape.
calorie_data_full <- readRDS("./output/final_calorie_distributions.rds")

cat("--- Protein and Calorie worlds have been summoned. ---\n")


# --- Step 2: Create the Final, Unified Analytical Dataset ---
# We join the two worlds together on their common keys (iso3, sex, age_group).
# We add suffixes to prevent any confusion with duplicate column names.

final_analysis_data <- protein_data_full %>%
  left_join(
    calorie_data_full, 
    by = c("iso3", "sex", "age_group"),
    suffix = c("_protein", "_calorie") # e.g., 'cv_protein', 'cv_calorie'
  )

# --- Final Verification ---
# Glimpse the data and check for any NAs in the key columns we are about to use.
cat("\n--- The final, unified dataset is ready for analysis: ---\n")
glimpse(final_analysis_data)

key_cols <- c("kcal_mean", "protein_kcal_share_mean", "cv_protein", "ear_mean_g_day", "cv_calorie")
na_check <- sapply(final_analysis_data[key_cols], function(x) sum(is.na(x)))

if(any(na_check > 0)) {
  cat("\n⚠️ WARNING: NAs detected in key columns after the join. This should not happen.\n")
  print(na_check[na_check > 0])
} else {
  cat("\n✅ SUCCESS: The final dataset is complete and key columns have no missing values.\n")
}


# ==========================================================================
# FINAL CALCULATIONS
# ==========================================================================

# --- Step 3: Calculate "True" Absolute Protein Intake ---
# We de-standardize the GDD protein intake data using our final calorie estimates.

final_results <- final_analysis_data %>%
  mutate(
    # Protein provides 4 kcal/gram
    protein_grams_true = (kcal_mean * protein_kcal_share_mean) / 4
  )

# --- Step 4: Calculate Final Protein & Calorie Inadequacy ---
# First, define our trusted inadequacy calculation function.
calculate_inadequacy <- function(mean_intake, cv_intake, distribution_type, requirement) {
  if (is.na(requirement) || requirement <= 0 || is.na(mean_intake) || is.na(cv_intake) || mean_intake <= 0 || cv_intake <= 0) return(NA_real_)
  if (distribution_type == "gamma") {
    shape_k <- 1 / (cv_intake^2); if (is.infinite(shape_k)) return(NA_real_)
    scale_theta <- mean_intake * cv_intake^2
    return(pgamma(requirement, shape = shape_k, scale = scale_theta))
  } else if (distribution_type == "log-normal") {
    meanlog <- log(mean_intake) - 0.5 * log(1 + cv_intake^2); sdlog <- sqrt(log(1 + cv_intake^2))
    if(is.na(sdlog) || sdlog <= 0) return(NA_real_)
    return(plnorm(requirement, meanlog = meanlog, sdlog = sdlog))
  } else { return(NA_real_) }
}

# Now, apply the function for both protein and calories.
final_results <- final_results %>%
  rowwise() %>%
  mutate(
    # A) Calculate PROTEIN inadequacy
    prevalence_protein_inadequate = calculate_inadequacy(
      mean_intake = protein_grams_true,      # Our new "true" protein mean
      cv_intake = cv_protein,                # The original PROTEIN CV
      distribution_type = best_dist_protein,   # The original PROTEIN distribution
      requirement = ear_mean_g_day           # The absolute PROTEIN requirement (EAR)
    ),
    
    # B) Calculate CALORIE inadequacy
    prevalence_calorie_inadequate = calculate_inadequacy(
      mean_intake = kcal_mean,               # Our final CALORIE mean
      cv_intake = cv_calorie,                # The CALORIE CV we matched
      distribution_type = best_dist_calorie, # The CALORIE distribution we matched
      requirement = eer_kcal_marco_mean    # The calorie requirement is the EER
    )
  ) %>%
  ungroup()


# --- Step 5: Report The Grand Finale ---
# Calculate the final, population-weighted global averages.

global_summary_final <- final_results %>%
  summarise(
    global_protein_inadequacy = weighted.mean(prevalence_protein_inadequate, w = population, na.rm = TRUE),
    global_calorie_inadequacy = weighted.mean(prevalence_calorie_inadequate, w = population, na.rm = TRUE)
  )

cat("\n\n======================================================\n")
cat("          THE GRAND FINALE: FINAL RESULTS\n")
cat("======================================================\n\n")
cat("Global Prevalence of Protein Inadequacy (True Intake):\n")
cat(">>> ", percent(global_summary_final$global_protein_inadequacy, accuracy = 0.1), "\n\n")
cat("Global Prevalence of Caloric Inadequacy (True Intake):\n")
cat(">>> ", percent(global_summary_final$global_calorie_inadequacy, accuracy = 0.1), "\n\n")
cat("======================================================\n")


# --- (Optional) Save the final, most complete dataset ---
saveRDS(final_results, file = "./output/final_analysis_with_all_inadequacy.rds")
cat("\n--- Final results object saved successfully. The work is done. ---\n")
