##load needed elements from script 1,2, 3
#calculate protein share of caloric intake by de-standardizing GDD data
#calculate absolute protein levels by multiplying protein share for actual calories
#calculate absolute protein inadequacy caloric inadequacy and check overlap

#new approach: aim for having prevalence of undernutrition around 8.2%
# this is the current FAO estimate. we can tune our waste on that, and
# see what's the prevalence of protein inadequacy at these level. This is 
# how we will tune our "med" waste scenario. 

#to estimate prevalence of undernutrition in our data, we will take
#FAO


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
calorie_data_full <- readRDS("./output/final_calorie_distributions_wide.rds")

cat("--- Protein and Calorie worlds have been summoned. ---\n")

# --- Force Script 3 calorie columns to have unique names before joining ---
calorie_data_full_s3 <- calorie_data_full %>%
  rename(
    kcal_mean_low_s3    = kcal_mean_low,
    kcal_mean_medium_s3 = kcal_mean_medium,
    kcal_mean_high_s3   = kcal_mean_high,
    best_dist_calorie_s3 = best_dist,
    cv_calorie_s3        = cv
  ) %>%
  select(iso3, sex, age_group,
         kcal_mean_low_s3, kcal_mean_medium_s3, kcal_mean_high_s3,
         best_dist_calorie_s3, cv_calorie_s3)

# --- Join protein world + Script 3 calorie world (renamed) ---
final_analysis_data <- protein_data_full %>%
  left_join(calorie_data_full_s3, by = c("iso3", "sex", "age_group"))

# --- Make the canonical calorie columns come from Script 3 ---
final_analysis_data <- final_analysis_data %>%
  mutate(
    kcal_mean_low    = kcal_mean_low_s3,
    kcal_mean_medium = kcal_mean_medium_s3,
    kcal_mean_high   = kcal_mean_high_s3,
    best_dist_calorie = best_dist_calorie_s3,
    cv_calorie        = cv_calorie_s3
  )



stopifnot(!any(is.na(final_analysis_data$kcal_mean_low)))
stopifnot(!any(is.na(final_analysis_data$kcal_mean_high)))
stopifnot(!any(is.na(final_analysis_data$best_dist_calorie)))
stopifnot(!any(is.na(final_analysis_data$cv_calorie)))


# --- Step 2: Create the Final, Unified Analytical Dataset ---
# We join the two worlds together on their common keys (iso3, sex, age_group).
# We add suffixes to prevent any confusion with duplicate column names.

# --- Make canonical PROTEIN columns (your protein side still uses generic names) ---
# In script2_final_results.rds, protein distribution columns are named: best_dist, cv
# Your downstream code expects: best_dist_protein, cv_protein
final_analysis_data <- final_analysis_data %>%
  mutate(
    best_dist_protein = best_dist,
    cv_protein        = cv
  )

stopifnot(!any(is.na(final_analysis_data$best_dist_protein)))
stopifnot(!any(is.na(final_analysis_data$cv_protein)))




# --- Sanity: make sure Stu's optimal thresholds are present (g/day) ---
# They come from Script 1 -> Script 2 and should now be in the protein side.
# Depending on column collisions, they might be either bare or suffixed with _protein.
if (!("opt_mean_g_day" %in% names(final_analysis_data))) {
  if ("opt_mean_g_day_protein" %in% names(final_analysis_data)) {
    final_analysis_data <- final_analysis_data %>%
      mutate(opt_mean_g_day = opt_mean_g_day_protein)
  }
}

if (!("ear_mean_g_day" %in% names(final_analysis_data))) {
  if ("ear_mean_g_day_protein" %in% names(final_analysis_data)) {
    final_analysis_data <- final_analysis_data %>%
      mutate(ear_mean_g_day = ear_mean_g_day_protein)
  }
}

# Minimal assert so we fail loudly if something’s off
stopifnot("opt_mean_g_day" %in% names(final_analysis_data))
stopifnot("ear_mean_g_day" %in% names(final_analysis_data))


# ==========================================================================
# PART 3: LOAD AND PREPARE MINIMUM DIETARY ENERGY REQUIREMENT (MDER) DATA
# ==========================================================================

cat("\n--- Loading FAO Minimum Dietary Energy Requirement (MDER) data... ---\n")

# --- Step 3a: Load the raw data from OWID ---
owid_mder_raw <- read_csv("./data/OWID_minimum-requirement-calories.csv", show_col_types = FALSE)

# --- Step 3b: Clean and select the 2018 data ---
owid_mder <- owid_mder_raw %>%
  filter(Year == 2018) %>%
  # Select and rename the column for clarity
  select(
    iso3 = Code, 
    mder_kcal = `Minimum dietary energy requirement  (kcal/cap/day) | 00021056 || Value | 006128 || kilocalories per person per day`
  ) %>%
  filter(!is.na(iso3) & nchar(iso3) == 3)

# --- Step 3c: Identify which countries need imputation ---
# THE FIX: First, create the 'all_countries_needed' object from our master dataset.
all_countries_needed <- final_analysis_data %>% 
  distinct(iso3)

# Now, use this list to check for missing countries in the MDER data.
countries_missing_mder <- all_countries_needed %>%
  anti_join(owid_mder, by = "iso3")

# --- Step 3d: Report the findings ---
if (nrow(countries_missing_mder) > 0) {
  cat("\n⚠️ WARNING:", nrow(countries_missing_mder), "countries in our analysis are MISSING from the MDER dataset and will require imputation:\n")
  print(countries_missing_mder)
} else {
  cat("\n✅ SUCCESS! All countries required for our analysis have MDER data.\n")
}

# ==========================================================================
# PART 3 (Continued): IMPUTE AND FINALIZE MDER DATA
# ==========================================================================

# --- Step 3e: Impute the single missing country (MHL) using its peer (FJI) ---

cat("\n--- Imputing MDER for the Marshall Islands (MHL) using Fiji (FJI) as a proxy... ---\n")

# First, get the MDER value for our proxy country, Fiji
fiji_mder_value <- owid_mder %>%
  filter(iso3 == "FJI") %>%
  pull(mder_kcal) # pull() extracts the single value

# Create a small tibble for the imputed row
mhl_imputed_row <- tibble(
  iso3 = "MHL",
  mder_kcal = fiji_mder_value
)

# Bind this imputed row to the main MDER dataset
mder_final_complete <- bind_rows(owid_mder, mhl_imputed_row)


# --- Step 3g: Join the MDER data to our main analytical dataset ---
# Now we add the 'mder_kcal' column to our main data frame.
final_analysis_data <- final_analysis_data %>%
  left_join(mder_final_complete, by = "iso3")

cat("\n--- MDER data successfully joined to the final analysis dataset. ---\n")
glimpse(final_analysis_data)



# ==========================================================================
# PART 4: DISAGGREGATE NATIONAL MDER TO STRATUM-SPECIFIC VALUES
# ==========================================================================

# --- Step 4a: Create a temporary data frame for the disaggregation ---
# This joins the national MDER to the master data which contains the EERs and population.
# The 'master_data' object was loaded at the beginning of this script.
mder_disaggregation_data <- final_analysis_data %>%
  # We only need a few columns for this calculation
  select(iso3, sex, age_group, population, eer_kcal_marco_mean) %>%
  # Join the national-level MDER value
  left_join(mder_final_complete, by = "iso3")

# --- Step 4b: Perform the MDER Disaggregation ---
# The logic is identical to how we disaggregated caloric supply.
cat("\n--- Disaggregating national MDER to stratum-specific values... ---\n")

mder_stratum_specific <- mder_disaggregation_data %>%
  # For each country...
  group_by(iso3) %>%
  
  # Calculate the population-weighted average EER for that country.
  # This serves as our reference point for the "average" person's need.
  mutate(
    avg_country_eer = weighted.mean(eer_kcal_marco_mean, w = population, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  
  # Now, calculate the adjustment factor for each stratum based on its EER.
  mutate(
    adjustment_factor = eer_kcal_marco_mean / avg_country_eer,
    
    # Finally, apply this factor to the national MDER to get the
    # final, stratum-specific minimum energy requirement.
    mder_stratum_kcal = mder_kcal * adjustment_factor
  ) %>%
  
  # Select only the columns we need to join back to our main data
  select(iso3, sex, age_group, mder_stratum_kcal)

cat("--- MDER disaggregation complete. ---\n")

# after you create mder_stratum_specific
mder_stratum_specific <- mder_stratum_specific %>%
  distinct(iso3, sex, age_group, .keep_all = TRUE)   # add year/scenario if they’re real keys

# re-run the reconstruction check with RIGHT side distinct
check_mder <- mder_disaggregation_data %>%
  left_join(mder_stratum_specific, by = c("iso3","sex","age_group")) %>%
  filter(!is.na(mder_stratum_kcal), population > 0) %>%
  group_by(iso3) %>%  # add year if applicable
  summarise(
    mder_nat   = first(mder_kcal),
    mder_recon = sum(mder_stratum_kcal * population) / sum(population),
    diff       = mder_recon - mder_nat,
    .groups = "drop"
  )
stopifnot(max(abs(check_mder$diff), na.rm = TRUE) < 1e-6)

# now the join back won’t warn
final_analysis_data <- final_analysis_data %>%
  left_join(mder_stratum_specific, by = c("iso3","sex","age_group"))

#####REMOVE 0-1 age group from the data
final_analysis_data <- final_analysis_data %>% dplyr::filter(age_group != "0-0.99")

# --- Verification ---
cat("\n--- Verifying the disaggregated MDER values for the USA ---\n")
final_analysis_data %>%
  filter(iso3 == "USA") %>%
  select(sex, age_group, `Stratum MDER` = mder_stratum_kcal) %>%
  arrange(desc(`Stratum MDER`)) %>%
  print(n=100)

glimpse(final_analysis_data)


# ==========================================================================
# FINAL CALCULATIONS
# ==========================================================================

##OPTIONAL: cap CVs at some values
# final_analysis_data <- final_analysis_data %>%
#   mutate(cv_calorie = pmin(cv_calorie, 0.35))   # cap extreme outliers
# final_analysis_data <- final_analysis_data %>%
#   mutate(cv_protein = pmin(cv_calorie, 0.45))   # cap extreme outliers



# --- Step 3: Calculate "True" Absolute Protein Intake ---
final_results <- final_analysis_data %>%
  mutate(
    # Protein provides 4 kcal/gram
    protein_grams_true = (kcal_mean_medium * protein_kcal_share_mean) / 4
  )

# --- Step 4: Calculate Final Protein & Calorie Inadequacy (EAR + OPT) ---
calculate_inadequacy <- function(mean_intake, cv_intake, distribution_type, requirement) {
  if (is.na(requirement) || requirement <= 0 || is.na(mean_intake) || is.na(cv_intake) || mean_intake <= 0 || cv_intake <= 0) return(NA_real_)
  if (distribution_type == "gamma") {
    shape_k <- 1 / (cv_intake^2); if (is.infinite(shape_k)) return(NA_real_)
    scale_theta <- mean_intake * cv_intake^2
    return(pgamma(requirement, shape = shape_k, scale = scale_theta))
  } else if (distribution_type == "log-normal") {
    meanlog <- log(mean_intake) - 0.5 * log(1 + cv_intake^2); sdlog <- sqrt(log(1 + cv_intake^2))
    if (is.na(sdlog) || sdlog <= 0) return(NA_real_)
    return(plnorm(requirement, meanlog = meanlog, sdlog = sdlog))
  } else {
    return(NA_real_)
  }
}

final_results <- final_analysis_data %>%
  mutate(
    # Protein provides 4 kcal/gram
    protein_grams_true = (kcal_mean_medium * protein_kcal_share_mean) / 4
  ) %>%
  rowwise() %>%
  mutate(
    # A) PROTEIN inadequacy — EAR (original)
    prevalence_protein_inadequate_EAR = calculate_inadequacy(
      mean_intake = protein_grams_true,
      cv_intake = cv_protein,
      distribution_type = best_dist_protein,
      requirement = ear_mean_g_day
    ),
    # B) PROTEIN inadequacy — OPT (Stu thresholds)
    prevalence_protein_inadequate_OPT = calculate_inadequacy(
      mean_intake = protein_grams_true,
      cv_intake = cv_protein,
      distribution_type = best_dist_protein,
      requirement = opt_mean_g_day
    ),
    # C) CALORIE inadequacy — MDER (realistic scenario)
    prevalence_calorie_inadequate = calculate_inadequacy(
      mean_intake = kcal_mean_medium,
      cv_intake = cv_calorie,
      distribution_type = best_dist_calorie,
      requirement = mder_stratum_kcal
    ),
    # D) CALORIE inadequacy — MDER under "energy-adequate" mean intake (EER scenario)
    prevalence_calorie_inadequate_EER = calculate_inadequacy(
      mean_intake = eer_kcal_marco_mean,
      cv_intake = cv_calorie,
      distribution_type = best_dist_calorie,
      requirement = mder_stratum_kcal
    )
  ) %>%
  ungroup()




# --- Step 5: Report The Grand Finale ---
# Calculate the final, population-weighted global averages.

global_summary_final <- final_results %>%
  summarise(
    global_protein_inadequacy_EAR = weighted.mean(prevalence_protein_inadequate_EAR, w = population, na.rm = TRUE),
    global_protein_inadequacy_OPT = weighted.mean(prevalence_protein_inadequate_OPT, w = population, na.rm = TRUE),
    global_calorie_inadequacy     = weighted.mean(prevalence_calorie_inadequate,     w = population, na.rm = TRUE),
    global_calorie_inad_EER       = weighted.mean(prevalence_calorie_inadequate_EER, w = population, na.rm = TRUE)
  )


# --- NEW: Global population-weighted averages by sex ---
global_summary_by_sex <- final_results %>%
  group_by(sex) %>%
  summarise(
    global_protein_inadequacy_EAR = weighted.mean(prevalence_protein_inadequate_EAR, w = population, na.rm = TRUE),
    global_protein_inadequacy_OPT = weighted.mean(prevalence_protein_inadequate_OPT, w = population, na.rm = TRUE),
    global_calorie_inadequacy     = weighted.mean(prevalence_calorie_inadequate,     w = population, na.rm = TRUE),
    global_calorie_inad_EER       = weighted.mean(prevalence_calorie_inadequate_EER, w = population, na.rm = TRUE),
    .groups = "drop"
  )



cat("\n\n======================================================\n")
cat("          THE GRAND FINALE: FINAL RESULTS\n")
cat("======================================================\n\n")
cat("Global Prevalence of Protein Inadequacy (True Intake):\n")
cat(">>>  EAR-based: ", percent(global_summary_final$global_protein_inadequacy_EAR, accuracy = 0.1), "\n")
cat(">>>  OPT-based: ", percent(global_summary_final$global_protein_inadequacy_OPT, accuracy = 0.1), "\n\n")
cat("Global Prevalence of Caloric Inadequacy:\n")
cat(">>>  Realistic (modeled intake): ", percent(global_summary_final$global_calorie_inadequacy, accuracy = 0.1), "\n")
cat(">>>  EER scenario (mean intake = EER): ", percent(global_summary_final$global_calorie_inad_EER, accuracy = 0.1), "\n\n")



# --- NEW: Print by-sex results (including EER-calorie counterfactual) ---
cat("\nBy-sex global prevalences (population-weighted):\n")

global_summary_by_sex_fmt <- global_summary_by_sex %>%
  mutate(
    global_protein_inadequacy_EAR = scales::percent(global_protein_inadequacy_EAR, accuracy = 0.1),
    global_protein_inadequacy_OPT = scales::percent(global_protein_inadequacy_OPT, accuracy = 0.1),
    global_calorie_inadequacy     = scales::percent(global_calorie_inadequacy,     accuracy = 0.1),
    global_calorie_inad_EER       = scales::percent(global_calorie_inad_EER,       accuracy = 0.1)
  ) %>%
  arrange(sex)

for (i in seq_len(nrow(global_summary_by_sex_fmt))) {
  cat("\nSex:", global_summary_by_sex_fmt$sex[i], "\n")
  cat("  Protein inadequacy (EAR): ", global_summary_by_sex_fmt$global_protein_inadequacy_EAR[i], "\n")
  cat("  Protein inadequacy (OPT): ", global_summary_by_sex_fmt$global_protein_inadequacy_OPT[i], "\n")
  cat("  Calorie inadequacy (realistic intake): ", global_summary_by_sex_fmt$global_calorie_inadequacy[i], "\n")
  cat("  Calorie inadequacy (EER scenario mean): ", global_summary_by_sex_fmt$global_calorie_inad_EER[i], "\n")
}
cat("\n")

# (Optional) Save a tidy CSV for later use
# write_csv(global_summary_by_sex, "./output/global_summary_by_sex.csv")


# --- (Optional) Save the final, most complete dataset ---
saveRDS(final_results, file = "./output/final_analysis_with_all_inadequacy.rds")
cat("\n--- Final results object saved successfully. The work is done. ---\n")


######


print(final_results %>%
  group_by(iso3) %>%
  summarise(
    mean_eer = mean(eer_kcal_marco_mean, na.rm = TRUE),
    mean_obs = mean(kcal_mean_medium, na.rm = TRUE),
    diff = mean_eer - mean_obs
  ) %>%
  filter(mean_eer > mean_obs) %>%
  arrange(desc(diff)), n = 100)


# --- NEW: EER vs observed by sex ---
# print(
#   final_results %>%
#     group_by(iso3, sex) %>%
#     summarise(
#       mean_eer = mean(eer_kcal_marco_mean, na.rm = TRUE),
#       mean_obs = mean(kcal_mean_medium,    na.rm = TRUE),
#       diff     = mean_eer - mean_obs,
#       .groups  = "drop"
#     ) %>%
#     filter(mean_eer > mean_obs) %>%
#     arrange(desc(diff)),
#   n = 100
# )


