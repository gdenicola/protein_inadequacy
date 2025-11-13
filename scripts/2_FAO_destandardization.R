
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


# --- Load Marco's Country-Specific EER Data ---
cat("\n--- Loading Marco's EER data from data_extract_040725.xlsx... ---\n")

# Load the data from the specified sheet
library(readxl)
marco_eer_raw <- read_excel("./data/data_extract_040725.xlsx", sheet = 1)
glimpse(marco_eer_raw)
summary(marco_eer_raw)
unique(marco_eer_raw$Region)
unique(marco_eer_raw$Sex)
unique(marco_eer_raw$Age)
unique(marco_eer_raw$Year)
unique(marco_eer_raw$Stats)
unique(marco_eer_raw$Value)
# --------------------------------------------------------------------------
# Using the Our World in Data (OWID) Caloric Supply File
# --------------------------------------------------------------------------

# --- Step 1: Load the OWID data file ---
owid_raw <- read_csv("./data/OWID_daily-per-capita-caloric-supply.csv", show_col_types = FALSE)
#from https://ourworldindata.org/grapher/daily-per-capita-caloric-supply
owid_protein_raw <- read_csv("./data/OWID_daily-per-capita-protein-supply.csv", show_col_types = FALSE)
#from https://ourworldindata.org/grapher/daily-per-capita-protein-supply
glimpse(owid_protein_raw)


# --- Step 2: Filter and clean the OWID data ---
# Based on the glimpse, the column names are simple and direct.
owid_calories <- owid_raw %>%
  filter(Year == 2018) %>%
  # Select the 'Code' (for iso3) and the calorie column, giving it a standard name
  select(
    iso3 = Code,
    fao_des_kcal = `Daily calorie supply per person` # Use this name to match later code if needed
  ) %>%
  # OWID files sometimes include continents or income groups which lack an iso3 'Code'
  # or have non-standard codes. This filter removes them.
  filter(!is.na(iso3) & nchar(iso3) == 3)


# --- Step 3: Final Validation ---
# The moment of truth: check this new, clean data against our list of required countries.

# Check for any countries in our main dataset that are missing from the OWID data
missing_from_owid <- prepped_data %>%
  distinct(iso3) %>%
  anti_join(owid_calories, by = "iso3")

cat("\n--- Data Join Validation (OWID vs. GDD) ---\n")
if (nrow(missing_from_owid) > 0) {
  cat("\n⚠️ WARNING: The following", nrow(missing_from_owid), "countries from your main dataset are STILL missing from the OWID data:\n")
  print(missing_from_owid$iso3)
} else {
  cat("\n✅✅✅ SUCCESS! All countries in your main dataset were successfully found in the OWID data.\n")
}

# --------------------------------------------------------------------------
# FINAL STEP: Targeted Imputation for the 12 Remaining Countries
# --------------------------------------------------------------------------

# We have high-quality OWID data for most countries.
# For the last 12, we will use a transparent imputation strategy.

# Step 1: Define the proxy relationships for the remaining countries
imputation_map_final <- tibble::tribble(
  ~iso3, ~proxy_iso3, ~reason,
  "BHR", "SAU",       "Regional peer (GCC)",
  "QAT", "SAU",       "Regional peer (GCC)",
  "SGP", "KOR",       "High-income Asian peer",
  "BRN", "MYS",       "Regional peer",
  "PSE", "JOR",       "Regional peer",
  "ERI", "ETH",       "Regional peer",
  "SSD", "SDN",       "Direct neighbor (if SDN is now available, otherwise use ETH)",
  "GNQ", "GAB",       "Regional peer",
  "BTN", "NPL",       "Regional peer (Himalayan)",
  "TON", "FJI",       "Regional peer (Pacific)",
  "FSM", "FJI",       "Regional peer (Pacific)",
  "MHL", "FJI",       "Regional peer (Pacific)"
)

# Step 2: Get the caloric values for the proxy countries from our clean OWID data
imputation_values <- owid_calories %>%
  filter(iso3 %in% imputation_map_final$proxy_iso3) %>%
  select(proxy_iso3 = iso3, fao_des_kcal)

# Step 3: Create the final imputation table
imputation_table_final <- imputation_map_final %>%
  left_join(imputation_values, by = "proxy_iso3", relationship = "many-to-one") %>%
  # Keep only the rows for countries that are ACTUALLY missing
  filter(iso3 %in% missing_from_owid$iso3) %>%
  select(iso3, fao_des_kcal) # We only need these two columns to bind

# Step 4: Bind the imputed rows to the main OWID data to create the complete dataset
calories_final_complete <- bind_rows(owid_calories, imputation_table_final)


# We will now use 'calories_final_complete' as our definitive source.
cat("\n--- Final Complete Caloric Supply Data (Head) ---\n")
print(head(calories_final_complete))


# --------------------------------------------------------------------------
# PART 3, STEP 1: Define the EER Reference Table and Its Derivation
# --------------------------------------------------------------------------

# The following 'eer_table' provides pre-calculated Estimated Energy Requirements
# (EERs) in kcal/day for each age-sex stratum. These values are used as
# relative weights to disaggregate national per-capita caloric supply into
# stratum-specific estimates.

# ---
# SOURCE & DERIVATION FOR PUBLICATION:
#
# The EER values are derived from the predictive equations published by the
# Institute of Medicine (IOM) of the National Academies.
#
# Primary Source:
# Institute of Medicine. (2005). Dietary Reference Intakes for Energy,
# Carbohydrate, Fiber, Fat, Fatty Acids, Cholesterol, Protein, and Amino Acids.
# Washington, DC: The National Academies Press. https://doi.org/10.17226/10490.
#
# Methodology:
# The values were calculated using the equations specified in Chapter 5,
# Table 5-2 ("Estimated Energy Requirement (EER) Equations for All Life-Stage
# Groups"). A consistent Physical Activity Level (PAL) of "Low Active" was
# assumed for all groups ages 3 and older, using the Physical Activity (PA)
# coefficients from Table 5-3 of the same source. Standard reference body
# weights and heights for each age-sex group were used as inputs for the
# equations, and the midpoint age for each stratum was used for the 'age'
# variable in the formulas.
# ---

eer_table <- tibble::tribble(
  ~age_group, ~sex,     ~eer_kcal,
  # --- Infants & Children ---
  "0-0.99",   "Males",   850,
  "0-0.99",   "Females", 850,
  "1-4",      "Males",   1300,
  "1-4",      "Females", 1250,
  "5-9",      "Males",   1750,
  "5-9",      "Females", 1650,
  # --- Adolescents ---
  "10-14",    "Males",   2300,
  "10-14",    "Females", 2000,
  "15-19",    "Males",   2900,
  "15-19",    "Females", 2200,
  # --- Adults 19+ ---
  "20-24",    "Males",   2800, "20-24",    "Females", 2150,
  "25-29",    "Males",   2800, "25-29",    "Females", 2150,
  "30-34",    "Males",   2800, "30-34",    "Females", 2150,
  "35-39",    "Males",   2800, "35-39",    "Females", 2150,
  "40-44",    "Males",   2700, "40-44",    "Females", 2050,
  "45-49",    "Males",   2700, "45-49",    "Females", 2050,
  "50-54",    "Males",   2600, "50-54",    "Females", 1950,
  "55-59",    "Males",   2600, "55-59",    "Females", 1950,
  "60-64",    "Males",   2400, "60-64",    "Females", 1850,
  "65-69",    "Males",   2400, "65-69",    "Females", 1850,
  "70-74",    "Males",   2200, "70-74",    "Females", 1750,
  "75-79",    "Males",   2200, "75-79",    "Females", 1750,
  "80-84",    "Males",   2000, "80-84",    "Females", 1650,
  "85-89",    "Males",   2000, "85-89",    "Females", 1650,
  "90-94",    "Males",   1800, "90-94",    "Females", 1550,
  "95-99",    "Males",   1800, "95-99",    "Females", 1550
)

# --------------------------------------------------------------------------
# PART 3, STEP 2: Create the Master Analytical Data Frame
# --------------------------------------------------------------------------

# As requested, this step combines our four primary data sources into a single,
# comprehensive data frame named 'master_data'.

cat("\n--- Joining all data sources into the master data frame... ---\n")

master_data <- prepped_data %>%
  # Join stratum-specific population counts
  left_join(population_data, by = c("iso3", "sex", "age_group")) %>%
  # Join country-level caloric supply
  left_join(calories_final_complete, by = "iso3") %>%
  # Join stratum-specific energy requirements (EER)
  left_join(eer_table, by = c("sex", "age_group"))


# --------------------------------------------------------------------------
# PART 3, STEP 3: Disaggregation and Final Calculation
# --------------------------------------------------------------------------

# We will now use the complete 'master_data' frame to perform the disaggregation
# and calculate the final, de-standardized protein intake estimates.

# --- Step 3a: Disaggregate National Caloric Supply ---
# This block calculates the weighted average EER for each country and then
# determines the stratum-specific caloric supply.

disaggregated_data <- master_data %>%
  # For each country...
  group_by(iso3) %>%
  # Calculate the population-weighted average EER for that country
  mutate(
    avg_country_eer = weighted.mean(eer_kcal, w = population, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  # Now, calculate the stratum-specific caloric intake based on the adjustment factor
  mutate(
    adjustment_factor = eer_kcal / avg_country_eer,
    # This is our "High" scenario (Supply = Consumption, 0% waste)
    disagg_kcal_supply = fao_des_kcal * adjustment_factor
  )

# the disaggregation was here made by using the EER table inputted above with population
# by stratum as weights, such that the sum of all still equals the same as by
# multiplying the countrywide average by the country population, thus preserving
# total claoric supply.

# --------------------------------------------------------------------------
# --- Step 3b (REVISED & RENAMED): Calculate "True" Protein Intake for all 9 Scenarios ---
# --------------------------------------------------------------------------
# This block calculates the de-standardized protein intake (in grams/day) for
# all 9 combinations of consumption scenarios (Low/Med/High waste) and
# dietary composition scenarios (Low/Med/High protein share).

final_data <- disaggregated_data %>%
  mutate(
    # --- First, create the 3 caloric consumption scenarios (adjusting for waste) ---
    est_kcal_consumed_high = disagg_kcal_supply * 0.90, # 10% waste
    est_kcal_consumed_med  = disagg_kcal_supply * 0.75, # 25% waste (Median/Central)
    est_kcal_consumed_low  = disagg_kcal_supply * 0.60, # 40% waste
    
    # --- Now, calculate all 9 protein intake scenarios (in grams/day) ---
    # Using the clearer "low, med, high" naming convention.
    
    # Central Estimate (Median Consumption x Median Protein Share)
    protein_g_med_med = (est_kcal_consumed_med * protein_kcal_share_mean) / 4, # Note: protein_kcal_share_mean is our "median" share
    
    # Sensitivity Scenarios
    protein_g_low_low   = (est_kcal_consumed_low  * protein_kcal_share_lower) / 4,
    protein_g_low_med   = (est_kcal_consumed_low  * protein_kcal_share_mean)  / 4,
    protein_g_low_high  = (est_kcal_consumed_low  * protein_kcal_share_upper) / 4,
    
    protein_g_med_low   = (est_kcal_consumed_med  * protein_kcal_share_lower) / 4,
    # protein_g_med_med is already defined above
    protein_g_med_high  = (est_kcal_consumed_med  * protein_kcal_share_upper) / 4,
    
    protein_g_high_low  = (est_kcal_consumed_high * protein_kcal_share_lower) / 4,
    protein_g_high_med  = (est_kcal_consumed_high * protein_kcal_share_mean)  / 4,
    protein_g_high_high = (est_kcal_consumed_high * protein_kcal_share_upper) / 4
  )

# --- Checkpoint: Inspect the new 'final_data' frame ---
cat("\n--- Data frame now contains all 9 protein intake scenarios (using 'med' naming) ---\n")
final_data %>%
  # Select a few key columns to see the results
  select(iso3, sex, age_group, protein_g_low_low, protein_g_med_med, protein_g_high_high) %>%
  head()


# --------------------------------------------------------------------------
# PART 5: FINAL INADEQUACY CALCULATION AND SENSITIVITY ANALYSIS (CORRECTED)
# --------------------------------------------------------------------------

# --- Step 5a: Define the inadequacy calculation function ---
calculate_inadequacy <- function(mean_intake, cv_intake, distribution_type, requirement) {
  if (is.na(requirement) || requirement <= 0 || is.na(mean_intake) || is.na(cv_intake) || mean_intake <= 0 || cv_intake <= 0) return(NA_real_)
  if (distribution_type == "gamma") {
    shape_k <- 1 / (cv_intake^2)
    if (shape_k == Inf) return(NA_real_) # Added check for CV near zero
    scale_theta <- mean_intake * cv_intake^2
    pgamma(requirement, shape = shape_k, scale = scale_theta)
  } else if (distribution_type == "log-normal") {
    if (1 + cv_intake^2 <= 0) return(NA_real_)
    meanlog <- log(mean_intake) - 0.5 * log(1 + cv_intake^2)
    sdlog <- sqrt(log(1 + cv_intake^2))
    if(is.na(sdlog) || sdlog <= 0) return(NA_real_)
    plnorm(requirement, meanlog = meanlog, sdlog = sdlog)
  } else { NA_real_ }
}

# --- Step 5b: Calculate inadequacy for all 9 scenarios ---
# This is a cleaner approach that creates 9 new columns directly.

cat("\n--- Calculating inadequacy for all 9 scenarios... ---\n")

inadequacy_results <- final_data %>%
  rowwise() %>%
  mutate(
    # Central estimate
    prev_med_med = calculate_inadequacy(protein_g_med_med, cv, best_dist, ear_mean_g_day),
    
    # Full sensitivity matrix
    prev_low_low = calculate_inadequacy(protein_g_low_low, cv, best_dist, ear_low_g_day),
    prev_low_med = calculate_inadequacy(protein_g_low_med, cv, best_dist, ear_mean_g_day),
    prev_low_high = calculate_inadequacy(protein_g_low_high, cv, best_dist, ear_high_g_day),
    
    prev_med_low = calculate_inadequacy(protein_g_med_low, cv, best_dist, ear_low_g_day),
    prev_med_high = calculate_inadequacy(protein_g_med_high, cv, best_dist, ear_high_g_day),
    
    prev_high_low = calculate_inadequacy(protein_g_high_low, cv, best_dist, ear_low_g_day),
    prev_high_med = calculate_inadequacy(protein_g_high_med, cv, best_dist, ear_mean_g_day),
    prev_high_high = calculate_inadequacy(protein_g_high_high, cv, best_dist, ear_high_g_day)
  ) %>%
  ungroup()

# --- Step 5c: Summarize the results ---

# Create a long data frame for easy summarization
prevalence_long <- inadequacy_results %>%
  select(iso3, sex, age_group, population, starts_with("prev_")) %>%
  pivot_longer(
    cols = starts_with("prev_"),
    names_to = "scenario_label",
    names_prefix = "prev_",
    values_to = "prevalence"
  )

# Calculate global prevalence for each of the 9 scenarios
global_prevalence_per_scenario <- prevalence_long %>%
  group_by(scenario_label) %>%
  summarise(
    global_prevalence = weighted.mean(prevalence, w = population, na.rm = TRUE)
  )

cat("\n--- Global Prevalence for Each of the 9 Scenarios ---\n")
print(global_prevalence_per_scenario)

# Now, find the overall central, min, and max estimates
final_summary <- global_prevalence_per_scenario %>%
  summarise(
    # Use a direct filter to get the single value for the central estimate
    central_estimate_prevalence = global_prevalence[scenario_label == "med_med"],
    min_estimate_prevalence = min(global_prevalence, na.rm = TRUE),
    max_estimate_prevalence = max(global_prevalence, na.rm = TRUE)
  )

cat("\n--- FINAL GLOBAL SUMMARY (CENTRAL, MIN, MAX) ---\n")
final_summary %>%
  mutate(across(everything(), ~ . * 100)) %>%
  rename_with(~ paste0(., " (%)")) %>%
  print()


# --- Step 5d: Summarize by Sex ---

final_summary_by_sex <- prevalence_long %>%
  group_by(sex, scenario_label) %>%
  summarise(
    sex_prevalence = weighted.mean(prevalence, w = population, na.rm = TRUE),
    .groups = "drop" # Drop grouping by scenario_label
  ) %>%
  group_by(sex) %>% # Re-group just by sex
  summarise(
    central_estimate_prevalence = sex_prevalence[scenario_label == "med_med"],
    min_estimate_prevalence = min(sex_prevalence, na.rm = TRUE),
    max_estimate_prevalence = max(sex_prevalence, na.rm = TRUE)
  )

cat("\n--- FINAL SUMMARY BY SEX (CENTRAL, MIN, MAX) ---\n")
final_summary_by_sex %>%
  mutate(across(-sex, ~ . * 100)) %>%
  rename_with(~ ifelse(. == "sex", ., paste0(., " (%)"))) %>%
  print()


# --------------------------------------------------------------------------
# PART 6: ANALYZE RESULTS BY AGE GROUP
# --------------------------------------------------------------------------

# We will now summarize the prevalence of inadequacy for each age group to see
# how risk varies across the lifespan.

# --- Step 6a: Calculate prevalence by age group for each scenario ---
# We use the 'prevalence_long' data frame we already created.

summary_by_age <- prevalence_long %>%
  # Group by age group and scenario
  group_by(age_group, scenario_label) %>%
  # Calculate the population-weighted average prevalence for that global age group
  summarise(
    age_group_prevalence = weighted.mean(prevalence, w = population, na.rm = TRUE),
    .groups = "drop" # Drop the scenario grouping
  )

# --- Step 6b: Find the Central, Min, and Max for each age group ---

final_summary_by_age <- summary_by_age %>%
  # Now group just by age group
  group_by(age_group) %>%
  summarise(
    # Extract the central "med_med" estimate
    central_estimate_prevalence = age_group_prevalence[scenario_label == "med_med"],
    # Find the absolute min and max across all 9 scenarios
    min_estimate_prevalence = min(age_group_prevalence, na.rm = TRUE),
    max_estimate_prevalence = max(age_group_prevalence, na.rm = TRUE)
  )

cat("\n--- FINAL SUMMARY BY AGE GROUP (CENTRAL, MIN, MAX) ---\n")
final_summary_by_age %>%
  mutate(across(-age_group, ~ . * 100)) %>%
  rename_with(~ paste0(., " (%)")) %>%
  print(n = 30)


# --- Step 6c: Visualize the results by age group ---

# Order the age groups for the plot's x-axis
age_group_levels_ordered <- c("0-0.99", "1-4", "5-9", "10-14", "15-19",
                              "20-24", "25-29", "30-34", "35-39",
                              "40-44", "45-49", "50-54", "55-59",
                              "60-64", "65-69", "70-74", "75-79",
                              "80-84", "85-89", "90-94", "95-99")

plot_inadequacy_by_age <- final_summary_by_age %>%
  mutate(age_group_f = factor(age_group, levels = age_group_levels_ordered)) %>%
  ggplot(aes(x = age_group_f, y = central_estimate_prevalence, group = 1)) +
  # Add a ribbon to show the full uncertainty range (min to max)
  geom_ribbon(aes(ymin = min_estimate_prevalence, ymax = max_estimate_prevalence),
              fill = "skyblue", alpha = 0.5) +
  # Plot the central estimate line on top
  geom_line(aes(y = central_estimate_prevalence), color = "dodgerblue4", linewidth = 1.2) +
  geom_point(aes(y = central_estimate_prevalence), color = "dodgerblue4", size = 3) +
  # Format the y-axis as percentages
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Global Prevalence of Inadequate Protein Intake by Age Group",
    subtitle = "Central estimate with full uncertainty range from 9 sensitivity scenarios",
    x = "Age Group",
    y = "Prevalence of Inadequacy"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_inadequacy_by_age)

#Over the years, attention to food waste has increased because
#one-third of food produced for human consumption is
#lost or wasted globally. From Gustavsson J, Cederberg C, Sonesson U, et al. 
#Global Food Losses and Food Waste. Extent, Causes and Prevention.
# 2011. Available at: http://www.fao.org/3/a-i2697e.pdf.


#low-med is prolly king
#next step: calculate the percentage of caloric inadequacy under the low-med scenario
#see if it's more or less than 20% (that's protein inadequacy)












####viz of relative protein share by sex, age, countries
#the 0-0.99 group looks funky for some countries

# library(ggplot2)
# library(viridis) # For nice color palettes
# 
# # Plot 1: Overall Distribution of Protein Calorie Share
# ggplot(prepped_data, aes(x = protein_kcal_share_mean)) +
#   geom_histogram(aes(y = ..density..), binwidth = 0.01, fill = "skyblue", color = "white", alpha = 0.8) +
#   geom_density(color = "dodgerblue4", linewidth = 1.2) +
#   scale_x_continuous(labels = scales::percent_format()) +
#   labs(
#     title = "Distribution of Protein's Share of Calories (from GDD data)",
#     subtitle = "Across all 7,770 country-sex-age strata",
#     x = "Protein Share of Standardized Calories",
#     y = "Density"
#   ) +
#   theme_minimal()
# 
# 
# # Define the age group order for plotting
# age_group_levels_ordered <- c("0-0.99", "1-4", "5-9", "10-14", "15-19", "20-24", "25-29",
#                               "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
#                               "60-64", "65-69", "70-74", "75-79", "80-84", "85-89",
#                               "90-94", "95-99")
# 
# prepped_data_viz <- prepped_data %>%
#   mutate(age_group_f = factor(age_group, levels = age_group_levels_ordered))
# 
# # Plot 2: Protein Share by Age and Sex (Boxplots)
# ggplot(prepped_data_viz, aes(x = age_group_f, y = protein_kcal_share_mean, fill = sex)) +
#   geom_boxplot(alpha = 0.7, outlier.shape = NA) + # Hide outliers for now to see the box
#   geom_jitter(aes(color=sex), width=0.2, alpha=0.1, size=0.5) + # Add points to see density
#   scale_y_continuous(labels = scales::percent_format(), limits = c(0, 0.40)) + # Zoom in, ignore the big outlier
#   scale_fill_brewer(palette = "Set1") +
#   scale_color_brewer(palette = "Set1") +
#   labs(
#     title = "Protein's Share of Calories by Age Group and Sex",
#     subtitle = "Each point is a country. Y-axis is zoomed to 0-40% to show main trends.",
#     x = "Age Group",
#     y = "Protein Share of Standardized Calories",
#     fill = "Sex",
#     color = "Sex"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.position = "top")
# 
# 
# # Plot 3: Facet by a sample of countries
# # Let's pick a sample of countries from different regions/income levels
# sample_iso3 <- c("USA", "JPN", "NGA", "BRA", "IND", "CHN", "FRA", "ETH", "MDG")
# 
# prepped_data_viz %>%
#   filter(iso3 %in% sample_iso3) %>%
#   ggplot(aes(x = age_group_f, y = protein_kcal_share_mean, color = sex, group = sex)) +
#   geom_line(linewidth = 1) +
#   geom_point(size = 2) +
#   facet_wrap(~iso3, scales = "free_y") + # Use free_y to see the range in each country
#   scale_y_continuous(labels = scales::percent_format()) +
#   scale_color_brewer(palette = "Set1") +
#   labs(
#     title = "Protein's Share of Calories Across the Lifespan",
#     subtitle = "For a diverse sample of countries",
#     x = "Age Group",
#     y = "Protein Share of Standardized Calories"
#   ) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8))

