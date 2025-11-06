
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
gdd_kcal_standards <- tibble::tribble(
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

#question: is it really not sex dependent? seems strange

# Correcting the 1-4 group based on the components:
# GDD: 1–<2 years = 1000 kcal; 2–5 years = 1300 kcal.
# Your 1-4 group is a 1-year band and a 3-year band.
# Weighted avg: (1 * 1000 + 3 * 1300) / 4 = 1225 kcal.
gdd_kcal_standards$standard_kcal[gdd_kcal_standards$age_group == "1-4"] <- 1225


# Join with standardization values and calculate protein's share of calories
prepped_data <- prepped_data %>%
  left_join(gdd_kcal_standards, by = "age_group") %>%
  mutate(
    # Protein provides 4 kcal/gram
    protein_kcal_share_mean  = (gdd_mean * 4) / standard_kcal,
    protein_kcal_share_lower = (gdd_lower * 4) / standard_kcal,
    protein_kcal_share_upper = (gdd_upper * 4) / standard_kcal
  )


# Find the row(s) with the highest protein share
prepped_data %>%
  arrange(desc(protein_kcal_share_mean)) %>%
  select(iso3, sex, age_group, gdd_mean, standard_kcal, protein_kcal_share_mean) %>%
  head()



# --------------------------------------------------------------------------
# NEW STRATEGY: Using the Our World in Data (OWID) Caloric Supply File
# --------------------------------------------------------------------------

# --- Step 1: Load the OWID data file ---
owid_raw <- read_csv("./data/OWID_daily-per-capita-caloric-supply.csv", show_col_types = FALSE)
#from https://ourworldindata.org/grapher/daily-per-capita-caloric-supply

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

# Explicitly check for Japan to be sure
japan_owid_check <- owid_calories %>% filter(iso3 == "JPN")

cat("\n--- Checking for Japan in OWID data ---\n")
if (nrow(japan_owid_check) > 0) {
  cat("✅ Found Japan in the OWID data!\n")
  print(japan_owid_check)
} else {
  cat("❌ Japan is STILL missing, even from OWID.\n")
}


# --- Step 4: Final Inspection of the Cleaned OWID Data ---
cat("\n--- Final, Cleaned OWID Caloric Supply Data (2018) ---\n")
print(head(owid_calories))
summary(owid_calories$fao_des_kcal)

# From now on, we will use 'owid_calories' as our source for national caloric supply.

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

# A special check for South Sudan's proxy. If Sudan (SDN) is now in our data, use it. Otherwise, fall back to Ethiopia.
if ("SDN" %in% owid_calories$iso3) {
  imputation_map_final$proxy_iso3[imputation_map_final$iso3 == "SSD"] <- "SDN"
} else {
  imputation_map_final$proxy_iso3[imputation_map_final$iso3 == "SSD"] <- "ETH"
}

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

cat("\n--- Final Imputation Table (for the 12 missing countries) ---\n")
print(imputation_table_final)

# Step 4: Bind the imputed rows to the main OWID data to create the complete dataset
calories_final_complete <- bind_rows(owid_calories, imputation_table_final)

# --- FINAL VALIDATION (This MUST pass) ---
missing_last_check <- prepped_data %>%
  distinct(iso3) %>%
  anti_join(calories_final_complete, by = "iso3")

cat("\n--- Final Check After Imputation ---\n")
if (nrow(missing_last_check) > 0) {
  cat("❌ The following", nrow(missing_last_check), "countries are STILL missing somehow:\n")
  print(missing_last_check$iso3)
} else {
  cat("✅✅✅ COMPLETE SUCCESS! All", nrow(prepped_data %>% distinct(iso3)), "countries now have a caloric supply value.\n")
}

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

# You can now proceed with the rest of the analysis using this 'eer_table' object.

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

# --- CHECKPOINT 1: Inspect the Master Data Frame ---
# This check is crucial to ensure that the joins worked as expected and did not
# introduce unexpected missing values. We expect zero NAs in the newly joined columns.

cat("\n--- CHECKPOINT 1: Validating the joined master_data ---\n")

# Summarize NA counts for the newly added columns
na_check <- master_data %>%
  summarise(
    NAs_in_population = sum(is.na(population)),
    NAs_in_calories = sum(is.na(fao_des_kcal)),
    NAs_in_eer = sum(is.na(eer_kcal))
  )

# Print the validation results
if (all(na_check == 0)) {
  cat("✅ SUCCESS: All joins were successful. No missing values introduced.\n")
} else {
  cat("⚠️ WARNING: Joins introduced NA values. Please inspect the summary below.\n")
}
print(na_check)

# As a robust measure, we will filter out any rows that might have NAs,
# although we do not expect any at this point.
master_data <- master_data %>%
  filter(!is.na(population), !is.na(fao_des_kcal), !is.na(eer_kcal))

# Let's glimpse the final structure of our "big table"
cat("\n--- Structure of the final master_data frame ---\n")
glimpse(master_data)


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

# --- CHECKPOINT 2: Validate the Disaggregation Math ---
# This check ensures our disaggregation logic is sound. We will average the
# new stratum-specific values back up and confirm it equals the original national value.
cat("\n--- CHECKPOINT 2: Validating the disaggregation for one country (e.g., USA) ---\n")

validation_check <- disaggregated_data %>%
  filter(iso3 == "USA") %>% # Pick any large, diverse country
  summarise(
    original_supply = first(fao_des_kcal),
    recalculated_average_supply = weighted.mean(disagg_kcal_supply, w = population, na.rm = TRUE)
  ) %>%
  mutate(difference = original_supply - recalculated_average_supply)

print(validation_check)
# We expect the 'difference' to be a very small number (e.g., 1e-12), proving our math is correct.


# --- Step 3b: Create Consumption Scenarios & Calculate "True" Protein Intake ---
# Finally, we apply waste scenarios and calculate the final protein grams.
final_data <- disaggregated_data %>%
  mutate(
    # --- Consumption Scenarios (adjusting for waste) ---
    # High Scenario: Assumes 0% waste (Consumption = Supply)
    est_kcal_consumed_high = disagg_kcal_supply * 1.00,
    # Mean Scenario: Assumes 15% waste
    est_kcal_consumed_mean = disagg_kcal_supply * 0.85,
    # Low Scenario: Assumes 30% waste
    est_kcal_consumed_low  = disagg_kcal_supply * 0.70,
    
    # --- "True" Protein Intake in Grams/Day ---
    # This combines the consumption scenarios with the protein share blueprints.
    # Central estimate (Mean consumption, Mean protein share)
    true_protein_g_mean = (est_kcal_consumed_mean * protein_kcal_share_mean) / 4,
    
    # Low estimate for sensitivity (Low consumption, Low protein share)
    true_protein_g_low = (est_kcal_consumed_low * protein_kcal_share_lower) / 4,
    
    # High estimate for sensitivity (High consumption, High protein share)
    true_protein_g_high = (est_kcal_consumed_high * protein_kcal_share_upper) / 4
  )

# --- CHECKPOINT 3: Inspect the Final Analytical Data ---
cat("\n--- CHECKPOINT 3: Final data with 'TRUE' protein intake estimates ---\n")
final_data %>%
  # Select key new columns to review
  select(iso3, sex, age_group, true_protein_g_low, true_protein_g_mean, true_protein_g_high, ear_mean_g_day) %>%
  head()

cat("\nSummary of the central estimate for 'true' protein intake (grams/day):\n")
summary(final_data$true_protein_g_mean)















# --- Step 1: Load FAO Data ---
# Load the raw data from the 'data' subfolder
fao_raw <- read_csv("./data/fao_data_full.csv")

cat("✅ FAO raw data loaded. Glimpse of the data:\n")
glimpse(fao_raw)


# --- FINAL & DEFINITIVE FAO DATA PROCESSING (with whitespace trim) ---

# This code block cleans whitespace, applies robust matching, and will produce
# the TRUE list of genuinely missing countries.

# Filter for the specific rows we need (no change here)
fao_calories_filtered <- fao_raw %>%
  filter(
    Item == "Grand Total",
    Element == "Food supply (kcal/capita/day)"
  )

# Select the identifier columns and the 2018 data column (no change here)
fao_calories_2018 <- fao_calories_filtered %>%
  select(
    `Area Code (M49)`,
    Area,
    fao_des_kcal = Y2018
  )

# Create the final data frame with the added str_trim() fix
fao_calories <- fao_calories_2018 %>%
  mutate(
    # --- THE CRUCIAL FIX IS HERE ---
    # First, trim leading/trailing whitespace from all Area names
    Area = stringr::str_trim(Area),
    
    # Now, the rest of the matching logic will work on the CLEANED names
    iso3_guess = countrycode(Area, origin = "country.name", destination = "iso3c", nomatch = NA),
    
    iso3 = case_when(
      Area == "Bolivia (Plurinational State of)" ~ "BOL",
      Area == "China, mainland" ~ "CHN",
      Area == "China, Taiwan Province of" ~ "TWN",
      Area == "Côte d'Ivoire" ~ "CIV",
      Area == "Czechia" ~ "CZE",
      Area == "Iran (Islamic Republic of)" ~ "IRN",
      Area == "Micronesia (Federated States of)" ~ "FSM",
      Area == "Netherlands (Kingdom of the)" ~ "NLD",
      Area == "Republic of Korea" ~ "KOR",
      Area == "Republic of Moldova" ~ "MDA",
      Area == "Russian Federation" ~ "RUS",
      Area == "Türkiye" ~ "TUR",
      Area == "United Kingdom of Great Britain and Northern Ireland" ~ "GBR",
      Area == "United Republic of Tanzania" ~ "TZA",
      Area == "United States of America" ~ "USA",
      Area == "Venezuela (Bolivarian Republic of)" ~ "VEN",
      Area == "Viet Nam" ~ "VNM",
      TRUE ~ iso3_guess
    )
  ) %>%
  # The rest of the pipeline is the same
  select(iso3, fao_des_kcal) %>%
  filter(!is.na(iso3), !is.na(fao_des_kcal)) %>%
  group_by(iso3) %>%
  summarise(fao_des_kcal = mean(fao_des_kcal, na.rm = TRUE), .groups = "drop")


# --- FINAL VALIDATION (RE-RUN) ---
missing_from_fao_final <- prepped_data %>%
  distinct(iso3) %>%
  anti_join(fao_calories, by = "iso3")

cat("\n--- FINAL LIST of countries missing from the FAO file (after trim) ---\n")
if (nrow(missing_from_fao_final) > 0) {
  cat("The following", nrow(missing_from_fao_final), "countries are GENUINELY MISSING from the source CSV file:\n")
  print(missing_from_fao_final$iso3)
} else {
  cat("✅ Miraculously, all countries were matched!\n")
}









#function to calculate protein inadequacy
calculate_inadequacy <- function(mean_intake, cv_intake, distribution_type, requirement) {
  # Ensure requirement is not NA or negative, handle if necessary
  if (is.na(requirement) || requirement < 0) return(NA_real_)
  # Ensure mean_intake and cv_intake are valid
  if (is.na(mean_intake) || is.na(cv_intake) || mean_intake <= 0 || cv_intake <= 0) return(NA_real_)
  
  if (distribution_type == "gamma") {
    if (cv_intake^2 == 0) return(NA_real_) # Avoid division by zero if CV is exactly 0
    shape_k <- 1 / (cv_intake^2)
    scale_theta <- mean_intake / shape_k # which is mean_intake * cv_intake^2
    p_inadequate <- pgamma(requirement, shape = shape_k, scale = scale_theta)
  } else if (distribution_type == "log-normal") {
    if (1 + cv_intake^2 <= 0) return(NA_real_) # log argument must be positive
    meanlog <- log(mean_intake) - 0.5 * log(1 + cv_intake^2)
    sdlog <- sqrt(log(1 + cv_intake^2))
    if (is.na(sdlog) || sdlog <=0) return(NA_real_) # sdlog must be positive
    p_inadequate <- plnorm(requirement, meanlog = meanlog, sdlog = sdlog)
  } else {
    # Handle unknown distribution types if any, or stop with an error
    warning(paste("Unknown distribution type:", distribution_type))
    p_inadequate <- NA_real_
  }
  return(p_inadequate)
}











####viz of relative protein share by sex, age, countries
#the 0-0.99 group looks funky for some countries

library(ggplot2)
library(viridis) # For nice color palettes

# Plot 1: Overall Distribution of Protein Calorie Share
ggplot(prepped_data, aes(x = protein_kcal_share_mean)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.01, fill = "skyblue", color = "white", alpha = 0.8) +
  geom_density(color = "dodgerblue4", linewidth = 1.2) +
  scale_x_continuous(labels = scales::percent_format()) +
  labs(
    title = "Distribution of Protein's Share of Calories (from GDD data)",
    subtitle = "Across all 7,770 country-sex-age strata",
    x = "Protein Share of Standardized Calories",
    y = "Density"
  ) +
  theme_minimal()


# Define the age group order for plotting
age_group_levels_ordered <- c("0-0.99", "1-4", "5-9", "10-14", "15-19", "20-24", "25-29",
                              "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
                              "60-64", "65-69", "70-74", "75-79", "80-84", "85-89",
                              "90-94", "95-99")

prepped_data_viz <- prepped_data %>%
  mutate(age_group_f = factor(age_group, levels = age_group_levels_ordered))

# Plot 2: Protein Share by Age and Sex (Boxplots)
ggplot(prepped_data_viz, aes(x = age_group_f, y = protein_kcal_share_mean, fill = sex)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) + # Hide outliers for now to see the box
  geom_jitter(aes(color=sex), width=0.2, alpha=0.1, size=0.5) + # Add points to see density
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 0.40)) + # Zoom in, ignore the big outlier
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Protein's Share of Calories by Age Group and Sex",
    subtitle = "Each point is a country. Y-axis is zoomed to 0-40% to show main trends.",
    x = "Age Group",
    y = "Protein Share of Standardized Calories",
    fill = "Sex",
    color = "Sex"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")


# Plot 3: Facet by a sample of countries
# Let's pick a sample of countries from different regions/income levels
sample_iso3 <- c("USA", "JPN", "NGA", "BRA", "IND", "CHN", "FRA", "ETH", "MDG")

prepped_data_viz %>%
  filter(iso3 %in% sample_iso3) %>%
  ggplot(aes(x = age_group_f, y = protein_kcal_share_mean, color = sex, group = sex)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~iso3, scales = "free_y") + # Use free_y to see the range in each country
  scale_y_continuous(labels = scales::percent_format()) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Protein's Share of Calories Across the Lifespan",
    subtitle = "For a diverse sample of countries",
    x = "Age Group",
    y = "Protein Share of Standardized Calories"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8))
