# Load required packages
library(dplyr)
library(readr)
library(ggplot2)

# Set working directory for your R script and clear the environment
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
options(scipen=999)
rm(list = ls())

# --- Step 1: Load GDD 2018 macronutrient data ---
# I've used the file numbers from your code.
protein_gdd <- read_csv("./data/GDD_FinalEstimates_01102022/Country-level estimates/v23_cnty.csv") # Protein (g/day)
carbs_gdd <- read_csv("./data/GDD_FinalEstimates_01102022/Country-level estimates/v22_cnty.csv") # Carbs (% kcal)
satfat_gdd <- read_csv("./data/GDD_FinalEstimates_01102022/Country-level estimates/v27_cnty.csv") # Sat Fat (% kcal)
mufa_gdd <- read_csv("./data/GDD_FinalEstimates_01102022/Country-level estimates/v28_cnty.csv") # MUFA (% kcal)
pufa_gdd <- read_csv("./data/GDD_FinalEstimates_01102022/Country-level estimates/v29_cnty.csv") # Omega-6 PUFA (% kcal)
# Note: v30 (Seafood Omega-3) and v31 (Plant Omega-3) are in mg/day.
# For this verification, their caloric contribution is negligible and would complicate the formula.
# We will exclude them for clarity, as this check won't be affected.


# --- Step 2: Create a function to clean and filter each dataset ---
# This avoids repeating the same code block multiple times.
preprocess_gdd_data <- function(df, value_col_name) {
  df %>%
    filter(year == 2018, edu == 999, urban == 999, age != 999, female != 999) %>%
    mutate(sex = ifelse(female == 1, "Females", "Males")) %>%
    select(iso3, sex, age, !!value_col_name := median) # Use !! and := to set column name
}

# Apply the function to each dataset
protein_clean <- preprocess_gdd_data(protein_gdd, "protein_g_day")
carbs_clean <- preprocess_gdd_data(carbs_gdd, "carbs_pct_kcal")
satfat_clean <- preprocess_gdd_data(satfat_gdd, "satfat_pct_kcal")
mufa_clean <- preprocess_gdd_data(mufa_gdd, "mufa_pct_kcal")
pufa_clean <- preprocess_gdd_data(pufa_gdd, "pufa_pct_kcal")


# --- Step 3: Join all data into a single master dataframe ---
# This is the correct approach instead of separate aggregations.
# We join by iso3, sex, and the original 'age' column.
verification_df <- protein_clean %>%
  left_join(carbs_clean, by = c("iso3", "sex", "age")) %>%
  left_join(satfat_clean, by = c("iso3", "sex", "age")) %>%
  left_join(mufa_clean, by = c("iso3", "sex", "age")) %>%
  left_join(pufa_clean, by = c("iso3", "sex", "age")) %>%
  # Remove rows where any of the essential values are missing
  na.omit()


# --- Step 4: Apply the CORRECT algebraic formula for verification ---
verification_final <- verification_df %>%
  mutate(
    # First, sum the fat percentages to get the total fat contribution to energy
    total_fat_pct_kcal = satfat_pct_kcal + mufa_pct_kcal + pufa_pct_kcal,
    
    # Now, apply the algebraic formula to solve for Total Kilocalories
    # Total Energy = (Energy from Protein) / (1 - %Energy from Carbs - %Energy from Fat)
    calculated_total_kcal = (protein_g_day * 4) / (1 - (carbs_pct_kcal / 100) - (total_fat_pct_kcal / 100))
  )

# Create the age_group variable for analysis AFTER the calculation
verification_final <- verification_final %>%
  mutate(age_group = case_when(
    age < 1 ~ "0-0.99",
    age < 5 ~ "1-4", # Combining for simplicity in verification
    age < 11 ~ "5-10",
    age < 75 ~ "11-74",
    age >= 75 ~ "75+",
    TRUE ~ NA_character_
  ))

# --- Step 5: Analyze the results ---

# Get a summary of the calculated calories for the main adult group
summary_adults <- verification_final %>%
  filter(age_group == "11-74") %>%
  summarise(
    count = n(),
    mean_calories = mean(calculated_total_kcal, na.rm = TRUE),
    median_calories = median(calculated_total_kcal, na.rm = TRUE),
    min_calories = min(calculated_total_kcal, na.rm = TRUE),
    max_calories = max(calculated_total_kcal, na.rm = TRUE),
    sd_calories = sd(calculated_total_kcal, na.rm = TRUE)
  )

cat("--- Verification Summary for Adult Group (Expected: 2000 kcal) ---\n")
print(summary_adults)

# Get a summary for the oldest group
summary_seniors <- verification_final %>%
  filter(age_group == "75+") %>%
  summarise(
    count = n(),
    mean_calories = mean(calculated_total_kcal, na.rm = TRUE),
    median_calories = median(calculated_total_kcal, na.rm = TRUE),
    min_calories = min(calculated_total_kcal, na.rm = TRUE),
    max_calories = max(calculated_total_kcal, na.rm = TRUE),
    sd_calories = sd(calculated_total_kcal, na.rm = TRUE)
  )

cat("\n--- Verification Summary for Senior Group (Expected: 1700 kcal) ---\n")
print(summary_seniors)


# --- Step 6: Visualize the distribution to confirm ---
ggplot(verification_final %>% filter(age_group == "11-74"), aes(x = calculated_total_kcal)) +
  geom_histogram(binwidth = 30, fill = "dodgerblue", color = "black", alpha = 0.8) +
  geom_vline(xintercept = 2000, color = "red", linetype = "dashed", size = 1.5) +
  labs(
    title = "Distribution of Calculated Calories for Ages 11-74",
    subtitle = "If standardized, values should cluster tightly around the red line (2000 kcal)",
    x = "Calculated Total Calories (kcal)",
    y = "Count of Country-Sex-Age Strata"
  ) +
  theme_minimal() +
  # Zoom in on the area of interest to see the distribution clearly
  coord_cartesian(xlim = c(0, 5000))

# --- Step 7: Filter for Plausible Values and Rank Countries ---

# First, let's filter the verification_final dataframe to include only
# values that are biologically plausible. We'll choose a generous range.
# A value below 500 or above 5000 for a stratum average is highly suspect.
plausible_calories_df <- verification_final %>%
  filter(calculated_total_kcal > 500 & calculated_total_kcal < 5000)

# Now, let's aggregate to the country level by taking the mean of these plausible values
country_level_summary <- plausible_calories_df %>%
  group_by(iso3) %>%
  summarise(
    avg_calculated_kcal = mean(calculated_total_kcal, na.rm = TRUE),
    n_strata = n() # Count how many valid strata we have for each country
  ) %>%
  ungroup()

# --- Find the 10 Countries with the HIGHEST Average Calculated Calories ---
highest_cal_countries <- country_level_summary %>%
  arrange(desc(avg_calculated_kcal)) %>%
  head(10)

cat("\n--- Countries with HIGHEST Average Calculated Calories (within plausible range) ---\n")
print(highest_cal_countries)


# --- Find the 10 Countries with the LOWEST Average Calculated Calories ---
lowest_cal_countries <- country_level_summary %>%
  arrange(avg_calculated_kcal) %>%
  head(10)

cat("\n--- Countries with LOWEST Average Calculated Calories (within plausible range) ---\n")
print(lowest_cal_countries)


artifact_counts <- verification_final %>%
  summarise(
    total_strata = n(),
    count_negative = sum(calculated_total_kcal < 0, na.rm = TRUE),
    count_over_5k = sum(calculated_total_kcal > 5000, na.rm = TRUE)
  ) %>%
  mutate(
    total_artifact_count = count_negative + count_over_5k,
    percent_artifact = (total_artifact_count / total_strata) * 100
  )

cat("\n--- Counts of Impossible 'Artifact' Calorie Values ---\n")
print(artifact_counts)


# --- Step 9: Inspect the Strata with Impossible 'Artifact' Values ---

# We continue using the `verification_final` dataframe.

# Create a dataframe for all strata that resulted in NEGATIVE calories
negative_calorie_strata <- verification_final %>%
  filter(calculated_total_kcal < 0) %>%
  # Select key columns for easy inspection and sort by the most extreme values
  select(iso3, sex, age, age_group, protein_g_day, carbs_pct_kcal, total_fat_pct_kcal, calculated_total_kcal) %>%
  arrange(calculated_total_kcal) # Sort to see the most negative values first

cat("\n\n--- Strata with NEGATIVE Calculated Calories (Most Extreme First) ---\n")
print(head(negative_calorie_strata, 20)) # Print the top 20 most extreme negative cases


# Create a dataframe for all strata that resulted in EXTREMELY HIGH calories
high_calorie_strata <- verification_final %>%
  filter(calculated_total_kcal > 5000) %>%
  # Select the same key columns for easy inspection
  select(iso3, sex, age, age_group, protein_g_day, carbs_pct_kcal, total_fat_pct_kcal, calculated_total_kcal) %>%
  arrange(desc(calculated_total_kcal)) # Sort to see the highest values first

cat("\n\n--- Strata with EXTREMELY HIGH (>5000) Calculated Calories (Most Extreme First) ---\n")
print(head(high_calorie_strata, 20)) # Print the top 20 most extreme high cases

# --- Optional but helpful: Summarize by country ---
# Let's count how many artifact strata each country has.

negative_counts_by_country <- negative_calorie_strata %>%
  group_by(iso3) %>%
  summarise(n_negative = n()) %>%
  arrange(desc(n_negative))

cat("\n\n--- Countries with the MOST Strata Resulting in Negative Calories ---\n")
print(head(negative_counts_by_country, 15))


high_counts_by_country <- high_calorie_strata %>%
  group_by(iso3) %>%
  summarise(n_high = n()) %>%
  arrange(desc(n_high))

cat("\n\n--- Countries with the MOST Strata Resulting in High (>5000) Calories ---\n")
print(head(high_counts_by_country, 15))

