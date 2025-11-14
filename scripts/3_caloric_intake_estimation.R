#load needed elements from scripts 1,2
#load FAO FBS data and WDI per-capita GDP (PPP) data
#model retail and consumer waste as a function of GDP
#distribute adjusted caloric intake to strata based on EER (from Marco's data)
#match distribution and CVs for calories using NutriR just as done for protein
#final output: caloric intake estimation and distribution for each strata


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

library(WDI)

# Define the indicator and year
indicator_code <- "NY.GDP.PCAP.PP.CD"
start_year <- 2018
end_year <- 2018

# Download the data
gdp_data_raw <- WDI(
  indicator = indicator_code,
  start = start_year,
  end = end_year,
  extra = TRUE # Set to TRUE to get country metadata like region, income level etc.
)

# Clean up the data for use
gdp_data_2018 <- gdp_data_raw %>%
  as_tibble() %>%
  # Rename the column for clarity
  rename(gdp_per_cap_ppp = !!indicator_code) %>%
  # Select only the columns you need
  select(iso3 = iso3c, iso2 = iso2c, country, year, gdp_per_cap_ppp, income, region) %>%
  # Filter out non-country aggregates (like "World", "Arab World", etc.)
  # The 'income' column is a good way to do this, as aggregates have NA for income level.
  filter(!is.na(income), income != "Aggregates")

# Verify the result
print(head(gdp_data_2018))
cat("\nNumber of countries/territories loaded:", nrow(gdp_data_2018), "\n")

# --- DIAGNOSIS: Investigate NA values in GDP data ---

# Filter for rows where gdp_per_cap_ppp is NA
missing_gdp_countries <- gdp_data_2018 %>%
  filter(is.na(gdp_per_cap_ppp))

# Print the list of countries with missing GDP for 2018
cat("The following", nrow(missing_gdp_countries), "countries have NA for GDP per capita (PPP) in 2018:\n")
print(missing_gdp_countries %>% select(iso3, country, income, region))

# --- FOCUSED IMPUTATION for the 6 required countries ---

# Imputation table for ONLY the countries that are in our analysis but missing GDP.
# The object name and structure are identical to your original for easy replacement.
imputation_table_focused <- tibble::tribble(
  ~iso3, ~imputed_gdp, ~country_name_fix,
  
  # Source: CIA World Factbook, 2017 estimate.
  "CUB", 12300,        "Cuba",
  
  # Source: African Development Bank (2019) & regional estimates.
  "ERI", 1500,         "Eritrea",
  
  # Source: IMF World Economic Outlook (Oct 2021) & regional estimates.
  "SSD", 1600,         "South Sudan",
  
  # Source: IMF World Economic Outlook (Oct 2021), direct PPP data.
  "TWN", 53000,        "Taiwan",
  
  # Source: IMF World Economic Outlook (Apr 2023), direct 2018 PPP estimate.
  "VEN", 12500,        "Venezuela, RB",
  
  # Source: World Bank/UN consensus estimates for conflict period.
  "YEM", 2500,         "Yemen, Rep."
)

# Make sure gdp_data_2018 is loaded from your previous step.
# If not, re-run the WDI call.
if (!exists("gdp_data_2018")) {
  # (Code to load gdp_data_2018 from WDI)
  gdp_data_raw <- WDI(indicator = "NY.GDP.PCAP.PP.CD", start = 2018, end = 2018, extra = TRUE)
  gdp_data_2018 <- gdp_data_raw %>% as_tibble() %>%
    rename(gdp_per_cap_ppp = NY.GDP.PCAP.PP.CD) %>%
    select(iso3 = iso3c, country, gdp_per_cap_ppp) %>%
    filter(!is.na(iso3) & nchar(iso3) == 3)
}


# Join the imputation table and fill the gaps
gdp_data_imputed <- gdp_data_2018 %>%
  left_join(imputation_table_focused, by = "iso3") %>%
  mutate(
    # Add a flag for transparency
    gdp_imputed_flag = if_else(!is.na(imputed_gdp), TRUE, FALSE),
    
    # Fill the NA GDP values
    gdp_per_cap_ppp = coalesce(gdp_per_cap_ppp, imputed_gdp),
    
    # Fix the missing country name for Taiwan
    country = coalesce(country, country_name_fix)
  ) %>%
  # Remove the temporary helper columns
  select(iso3, country, gdp_per_cap_ppp, gdp_imputed_flag)


# --- VERIFICATION ---

cat("\n--- Review of Imputed Rows ---\n")
gdp_data_imputed %>%
  filter(gdp_imputed_flag == TRUE) %>%
  print()








# --- Step 1: Load the OWID FAO FBS data file ---
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