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

# Load the pre-processed data from the previous script
#prepped_data <- readRDS("./output/gdd_prepped_data_for_caloric_adjustment.rds")

# --- Step 1: Load the single, complete data object from Script 2 ---
# This object contains everything we need: GDD data, population, and Marco's EER.
master_data <- readRDS("./output/script2_final_results.rds")

cat("--- Complete data from Script 2 loaded successfully. ---\n")

# ==========================================================================
# PART 1: DEFINE MASTER COUNTRY LIST
# ==========================================================================
all_countries_needed <- master_data %>% 
  distinct(iso3)

cat("--- Master list defined. Total countries required for analysis:", nrow(all_countries_needed), "---\n")


# ==========================================================================
# PART 2: PREPARE COMPLETE GDP DATASET (MANUALLY ADDING TAIWAN)
# ==========================================================================
library(WDI)

# --- Download raw GDP data ---
gdp_data_raw <- WDI(indicator = "NY.GDP.PCAP.PP.CD", start = 2018, end = 2018, extra = TRUE)

# --- Clean the data ---
# This will result in a dataframe that is MISSING Taiwan, which is what we expect now.
gdp_data_2018_no_twn <- gdp_data_raw %>%
  as_tibble() %>%
  rename(gdp_per_cap_ppp = NY.GDP.PCAP.PP.CD) %>%
  filter(nchar(iso3c) == 3) %>%
  select(iso3 = iso3c, country, gdp_per_cap_ppp)

# --- Manually create the Taiwan row ---
# We know the iso3 code, the country name, and the imputed GDP value.
taiwan_row_to_add <- tibble(
  iso3 = "TWN",
  country = "Taiwan",
  gdp_per_cap_ppp = 53000  # The value from our imputation research
)

# --- Bind the manually created Taiwan row to the main dataset ---
gdp_data_with_twn <- bind_rows(gdp_data_2018_no_twn, taiwan_row_to_add)

# --- Now, proceed with imputing the REMAINING missing countries ---
# The list of missing countries will now be just the other 5.
countries_missing_gdp <- all_countries_needed %>%
  anti_join(gdp_data_with_twn %>% filter(!is.na(gdp_per_cap_ppp)), by = "iso3")

cat("--- Countries still needing GDP imputation (after adding Taiwan):", paste(countries_missing_gdp$iso3, collapse=", "), "---\n")

# Define the imputation table for the remaining countries
gdp_imputation_table <- tibble::tribble(
  ~iso3, ~imputed_gdp,
  "CUB", 12300,
  "ERI", 1500,
  "SSD", 1600,
  "VEN", 12500,
  "YEM", 2500
)
####Sources to be written

# --- Create the final, complete GDP dataset ---
gdp_final_complete <- gdp_data_with_twn %>%
  left_join(gdp_imputation_table, by = "iso3") %>%
  mutate(
    gdp_per_cap_ppp = coalesce(gdp_per_cap_ppp, imputed_gdp)
  ) %>%
  select(iso3, gdp_per_cap_ppp) %>%
  # Final safety check to keep only the countries we need
  semi_join(all_countries_needed, by = "iso3")

cat("\n--- Final GDP dataset created. Total countries:", nrow(gdp_final_complete), "---\n")

# ==========================================================================
# PART 3: PREPARE COMPLETE CALORIE SUPPLY DATASET
# ==========================================================================
owid_raw <- read_csv("./data/OWID_daily-per-capita-caloric-supply.csv", show_col_types = FALSE)

owid_calories <- owid_raw %>%
  filter(Year == 2018) %>%
  select(iso3 = Code, fao_des_kcal = `Daily calorie supply per person`) %>%
  filter(!is.na(iso3) & nchar(iso3) == 3)

# --- Identify calorie values needed for imputation ---
countries_missing_calories <- all_countries_needed %>%
  anti_join(owid_calories, by = "iso3")

# --- Define the calorie imputation table ---
calorie_imputation_table <- tibble::tribble(
  ~iso3, ~proxy_iso3,
  "BHR", "SAU", "BRN", "MYS", "BTN", "NPL", "ERI", "ETH", "FSM", "FJI",
  "GNQ", "GAB", "MHL", "FJI", "PSE", "JOR", "QAT", "SAU", "SGP", "KOR",
  "SSD", "SDN", "TON", "FJI", "TWN", "KOR" # Taiwan needs to be imputed here too
)

imputation_values <- owid_calories %>%
  filter(iso3 %in% calorie_imputation_table$proxy_iso3) %>%
  select(proxy_iso3 = iso3, fao_des_kcal)

calories_to_add <- calorie_imputation_table %>%
  filter(iso3 %in% countries_missing_calories$iso3) %>%
  left_join(imputation_values, by = "proxy_iso3") %>%
  select(iso3, fao_des_kcal)

# --- Create the final, complete calorie dataset ---
calories_final_complete <- bind_rows(owid_calories, calories_to_add) %>%
  semi_join(all_countries_needed, by = "iso3")

cat("--- Final calorie dataset created. Total countries:", nrow(calories_final_complete), "---\n")


# ==========================================================================
# PART 4: CREATE MASTER DATA & MODEL WASTE
# ==========================================================================

# Start with the master list and join the two complete datasets.
country_level_data <- all_countries_needed %>%
  left_join(gdp_final_complete, by = "iso3") %>%
  left_join(calories_final_complete, by = "iso3")

# --- Now, modeling waste as a function of GDP (PPP)---

# Define the waste function
calculate_waste_pct <- function(gdp) {
  min_waste <- 0.05; max_waste <- 0.35; gdp_midpoint <- 10000; k <- 0.0001
  waste_percentage <- min_waste + (max_waste - min_waste) / (1 + exp(-k * (gdp - gdp_midpoint)))
  return(waste_percentage)
}

# Apply waste model
country_level_data <- country_level_data %>%
  mutate(
    downstream_waste_pct = calculate_waste_pct(gdp_per_cap_ppp),
    est_kcal_consumed_per_capita = fao_des_kcal * (1 - downstream_waste_pct)
  )


# ==========================================================================
# PART 5: VISUALIZE AND VERIFY THE WASTE MODEL
# ==========================================================================
library(ggplot2)
# install.packages("ggrepel") # Run this line once if you don't have the ggrepel package
library(ggrepel)

# --- Step 1: Select a sample of countries to label on the plot ---
# We choose countries from across the GDP spectrum to illustrate the model's behavior.
countries_to_label <- country_level_data %>%
  filter(iso3 %in% c(
    # High Income / High Waste
    "CHE", "SGP", "USA", 
    # Upper-Middle Income (on the steep part of the curve)
    "CHN", "BRA", "ZAF",
    # Lower-Middle / Low Income
    "IND", "NGA", "BGD",
    # Very Low Income (on the floor of the curve)
    "SSD", "ERI"
  ))


# --- Step 2: Create the plot ---

waste_model_plot <- ggplot(country_level_data, aes(x = gdp_per_cap_ppp, y = downstream_waste_pct)) +
  
  # Plot all 185 countries as points to show the full distribution
  geom_point(alpha = 0.5, color = "gray60") +
  
  # Highlight the sample countries with colored points
  geom_point(data = countries_to_label, color = "dodgerblue", size = 3) +
  
  # Add labels to the sample countries using ggrepel to avoid overlap
  geom_text_repel(
    data = countries_to_label,
    aes(label = iso3),
    fontface = "bold",
    size = 3.5,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = 'grey50'
  ) +
  
  # Format the axes for readability
  scale_x_continuous(
    name = "GDP per capita, PPP (2018)",
    labels = scales::dollar_format(prefix = "$", scale = 1/1000, suffix = "k") # e.g., $20k
  ) +
  scale_y_continuous(
    name = "Estimated Downstream Waste %",
    labels = scales::percent_format(accuracy = 1)
  ) +
  
  # Add informative titles
  labs(
    title = "Modeled Relationship Between GDP and Downstream Food Waste",
    subtitle = "Based on a logistic function with a 2% floor and 30% ceiling"
  ) +
  
  # Use a clean theme
  theme_minimal(base_size = 14)

# --- Step 3: Print the plot ---
print(waste_model_plot)


# ==========================================================================
# PART 6: DISAGGREGATE NATIONAL CALORIC CONSUMPTION TO STRATA
# ==========================================================================

# Our goal is to take the country-level 'est_kcal_consumed_per_capita' and
# distribute it among the various age-sex strata within each country.

# --- Step 1: Create the master analytical dataframe ---
# This joins the country-level data (with our new consumption estimates)
# to the stratum-level data (which has EER and population).

# The 'master_data' object already contains everything we need:
# iso3, sex, age_group, population, and eer_kcal_marco_mean.
# The 'country_level_data' has iso3 and our new est_kcal_consumed_per_capita.

analysis_data_full <- master_data %>%
  # Join our newly calculated country-level consumption estimates
  left_join(
    country_level_data %>% select(iso3, est_kcal_consumed_per_capita),
    by = "iso3"
  )

# --- Step 2: Perform the disaggregation ---
# The logic is identical to the one we outlined before: use EERs as relative
# weights to distribute the national total.

disaggregated_data <- analysis_data_full %>%
  # For each country...
  group_by(iso3) %>%
  
  # Calculate the population-weighted average EER for that country.
  # This represents the "average requirement" for a person in that country,
  # given its specific demographic structure.
  mutate(
    avg_country_eer = weighted.mean(eer_kcal_marco_mean, w = population, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  
  # Now, calculate the adjustment factor for each stratum. A stratum with an
  # EER higher than the country average will get a factor > 1, and vice-versa.
  mutate(
    adjustment_factor = eer_kcal_marco_mean / avg_country_eer,
    
    # Finally, apply this factor to the national per-capita CONSUMPTION estimate
    # to get the final estimated caloric consumption for that specific stratum.
    kcal_consumed_stratum = est_kcal_consumed_per_capita * adjustment_factor
  )

# ==========================================================================
# PART 7: VISUALIZE FINAL DISAGGREGATED CALORIC CONSUMPTION
# ==========================================================================
library(ggplot2)

# --- Step 1: Define the sample of countries and the age group order ---
countries_to_plot <- c("MDG", "BRA", "ECU", "CHN", "ITA", "DEU", "USA")

age_group_levels_ordered <- c("0-0.99", "1-4", "5-9", "10-14", "15-19",
                              "20-24", "25-29", "30-34", "35-39",
                              "40-44", "45-49", "50-54", "55-59",
                              "60-64", "65-69", "70-74", "75-79",
                              "80-84", "85-89", "90-94", "95-99")

# --- Step 2: Prepare the data for plotting ---
# Filter for our sample countries and create an ordered factor for the x-axis.
plot_data_calories <- disaggregated_data %>%
  filter(iso3 %in% countries_to_plot) %>%
  mutate(
    # Create an ordered factor for the x-axis
    age_group_f = factor(age_group, levels = age_group_levels_ordered),
    # Create a new column that combines country and sex for distinct lines
    country_sex = paste(iso3, sex, sep = " - ")
  )

# --- Step 3: Create the plot ---
# We will use facets to create a separate panel for each country.

disaggregated_calories_plot <- ggplot(plot_data_calories, 
                                      aes(x = age_group_f, 
                                          y = kcal_consumed_stratum, 
                                          group = sex, 
                                          color = sex)) +
  
  # Draw the lines connecting the points for each sex
  geom_line(linewidth = 1.2, alpha = 0.8) +
  
  # Add points at each age group to make the values clear
  geom_point(size = 2.5) +
  
  # Create a separate plot panel for each country
  facet_wrap(~iso3, scales = "free_y", ncol = 4) +
  
  # Use a clear color scheme
  scale_color_brewer(palette = "Set1") +
  
  # Add informative labels and a title
  labs(
    title = "Estimated Caloric Consumption Across the Lifespan",
    subtitle = "For a diverse sample of countries, disaggregated from national supply",
    x = "Age Group",
    y = "Estimated Caloric Consumption (kcal/day)",
    color = "Sex"
  ) +
  
  # Use a clean theme and improve readability
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9), # Smaller text for crowded facets
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 12) # Make country labels in facets bold
  )

# --- Step 4: Print the plot ---
print(disaggregated_calories_plot)
