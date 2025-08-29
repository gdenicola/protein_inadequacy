# Load required packages
library(dplyr)
library(readr)
library(nutriR)
library(ggplot2)

# Set working directory for your R script and clear the environment
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
options(scipen=999)
rm(list = ls())

# Step 1: Load nutriR fat and carbs distribution data
protein_nutriR <- nutriR::get_dists(nutrients = "Protein")
carbs_nutriR <- nutriR::get_dists(nutrients = "Carbohydrates")
fats_nutriR <- nutriR::get_dists(nutrients = "Fat")

# Step 2: Load GDD 2018 macronutrient data
protein_gdd <- read_csv("./data/GDD_FinalEstimates_01102022/Country-level estimates/v23_cnty.csv")
carbs_gdd <- read_csv("./data/GDD_FinalEstimates_01102022/Country-level estimates/v22_cnty.csv")
fat1 <- read_csv("./data/GDD_FinalEstimates_01102022/Country-level estimates/v27_cnty.csv")
fat2 <- read_csv("./data/GDD_FinalEstimates_01102022/Country-level estimates/v28_cnty.csv")
fat3 <- read_csv("./data/GDD_FinalEstimates_01102022/Country-level estimates/v29_cnty.csv")
fat4 <- read_csv("./data/GDD_FinalEstimates_01102022/Country-level estimates/v30_cnty.csv")
fat5 <- read_csv("./data/GDD_FinalEstimates_01102022/Country-level estimates/v31_cnty.csv")

# Step 3: Preprocess and filter for 2018 data for all macronutrients
protein_gdd_natl <- protein_gdd %>%
  filter(year == 2018, edu == 999, urban == 999, age != 999, female != 999) %>%
  mutate(sex = ifelse(female == 1, "Females", "Males")) %>%
  select(iso3, sex, age, median, lowerci_95, upperci_95) # relevant protein columns

carbs_gdd_natl <- carbs_gdd %>%
  filter(year == 2018, edu == 999, urban == 999, age != 999, female != 999) %>%
  mutate(sex = ifelse(female == 1, "Females", "Males")) %>%
  select(iso3, sex, age, median, lowerci_95, upperci_95) # relevant carbs columns

# Combine all fat datasets into one
fat_gdd_natl <- bind_rows(
  fat1 %>% filter(year == 2018, edu == 999, urban == 999, female != 999) %>%
    mutate(sex = ifelse(female == 1, "Females", "Males")) %>%
    select(iso3, sex, age, median, lowerci_95, upperci_95),
  
  fat2 %>% filter(year == 2018, edu == 999, urban == 999, female != 999) %>%
    mutate(sex = ifelse(female == 1, "Females", "Males")) %>%
    select(iso3, sex, age, median, lowerci_95, upperci_95),
  
  fat3 %>% filter(year == 2018, edu == 999, urban == 999, female != 999) %>%
    mutate(sex = ifelse(female == 1, "Females", "Males")) %>%
    select(iso3, sex, age, median, lowerci_95, upperci_95),
  
  fat4 %>% filter(year == 2018, edu == 999, urban == 999, female != 999) %>%
    mutate(sex = ifelse(female == 1, "Females", "Males")) %>%
    select(iso3, sex, age, median, lowerci_95, upperci_95),
  
  fat5 %>% filter(year == 2018, edu == 999, urban == 999, female != 999) %>%
    mutate(sex = ifelse(female == 1, "Females", "Males")) %>%
    select(iso3, sex, age, median, lowerci_95, upperci_95)
)

# Step 4: Create age_group variable for all datasets
create_age_group <- function(df) {
  df %>% mutate(age_group = case_when(
    age < 1 ~ "0-0.99",
    age < 2 ~ "1-1.99",
    age < 5 ~ "2-4",
    age < 10 ~ "5-9",
    age < 15 ~ "10-14",
    age < 20 ~ "15-19",
    age < 25 ~ "20-24",
    age < 30 ~ "25-29",
    age < 35 ~ "30-34",
    age < 40 ~ "35-39",
    age < 45 ~ "40-44",
    age < 50 ~ "45-49",
    age < 55 ~ "50-54",
    age < 60 ~ "55-59",
    age < 65 ~ "60-64",
    age < 70 ~ "65-69",
    age < 75 ~ "70-74",
    age < 80 ~ "75-79",
    age < 85 ~ "80-84",
    age < 90 ~ "85-89",
    age < 95 ~ "90-94",
    age < 100 ~ "95-99",
    TRUE ~ NA_character_
  ))
}

# Apply the create_age_group function
protein_gdd_natl <- create_age_group(protein_gdd_natl)
carbs_gdd_natl <- create_age_group(carbs_gdd_natl)
fat_gdd_natl <- create_age_group(fat_gdd_natl)

# Step 5: Aggregate means for all macronutrients by iso3, sex, and age_group
agg_protein <- protein_gdd_natl %>%
  group_by(iso3, sex, age_group) %>%
  summarise(
    gdd_mean_protein = mean(median, na.rm = TRUE),
    gdd_lower_protein = mean(lowerci_95, na.rm = TRUE),
    gdd_upper_protein = mean(upperci_95, na.rm = TRUE),
    .groups = "drop"
  )

agg_carbs <- carbs_gdd_natl %>%
  group_by(iso3, sex, age_group) %>%
  summarise(
    gdd_mean_carbs = mean(median, na.rm = TRUE),
    gdd_lower_carbs = mean(lowerci_95, na.rm = TRUE),
    gdd_upper_carbs = mean(upperci_95, na.rm = TRUE),
    .groups = "drop"
  )

agg_fats <- fat_gdd_natl %>%
  group_by(iso3, sex, age_group) %>%
  summarise(
    gdd_mean_fats = mean(median, na.rm = TRUE),
    gdd_lower_fats = mean(lowerci_95, na.rm = TRUE),
    gdd_upper_fats = mean(upperci_95, na.rm = TRUE),
    .groups = "drop"
  )

# Step 6: Merge aggregated means into one dataframe
caloric_intake_agg <- agg_protein %>%
  left_join(agg_carbs, by = c("iso3", "sex", "age_group")) %>%
  left_join(agg_fats, by = c("iso3", "sex", "age_group"))

# Step 7: Calculate Total Calories using Atwater factors
caloric_intake_final <- caloric_intake_agg %>%
  mutate(
    total_calories_mean = (gdd_mean_protein * 4) + 
      (gdd_mean_carbs * 4) + 
      (gdd_mean_fats * 9), # Atwater factors
    total_calories_lower = (gdd_lower_protein * 4) + 
      (gdd_lower_carbs * 4) + 
      (gdd_lower_fats * 9),
    total_calories_upper = (gdd_upper_protein * 4) + 
      (gdd_upper_carbs * 4) + 
      (gdd_upper_fats * 9)
  )

# Step 8: Summary of Total Caloric Intake
summary_caloric_intake <- caloric_intake_final %>%
  summarise(
    count = n(),
    mean_calories = mean(total_calories_mean, na.rm = TRUE),
    median_calories = median(total_calories_mean, na.rm = TRUE),
    min_calories = min(total_calories_mean, na.rm = TRUE),
    max_calories = max(total_calories_mean, na.rm = TRUE),
    sd_calories = sd(total_calories_mean, na.rm = TRUE)
  )

# Print the summary of total caloric intake across all groups
print(summary_caloric_intake)

# Step 9: Visualize the total caloric intake distribution
boxplot_caloric_intake <- ggplot(caloric_intake_final, aes(x = sex, y = total_calories_mean)) +
  geom_boxplot() +
  labs(title = "Distribution of Total Caloric Intake by Sex",
       x = "Sex",
       y = "Mean Total Calories (kcal/day)") +
  theme_minimal()

print(boxplot_caloric_intake)

# Optional: Highlight some richest and poorest countries by total caloric intake
richest_countries <- caloric_intake_final %>%
  arrange(desc(total_calories_mean)) %>%
  slice(1:10) # Top 10 richest based on total calories

poorest_countries <- caloric_intake_final %>%
  arrange(total_calories_mean) %>%
  slice(1:10) # Bottom 10 poorest based on total calories

cat("Richest Countries (Top 10 by Total Caloric Intake):\n")
print(richest_countries)

cat("\nPoorest Countries (Bottom 10 by Total Caloric Intake):\n")
print(poorest_countries)


# Define the countries to filter (as their ISO3 codes)
#selected_countries <- c("USA", "ITA", "GRC", "AFG", "MDG", "ETH", "BWA", "LAO", "GBR")
#selected_countries <- c("YEM", "SSD", "CAF", "TCD", "NER", "SDN", "AFG", "HTI", "BFA", "ETH")
selected_countries <- c("BDI", "TLS", "NER", "MDG", "NPL", "USA", "GBR", "AUS", "WSM", "KIR")

# Step to Filter the caloric intake data for the selected countries
filtered_caloric_intake <- caloric_intake_final %>%
  filter(iso3 %in% selected_countries)

# Ensure the order of countries is based on `selected_countries`
filtered_caloric_intake$iso3 <- factor(filtered_caloric_intake$iso3, levels = selected_countries)

# Check the results for the specified countries
print(filtered_caloric_intake)

# (Optional) Summary statistics for the selected countries
summary_filtered_caloric_intake <- filtered_caloric_intake %>%
  summarise(
    count = n(),
    mean_calories = mean(total_calories_mean, na.rm = TRUE),
    median_calories = median(total_calories_mean, na.rm = TRUE),
    min_calories = min(total_calories_mean, na.rm = TRUE),
    max_calories = max(total_calories_mean, na.rm = TRUE),
    sd_calories = sd(total_calories_mean, na.rm = TRUE)
  )

# Print summary statistics
print(summary_filtered_caloric_intake)

# (Optional) Visualization for selected countries
ggplot(filtered_caloric_intake, aes(x = iso3, y = total_calories_mean, fill = sex)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Total Caloric Intake for Selected Countries",
       x = "Country (ISO3)",
       y = "Mean Total Calories (kcal/day)") +
  theme_minimal() +
  scale_x_discrete(limits = selected_countries) # Maintain the specified order in the x-axis

library(tidyr)

# Filter for Burundi and calculate mean total calories by sex and age group
burundi_caloric_analysis_wide <- filtered_caloric_intake %>%
  filter(iso3 == "USA") %>% # Filter for Burundi
  group_by(age_group, sex) %>%
  summarise(
    average_calories = mean(total_calories_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = sex, values_from = average_calories, 
              values_fill = list(average_calories = NA)) # Fill missing values with NAs

# Print the wide version of the average caloric intake for Burundi
View(burundi_caloric_analysis_wide)



# Load required packages
library(dplyr)
library(readr)
library(ggplot2)

# Set working directory for your R script and clear the environment
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
options(scipen=999)
rm(list = ls())

# Step 1: Load GDD 2018 macronutrient protein data only
protein_gdd <- read_csv("./data/GDD_FinalEstimates_01102022/Country-level estimates/v23_cnty.csv")

# Step 2: Preprocess and filter for 2018 data for protein
protein_gdd_natl <- protein_gdd %>%
  filter(year == 2018, edu == 999, urban == 999, age != 999, female != 999) %>%
  mutate(sex = ifelse(female == 1, "Females", "Males")) %>%
  select(iso3, sex, age, median, lowerci_95, upperci_95) # relevant protein columns

# Step 3: Create age_group variable for protein data
create_age_group <- function(df) {
  df %>% mutate(age_group = case_when(
    age < 1 ~ "0-0.99",
    age < 2 ~ "1-1.99",
    age < 5 ~ "2-4",
    age < 10 ~ "5-9",
    age < 15 ~ "10-14",
    age < 20 ~ "15-19",
    age < 25 ~ "20-24",
    age < 30 ~ "25-29",
    age < 35 ~ "30-34",
    age < 40 ~ "35-39",
    age < 45 ~ "40-44",
    age < 50 ~ "45-49",
    age < 55 ~ "50-54",
    age < 60 ~ "55-59",
    age < 65 ~ "60-64",
    age < 70 ~ "65-69",
    age < 75 ~ "70-74",
    age < 80 ~ "75-79",
    age < 85 ~ "80-84",
    age < 90 ~ "85-89",
    age < 95 ~ "90-94",
    age < 100 ~ "95-99",
    TRUE ~ NA_character_
  ))
}

# Apply the create_age_group function to protein data
protein_gdd_natl <- create_age_group(protein_gdd_natl)

# Step 4: Aggregate means for protein data by iso3, sex, and age_group
agg_protein <- protein_gdd_natl %>%
  group_by(iso3, sex, age_group) %>%
  summarise(
    gdd_mean_protein = mean(median, na.rm = TRUE),
    gdd_lower_protein = mean(lowerci_95, na.rm = TRUE),
    gdd_upper_protein = mean(upperci_95, na.rm = TRUE),
    .groups = "drop"
  )

# Define the countries to filter (as their ISO3 codes)
selected_countries <- c("BDI", "TLS", "NER", "MDG", "NPL", "USA", "GBR", "AUS", "WSM", "KIR")

# Step to Filter the protein intake data for the selected countries
filtered_protein_intake <- agg_protein %>%
  filter(iso3 %in% selected_countries)

# Ensure the order of countries is based on `selected_countries`
filtered_protein_intake$iso3 <- factor(filtered_protein_intake$iso3, levels = selected_countries)

# Check the results for the specified countries
print(filtered_protein_intake)

# (Optional) Summary statistics for the selected countries
summary_filtered_protein_intake <- filtered_protein_intake %>%
  summarise(
    count = n(),
    mean_protein = mean(gdd_mean_protein, na.rm = TRUE),
    median_protein = median(gdd_mean_protein, na.rm = TRUE),
    min_protein = min(gdd_mean_protein, na.rm = TRUE),
    max_protein = max(gdd_mean_protein, na.rm = TRUE),
    sd_protein = sd(gdd_mean_protein, na.rm = TRUE)
  )

# Print summary statistics for protein intake
print(summary_filtered_protein_intake)

# (Optional) Visualization for selected countries
ggplot(filtered_protein_intake, aes(x = iso3, y = gdd_mean_protein, fill = sex)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Average Protein Intake for Selected Countries",
       x = "Country (ISO3)",
       y = "Mean Protein Intake (g/day)") +
  theme_minimal() +
  scale_x_discrete(limits = selected_countries) # Maintain the specified order in the x-axis

