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
  min_waste <- 0.00; max_waste <- 0.35; gdp_midpoint <- 8000; k <- 0.00015
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
    "SSD", "ERI", "MDG"
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
countries_to_plot <- c("MDG", "BRA", "ECU", "CHN", "ITA", "DEU", "USA", "MEX")

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


















































# Step 1: Load nutriR protein distribution data
protein_nutriR <- nutriR::get_dists(nutrients = "Protein")

# Convert age_group to character to avoid factor issues
protein_nutriR <- protein_nutriR %>%
  mutate(age_group = as.character(age_group))

# Duplicate "0-4" nutriR distributions to match GDD's under-5 splits
under5_expansion <- protein_nutriR %>%
  filter(age_group == "0-4") %>%
  mutate(age_group = list(c("0-0.99", "1-1.99", "2-4"))) %>%
  unnest(age_group)

# Bind back to original nutriR
protein_nutriR <- bind_rows(
  protein_nutriR %>% filter(age_group != "0-4"),
  under5_expansion
)

# Step 2: Load GDD 2018 protein data (v23 = total protein intake)
protein_gdd <- read_csv("./data/GDD_FinalEstimates_01102022/Country-level estimates/v23_cnty.csv")

# Step 3: Filter to national-level rows and add sex label
protein_gdd_natl <- protein_gdd %>%
  filter(year == 2018, edu == 999, urban == 999, age != 999, female != 999) %>%
  mutate(sex = ifelse(female == 1, "Females", "Males"))

# Step 4: Define age groups
age_breaks <- c(0, 0.99, 1.99, 4.99, 9.99, 14.99, 19.99, 24.99, 29.99, 34.99,
                39.99, 44.99, 49.99, 54.99, 59.99, 64.99, 69.99, 74.99,
                79.99, 84.99, 89.99, 94.99, 99.99)

age_labels <- c("0-0.99","1-1.99","2-4", "5-9", "10-14", "15-19", "20-24", "25-29",
                "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
                "60-64", "65-69", "70-74", "75-79", "80-84", "85-89",
                "90-94", "95-99")

protein_gdd_natl <- protein_gdd_natl %>%
  mutate(age_group = cut(age, breaks = age_breaks, labels = age_labels, right = TRUE))

# Step 5: Aggregate to one row per iso3–sex–age_group
# gdd_0_4 <- protein_gdd_natl %>%
#   filter(age_group == "0-4", age %in% c(0.5, 1.5, 3.5)) %>%
#   mutate(weight = case_when(age == 0.5 ~ 1, age == 1.5 ~ 1, age == 3.5 ~ 3)) %>%
#   group_by(iso3, sex, age_group) %>%
#   summarise(
#     gdd_mean = weighted.mean(median, weight, na.rm = TRUE),
#     gdd_lower = weighted.mean(lowerci_95, weight, na.rm = TRUE),
#     gdd_upper = weighted.mean(upperci_95, weight, na.rm = TRUE),
#     .groups = "drop"
#   )

gdd_clean <- protein_gdd_natl %>%
  group_by(iso3, sex, age_group) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    iso3, sex, age_group,
    gdd_mean = median,
    gdd_lower = lowerci_95,
    gdd_upper = upperci_95
  )

protein_gdd_agg <- gdd_clean %>%
  arrange(iso3, sex, age_group)

#convert age to character to prevent later matching issues
protein_gdd_agg <- protein_gdd_agg %>%
  mutate(age_group = as.character(age_group))

# -----------------------------
# MATCHING STEPS 1–3
# -----------------------------
gdd_keys <- protein_gdd_agg %>% select(iso3, sex, age_group) %>% distinct()
nutrir_keys <- protein_nutriR %>% select(iso3, sex, age_group) %>% distinct()
nutrir_keys <- protein_nutriR %>%
  mutate(age_group = as.character(age_group)) %>%
  select(iso3, sex, age_group) %>%
  distinct()


# Step 1: Exact match
matched_1_exact <- inner_join(gdd_keys, nutrir_keys, by = c("iso3", "sex", "age_group")) %>%
  mutate(match_iso3 = iso3, match_sex = sex, match_age_group = age_group, source = "exact")

# Step 2: Nearest age in same iso3/sex
find_nearest_age <- function(target, pool) {
  target_mid <- as.numeric(sub("([0-9]+)-.*", "\\1", target))
  pool_mids <- as.numeric(sub("([0-9]+)-.*", "\\1", pool))
  pool[which.min(abs(pool_mids - target_mid))]
}

unmatched_1 <- anti_join(gdd_keys, matched_1_exact, by = c("iso3", "sex", "age_group"))

matched_2_age <- unmatched_1 %>%
  mutate(match = pmap(list(iso3, sex, age_group), function(iso, s, ag) {
    pool <- nutrir_keys %>% filter(iso3 == iso, sex == s)
    if (nrow(pool) == 0) return(NULL)
    nearest_ag <- find_nearest_age(ag, pool$age_group)
    pool %>% filter(age_group == nearest_ag)
  })) %>%
  unnest(match, names_sep = "_") %>%
  mutate(source = "nearest_age")

# Step 3: Opposite sex, same iso3 + age group
unmatched_2 <- anti_join(unmatched_1, matched_2_age, by = c("iso3", "sex", "age_group"))

matched_3_sex <- unmatched_2 %>%
  mutate(match = pmap(list(iso3, sex, age_group), function(iso, s, ag) {
    opposite <- ifelse(s == "Males", "Females", "Males")
    nutrir_keys %>%
      filter(iso3 == iso, sex == opposite, as.character(age_group) == as.character(ag))
  })) %>%
  unnest(match, names_sep = "_") %>%
  mutate(source = "opposite_sex")

all_matches_1_3 <- bind_rows(
  matched_1_exact %>% select(iso3, sex, age_group, match_iso3, match_sex, match_age_group, source),
  matched_2_age,
  matched_3_sex
)

# -----------------------------
# STEP 4: Fallback to nearest country with known shape
# -----------------------------

# Re-define helper in case not scoped
find_nearest_age <- function(target, pool) {
  target_mid <- as.numeric(sub("([0-9]+)-.*", "\\1", target))
  pool_mids <- as.numeric(sub("([0-9]+)-.*", "\\1", pool))
  pool[which.min(abs(pool_mids - target_mid))]
}

# Identify unmatched from Steps 1–3
unmatched_3 <- anti_join(gdd_keys, all_matches_1_3, by = c("iso3", "sex", "age_group"))

# Build GDD country means
gdd_country_means <- protein_gdd_agg %>%
  group_by(iso3) %>%
  summarise(protein_mean = mean(gdd_mean, na.rm = TRUE), .groups = "drop")

# Euclidean distance matrix between countries
dist_mat <- as.matrix(dist(column_to_rownames(gdd_country_means, "iso3")))

# Identify fallback iso3 (closest country w/ matched distributions)
shape_donor_iso3s <- unique(all_matches_1_3$match_iso3)
iso_match_key <- tibble(iso3 = unique(unmatched_3$iso3)) %>%
  rowwise() %>%
  mutate(fallback_iso3 = names(sort(dist_mat[iso3, shape_donor_iso3s]))[1]) %>%
  ungroup()

# Perform fallback match with 4 levels of logic
step4_fallback <- unmatched_3 %>%
  left_join(iso_match_key, by = "iso3") %>%
  rowwise() %>%
  mutate(
    match = list({
      fallback <- fallback_iso3
      target_sex <- sex
      target_age <- age_group
      opposite_sex <- ifelse(sex == "Males", "Females", "Males")
      pool <- all_matches_1_3
      
      # 1. Exact sex + age
      m1 <- pool %>%
        filter(
          match_iso3 == fallback,
          match_sex == target_sex,
          as.character(match_age_group) == as.character(target_age)
        )
      if (nrow(m1) > 0) {
        m1 <- m1 %>% slice(1) %>% mutate(source = "nearest_country_exact")
        return_value <- m1
      } else {
        # 2. Nearest age, same sex
        pool_same_sex <- pool %>%
          filter(match_iso3 == fallback, match_sex == target_sex)
        if (nrow(pool_same_sex) > 0) {
          nearest_age <- find_nearest_age(target_age, pool_same_sex$match_age_group)
          m2 <- pool_same_sex %>% filter(match_age_group == nearest_age)
          if (nrow(m2) > 0) {
            m2 <- m2 %>% slice(1) %>% mutate(source = "nearest_country_nearest_age")
            return_value <- m2
          } else {
            # 3. Same age, opposite sex
            m3 <- pool %>%
              filter(match_iso3 == fallback, match_sex == opposite_sex, match_age_group == target_age)
            if (nrow(m3) > 0) {
              m3 <- m3 %>% slice(1) %>% mutate(source = "nearest_country_opposite_sex")
              return_value <- m3
            } else {
              # 4. Nearest age + opposite sex
              pool_opp <- pool %>%
                filter(match_iso3 == fallback, match_sex == opposite_sex)
              if (nrow(pool_opp) > 0) {
                nearest_age_opp <- find_nearest_age(target_age, pool_opp$match_age_group)
                m4 <- pool_opp %>% filter(match_age_group == nearest_age_opp)
                if (nrow(m4) > 0) {
                  m4 <- m4 %>% slice(1) %>% mutate(source = "ultimate_fallback_opposite_sex_nearest_age")
                  return_value <- m4
                } else {
                  return_value <- NULL
                }
              } else {
                return_value <- NULL
              }
            }
          }
        } else {
          # Fall back directly to opposite sex
          m3 <- pool %>%
            filter(match_iso3 == fallback, match_sex == opposite_sex, match_age_group == target_age)
          if (nrow(m3) > 0) {
            m3 <- m3 %>% slice(1) %>% mutate(source = "nearest_country_opposite_sex")
            return_value <- m3
          } else {
            pool_opp <- pool %>%
              filter(match_iso3 == fallback, match_sex == opposite_sex)
            if (nrow(pool_opp) > 0) {
              nearest_age_opp <- find_nearest_age(target_age, pool_opp$match_age_group)
              m4 <- pool_opp %>% filter(match_age_group == nearest_age_opp)
              if (nrow(m4) > 0) {
                m4 <- m4 %>% slice(1) %>% mutate(source = "ultimate_fallback_opposite_sex_nearest_age")
                return_value <- m4
              } else {
                return_value <- NULL
              }
            } else {
              return_value <- NULL
            }
          }
        }
      }
      
      return_value
    })
  ) %>%
  unnest(match, names_sep = "_", names_repair = "unique") %>%
  ungroup()

all_matches_final <- bind_rows(all_matches_1_3, step4_fallback)
final_unmatched <- anti_join(gdd_keys, all_matches_final, by = c("iso3", "sex", "age_group"))

cat("✅ Total matched:", nrow(all_matches_final), "\n")
cat("❌ Still unmatched after Step 4:", nrow(final_unmatched), "\n")

all_matches_final <- all_matches_final %>%
  mutate(source_final = coalesce(source, match_source)) %>%
  select(-source, -match_source)  # Optionally drop old columns

#drop more redundant columns
all_matches_final <- all_matches_final %>%
  select(-match_match_iso3, -match_match_sex, -match_match_age_group)


####NOW WE HAVE FULL MATCHING DICTIONARY#########
##########################################################################


# Merge GDD mean protein intakes with matched nutriR distribution shapes.
# Each row represents a country–sex–age_group group from GDD, matched to a 
# distribution from nutriR (via exact or 4-fold fallback logic) to allow reconstruction
# of full intake distributions and downstream estimation of inadequacy.
gdd_distributions <- all_matches_final %>%
  left_join(protein_gdd_agg, by = c("iso3", "sex", "age_group")) %>%
  left_join(protein_nutriR, 
            by = c("match_iso3" = "iso3", "match_sex" = "sex", "match_age_group" = "age_group"))

#create leaner version by keeping only needed variables
#only keep variables needed to get full distributions
gdd_distributions_lean <- gdd_distributions %>%
  select(
    iso3, sex, age_group,      # group identifiers
    gdd_mean, gdd_lower, gdd_upper,  # mean intake and uncertainty
    best_dist,                 # distribution type
    cv                         # matched coefficient of variation
  )

#we now have fully specified distributions for every subgroup!