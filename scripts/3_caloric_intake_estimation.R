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
# PART 4: MODELING THREE INDEPENDENT WASTE SCENARIOS
# ==========================================================================

# Start with the master country-level data (as before)
country_level_data <- all_countries_needed %>%
  left_join(gdp_final_complete, by = "iso3") %>%
  left_join(calories_final_complete, by = "iso3")

# --- Step 4a: Define THREE independent waste functions ---
# This gives you full control to tune each scenario's parameters separately.

# 1. LOW Waste Scenario Function
calculate_waste_pct_low <- function(gdp) {
  min_waste    <- 0.03
  max_waste    <- 0.25
  gdp_midpoint <- 12000 # Rises a bit later
  k            <- 0.0001
  waste_percentage <- min_waste + (max_waste - min_waste) / (1 + exp(-k * (gdp - gdp_midpoint)))
  return(waste_percentage)
}

# 2. MEDIUM Waste Scenario Function (Your Current/Central Model)
calculate_waste_pct_medium <- function(gdp) {
  min_waste    <- 0.05
  max_waste    <- 0.30
  gdp_midpoint <- 10000 # Rises earlier
  k            <- 0.0001
  waste_percentage <- min_waste + (max_waste - min_waste) / (1 + exp(-k * (gdp - gdp_midpoint)))
  return(waste_percentage)
}

# 3. HIGH Waste Scenario Function
calculate_waste_pct_high <- function(gdp) {
  min_waste    <- 0.05
  max_waste    <- 0.35
  gdp_midpoint <- 8000 # Rises very early
  k            <- 0.00015 # Rises a bit steeper
  waste_percentage <- min_waste + (max_waste - min_waste) / (1 + exp(-k * (gdp - gdp_midpoint)))
  return(waste_percentage)
}


# --- Step 4b: Apply all three functions to create three consumption estimates ---
cat("\n--- Applying Low, Medium, and High waste models... ---\n")

country_level_data_scenarios <- country_level_data %>%
  mutate(
    # Calculate the resulting national consumption for each scenario
    # We apply each specific function to the GDP column.
    est_kcal_consumed_low    = fao_des_kcal * (1 - calculate_waste_pct_low(gdp_per_cap_ppp)),
    est_kcal_consumed_medium = fao_des_kcal * (1 - calculate_waste_pct_medium(gdp_per_cap_ppp)),
    est_kcal_consumed_high   = fao_des_kcal * (1 - calculate_waste_pct_high(gdp_per_cap_ppp))
  )

# --- Step 4c: Visualize the three scenarios to verify ---
# (This visualization code remains the same as my previous response, it will now
#  just use the newly calculated columns to generate the plot)
waste_scenarios_long <- country_level_data_scenarios %>%
  select(iso3, gdp_per_cap_ppp, starts_with("est_kcal_consumed_")) %>%
  pivot_longer(
    cols = starts_with("est_kcal_consumed_"),
    names_to = "scenario",
    names_prefix = "est_kcal_consumed_",
    values_to = "kcal_consumed"
  ) %>%
  # We can also calculate the waste % for plotting
  left_join(country_level_data %>% select(iso3, fao_des_kcal), by = "iso3") %>%
  mutate(waste_pct = 1 - (kcal_consumed / fao_des_kcal)) %>%
  mutate(scenario = factor(scenario, levels = c("low", "medium", "high")))

waste_scenarios_plot <- ggplot(waste_scenarios_long, aes(x = gdp_per_cap_ppp, y = waste_pct, color = scenario)) +
  geom_point(alpha = 0.4) +
  scale_x_continuous(name = "GDP per capita, PPP (2018)", labels = scales::dollar_format(prefix = "$", scale = 1/1000, suffix = "k")) +
  scale_y_continuous(name = "Estimated Downstream Waste %", labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(
    name = "Waste Scenario",
    values = c("low" = "darkgreen", "medium" = "darkorange", "high" = "darkred"),
    labels = c("Low (25% max)", "Medium (30% max)", "High (35% max)")
  ) +
  labs(
    title = "Three Modeled Scenarios for Downstream Food Waste",
    subtitle = "Each scenario uses an independently tuned logistic function"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(waste_scenarios_plot)

cat("--- Three consumption scenarios created successfully. ---\n")

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

# ==========================================================================
# PART 6 (REVISED): DISAGGREGATE ALL THREE CALORIC SCENARIOS
# ==========================================================================

# --- Step 6a: Create the master analytical dataframe ---
# Join the 3 country-level consumption estimates to the main stratum-level data.
cat("\n--- Joining 3 consumption scenarios to master data... ---\n")
analysis_data_full <- master_data %>%
  left_join(
    # Select the iso3 and the three new consumption estimate columns
    country_level_data_scenarios %>% select(iso3, est_kcal_consumed_low, est_kcal_consumed_medium, est_kcal_consumed_high),
    by = "iso3"
  )

# --- Step 6b: Perform disaggregation for EACH scenario separately ---
# This creates three new columns for the final stratum-level mean consumption.
cat("--- Disaggregating all 3 scenarios to the stratum level... ---\n")
disaggregated_data <- analysis_data_full %>%
  # For each country...
  group_by(iso3) %>%
  
  # Calculate the single adjustment factor (this is the same for all scenarios)
  mutate(
    avg_country_eer = weighted.mean(eer_kcal_marco_mean, w = population, na.rm = TRUE),
    adjustment_factor = eer_kcal_marco_mean / avg_country_eer
  ) %>%
  ungroup() %>%
  
  # Now, apply the factor to each of the three national consumption estimates
  mutate(
    kcal_consumed_stratum_low    = est_kcal_consumed_low * adjustment_factor,
    kcal_consumed_stratum_medium = est_kcal_consumed_medium * adjustment_factor,
    kcal_consumed_stratum_high   = est_kcal_consumed_high * adjustment_factor
  )

# --- Verification ---
cat("--- Disaggregation for 3 scenarios complete. Glimpsing the result: ---\n")
glimpse(disaggregated_data)

# ==========================================================================
# PART 7 (REVISED): VISUALIZE ALL THREE DISAGGREGATED SCENARIOS
# ==========================================================================
library(ggplot2)

# --- Step 7a: Prepare the data for plotting ---
# We need to reshape the data from wide to long to plot all three scenarios.

cat("\n--- Preparing data to visualize the 3 disaggregated scenarios... ---\n")
plot_data_long <- disaggregated_data %>%
  # Select the key columns, including the three new stratum consumption columns
  select(iso3, sex, age_group, starts_with("kcal_consumed_stratum_")) %>%
  # Pivot these three columns into a long format
  pivot_longer(
    cols = starts_with("kcal_consumed_stratum_"),
    names_to = "waste_scenario",
    names_prefix = "kcal_consumed_stratum_",
    values_to = "kcal_consumed"
  ) %>%
  # Make the scenario names a clean, ordered factor for the plot legend
  mutate(
    waste_scenario = factor(waste_scenario, levels = c("low", "medium", "high"))
  )

# --- Step 7b: Create the multi-scenario plot ---
# Define the countries and age order for the plot
countries_to_plot <- c("MDG", "BRA", "ECU", "CHN", "ITA", "DEU", "USA", "NGA")
age_group_levels_ordered <- c("0-0.99", "1-4", "5-9", "10-14", "15-19",
                              "20-24", "25-29", "30-34", "35-39",
                              "40-44", "45-49", "50-54", "55-59",
                              "60-64", "65-69", "70-74", "75-79",
                              "80-84", "85-89", "90-94", "95-99")

# Filter the long data and create the final plot object
disaggregated_scenarios_plot <- plot_data_long %>%
  filter(iso3 %in% countries_to_plot) %>%
  mutate(age_group_f = factor(age_group, levels = age_group_levels_ordered)) %>%
  
  ggplot(aes(x = age_group_f, y = kcal_consumed, group = interaction(sex, waste_scenario), color = sex, linetype = waste_scenario)) +
  
  # Draw the lines for each scenario
  geom_line(linewidth = 1.1, alpha = 0.9) +
  
  # Use facets to create a panel for each country
  facet_wrap(~iso3, scales = "free_y", ncol = 4) +
  
  # Define colors for sex and linetypes for scenario
  scale_color_brewer(palette = "Set1") +
  scale_linetype_manual(
    name = "Waste Scenario",
    values = c("low" = "dotted", "medium" = "solid", "high" = "dashed"),
    labels = c("Low Waste", "Medium Waste", "High Waste")
  ) +
  
  # Add informative labels
  labs(
    title = "Estimated Caloric Consumption Across the Lifespan (3 Scenarios)",
    subtitle = "For a diverse sample of countries, showing low, medium, and high waste models",
    x = "Age Group",
    y = "Estimated Caloric Consumption (kcal/day)",
    color = "Sex"
  ) +
  
  # Use a clean theme
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 12)
  )

# --- Step 7c: Print the plot ---
cat("--- Generating the final visualization... ---\n")
print(disaggregated_scenarios_plot)


# ============================
# PART 8, CALORIES (REVISED): Match means to NutriR Energy shapes for 3 scenarios
# ============================

# ============================================================
# ENERGY (kcal) — run the same matching logic 3 times (low/medium/high)
# with pre-matching CV repair at the source (scenario-invariant)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(rlang)
})

# --- Normalizers (same style as earlier) ---
norm_sex <- function(x) {
  dplyr::case_when(
    x %in% c("Males","Male","M") ~ "Males",
    x %in% c("Females","Female","F") ~ "Females",
    TRUE ~ x
  )
}
norm_age <- function(a) {
  a <- as.character(a)
  a <- str_replace_all(a, "\\s", "")
  a <- str_replace_all(a, "–", "-")
  a
}

# nearest-age helper (used multiple times)
find_nearest_age <- function(target, pool) {
  target_mid <- suppressWarnings(as.numeric(sub("([0-9]+)-.*", "\\1", target)))
  pool_mids  <- suppressWarnings(as.numeric(sub("([0-9]+)-.*", "\\1", pool)))
  if (length(pool) == 0 || all(is.na(pool_mids))) return(NA_character_)
  pool[which.min(abs(pool_mids - target_mid))]
}

# ------------------------------------------------------------
# A) Build NutriR Energy shapes once; expand 0-4; synthesize 1-4
# ------------------------------------------------------------
energy_nutriR_raw <- nutriR::get_dists(nutrients = "Energy") %>%
  mutate(
    sex = norm_sex(sex),
    age_group = as.character(age_group)
  )

# Expand 0-4 into fine bins
u5_expand <- energy_nutriR_raw %>%
  filter(age_group == "0-4") %>%
  mutate(age_group = list(c("0-0.99","1-1.99","2-4"))) %>%
  unnest(age_group)

energy_nutriR_fine <- bind_rows(
  energy_nutriR_raw %>% filter(age_group != "0-4"),
  u5_expand
)

# Create synthetic 1-4 by copying shapes from 2-4
energy_nutriR_u14 <- energy_nutriR_fine %>%
  filter(age_group == "2-4") %>%
  mutate(age_group = "1-4")

energy_shapes_base <- bind_rows(energy_nutriR_fine, energy_nutriR_u14) %>%
  mutate(age_group = norm_age(age_group)) %>%
  group_by(iso3, sex, age_group) %>%
  slice(1) %>%
  ungroup()

# ------------------------------------------------------------
# B) PRE-MATCH CV REPAIR (scenario-invariant; fixes tiny CVs at the source)
# ------------------------------------------------------------
cv_floor <- 0.05  # minimum acceptable CV (tweakable)

age_levels_for_order <- c("0-0.99","1-1.99","1-4","2-4","5-9","10-14","15-19","20-24","25-29","30-34",
                          "35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79",
                          "80-84","85-89","90-94","95-99")

energy_shapes <- (function(es, floor_cv) {
  es <- es %>%
    mutate(age_group_chr = as.character(age_group),
           age_group_ord = factor(age_group_chr, levels = age_levels_for_order))
  
  to_fix <- es %>% filter(!is.na(cv) & cv < floor_cv)
  if (nrow(to_fix) == 0) {
    message(sprintf("🔎 CV repair: nothing to fix (no cv < %.3f).", floor_cv))
    return(es %>% select(-age_group_chr, -age_group_ord))
  }
  
  find_replacement_row <- function(iso3_i, sex_i, age_i) {
    pool_ok <- es %>% filter(iso3 == iso3_i, sex == sex_i, !is.na(cv) & cv >= floor_cv)
    if (nrow(pool_ok) > 0) {
      nearest_age <- find_nearest_age(age_i, as.character(pool_ok$age_group))
      pool_ok %>% filter(as.character(age_group) == nearest_age) %>% slice(1)
    } else {
      pool_any <- es %>% filter(iso3 == iso3_i, sex == sex_i)
      if (nrow(pool_any) > 0) {
        nearest_age <- find_nearest_age(age_i, as.character(pool_any$age_group))
        pool_any %>% filter(as.character(age_group) == nearest_age) %>% slice(1)
      } else NULL
    }
  }
  
  replacements <- to_fix %>%
    rowwise() %>%
    mutate(.rep = list(find_replacement_row(iso3, sex, as.character(age_group)))) %>%
    ungroup()
  
  es_fixed <- es
  for (i in seq_len(nrow(replacements))) {
    row_i <- replacements[i,]
    rep_i <- row_i$.rep[[1]]
    if (!is.null(rep_i) && nrow(rep_i) == 1) {
      es_fixed <- es_fixed %>%
        mutate(
          best_dist = if_else(iso3 == row_i$iso3 & sex == row_i$sex &
                                as.character(age_group) == as.character(row_i$age_group),
                              rep_i$best_dist, best_dist),
          cv = if_else(iso3 == row_i$iso3 & sex == row_i$sex &
                         as.character(age_group) == as.character(row_i$age_group),
                       rep_i$cv, cv)
        )
    } else {
      es_fixed <- es_fixed %>%
        mutate(cv = if_else(iso3 == row_i$iso3 & sex == row_i$sex &
                              as.character(age_group) == as.character(row_i$age_group) & cv < floor_cv,
                            floor_cv, cv))
    }
  }
  es_fixed %>% select(-age_group_chr, -age_group_ord)
})(energy_shapes_base, cv_floor)

nutrir_keys <- energy_shapes %>%
  select(iso3, sex, age_group) %>% distinct()

# ------------------------------------------------------------
# C) Scenario-matching function (runs full 4-step matching + post-join repairs)
# ------------------------------------------------------------
match_shapes_for_scenario <- function(stratum_col, country_col, scenario_name) {
  
  # (1) Build scenario-specific means
  calories_means <- disaggregated_data %>%
    transmute(
      iso3,
      sex = norm_sex(sex),
      age_group = norm_age(age_group),
      kcal_mean = .data[[stratum_col]]
    ) %>%
    arrange(iso3, sex, age_group)
  
  cal_keys <- calories_means %>%
    select(iso3, sex, age_group) %>% distinct()
  
  # (2) Steps 1–3 (exact; nearest age same sex; opposite sex same age)
  matched_1_exact <- inner_join(cal_keys, nutrir_keys, by = c("iso3","sex","age_group")) %>%
    mutate(match_iso3 = iso3, match_sex = sex, match_age_group = age_group, source = "exact")
  
  unmatched_1 <- anti_join(cal_keys, matched_1_exact, by = c("iso3","sex","age_group"))
  
  matched_2_age <- unmatched_1 %>%
    mutate(match = pmap(list(iso3, sex, age_group), function(iso, s, ag) {
      pool <- nutrir_keys %>% filter(iso3 == iso, sex == s)
      if (nrow(pool) == 0) return(NULL)
      nearest_ag <- find_nearest_age(ag, pool$age_group)
      pool %>% filter(age_group == nearest_ag)
    })) %>%
    unnest(match, names_sep = "_") %>%
    mutate(source = "nearest_age")
  
  unmatched_2 <- anti_join(unmatched_1, matched_2_age, by = c("iso3","sex","age_group"))
  
  matched_3_sex <- unmatched_2 %>%
    mutate(match = pmap(list(iso3, sex, age_group), function(iso, s, ag) {
      opposite <- ifelse(s == "Males", "Females", "Males")
      nutrir_keys %>% filter(iso3 == iso, sex == opposite, as.character(age_group) == as.character(ag))
    })) %>%
    unnest(match, names_sep = "_") %>%
    mutate(source = "opposite_sex")
  
  all_matches_1_3 <- bind_rows(
    matched_1_exact %>% select(iso3, sex, age_group, match_iso3, match_sex, match_age_group, source),
    matched_2_age,
    matched_3_sex
  )
  
  # (3) Step 4 (nearest country) — distance in scenario's national kcal
  ck <- country_level_data_scenarios %>%
    select(iso3, kcal = all_of(country_col))
  stopifnot(!any(is.na(ck$kcal)))
  
  m <- as.matrix(ck$kcal)
  rownames(m) <- ck$iso3
  dist_mat <- as.matrix(stats::dist(m))
  
  unmatched_3 <- anti_join(cal_keys, all_matches_1_3, by = c("iso3","sex","age_group"))
  shape_donor_iso3s <- unique(all_matches_1_3$match_iso3)
  
  if (length(shape_donor_iso3s) == 0 || nrow(unmatched_3) == 0) {
    step4_fallback <- tibble(
      iso3 = character(), sex = character(), age_group = character(),
      match_iso3 = character(), match_sex = character(), match_age_group = character(),
      source = character()
    )
  } else {
    # for each unmatched iso3, choose nearest donor iso3 by dist on kcal
    iso_match_key <- tibble(iso3 = unique(unmatched_3$iso3)) %>%
      rowwise() %>%
      mutate(
        fallback_iso3 = {
          neigh <- shape_donor_iso3s[shape_donor_iso3s %in% rownames(dist_mat)]
          # handle any missing rows defensively
          if (iso3 %in% rownames(dist_mat) && length(neigh) > 0) {
            names(sort(dist_mat[iso3, neigh]))[1]
          } else NA_character_
        }
      ) %>%
      ungroup()
    
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
          return_value <- NULL
          
          if (!is.na(fallback)) {
            m1 <- pool %>% filter(match_iso3 == fallback,
                                  match_sex == target_sex,
                                  as.character(match_age_group) == as.character(target_age))
            if (nrow(m1) > 0) {
              return_value <- m1 %>% slice(1) %>% mutate(source = "nearest_country_exact")
            } else {
              pool_same <- pool %>% filter(match_iso3 == fallback, match_sex == target_sex)
              if (nrow(pool_same) > 0) {
                nearest_age <- find_nearest_age(target_age, pool_same$match_age_group)
                m2 <- pool_same %>% filter(match_age_group == nearest_age)
                if (nrow(m2) > 0) {
                  return_value <- m2 %>% slice(1) %>% mutate(source = "nearest_country_nearest_age")
                }
              }
              if (is.null(return_value)) {
                m3 <- pool %>% filter(match_iso3 == fallback,
                                      match_sex == opposite_sex,
                                      as.character(match_age_group) == as.character(target_age))
                if (nrow(m3) > 0) {
                  return_value <- m3 %>% slice(1) %>% mutate(source = "nearest_country_opposite_sex")
                }
              }
              if (is.null(return_value)) {
                pool_opp <- pool %>% filter(match_iso3 == fallback, match_sex == opposite_sex)
                if (nrow(pool_opp) > 0) {
                  nearest_age_opp <- find_nearest_age(target_age, pool_opp$match_age_group)
                  m4 <- pool_opp %>% filter(match_age_group == nearest_age_opp)
                  if (nrow(m4) > 0) {
                    return_value <- m4 %>% slice(1) %>% mutate(source = "ultimate_fallback_opposite_sex_nearest_age")
                  }
                }
              }
            }
          }
          return_value
        })
      ) %>%
      unnest(match, names_sep = "_", names_repair = "unique") %>%
      ungroup()
  }
  
  cal_matches_final <- bind_rows(all_matches_1_3, step4_fallback) %>%
    mutate(source_final = coalesce(source, match_source)) %>%
    select(-source, -match_source, -match_match_iso3, -match_match_sex, -match_match_age_group)
  
  final_unmatched_cal <- anti_join(cal_keys, cal_matches_final, by = c("iso3","sex","age_group"))
  cat(sprintf("✅ [%s] Total matched: %d\n", scenario_name, nrow(cal_matches_final)))
  cat(sprintf("❌ [%s] Still unmatched after Step 4: %d\n", scenario_name, nrow(final_unmatched_cal)))
  
  # (4) Join means + repaired shapes/CV
  calories_distributions <- cal_matches_final %>%
    left_join(calories_means, by = c("iso3","sex","age_group")) %>%
    left_join(energy_shapes,
              by = c("match_iso3"="iso3","match_sex"="sex","match_age_group"="age_group"))
  
  # (5) Post-join repair (nearest age same sex within donor, then opposite sex + nearest age)
  missing_shape_n <- sum(is.na(calories_distributions$best_dist))
  cat(sprintf("🔎 [%s] Missing shapes after first join: %d\n", scenario_name, missing_shape_n))
  
  if (missing_shape_n > 0) {
    avail_by_iso_sex <- energy_shapes %>%
      distinct(iso3, sex, age_group) %>%
      group_by(iso3, sex) %>%
      summarise(avail_ages = list(sort(unique(as.character(age_group)))), .groups = "drop")
    
    # A) Same sex nearest age
    missing_keys <- calories_distributions %>%
      filter(is.na(best_dist)) %>%
      distinct(match_iso3, match_sex, match_age_group)
    
    remap_same <- missing_keys %>%
      left_join(avail_by_iso_sex, by = c("match_iso3"="iso3","match_sex"="sex")) %>%
      mutate(new_age = ifelse(lengths(avail_ages) > 0,
                              map2_chr(match_age_group, avail_ages, find_nearest_age),
                              NA_character_)) %>%
      select(match_iso3, match_sex, match_age_group, new_age)
    
    cal_matches_final_v2 <- cal_matches_final %>%
      left_join(remap_same, by = c("match_iso3","match_sex","match_age_group")) %>%
      mutate(match_age_group = if_else(!is.na(new_age), new_age, match_age_group)) %>%
      select(-new_age)
    
    calories_distributions_v2 <- cal_matches_final_v2 %>%
      left_join(calories_means, by = c("iso3","sex","age_group")) %>%
      left_join(energy_shapes,
                by = c("match_iso3"="iso3","match_sex"="sex","match_age_group"="age_group"))
    
    na_after_A <- sum(is.na(calories_distributions_v2$best_dist))
    cat(sprintf("🛠  [%s] Repair A (nearest age same sex) — remaining missing: %d\n", scenario_name, na_after_A))
    
    # B) Opposite sex nearest age (only those still missing)
    if (na_after_A > 0) {
      opp <- function(s) ifelse(s == "Males", "Females", "Males")
      
      still_missing <- calories_distributions_v2 %>%
        filter(is.na(best_dist)) %>%
        distinct(match_iso3, match_sex, match_age_group) %>%
        mutate(match_sex_opp = opp(match_sex))
      
      remap_opp <- still_missing %>%
        left_join(avail_by_iso_sex, by = c("match_iso3"="iso3","match_sex_opp"="sex")) %>%
        mutate(new_age_opp = ifelse(lengths(avail_ages) > 0,
                                    map2_chr(match_age_group, avail_ages, find_nearest_age),
                                    NA_character_)) %>%
        transmute(
          match_iso3,
          old_match_sex = match_sex,
          old_match_age_group = match_age_group,
          match_sex = match_sex_opp,
          match_age_group = new_age_opp
        )
      
      cal_matches_final_v3 <- cal_matches_final_v2 %>%
        left_join(
          remap_opp,
          by = c("match_iso3",
                 "match_sex" = "old_match_sex",
                 "match_age_group" = "old_match_age_group")
        ) %>%
        mutate(
          match_sex       = if_else(!is.na(match_sex.y), match_sex.y, match_sex.x),
          match_age_group = if_else(!is.na(match_age_group.y), match_age_group.y, match_age_group.x)
        ) %>%
        select(-ends_with(".x"), -ends_with(".y"))
      
      calories_distributions_v3 <- cal_matches_final_v3 %>%
        left_join(calories_means, by = c("iso3","sex","age_group")) %>%
        left_join(energy_shapes,
                  by = c("match_iso3"="iso3","match_sex"="sex","match_age_group"="age_group"))
      
      na_after_B <- sum(is.na(calories_distributions_v3$best_dist))
      cat(sprintf("🛠  [%s] Repair B (nearest age opposite sex) — remaining missing: %d\n", scenario_name, na_after_B))
      
      if (na_after_B < na_after_A) {
        calories_distributions <- calories_distributions_v3
      } else {
        calories_distributions <- calories_distributions_v2
      }
    } else {
      calories_distributions <- calories_distributions_v2
    }
  }
  
  # (6) Lean output + assertions
  calories_distributions_lean <- calories_distributions %>%
    select(
      iso3, sex, age_group,      # identifiers
      kcal_mean,                 # scenario-specific mean
      best_dist,                 # distribution type
      cv,                        # coefficient of variation (repaired at source)
      source_final               # diagnostic of where the shape came from
    ) %>%
    mutate(waste_scenario = scenario_name)
  
  stopifnot(!any(is.na(calories_distributions_lean$kcal_mean)))
  stopifnot(!any(is.na(calories_distributions_lean$best_dist)))
  stopifnot(!any(is.na(calories_distributions_lean$cv)))
  
  calories_distributions_lean
}

# ------------------------------------------------------------
# D) Run the matcher for LOW, MEDIUM, HIGH and save results
# ------------------------------------------------------------
res_low    <- match_shapes_for_scenario(
  stratum_col = "kcal_consumed_stratum_low",
  country_col = "est_kcal_consumed_low",
  scenario_name = "low"
)
res_med    <- match_shapes_for_scenario(
  stratum_col = "kcal_consumed_stratum_medium",
  country_col = "est_kcal_consumed_medium",
  scenario_name = "medium"
)
res_high   <- match_shapes_for_scenario(
  stratum_col = "kcal_consumed_stratum_high",
  country_col = "est_kcal_consumed_high",
  scenario_name = "high"
)

calories_distributions_all <- bind_rows(res_low, res_med, res_high)

# Build a single wide table: one row per (iso3, sex, age_group),
# scenario-invariant shape columns, plus 3 kcal columns (low/medium/high).

# (1) Wide kcal means
wide_kcal <- calories_distributions_all %>%
  select(iso3, sex, age_group, waste_scenario, kcal_mean) %>%
  tidyr::pivot_wider(
    names_from  = waste_scenario,
    values_from = kcal_mean,
    names_prefix = "kcal_mean_"
  )

# (2) Scenario-invariant shape/CV (repaired at source) — keep one copy
shape_cv_unique <- calories_distributions_all %>%
  select(iso3, sex, age_group, best_dist, cv) %>%
  distinct()

# (3) Join into final wide output
final_calorie_distributions_wide <- shape_cv_unique %>%
  left_join(wide_kcal, by = c("iso3","sex","age_group"))

# --- Flip kcal_mean_low and kcal_mean_high column names safely ---

# Ensure the columns exist (create if missing so the swap doesn't error)
for (nm in c("kcal_mean_low", "kcal_mean_high")) {
  if (!nm %in% names(final_calorie_distributions_wide)) {
    final_calorie_distributions_wide[[nm]] <- NA_real_
  }
}

# flip low and high scenarios to signify caloric scenario instead of 
# waste scenario with naming
final_calorie_distributions_wide_flipped <- final_calorie_distributions_wide %>%
  mutate(
    .kcal_low_tmp  = .data$kcal_mean_low,
    kcal_mean_low  = .data$kcal_mean_high,
    kcal_mean_high = .kcal_low_tmp
  ) %>%
  select(-.kcal_low_tmp)

# Overwrite the saved file with flipped names
saveRDS(final_calorie_distributions_wide_flipped, file = "./output/final_calorie_distributions_wide.rds")

cat("✅ Saved file with flipped column names\n")
cat("  - kcal_mean_low ↔ kcal_mean_high swapped in final_calorie_distributions_wide.rds\n")
