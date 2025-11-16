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




# ============================
# PART 8, CALORIES: Match in-memory means to NutriR Energy shapes (Option A)
# ============================

# ============================================================
# ENERGY (kcal) — match NutriR shapes/CVs via same 4-step algo
# with pre-matching CV repair (fix tiny BFA-F 20–49 CVs at source)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
})

# --- Normalizers (same style as protein script) ---
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

# nearest-age helper used multiple times
find_nearest_age <- function(target, pool) {
  target_mid <- suppressWarnings(as.numeric(sub("([0-9]+)-.*", "\\1", target)))
  pool_mids  <- suppressWarnings(as.numeric(sub("([0-9]+)-.*", "\\1", pool)))
  if (length(pool) == 0 || all(is.na(pool_mids))) return(NA_character_)
  pool[which.min(abs(pool_mids - target_mid))]
}

# ------------------------------------------------------------
# 0) Build Energy means at iso3/sex/age_group (already 0-0.99, 1-4, 5-9, ...)
# ------------------------------------------------------------
calories_means <- disaggregated_data %>%
  transmute(
    iso3,
    sex = norm_sex(sex),
    age_group = norm_age(age_group),       # includes "1-4" in your in-memory data
    kcal_mean = kcal_consumed_stratum
  ) %>%
  arrange(iso3, sex, age_group)

cal_keys <- calories_means %>% select(iso3, sex, age_group) %>% distinct()

# ------------------------------------------------------------
# 1) NutriR Energy shapes; expand 0-4 to 0-0.99/1-1.99/2-4; synthesize "1-4" from "2-4"
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

# Create synthetic 1-4 by copying shapes from 2-4 (same convention as protein script)
energy_nutriR_u14 <- energy_nutriR_fine %>%
  filter(age_group == "2-4") %>%
  mutate(age_group = "1-4")

energy_shapes_base <- bind_rows(energy_nutriR_fine, energy_nutriR_u14) %>%
  mutate(age_group = norm_age(age_group)) %>%
  group_by(iso3, sex, age_group) %>%
  slice(1) %>%     # de-dup
  ungroup()

# ------------------------------------------------------------
# 1B) PRE-MATCH CV REPAIR (fix tiny CVs by re-borrowing within the same donor)
# ------------------------------------------------------------
cv_floor <- 0.05  # change if you want a different minimum acceptable CV

age_levels_for_order <- c("0-0.99","1-1.99","1-4","2-4","5-9","10-14","15-19","20-24","25-29","30-34",
                          "35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79",
                          "80-84","85-89","90-94","95-99")

energy_shapes_fix <- (function(es, floor_cv) {
  es <- es %>%
    mutate(age_group_chr = as.character(age_group),
           age_group_ord = factor(age_group_chr, levels = age_levels_for_order))
  
  # rows that need fixing
  to_fix <- es %>% filter(!is.na(cv) & cv < floor_cv)
  
  if (nrow(to_fix) == 0) {
    message(sprintf("🔎 CV repair: nothing to fix (no cv < %.3f).", floor_cv))
    return(es %>% select(-age_group_chr, -age_group_ord))
  }
  
  # helper: nearest acceptable age within same (iso3,sex)
  find_replacement_row <- function(iso3_i, sex_i, age_i) {
    pool_ok <- es %>% filter(iso3 == iso3_i, sex == sex_i, !is.na(cv) & cv >= floor_cv)
    if (nrow(pool_ok) > 0) {
      # nearest by age midpoint
      nearest_age <- find_nearest_age(age_i, as.character(pool_ok$age_group))
      pool_ok %>% filter(as.character(age_group) == nearest_age) %>% slice(1)
    } else {
      # fallback: any nearest age within donor, even if CV < floor (then we’ll floor it)
      pool_any <- es %>% filter(iso3 == iso3_i, sex == sex_i)
      if (nrow(pool_any) > 0) {
        nearest_age <- find_nearest_age(age_i, as.character(pool_any$age_group))
        pool_any %>% filter(as.character(age_group) == nearest_age) %>% slice(1)
      } else {
        NULL
      }
    }
  }
  
  # build a table of replacements
  replacements <- to_fix %>%
    rowwise() %>%
    mutate(.rep = list(find_replacement_row(iso3, sex, as.character(age_group)))) %>%
    ungroup()
  
  # apply replacements
  es_fixed <- es
  n_replaced <- 0L
  for (i in seq_len(nrow(replacements))) {
    row_i <- replacements[i,]
    rep_i <- row_i$.rep[[1]]
    key <- with(row_i, paste(iso3, sex, as.character(age_group)))
    if (!is.null(rep_i) && nrow(rep_i) == 1) {
      # Copy both best_dist and cv from replacement
      es_fixed <- es_fixed %>%
        mutate(
          best_dist = if_else(iso3 == row_i$iso3 & sex == row_i$sex &
                                as.character(age_group) == as.character(row_i$age_group),
                              rep_i$best_dist, best_dist),
          cv = if_else(iso3 == row_i$iso3 & sex == row_i$sex &
                         as.character(age_group) == as.character(row_i$age_group),
                       rep_i$cv, cv)
        )
      n_replaced <- n_replaced + 1L
    } else {
      # Last resort: just floor CV (keep best_dist)
      es_fixed <- es_fixed %>%
        mutate(
          cv = if_else(iso3 == row_i$iso3 & sex == row_i$sex &
                         as.character(age_group) == as.character(row_i$age_group) & cv < floor_cv,
                       floor_cv, cv)
        )
      n_replaced <- n_replaced + 1L
    }
  }
  
  message(sprintf("🩹 CV repair: adjusted %d rows with cv < %.3f", n_replaced, floor_cv))
  
  es_fixed %>% select(-age_group_chr, -age_group_ord)
})(energy_shapes_base, cv_floor)

# Keep the repaired shapes under the conventional name
energy_shapes <- energy_shapes_fix

nutrir_keys <- energy_shapes %>%
  select(iso3, sex, age_group) %>% distinct()

# ------------------------------------------------------------
# 2) Steps 1–3 (exact; nearest age same sex; opposite sex same age)
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# 3) Step 4 (nearest country) — distance in country kcal (like your protein mean distance)
# ------------------------------------------------------------
country_kcal <- country_level_data %>%
  select(iso3, est_kcal_consumed_per_capita)

stopifnot(!any(is.na(country_kcal$est_kcal_consumed_per_capita)))
dist_mat <- as.matrix(dist(column_to_rownames(country_kcal, "iso3")))

unmatched_3 <- anti_join(cal_keys, all_matches_1_3, by = c("iso3","sex","age_group"))
shape_donor_iso3s <- unique(all_matches_1_3$match_iso3)

if (length(shape_donor_iso3s) == 0 || nrow(unmatched_3) == 0) {
  step4_fallback <- tibble(
    iso3 = character(), sex = character(), age_group = character(),
    match_iso3 = character(), match_sex = character(), match_age_group = character(),
    source = character()
  )
} else {
  iso_match_key <- tibble(iso3 = unique(unmatched_3$iso3)) %>%
    rowwise() %>%
    mutate(fallback_iso3 = names(sort(dist_mat[iso3, shape_donor_iso3s]))[1]) %>%
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
        
        # 4a exact sex+age in fallback
        m1 <- pool %>% filter(match_iso3 == fallback,
                              match_sex == target_sex,
                              as.character(match_age_group) == as.character(target_age))
        if (nrow(m1) > 0) {
          return_value <- m1 %>% slice(1) %>% mutate(source = "nearest_country_exact")
        } else {
          # 4b nearest age same sex
          pool_same <- pool %>% filter(match_iso3 == fallback, match_sex == target_sex)
          if (nrow(pool_same) > 0) {
            nearest_age <- find_nearest_age(target_age, pool_same$match_age_group)
            m2 <- pool_same %>% filter(match_age_group == nearest_age)
            if (nrow(m2) > 0) {
              return_value <- m2 %>% slice(1) %>% mutate(source = "nearest_country_nearest_age")
            }
          }
          # 4c same age opposite sex
          if (is.null(return_value)) {
            m3 <- pool %>% filter(match_iso3 == fallback,
                                  match_sex == opposite_sex,
                                  as.character(match_age_group) == as.character(target_age))
            if (nrow(m3) > 0) {
              return_value <- m3 %>% slice(1) %>% mutate(source = "nearest_country_opposite_sex")
            }
          }
          # 4d nearest age opposite sex
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
cat("✅ [Calories] Total matched:", nrow(cal_matches_final), "\n")
cat("❌ [Calories] Still unmatched after Step 4:", nrow(final_unmatched_cal), "\n")

# ------------------------------------------------------------
# 4) Join means + (repaired) shapes/CV
# ------------------------------------------------------------
calories_distributions <- cal_matches_final %>%
  left_join(calories_means, by = c("iso3","sex","age_group")) %>%
  left_join(energy_shapes,
            by = c("match_iso3"="iso3","match_sex"="sex","match_age_group"="age_group"))

missing_shape_n <- sum(is.na(calories_distributions$best_dist))
cat("🔎 [Calories] Missing shapes after first join:", missing_shape_n, "\n")

# ------------------------------------------------------------
# 5) Repair after join (nearest age same sex within donor, then opposite sex + nearest age)
# ------------------------------------------------------------
if (missing_shape_n > 0) {
  # Map available ages per (iso3, sex) from *repaired* energy_shapes
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
  cat("🛠  Repair A (nearest age same sex) — remaining missing:", na_after_A, "\n")
  
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
    cat("🛠  Repair B (nearest age opposite sex) — remaining missing:", na_after_B, "\n")
    
    # pick best version
    if (na_after_B < na_after_A) {
      cal_matches_final      <- cal_matches_final_v3
      calories_distributions <- calories_distributions_v3
    } else {
      cal_matches_final      <- cal_matches_final_v2
      calories_distributions <- calories_distributions_v2
    }
  } else {
    cal_matches_final      <- cal_matches_final_v2
    calories_distributions <- calories_distributions_v2
  }
}

# ------------------------------------------------------------
# 6) Lean output + assertions (mirrors protein)
# ------------------------------------------------------------
calories_distributions_lean <- calories_distributions %>%
  select(
    iso3, sex, age_group,      # identifiers
    kcal_mean,                 # in-memory mean (already harmonized bins)
    best_dist,                 # distribution type
    cv,                        # coefficient of variation (now repaired at source)
    source_final               # diagnostic of where the shape came from
  )


# This is the last command in your Script 3
saveRDS(calories_distributions_lean, file = "./output/final_calorie_distributions.rds")

