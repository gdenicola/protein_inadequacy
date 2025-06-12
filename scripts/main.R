
# Load required packages
library(dplyr)
library(readr)
library(nutriR)
library(ggplot2)
library(tidyr)
library(purrr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
options(scipen=999)
rm(list = ls())

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
#####DO HERE TRIPLICATION OF NUTRIR AGE GROUPS TO ENSURE FULL MATCHING


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


#code to aggregate GDD classes 1-2 and 2-4, for martching with weight data later

# Step 1: Weighted average from 1-1.99 and 2-4 (weights: 1 year and 3 years)
weights <- tibble(
  age_group = c("1-1.99", "2-4"),
  weight = c(1, 3)
)

# Step 2: Weighted aggregation
agg_values <- gdd_distributions_lean %>%
  filter(age_group %in% c("1-1.99", "2-4")) %>%
  left_join(weights, by = "age_group") %>%
  group_by(iso3, sex) %>%
  summarise(
    age_group = "1-4",
    gdd_mean  = weighted.mean(gdd_mean,  weight, na.rm = TRUE),
    gdd_lower = weighted.mean(gdd_lower, weight, na.rm = TRUE),
    gdd_upper = weighted.mean(gdd_upper, weight, na.rm = TRUE),
    .groups   = "drop"
  )

# Step 3: Just get best_dist and cv from 2-4
shape_info <- gdd_distributions_lean %>%
  filter(age_group == "2-4") %>%
  select(iso3, sex, best_dist, cv)

# Step 4: Merge them
agg_1_4 <- agg_values %>%
  left_join(shape_info, by = c("iso3", "sex"))

# Step 5: Bind back to dataset, drop old 1-1.99 and 2-4
gdd_distributions_agg <- gdd_distributions_lean %>%
  filter(!age_group %in% c("1-1.99", "2-4")) %>%
  bind_rows(agg_1_4) %>%
  arrange(iso3, sex, age_group)

library(tibble)

#official protein RDA and EAR values
#source: https://nap.nationalacademies.org/read/10490/chapter/12#633

protein_rda <- tibble::tibble(
  age_range = c("0-0.5", "0.6-0.99", "1-3", "4-8", "9-13", "14-18", "19+"),
  age_lower = c(0.0, 0.6, 1.0, 4.0, 9.0, 14.0, 19.0),
  age_upper = c(0.5, 1.0, 3.0, 8.0, 13.0, 18.0, Inf),
  rda_g_per_kg = c(1.52, 1.20, 1.05, 0.95, 0.95, 0.85, 0.80)
)

protein_ear <- tibble::tibble(
  age_range = c("0-0.5", "0.6-0.99", "1-3", "4-8", "9-13", "14-18", "19+"),
  age_lower = c(0.0, 0.6, 1.0, 4.0, 9.0, 14.0, 19.0),
  age_upper = c(0.5, 1.0, 3.0, 8.0, 13.0, 18.0, Inf),
  ear_g_per_kg = c(1.52, 1.00, 0.87, 0.76, 0.76, 0.72, 0.66)
)

#now the same RDS values, but converted to the GDD age groups 
#through weighted averages:

protein_rda_gdd <- tibble::tibble(
  age_group = c("0-0.99", "1-4", "5-9", "10-14", "15-19",
                "20-24", "25-29", "30-34", "35-39", "40-44", "45-49",
                "50-54", "55-59", "60-64", "65-69", "70-74", "75-79",
                "80-84", "85-89", "90-94", "95-99"),
  rda_g_per_kg = c(
    1.36,   # 0-0.99
    1.025,   # 1-4
    0.95,   # 5-9
    0.93,   # 10-14
    0.84,   # 15-19
    rep(0.80, 16) # 20+ (20-24 to 95-99)
  )
)

#now the same EAR values, but converted to the GDD age groups 
#through weighted averages:

protein_ear_gdd <- tibble::tibble(
  age_group = c("0-0.99", "1-4", "5-9", "10-14", "15-19",
                "20-24", "25-29", "30-34", "35-39", "40-44", "45-49",
                "50-54", "55-59", "60-64", "65-69", "70-74", "75-79",
                "80-84", "85-89", "90-94", "95-99"),
  ear_g_per_kg = c(
    1.26,   # 0-0.99
    0.8425,   # 1-4
    0.76,   # 5-9
    0.752,   # 10-14
    0.706,   # 15-19
    rep(0.66, 16) # 20+ (20-24 to 95-99)
  )
)



####now, incorporate weight data and match it 
library(readxl)
weight_data <- read_excel("./data/data_extract_040725.xlsx", sheet = 2)
# --- Step 1: Pivot weight_data to wide format ---
weight_data_wide <- weight_data %>%
  pivot_wider(
    names_from = Stats,
    values_from = Value,
    names_prefix = "weight_",
    values_fill = list(Value = NA_real_)
  )

# --- Step 2, 3, 4, 5: Filter for MLE/FML, Map Sex, Harmonize Age, Select ---
male_identifier <- "MLE"
female_identifier <- "FML"
single_sex_identifiers <- c(male_identifier, female_identifier)

weight_data_processed <- weight_data_wide %>%
  # Filter for relevant Year (e.g., 2018)
  filter(Year == 2018) %>%
  # Filter for single-sex data ("MLE", "FML"). Rows with Sex=="BTH" will be excluded.
  filter(Sex %in% single_sex_identifiers) %>%
  # Remove aggregate/unwanted age groups FROM THE ORIGINAL 'Age' COLUMN
  filter(!Age %in% c("all-a", "20+")) %>% # Add other specific strings from Age column to remove
  
  # Rename Region to iso3 (assuming Region contains ISO3 codes)
  rename(iso3 = Region) %>%
  
  # Map original Sex column ("MLE"/"FML") to GDD "Males"/"Females"
  mutate(
    sex_gdd = case_when(
      Sex == male_identifier   ~ "Males",
      Sex == female_identifier ~ "Females",
      TRUE ~ NA_character_ # Should not happen if filter above is correct
    )
  ) %>%
  # Filter out any rows that didn't map to Males/Females (e.g. if 'Sex' was something unexpected)
  filter(!is.na(sex_gdd)) %>%
  select(-Sex) %>% # Remove original "Sex" column
  rename(sex = sex_gdd) %>% # Rename to 'sex' to match GDD data
  
  # Harmonize Age Groups from weight_data$Age to gdd_distributions_agg$age_group
  mutate(
    age_group_mapped = case_when(
      Age == "<1" ~ "0-0.99",
      Age == "95+" ~ "95-99", # Mapping for 95+
      # Add other explicit mappings if weight_data$Age strings differ from GDD age_group strings
      # e.g., Age == "01-04" ~ "1-4",
      TRUE ~ Age # Default: if Age string matches a GDD age_group string, use it
    )
  ) %>%
  # Keep only weights for age groups that are present in your target gdd_distributions_agg
  filter(age_group_mapped %in% unique(gdd_distributions_agg$age_group)) %>%
  
  # Select final columns (weights are assumed to be in kg)
  select(
    iso3,
    sex, # This is now "Males" or "Females"
    age_group = age_group_mapped,
    year = Year,
    weight_mean_kg = weight_mean,
    weight_low_kg = weight_low,
    weight_high_kg = weight_high
  ) %>%
  distinct() # Ensure no accidental duplicates


# --- Step 6: Merge with gdd_distributions_agg ---
gdd_data_with_weights <- gdd_distributions_agg %>%
  left_join(weight_data_processed, by = c("iso3", "sex", "age_group"))


#Impute the bodyweight data for South Sudan with Sudan
# 1. Impute SSD (South Sudan) using SDN (Sudan) data
# First, prepare the Sudan data that will be used for imputation
sudan_weights_for_imputation <- weight_data_processed %>%
  filter(iso3 == "SDN") %>% 
  select(sex, age_group, 
         impute_weight_mean_kg = weight_mean_kg,
         impute_weight_low_kg = weight_low_kg,
         impute_weight_high_kg = weight_high_kg) %>%
  distinct() 

if(nrow(sudan_weights_for_imputation) == 0) {
  cat("WARNING: No weight data found for SDN (Sudan) to use for imputing SSD. Skipping SSD imputation.\n")
} else {
  cat("Imputing missing weights for SSD using SDN data...\n")
  
  # Ensure 'weight_imputed_source' column exists before trying to update it
  if (!"weight_imputed_source" %in% names(gdd_data_with_weights)) {
    gdd_data_with_weights$weight_imputed_source <- NA_character_
  }
  
  gdd_data_with_weights <- gdd_data_with_weights %>%
    left_join(
      sudan_weights_for_imputation,
      by = c("sex", "age_group") 
    ) %>%
    mutate(
      # Store original weight_mean_kg to check if it was NA *before* this imputation step's join
      # This is important for the weight_imputed_source logic
      original_na_flag_mean = is.na(weight_mean_kg), 
      original_na_flag_low = is.na(weight_low_kg),
      original_na_flag_high = is.na(weight_high_kg),
      
      weight_mean_kg = ifelse(iso3 == "SSD" & original_na_flag_mean & !is.na(impute_weight_mean_kg), 
                              impute_weight_mean_kg, 
                              weight_mean_kg),
      weight_low_kg  = ifelse(iso3 == "SSD" & original_na_flag_low & !is.na(impute_weight_low_kg), 
                              impute_weight_low_kg,  
                              weight_low_kg),
      weight_high_kg = ifelse(iso3 == "SSD" & original_na_flag_high & !is.na(impute_weight_high_kg), 
                              impute_weight_high_kg, 
                              weight_high_kg),
      
      # Updated logic for weight_imputed_source
      weight_imputed_source = ifelse(
        iso3 == "SSD" & original_na_flag_mean & !is.na(impute_weight_mean_kg), # Condition for this imputation
        "Imputed from SDN",                                                    # Value if imputed now
        weight_imputed_source                                                  # Else, keep existing value (NA or from previous imputation)
      )
    ) %>%
    select(-starts_with("impute_weight_"), -starts_with("original_na_flag_")) # remove temp columns
  
  ssd_check_after_impute <- gdd_data_with_weights %>% filter(iso3 == "SSD")
  cat("SSD rows after attempting imputation from SDN:\n")
  print(head(ssd_check_after_impute %>% select(iso3, sex, age_group, weight_mean_kg, weight_imputed_source)))
  cat("Number of remaining NAs for SSD weight_mean_kg:", sum(is.na(ssd_check_after_impute$weight_mean_kg)), "\n")
}

# 2. Impute KIR & MHL for age_group "95-99" using their "90-94" data
# Create a dataset of weights from the "90-94" age group for KIR and MHL
# We take these from gdd_data_with_weights as it contains the most up-to-date weights
impute_90_94_lookup <- gdd_data_with_weights %>%
  filter(iso3 %in% c("KIR", "MHL") & age_group == "90-94") %>%
  select(iso3, sex, 
         # Use distinct names for these lookup columns to avoid any clashes
         lookup_9094_weight_mean_kg = weight_mean_kg,
         lookup_9094_weight_low_kg  = weight_low_kg,
         lookup_9094_weight_high_kg = weight_high_kg) %>%
  # Ensure we only use rows where the 90-94 data actually exists and is not NA
  filter(!is.na(lookup_9094_weight_mean_kg)) 

if(nrow(impute_90_94_lookup) == 0) {
  cat("WARNING: No non-NA weight data found for KIR/MHL in the 90-94 age group to use for imputation. Skipping this imputation.\n")
} else {
  # Ensure 'weight_imputed_source' column exists (it should from SSD step or initialization)
  if (!"weight_imputed_source" %in% names(gdd_data_with_weights)) {
    gdd_data_with_weights$weight_imputed_source <- NA_character_
  }
  
  gdd_data_with_weights <- gdd_data_with_weights %>%
    left_join(
      impute_90_94_lookup, # Joining with the prepared lookup table
      by = c("iso3", "sex") 
    ) %>%
    mutate(
      # Capture NA status *before* this specific imputation for relevant rows
      # This flag checks if the target row (KIR/MHL, 95-99) currently has an NA weight
      original_na_flag_mean_kirmhl9599 = is.na(weight_mean_kg) & (iso3 %in% c("KIR", "MHL") & age_group == "95-99"),
      original_na_flag_low_kirmhl9599  = is.na(weight_low_kg)  & (iso3 %in% c("KIR", "MHL") & age_group == "95-99"),
      original_na_flag_high_kirmhl9599 = is.na(weight_high_kg) & (iso3 %in% c("KIR", "MHL") & age_group == "95-99"),
      
      # Impute weights using the lookup values
      weight_mean_kg = ifelse(original_na_flag_mean_kirmhl9599 & !is.na(lookup_9094_weight_mean_kg),
                              lookup_9094_weight_mean_kg,
                              weight_mean_kg),
      weight_low_kg  = ifelse(original_na_flag_low_kirmhl9599 & !is.na(lookup_9094_weight_low_kg),
                              lookup_9094_weight_low_kg,
                              weight_low_kg),
      weight_high_kg = ifelse(original_na_flag_high_kirmhl9599 & !is.na(lookup_9094_weight_high_kg),
                              lookup_9094_weight_high_kg,
                              weight_high_kg),
      
      # Set imputation flag correctly
      weight_imputed_source = case_when(
        # Condition: if this row was targeted, was NA, AND successfully received a value from lookup
        original_na_flag_mean_kirmhl9599 & !is.na(lookup_9094_weight_mean_kg) ~ "Imputed from 90-94",
        # Otherwise, keep the existing flag (could be NA, or "Imputed from SDN")
        TRUE ~ weight_imputed_source 
      )
    ) %>%
    # Remove the temporary lookup and flag columns
    select(
      -starts_with("lookup_9094_"), 
      -starts_with("original_na_flag_")
    )
  
  # Verify KIR/MHL imputation
  kirmhl_check_after_impute <- gdd_data_with_weights %>% 
    filter(iso3 %in% c("KIR", "MHL") & age_group == "95-99")
  cat("KIR/MHL 95-99 rows after attempting imputation from 90-94 (with corrected flag logic):\n")
  print(kirmhl_check_after_impute %>% select(iso3, sex, age_group, weight_mean_kg, weight_imputed_source))
  cat("Number of remaining NAs for KIR/MHL 95-99 weight_mean_kg:", sum(is.na(kirmhl_check_after_impute$weight_mean_kg)), "\n")
}

# Final clean-up: Remove the 'weight_imputed_source' and 'year' columns
gdd_data_with_weights <- gdd_data_with_weights %>%
  select(
    -weight_imputed_source, # Remove the imputation flag column
    -year                   # Remove the year column (if it came from weight data)
  )


#ok, now let's match RDA data to obtain absolute RDA values (in grams):

# Assume 'gdd_data_with_weights' and 'protein_rda_gdd' are in your environment.

# Ensure 'protein_rda_gdd' has the correct column names for joining
# (it should already, but a good check)
# Expected: protein_rda_gdd has 'age_group' and 'rda_g_per_kg'

# Merge the EAR per kg data with your main dataset
gdd_data_with_abs_ear <- gdd_data_with_weights %>%
  left_join(protein_ear_gdd, by = "age_group")

# Calculate absolute ear in g/day
# We'll create ear_mean_g_day.
# If you also have weight_low_kg and weight_high_kg and want ear ranges,
# you can calculate ear_low_g_day and ear_high_g_day similarly.
gdd_data_with_abs_ear <- gdd_data_with_abs_ear %>%
  mutate(
    ear_mean_g_day = ear_g_per_kg * weight_mean_kg,
    # Optional: Calculate ear based on lower and upper body weights if available and desired
    ear_low_g_day  = ear_g_per_kg * weight_low_kg,
    ear_high_g_day = ear_g_per_kg * weight_high_kg
  )

#now,let's calculate the proportion below ear:


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

# Apply this function row-wise
gdd_inadequacy_estimates <- gdd_data_with_abs_ear %>%
  rowwise() %>% # Process row by row
  mutate(
    prevalence_inadequate = calculate_inadequacy(
      gdd_mean,         # Corresponds to mean_intake
      cv,               # Corresponds to cv_intake
      best_dist,        # Corresponds to distribution_type
      ear_mean_g_day    # Corresponds to requirement
    )
  ) %>%
  ungroup() # Important to ungroup after rowwise operations

# Display a few rows with the new prevalence_inadequate column
# print(head(gdd_inadequacy_estimates %>%
#              select(iso3, sex, age_group, gdd_mean, cv, best_dist, ear_mean_g_day, prevalence_inadequate)))

# (Optional) Check for NAs in the prevalence_inadequate column
# NAs could arise if gdd_mean, cv, or ear_mean_g_day were NA, or if cv was 0 for gamma.
# cat("NAs in prevalence_inadequate:", sum(is.na(gdd_inadequacy_estimates$prevalence_inadequate)), "\n")

# (Optional) Summary of prevalence estimates
#summary(gdd_inadequacy_estimates$prevalence_inadequate)

########all combinations now

# Assume 'gdd_inadequacy_estimates' is your dataframe that already includes:
# gdd_mean, gdd_lower, gdd_upper
# weight_mean_kg, weight_low_kg, weight_high_kg
# ear_g_per_kg
# cv, best_dist
# (It might not have ear_mean_g_day, ear_low_g_day, ear_high_g_day from the *previous* step,
# or if it does, we will recalculate them scenario-specifically)

# We'll use the 'calculate_inadequacy' function from the previous step.

# Create a data frame defining the scenarios
# We need to handle cases where _low or _high might be NA (e.g. if only mean weight was available)
# For this example, I'll assume gdd_lower, gdd_upper, weight_low_kg, weight_high_kg exist
# and are populated. If they can be NA, the logic below needs more coalesce() or filtering.

# Select base columns that don't change per scenario for a stratum
base_data_for_sensitivity <- gdd_inadequacy_estimates %>%
  select(iso3, sex, age_group, 
         # GDD intake estimates
         gdd_L = gdd_lower, gdd_M = gdd_mean, gdd_U = gdd_upper,
         # Weight estimates
         wt_L = weight_low_kg, wt_M = weight_mean_kg, wt_H = weight_high_kg,
         # ear per kg (constant for an age_group)
         ear_g_per_kg,
         # Distribution shape (constant for a stratum in this sensitivity analysis)
         cv, best_dist,
         # Original prevalence if you want to compare
         # prevalence_inadequate_MM = prevalence_inadequate # From previous step
  )

# Reshape to long format to easily create all 9 combinations
sensitivity_df_long <- base_data_for_sensitivity %>%
  pivot_longer(
    cols = starts_with("gdd_"),
    names_to = "gdd_level_char",
    names_prefix = "gdd_",
    values_to = "scenario_gdd_intake"
  ) %>%
  pivot_longer(
    cols = starts_with("wt_"),
    names_to = "wt_level_char",
    names_prefix = "wt_",
    values_to = "scenario_body_weight"
  ) %>%
  mutate(
    scenario_label = paste0(gdd_level_char, wt_level_char) # e.g., "LM", "MH"
  )

# Now calculate scenario-specific ear and prevalence
sensitivity_results <- sensitivity_df_long %>%
  mutate(
    scenario_ear_g_day = ear_g_per_kg * scenario_body_weight
  ) %>%
  rowwise() %>% # Important for calculate_inadequacy
  mutate(
    scenario_prevalence_inadequate = calculate_inadequacy(
      mean_intake = scenario_gdd_intake,
      cv_intake = cv,
      distribution_type = best_dist,
      requirement = scenario_ear_g_day
    )
  ) %>%
  ungroup()

# You can view the results:
# print(head(sensitivity_results %>%
#              select(iso3, sex, age_group, scenario_label, scenario_gdd_intake, scenario_body_weight, 
#                     scenario_ear_g_day, scenario_prevalence_inadequate)))

# To get a summary for each original stratum (iso3, sex, age_group),
# you might want to find the min, mean, and max prevalence from these 9 scenarios.
summary_sensitivity_per_stratum <- sensitivity_results %>%
  group_by(iso3, sex, age_group) %>%
  summarise(
    min_prevalence_sensitivity = min(scenario_prevalence_inadequate, na.rm = TRUE),
    mean_prevalence_sensitivity = mean(scenario_prevalence_inadequate, na.rm = TRUE), # Mean of the 9 scenarios
    max_prevalence_sensitivity = max(scenario_prevalence_inadequate, na.rm = TRUE),
    # You can also pull out the original MM estimate if you didn't carry it
    # or recalculate it here for comparison.
    # For instance, find the scenario that was GDD_Mean and Weight_Mean
    prevalence_MM = scenario_prevalence_inadequate[scenario_label == "MM"], # Assuming MM label
    .groups = "drop"
  )

# print(head(summary_sensitivity_per_stratum))

# Optional: Merge this summary back to your original gdd_inadequacy_estimates if desired
gdd_inadequacy_estimates_with_sensitivity <- gdd_inadequacy_estimates %>%
  left_join(summary_sensitivity_per_stratum, by = c("iso3", "sex", "age_group"))

# print(head(gdd_inadequacy_estimates_with_sensitivity %>%
#              select(iso3, sex, age_group, prevalence_inadequate, 
#                     min_prevalence_sensitivity, max_prevalence_sensitivity)))


####let us now include popualtion dada
# library(devtools)
# options(timeout = 600)
# install_github("PPgp/wpp2024")
library(wpp2024)
data(popAge1dt)
pop2018_dt <- popAge1dt[year == 2018] #filter 2018 only

#next steps: match iso3 with country codes, then merge with protein data, then 
#calculate total population that is protein deficient per each stratum 
#(percentage*stratum population). the will also be able to calculate percentage
#of global popilation that is protein deficient (using column totals)

# Install countrycode if you haven't already
#install.packages("countrycode")

# Load necessary libraries
#library(data.table)
library(countrycode)
# Install countrycode if you haven't already
# install.packages("countrycode")


# --- Add iso3 column using countrycode ---
pop2018_dt[, iso3 := countrycode(sourcevar = country_code,
                                 origin    = "un",
                                 destination = "iso3c",
                                 nomatch   = NA_character_)]

# --- MANUAL FIX for Taiwan (UN M49 code 158) ---
# Check if country_code 158 exists and its current iso3 is NA
if (any(pop2018_dt$country_code == 158 & is.na(pop2018_dt$iso3))) {
  pop2018_dt[country_code == 158, iso3 := "TWN"]
  cat("Manually assigned iso3 = 'TWN' to country_code == 158 (Taiwan).\n")
} else {
  cat("Country_code 158 (Taiwan) not found or already had an iso3 mapping. Manual assignment skipped.\n")
}

# --- MANUAL FIX for Kosovo (UN M49 code 412) - Optional, if needed ---
# Check your GDD data to see if you need "XKX" for Kosovo.
# The UN M49 code for Kosovo under UNSCR 1244 is often 412.
# If you need it and it's NA:
 if (any(pop2018_dt$country_code == 412 & is.na(pop2018_dt$iso3))) {
   pop2018_dt[country_code == 412, iso3 := "XKX"] # XKX is a user-assigned code for Kosovo
   cat("Manually assigned iso3 = 'XKX' to country_code == 412 (Kosovo).\n")
 }

# --- Prepare the final dataset for your analysis (country-level) ---
population_data_iso3 <- pop2018_dt[!is.na(iso3), .(
  iso3,
  original_country_name = name,
  age,
  pop_male = popM,
  pop_female = popF
)]



###now, match population data with the gdd protein distributions

cat("\n--- Processing Population Data (dplyr/tibble focus) ---\n")

# Convert population_data_iso3 (data.table) to a tibble for dplyr operations
population_data_iso3_tbl <- as_tibble(population_data_iso3)

# --- Step 1 (Population): Reshape to long format using dplyr/tidyr ---
population_long_tbl <- population_data_iso3_tbl %>%
  pivot_longer(cols = c(pop_male, pop_female),
               names_to = "sex_raw",
               values_to = "population_thousands") %>%
  mutate(sex = case_when(
    sex_raw == "pop_male" ~ "Males",
    sex_raw == "pop_female" ~ "Females",
    TRUE ~ NA_character_
  )) %>%
  select(-sex_raw) # Remove the temporary sex_raw column

# --- Step 2 (Population): Create GDD Age Bands using dplyr ---
# Define your GDD age bands and their corresponding lower/upper single ages
age_breaks_pop <- c(0, 1, 2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 101) # Upper bound is exclusive in cut
age_labels_pop <- c("0-0.99", "1-1.99", "2-4", "5-9", "10-14", "15-19",
                    "20-24", "25-29", "30-34", "35-39", "40-44", "45-49",
                    "50-54", "55-59", "60-64", "65-69", "70-74", "75-79",
                    "80-84", "85-89", "90-94", "95-99")

population_long_tbl <- population_long_tbl %>%
  mutate(
    age_group_gdd_temp = cut(age,
                             breaks = age_breaks_pop,
                             labels = age_labels_pop,
                             right = FALSE,
                             include.lowest = TRUE),
    # Convert factor to character for re-labeling
    age_group_gdd_temp = as.character(age_group_gdd_temp),
    # Combine "1-1.99" and "2-4" into a "1-4" GDD target group
    age_group = ifelse(age_group_gdd_temp %in% c("1-1.99", "2-4"), "1-4", age_group_gdd_temp)
  ) %>%
  select(-age_group_gdd_temp) # Remove temporary column

# --- Step 3 (Population): Aggregate Population by iso3, sex, and GDD age_group ---
population_aggregated_tbl <- population_long_tbl %>%
  group_by(iso3, sex, age_group) %>%
  summarise(population_total_thousands = sum(population_thousands, na.rm = TRUE),
            .groups = "drop") # Summing populations

# --- Step 4 (Population): Merge with your GDD inadequacy estimates ---
# 'gdd_inadequacy_estimates_with_sensitivity' is your main GDD results tibble
gdd_data_with_population <- gdd_inadequacy_estimates_with_sensitivity %>%
  left_join(population_aggregated_tbl, by = c("iso3", "sex", "age_group"))

# Convert population to absolute numbers and remove the thousands column using dplyr
gdd_data_with_population <- gdd_data_with_population %>%
  mutate(
    population = population_total_thousands * 1000
  ) %>%
  select(-population_total_thousands)

# Your gdd_data_with_population is now ready!


####now, let's calculate percentages

# Load dplyr if not already loaded
library(dplyr)

# Assume 'gdd_data_with_population' is your dataframe from the previous step.
# It should contain:
# - population (absolute number of people in the stratum)
# - prevalence_inadequate (your central estimate using mean GDD, mean weight, EAR cut-point)
#   OR prevalence_MM (if you named the central estimate from sensitivity analysis this way)
# - min_prevalence_sensitivity
# - max_prevalence_sensitivity

# Make sure we have a consistent name for the central prevalence estimate.
# If 'prevalence_MM' exists and is what you want as central, and 'prevalence_inadequate' might be different,
# let's ensure we use 'prevalence_MM'.
# Otherwise, if 'prevalence_inadequate' is your main central estimate, that's fine.
# For this code, I'll assume 'prevalence_inadequate' is the primary central estimate
# and 'prevalence_MM' might be identical or also present. If you prefer to use
# 'prevalence_MM' explicitly as the central one, adjust the column name below.

# If 'prevalence_MM' is the definitive central estimate from the sensitivity analysis,
# and you want to ensure it's used:
if ("prevalence_MM" %in% names(gdd_data_with_population) && 
    !"prevalence_inadequate_central" %in% names(gdd_data_with_population)) {
  gdd_data_with_population <- gdd_data_with_population %>%
    rename(prevalence_inadequate_central = prevalence_MM)
} else if ("prevalence_inadequate" %in% names(gdd_data_with_population) &&
           !"prevalence_inadequate_central" %in% names(gdd_data_with_population)) {
  gdd_data_with_population <- gdd_data_with_population %>%
    rename(prevalence_inadequate_central = prevalence_inadequate)
} else if (!"prevalence_inadequate_central" %in% names(gdd_data_with_population)){
  stop("A central prevalence estimate column (e.g. 'prevalence_MM' or 'prevalence_inadequate') needs to be identified or renamed to 'prevalence_inadequate_central'.")
}


# --- Step 1: Calculate number of people with inadequate intake for each stratum ---
gdd_counts_inadequate <- gdd_data_with_population %>%
  mutate(
    # Using the central estimate (Mean GDD, Mean Weight, EAR cut-point)
    inadequate_count_central = prevalence_inadequate_central * population,
    
    # Using the minimum prevalence from sensitivity analysis
    inadequate_count_min_sensitivity = min_prevalence_sensitivity * population,
    
    # Using the maximum prevalence from sensitivity analysis
    inadequate_count_max_sensitivity = max_prevalence_sensitivity * population
  )

# Display a few rows with the new count columns
# print(head(gdd_counts_inadequate %>%
#              select(iso3, sex, age_group, population, prevalence_inadequate_central, inadequate_count_central,
#                     min_prevalence_sensitivity, inadequate_count_min_sensitivity,
#                     max_prevalence_sensitivity, inadequate_count_max_sensitivity)))


# --- Step 2: Calculate Global Totals ---
# Summing up population and inadequate counts across all strata

global_summary <- gdd_counts_inadequate %>%
  summarise(
    total_population_global = sum(population, na.rm = TRUE),
    
    total_inadequate_central_global = sum(inadequate_count_central, na.rm = TRUE),
    total_inadequate_min_sensitivity_global = sum(inadequate_count_min_sensitivity, na.rm = TRUE),
    total_inadequate_max_sensitivity_global = sum(inadequate_count_max_sensitivity, na.rm = TRUE)
  ) %>%
  mutate(
    # Calculate global percentage of inadequacy for each scenario
    percent_inadequate_central_global = (total_inadequate_central_global / total_population_global),
    percent_inadequate_min_sensitivity_global = (total_inadequate_min_sensitivity_global / total_population_global),
    percent_inadequate_max_sensitivity_global = (total_inadequate_max_sensitivity_global / total_population_global)
  )

cat("\n--- Global Summary of Protein Inadequacy (EAR-based, 2018) ---\n")
print(global_summary)

# Nicer print format for percentages
cat(sprintf("\nGlobal Total Population Analyzed: %.0f\n", global_summary$total_population_global))
cat(sprintf("Central Estimate - Total Inadequate: %.0f (%.2f%% of global total)\n", 
            global_summary$total_inadequate_central_global, 
            global_summary$percent_inadequate_central_global * 100))
cat(sprintf("Sensitivity Min - Total Inadequate: %.0f (%.2f%% of global total)\n", 
            global_summary$total_inadequate_min_sensitivity_global, 
            global_summary$percent_inadequate_min_sensitivity_global * 100))
cat(sprintf("Sensitivity Max - Total Inadequate: %.0f (%.2f%% of global total)\n", 
            global_summary$total_inadequate_max_sensitivity_global, 
            global_summary$percent_inadequate_max_sensitivity_global * 100))

###works til here
##note: population world in 2018 was 7,729,902,781 according to worldometer
temp_country_name_map_source <- as_tibble(pop2018_dt[!is.na(iso3)]) # Use the data.table that had iso3 added
country_name_map <- temp_country_name_map_source %>%
  distinct(iso3, original_country_name = name) # 'name' is the original UN country name


# --- Now proceed with summary_by_country ---
summary_by_country <- gdd_counts_inadequate %>%
  group_by(iso3) %>%
  summarise(
    country_total_population = sum(population, na.rm = TRUE),
    country_inadequate_central = sum(inadequate_count_central, na.rm = TRUE),
    country_inadequate_min = sum(inadequate_count_min_sensitivity, na.rm = TRUE),
    country_inadequate_max = sum(inadequate_count_max_sensitivity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    country_percent_inadequate_central = ifelse(country_total_population > 0, country_inadequate_central / country_total_population, 0),
    country_percent_inadequate_min = ifelse(country_total_population > 0, country_inadequate_min / country_total_population, 0),
    country_percent_inadequate_max = ifelse(country_total_population > 0, country_inadequate_max / country_total_population, 0)
  ) %>%
  # Join with the prepared country name map
  left_join(country_name_map, by = "iso3") %>%
  # Reorder columns to have name next to iso3 if desired
  select(iso3, original_country_name, everything()) %>% 
  arrange(desc(country_percent_inadequate_central)) # Sort by highest central prevalence


cat("\n--- Top 10 Countries by Central Estimate of Percent Inadequate ---\n")
if ("original_country_name" %in% names(summary_by_country)) {
  print(head(summary_by_country %>% 
               select(iso3, original_country_name, country_total_population, 
                      country_inadequate_central, country_percent_inadequate_central,
                      country_percent_inadequate_min, country_percent_inadequate_max), 10))
} else {
  print(head(summary_by_country %>% 
               select(iso3, country_total_population, 
                      country_inadequate_central, country_percent_inadequate_central,
                      country_percent_inadequate_min, country_percent_inadequate_max), 10))
  cat("Note: 'original_country_name' was not found in summary_by_country for printing.\n")
}


# The rest of your code for summary_by_age_global should be fine
# (Make sure age_group_levels_ordered is defined if you use age_group_f)
age_group_levels_ordered <- c("0-0.99", "1-4", "5-9", "10-14", "15-19",
                              "20-24", "25-29", "30-34", "35-39", "40-44", 
                              "45-49", "50-54", "55-59", "60-64", "65-69", 
                              "70-74", "75-79", "80-84", "85-89", "90-94", "95-99")

summary_by_age_global <- gdd_counts_inadequate %>%
  mutate(age_group_f = factor(age_group, levels = age_group_levels_ordered)) %>% # Ensure factor is created
  group_by(age_group_f) %>% 
  summarise(
    age_total_population = sum(population, na.rm = TRUE),
    age_inadequate_central = sum(inadequate_count_central, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    age_percent_inadequate_central = ifelse(age_total_population > 0, age_inadequate_central / age_total_population, 0)
  ) %>%
  arrange(age_group_f)

cat("\n--- Global Prevalence of Inadequacy by Age Group (Central Estimate) ---\n")
print(n = 21, summary_by_age_global %>% select(age_group_f, age_total_population, age_inadequate_central, age_percent_inadequate_central))


# --- (Your full script runs up to the point where 'gdd_counts_inadequate' is created) ---
# gdd_counts_inadequate contains:
# - iso3, sex, age_group
# - population
# - prevalence_inadequate_central (or your chosen central estimate column)
# - min_prevalence_sensitivity, max_prevalence_sensitivity
# - inadequate_count_central, inadequate_count_min_sensitivity, inadequate_count_max_sensitivity

# --- Calculate Global Totals by Sex ---
global_summary_by_sex <- gdd_counts_inadequate %>%
  group_by(sex) %>% # Group by sex
  summarise(
    total_population_by_sex = sum(population, na.rm = TRUE),
    
    total_inadequate_central_by_sex = sum(inadequate_count_central, na.rm = TRUE),
    total_inadequate_min_sensitivity_by_sex = sum(inadequate_count_min_sensitivity, na.rm = TRUE),
    total_inadequate_max_sensitivity_by_sex = sum(inadequate_count_max_sensitivity, na.rm = TRUE),
    .groups = "drop" # Ungroup after summarising
  ) %>%
  mutate(
    # Calculate global percentage of inadequacy for each scenario, BY SEX
    percent_inadequate_central_by_sex = ifelse(total_population_by_sex > 0, total_inadequate_central_by_sex / total_population_by_sex, 0),
    percent_inadequate_min_sensitivity_by_sex = ifelse(total_population_by_sex > 0, total_inadequate_min_sensitivity_by_sex / total_population_by_sex, 0),
    percent_inadequate_max_sensitivity_by_sex = ifelse(total_population_by_sex > 0, total_inadequate_max_sensitivity_by_sex / total_population_by_sex, 0)
  )

cat("\n\n--- Global Summary of Protein Inadequacy by Sex (EAR-based, 2018) ---\n")
print(global_summary_by_sex)

# Nicer print format for the by-sex summary
cat("\n--- Formatted Global Summary by Sex ---\n")
for (s_idx in 1:nrow(global_summary_by_sex)) {
  sex_label <- global_summary_by_sex$sex[s_idx]
  cat(sprintf("\nSex: %s\n", sex_label))
  cat(sprintf("  Total Population Analyzed (%s): %.0f\n", sex_label, global_summary_by_sex$total_population_by_sex[s_idx]))
  cat(sprintf("  Central Estimate - Total Inadequate (%s): %.0f (%.2f%% of %s total)\n", 
              sex_label,
              global_summary_by_sex$total_inadequate_central_by_sex[s_idx], 
              global_summary_by_sex$percent_inadequate_central_by_sex[s_idx] * 100,
              sex_label))
  cat(sprintf("  Sensitivity Min - Total Inadequate (%s): %.0f (%.2f%% of %s total)\n", 
              sex_label,
              global_summary_by_sex$total_inadequate_min_sensitivity_by_sex[s_idx], 
              global_summary_by_sex$percent_inadequate_min_sensitivity_by_sex[s_idx] * 100,
              sex_label))
  cat(sprintf("  Sensitivity Max - Total Inadequate (%s): %.0f (%.2f%% of %s total)\n", 
              sex_label,
              global_summary_by_sex$total_inadequate_max_sensitivity_by_sex[s_idx], 
              global_summary_by_sex$percent_inadequate_max_sensitivity_by_sex[s_idx] * 100,
              sex_label))
}



###############################################
##########VISUALIZATIONS ONLY FROM HERE#########
###############################################


# Load ggplot2 if not already loaded
library(ggplot2)
library(dplyr) # For data manipulation for plotting
library(scales) # For nice percentage labels

# Data needed: global_summary tibble
# --- Revised Plot 1: Global Prevalence of Inadequacy (Point with Error Bar) ---
plot_global_prevalence_revised <- ggplot(global_summary, aes(x = "Global Estimate")) +
  geom_point(aes(y = percent_inadequate_central_global), 
             color = "dodgerblue3", size = 4, shape = 18) + # Diamond shape for point
  geom_errorbar(aes(ymin = percent_inadequate_min_sensitivity_global, 
                    ymax = percent_inadequate_max_sensitivity_global),
                width = 0.1, linewidth = 0.8, color = "gray50") +
  geom_text(aes(y = percent_inadequate_central_global, 
                label = sprintf("%.1f%%", percent_inadequate_central_global * 100)),
            vjust = -1.5, size = 4.5, color = "black", fontface = "bold") +
  geom_text(aes(y = percent_inadequate_min_sensitivity_global, 
                label = sprintf("%.1f%%\n(Min)", percent_inadequate_min_sensitivity_global * 100)),
            vjust = 1.5, size = 3.5, color = "gray30", lineheight = 0.8) +
  geom_text(aes(y = percent_inadequate_max_sensitivity_global, 
                label = sprintf("%.1f%%\n(Max)", percent_inadequate_max_sensitivity_global * 100)),
            vjust = -0.8, size = 3.5, color = "gray30", lineheight = 0.8) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                     limits = c(0, max(global_summary$percent_inadequate_max_sensitivity_global, na.rm = TRUE) * 1.1), # Dynamic upper limit
                     expand = expansion(mult = c(0.05, 0.1))) + # Add a bit more space at top
  labs(
    title = "Global Prevalence of Inadequate Protein Intake (2018, EAR-based)",
    subtitle = "Central estimate with range from sensitivity analysis (GDD & Weight low/high)",
    x = "",
    y = "Prevalence of Inadequacy"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(size=12, face="bold"), 
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(), # Remove vertical grid line
        panel.grid.minor.x = element_blank())

print(plot_global_prevalence_revised)


# --- Revised Plot 2: Global Prevalence by Sex (Points with Error Bars) ---
plot_global_prevalence_sex_revised <- ggplot(global_summary_by_sex, 
                                             aes(x = sex, y = percent_inadequate_central_by_sex, color = sex)) +
  geom_point(size = 4, shape = 18, show.legend = FALSE) + # Diamond shape for points
  geom_errorbar(aes(ymin = percent_inadequate_min_sensitivity_by_sex, 
                    ymax = percent_inadequate_max_sensitivity_by_sex),
                width = 0.15, linewidth = 0.8, show.legend = FALSE) + # Error bars will take color from aes(color=sex)
  geom_text(aes(label = sprintf("%.1f%%", percent_inadequate_central_by_sex * 100)),
            vjust = -1.8, size = 4, color = "black", fontface="bold", show.legend = FALSE) +
  # Optional: Add text for min/max if not too cluttered
  # geom_text(aes(y = percent_inadequate_min_sensitivity_by_sex, 
  #               label = sprintf("%.1f%%", percent_inadequate_min_sensitivity_by_sex * 100)),
  #           vjust = 1.5, size = 3, show.legend = FALSE) +
  # geom_text(aes(y = percent_inadequate_max_sensitivity_by_sex, 
  #               label = sprintf("%.1f%%", percent_inadequate_max_sensitivity_by_sex * 100)),
  #           vjust = -0.5, size = 3, show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                     limits = c(0, max(global_summary_by_sex$percent_inadequate_max_sensitivity_by_sex, na.rm = TRUE) * 1.1),
                     expand = expansion(mult = c(0.05, 0.1))) +
  scale_color_brewer(palette = "Set1") + # Colors for points and error bars
  labs(
    title = "Global Prevalence of Inadequate Protein Intake by Sex (2018, EAR-based)",
    subtitle = "Central estimate with range from sensitivity analysis",
    x = "Sex",
    y = "Prevalence of Inadequacy"
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

print(plot_global_prevalence_sex_revised)


# --- Plot 3: Global Prevalence by Age Group (Line Plot with Uncertainty Band) ---
# Data needed: gdd_counts_inadequate (to calculate min/max by age) and summary_by_age_global (for central line)

# First, calculate min and max prevalence by age group globally
summary_by_age_sensitivity_global <- gdd_counts_inadequate %>%
  # Ensure age_group_f is created for proper ordering
  mutate(age_group_f = factor(age_group, levels = age_group_levels_ordered)) %>%
  group_by(age_group_f) %>%
  summarise(
    total_population = sum(population, na.rm = TRUE),
    total_inadequate_min = sum(inadequate_count_min_sensitivity, na.rm = TRUE),
    total_inadequate_max = sum(inadequate_count_max_sensitivity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    percent_inadequate_min = ifelse(total_population > 0, total_inadequate_min / total_population, 0),
    percent_inadequate_max = ifelse(total_population > 0, total_inadequate_max / total_population, 0)
  ) %>%
  # Join with the central estimate from summary_by_age_global
  left_join(summary_by_age_global %>% select(age_group_f, age_percent_inadequate_central), by = "age_group_f") %>%
  filter(!is.na(age_group_f)) # Remove any NA age_group_f that might have appeared


plot_age_prevalence_sensitivity <- ggplot(summary_by_age_sensitivity_global, aes(x = age_group_f)) +
  geom_ribbon(aes(ymin = percent_inadequate_min, ymax = percent_inadequate_max, group = 1), 
              fill = "skyblue", alpha = 0.3, linetype="dashed", color="grey70") +
  geom_line(aes(y = age_percent_inadequate_central, group = 1), color = "dodgerblue3", linewidth = 1) +
  geom_point(aes(y = age_percent_inadequate_central), color = "dodgerblue4", size = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Global Prevalence of Inadequate Protein Intake by Age Group (2018, EAR-based)",
    subtitle = "Central estimate (line) with sensitivity range (shaded band)",
    x = "Age Group",
    y = "Prevalence of Inadequacy"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_age_prevalence_sensitivity)

# --- Plot 4: Number of People Affected (Stacked Bar Chart - Central Estimate) ---
# We need broader age categories for a readable stacked bar.
# Let's define them:
broad_age_map <- tibble(
  age_group = age_group_levels_ordered, # Your fine GDD age groups
  broad_age_category = case_when(
    age_group %in% c("0-0.99", "1-4") ~ "0-4 Years",
    age_group %in% c("5-9", "10-14") ~ "5-14 Years",
    age_group %in% c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49") ~ "15-49 Years",
    age_group %in% c("50-54", "55-59", "60-64", "65-69") ~ "50-69 Years",
    TRUE ~ "70+ Years" # All remaining older groups
  )
) %>% 
  mutate(broad_age_category = factor(broad_age_category, 
                                     levels = c("0-4 Years", "5-14 Years", "15-49 Years", "50-69 Years", "70+ Years")))


# Join this map with gdd_counts_inadequate and sum up by broad_age_category and sex
counts_by_broad_age_sex <- gdd_counts_inadequate %>%
  left_join(broad_age_map, by = "age_group") %>%
  filter(!is.na(broad_age_category)) %>% # Ensure all age groups mapped
  group_by(broad_age_category, sex) %>%
  summarise(
    total_inadequate_central = sum(inadequate_count_central, na.rm = TRUE),
    .groups = "drop"
  )

plot_stacked_bar_counts <- ggplot(counts_by_broad_age_sex, 
                                  aes(x = sex, y = total_inadequate_central, fill = broad_age_category)) +
  geom_col(position = "stack", alpha = 0.8) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix = " M")) + # Show in Millions
  scale_fill_viridis_d(option = "plasma", name = "Age Category", direction = -1) +
  labs(
    title = "Global Number of People with Inadequate Protein Intake (Central Estimate)",
    subtitle = "By Sex and Broad Age Category (2018, EAR-based)",
    x = "Sex",
    y = "Number of People with Inadequate Intake"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

print(plot_stacked_bar_counts)
