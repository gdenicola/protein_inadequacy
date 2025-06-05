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

# Step 2: Load GDD 2018 protein data (v23 = total protein intake)
protein_gdd <- read_csv("./data/GDD_FinalEstimates_01102022/Country-level estimates/v23_cnty.csv")

# Step 3: Filter to national-level rows and add sex label
protein_gdd_natl <- protein_gdd %>%
  filter(year == 2018, edu == 999, urban == 999, age != 999) %>%
  mutate(sex = ifelse(female == 1, "Females", "Males"))

# Step 4: Define age groups
age_breaks <- c(-Inf, 4.99, 9.99, 14.99, 19.99, 24.99, 29.99, 34.99,
                39.99, 44.99, 49.99, 54.99, 59.99, 64.99, 69.99, 74.99,
                79.99, 84.99, 89.99, 94.99, 99.99)

age_labels <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
                "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
                "60-64", "65-69", "70-74", "75-79", "80-84", "85-89",
                "90-94", "95-99")

protein_gdd_natl <- protein_gdd_natl %>%
  mutate(age_group = cut(age, breaks = age_breaks, labels = age_labels, right = TRUE))

# Step 5: Aggregate to one row per iso3–sex–age_group
gdd_0_4 <- protein_gdd_natl %>%
  filter(age_group == "0-4", age %in% c(0.5, 1.5, 3.5)) %>%
  mutate(weight = case_when(age == 0.5 ~ 1, age == 1.5 ~ 1, age == 3.5 ~ 3)) %>%
  group_by(iso3, sex, age_group) %>%
  summarise(
    gdd_mean = weighted.mean(median, weight, na.rm = TRUE),
    gdd_lower = weighted.mean(lowerci_95, weight, na.rm = TRUE),
    gdd_upper = weighted.mean(upperci_95, weight, na.rm = TRUE),
    .groups = "drop"
  )

gdd_rest <- protein_gdd_natl %>%
  filter(age_group != "0-4") %>%
  group_by(iso3, sex, age_group) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    iso3, sex, age_group,
    gdd_mean = median,
    gdd_lower = lowerci_95,
    gdd_upper = upperci_95
  )

protein_gdd_agg <- bind_rows(gdd_0_4, gdd_rest) %>%
  arrange(iso3, sex, age_group)

# -----------------------------
# MATCHING STEPS 1–3
# -----------------------------
gdd_keys <- protein_gdd_agg %>% select(iso3, sex, age_group) %>% distinct()
nutrir_keys <- protein_nutriR %>% select(iso3, sex, age_group) %>% distinct()

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
    nutrir_keys %>% filter(iso3 == iso, sex == opposite, age_group == ag)
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
        filter(match_iso3 == fallback, match_sex == target_sex, match_age_group == target_age)
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


####NOW WE HAVE FULL MATCHING DICTIONARY#########

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

library(tibble)

#official protein RDA values
#source: https://nap.nationalacademies.org/read/10490/chapter/12#633

protein_rda <- tibble::tibble(
  age_range = c("0-0.5", "0.6-1", "1-3", "4-8", "9-13", "14-18", "19+"),
  age_lower = c(0.0, 0.6, 1.0, 4.0, 9.0, 14.0, 19.0),
  age_upper = c(0.5, 1.0, 3.0, 8.0, 13.0, 18.0, Inf),
  rda_g_per_kg = c(1.52, 1.20, 1.05, 0.95, 0.95, 0.85, 0.80)
)


# Define GDD-style age groups
gdd_age_bins <- tribble(
  ~age_group, ~ages,
  "0-4",      c(0.5, 2.5, 5.0),
  "5-9",      c(5.0, 10.0),
  "10-14",    c(10.0, 15.0),
  "15-19",    c(15.0, 22.5),
  "20-24",    c(22.5),
  "25-29",    c(22.5),
  "30-34",    c(45.0),
  "35-39",    c(45.0),
  "40-44",    c(45.0),
  "45-49",    c(45.0),
  "50-54",    c(45.0),
  "55-59",    c(45.0),
  "60-64",    c(65.0),
  "65-69",    c(65.0),
  "70-74",    c(65.0),
  "75-79",    c(65.0),
  "80-84",    c(65.0),
  "85-89",    c(65.0),
  "90-94",    c(65.0),
  "95-99",    c(65.0)
)

# Compute weighted means by GDD age group
protein_requirements_gdd <- gdd_age_bins %>%
  rowwise() %>%
  mutate(
    male_mean = mean(who_protein_raw$male_mean[who_protein_raw$age %in% ages]),
    male_sd   = mean(who_protein_raw$male_sd[who_protein_raw$age %in% ages]),
    female_mean = mean(who_protein_raw$female_mean[who_protein_raw$age %in% ages]),
    female_sd   = mean(who_protein_raw$female_sd[who_protein_raw$age %in% ages])
  ) %>%
  ungroup()
