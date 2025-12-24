# ==============================================================
# Script 7 — Global / sex / World Bank region summaries
#   Uses dat_quality.rds (output of Script 5)
#   Produces a tidy table of:
#     - Protein inadequacy (EAR), optimal calories
#     - Protein inadequacy (EAR), realistic MM calories
#     - Sub-optimal protein (OPT), optimal calories
#     - Sub-optimal protein (OPT), realistic MM calories
#     - Quality-adjusted protein inadequacy (EAR), optimal calories
#     - Quality-adjusted protein inadequacy (EAR), realistic MM calories
#   Each: global, by sex, by World Bank region, and by region × sex
# ==============================================================

# ---- Libraries ----
library(dplyr)
library(readr)
library(tidyr)
library(countrycode)
library(writexl)


# ---- Setup ----
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); setwd("..")
options(scipen = 999)

# ---- Shared helper for inadequacy (same as in Scripts 4/6) ----
calculate_inadequacy <- function(mean_intake, cv_intake, distribution_type, requirement) {
  if (is.na(requirement) || requirement <= 0 ||
      is.na(mean_intake) || is.na(cv_intake) ||
      mean_intake <= 0 || cv_intake <= 0) {
    return(NA_real_)
  }
  if (distribution_type == "gamma") {
    shape_k <- 1 / (cv_intake^2)
    if (!is.finite(shape_k)) return(NA_real_)
    scale_theta <- mean_intake * (cv_intake^2)
    return(pgamma(requirement, shape = shape_k, scale = scale_theta))
  } else if (distribution_type == "log-normal") {
    meanlog <- log(mean_intake) - 0.5 * log(1 + cv_intake^2)
    sdlog   <- sqrt(log(1 + cv_intake^2))
    if (is.na(sdlog) || sdlog <= 0) return(NA_real_)
    return(plnorm(requirement, meanlog = meanlog, sdlog = sdlog))
  } else {
    return(NA_real_)
  }
}

# ==============================================================
# 1. Load dat_quality and attach World Bank regions
# ==============================================================

dat <- readRDS("./output/dat_quality.rds") %>%
  distinct(iso3, sex, age_group, .keep_all = TRUE) %>%
  filter(age_group != "0-0.99")

dat <- readRDS("./output/dat_quality.rds")
dat %>%
  summarise(
    low_min = min(kcal_mean_low, na.rm=TRUE),
    low_max = max(kcal_mean_low, na.rm=TRUE),
    high_min = min(kcal_mean_high, na.rm=TRUE),
    high_max = max(kcal_mean_high, na.rm=TRUE)
  )

cal <- readRDS("./output/final_calorie_distributions_wide.rds")
cal %>%
  summarise(
    low_min = min(kcal_mean_low, na.rm=TRUE),
    low_max = max(kcal_mean_low, na.rm=TRUE),
    high_min = min(kcal_mean_high, na.rm=TRUE),
    high_max = max(kcal_mean_high, na.rm=TRUE)
  )

# World Bank regions (countrycode "region" is WB region)
region_lookup <- dat %>%
  distinct(iso3) %>%
  mutate(region_WB = countrycode(iso3, "iso3c", "region"))

dat <- dat %>%
  left_join(region_lookup, by = "iso3")

cat("\nSTART OF SCRIPT 7, just loaded dat_quality.rds:\n")
cat("Loaded from: ", normalizePath("./output/dat_quality.rds"), "\n")
summary(dat$kcal_mean_low)
summary(dat$kcal_mean_high)
summary(dat$kcal_mean_medium)

# ==============================================================
# 1b. Build long table for 3 calorie scenarios (low/medium/high)
# ==============================================================

scen3_long <- dat %>%
  transmute(
    iso3, sex, age_group, population, region_WB,
    # protein side
    best_dist_protein, cv_protein,
    ear_mean_g_day, ear_mean_g_day_adj, opt_mean_g_day,
    protein_kcal_share_mean,
    # calorie side
    best_dist_calorie, cv_calorie, mder_stratum_kcal,
    kcal_mean_low, kcal_mean_medium, kcal_mean_high
  ) %>%
  pivot_longer(
    cols = c(kcal_mean_low, kcal_mean_medium, kcal_mean_high),
    names_to = "calorie_scenario",
    values_to = "kcal_mean"
  ) %>%
  mutate(
    calorie_scenario = recode(
      calorie_scenario,
      kcal_mean_low = "low",
      kcal_mean_medium = "medium",
      kcal_mean_high = "high"
    )
  ) %>%
  rowwise() %>%
  mutate(
    # kcal -> protein grams using mean share
    protein_grams = (kcal_mean * protein_kcal_share_mean) / 4,
    
    # ---- Protein inadequacy prevalences (distribution-based) ----
    prot_inad_EAR = calculate_inadequacy(
      mean_intake = protein_grams, cv_intake = cv_protein,
      distribution_type = best_dist_protein, requirement = ear_mean_g_day
    ),
    prot_inad_OPT = calculate_inadequacy(
      mean_intake = protein_grams, cv_intake = cv_protein,
      distribution_type = best_dist_protein, requirement = opt_mean_g_day
    ),
    prot_inad_EAR_QA = calculate_inadequacy(
      mean_intake = protein_grams, cv_intake = cv_protein,
      distribution_type = best_dist_protein, requirement = ear_mean_g_day_adj
    ),
    
    # ---- Calorie inadequacy prevalence ----
    cal_inad = calculate_inadequacy(
      mean_intake = kcal_mean, cv_intake = cv_calorie,
      distribution_type = best_dist_calorie, requirement = mder_stratum_kcal
    )
  ) %>%
  ungroup()


# Optional sanity:
# table(dat$region_WB, useNA = "ifany")


# ==== DEBUG OUTPUT: Check if scenario cal_inad prevalence varies ====
cat("\nSCENARIO TEST: Are scenario means and prevalence changing?\n")
library(dplyr)
cat("\nWeighted means and range of kcal_mean/cal_inad by calorie_scenario:\n")
scen3_long %>%
  group_by(calorie_scenario) %>%
  summarise(
    min_kcal  = min(kcal_mean),
    max_kcal  = max(kcal_mean),
    mean_kcal = mean(kcal_mean),
    mean_prev = weighted.mean(cal_inad, w=population, na.rm=TRUE)
  ) %>%
  print()
cat("\n------------------------------------------------\n")

# ==============================================================
# Global caloric inadequacy for each calorie scenario (low/medium/high)
# ==============================================================

global_cal_inad_3scen <- scen3_long %>%
  group_by(calorie_scenario) %>%
  summarise(
    global_calorie_inadequacy = weighted.mean(cal_inad, w = population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(global_calorie_inadequacy_pct = 100 * global_calorie_inadequacy)

print(global_cal_inad_3scen)

#export
write_xlsx(
  list(
    global_calorie_inadequacy_3scenarios = global_cal_inad_3scen
  ),
  "./output/global_calorie_inadequacy_3scenarios.xlsx"
)



# ==============================================================
# 2. Optimal-calorie world (exact calories)
#    - Uses kcal_consumed_optimal if available,
#      otherwise EER (eer_kcal_marco_mean)
# ==============================================================

prot2_central <- dat %>%
  transmute(
    iso3,
    sex,
    age_group,
    population,
    region_WB,
    best_dist_protein,
    cv_protein,
    ear_mean_g_day,
    ear_mean_g_day_adj,     # quality-adjusted EAR
    opt_mean_g_day,
    protein_kcal_share_mean,
    exact_kcal = ifelse(!is.na(kcal_consumed_optimal),
                        kcal_consumed_optimal,
                        eer_kcal_marco_mean)
  ) %>%
  rowwise() %>%
  mutate(
    protein_grams_exact = (exact_kcal * protein_kcal_share_mean) / 4,
    # Standard EAR inadequacy
    # ---- Mean-only adequacy (no distribution) ----
    prot_mean_pct_EAR_exact    = 100 * (protein_grams_exact / ear_mean_g_day),
    prot_mean_pct_OPT_exact    = 100 * (protein_grams_exact / opt_mean_g_day),
    prot_mean_pct_EAR_exact_QA = 100 * (protein_grams_exact / ear_mean_g_day_adj),
    prot_inad_EAR_exact = calculate_inadequacy(
      mean_intake       = protein_grams_exact,
      cv_intake         = cv_protein,
      distribution_type = best_dist_protein,
      requirement       = ear_mean_g_day
    ),
    # OPT ("sub-optimal") inadequacy
    prot_inad_OPT_exact = calculate_inadequacy(
      mean_intake       = protein_grams_exact,
      cv_intake         = cv_protein,
      distribution_type = best_dist_protein,
      requirement       = opt_mean_g_day
    ),
    # Quality-adjusted EAR inadequacy
    prot_inad_EAR_exact_QA = calculate_inadequacy(
      mean_intake       = protein_grams_exact,
      cv_intake         = cv_protein,
      distribution_type = best_dist_protein,
      requirement       = ear_mean_g_day_adj
    )
  ) %>%
  ungroup()

# ==============================================================
# 3. Realistic-calorie world (MM scenario)
#    - Mean protein share, medium kcal scenario
# ==============================================================

scen_central <- dat %>%
  transmute(
    iso3,
    sex,
    age_group,
    population,
    region_WB,
    best_dist_protein,
    cv_protein,
    ear_mean_g_day,
    ear_mean_g_day_adj,   # quality-adjusted EAR
    opt_mean_g_day,
    best_dist_calorie,
    cv_calorie,
    mder_stratum_kcal,
    protein_kcal_share_mean,
    kcal_mean_medium
  ) %>%
  rowwise() %>%
  mutate(
    protein_grams_MM = (kcal_mean_medium * protein_kcal_share_mean) / 4,
    # ---- Mean-only adequacy (no distribution) ----
    prot_mean_pct_EAR_MM    = 100 * (protein_grams_MM / ear_mean_g_day),
    prot_mean_pct_OPT_MM    = 100 * (protein_grams_MM / opt_mean_g_day),
    prot_mean_pct_EAR_MM_QA = 100 * (protein_grams_MM / ear_mean_g_day_adj),
    # Standard EAR inadequacy
    prot_inad_EAR_MM = calculate_inadequacy(
      mean_intake       = protein_grams_MM,
      cv_intake         = cv_protein,
      distribution_type = best_dist_protein,
      requirement       = ear_mean_g_day
    ),
    # OPT ("sub-optimal") inadequacy
    prot_inad_OPT_MM = calculate_inadequacy(
      mean_intake       = protein_grams_MM,
      cv_intake         = cv_protein,
      distribution_type = best_dist_protein,
      requirement       = opt_mean_g_day
    ),
    # Quality-adjusted EAR inadequacy
    prot_inad_EAR_MM_QA = calculate_inadequacy(
      mean_intake       = protein_grams_MM,
      cv_intake         = cv_protein,
      distribution_type = best_dist_protein,
      requirement       = ear_mean_g_day_adj
    )
    # If you want calorie inadequacy summaries later, you can add:
    # cal_inad_MM = calculate_inadequacy(
    #   mean_intake       = kcal_mean_medium,
    #   cv_intake         = cv_calorie,
    #   distribution_type = best_dist_calorie,
    #   requirement       = mder_stratum_kcal
    # )
  ) %>%
  ungroup()

# ==============================================================
# 4. Helper to generate global / sex / region / region×sex summaries
# ==============================================================

make_summary_blocks <- function(df, value_col, scenario_label, metric_label) {
  
  # ---- GLOBAL ----
  global <- df %>%
    summarise(
      prevalence = weighted.mean(.data[[value_col]], w = population, na.rm = TRUE)
    ) %>%
    mutate(
      group_level = "global",
      sex        = NA_character_,
      region_WB  = NA_character_
    )
  
  # ---- BY SEX (world total, male vs female) ----
  by_sex <- df %>%
    group_by(sex) %>%
    summarise(
      prevalence = weighted.mean(.data[[value_col]], w = population, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    mutate(
      group_level = "sex",
      region_WB   = NA_character_
    )
  
  # ---- BY REGION (aggregated over both sexes) ----
  by_region <- df %>%
    group_by(region_WB) %>%
    summarise(
      prevalence = weighted.mean(.data[[value_col]], w = population, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    mutate(
      group_level = "region",
      sex         = NA_character_
    )
  
  # ---- BY REGION × SEX ----
  by_region_sex <- df %>%
    group_by(region_WB, sex) %>%
    summarise(
      prevalence = weighted.mean(.data[[value_col]], w = population, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    mutate(
      group_level = "region_sex"
    )
  
  bind_rows(global, by_sex, by_region, by_region_sex) %>%
    mutate(
      scenario = scenario_label,
      metric   = metric_label
    )
}

make_mean_blocks <- function(df, value_col, scenario_label, metric_label) {
  
  # ---- GLOBAL ----
  global <- df %>%
    summarise(
      value = weighted.mean(.data[[value_col]], w = population, na.rm = TRUE)
    ) %>%
    mutate(
      group_level = "global",
      sex        = NA_character_,
      region_WB  = NA_character_
    )
  
  # ---- BY SEX ----
  by_sex <- df %>%
    group_by(sex) %>%
    summarise(
      value = weighted.mean(.data[[value_col]], w = population, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      group_level = "sex",
      region_WB   = NA_character_
    )
  
  # ---- BY REGION ----
  by_region <- df %>%
    group_by(region_WB) %>%
    summarise(
      value = weighted.mean(.data[[value_col]], w = population, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      group_level = "region",
      sex         = NA_character_
    )
  
  # ---- BY REGION × SEX ----
  by_region_sex <- df %>%
    group_by(region_WB, sex) %>%
    summarise(
      value = weighted.mean(.data[[value_col]], w = population, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      group_level = "region_sex"
    )
  
  bind_rows(global, by_sex, by_region, by_region_sex) %>%
    mutate(
      scenario = scenario_label,
      metric   = metric_label
    )
}


# ==============================================================
# 5. Build all requested summaries
# ==============================================================

# 1) Global inadequacy (EAR), optimal calories
summary_opt_EAR <- make_summary_blocks(
  prot2_central,
  value_col      = "prot_inad_EAR_exact",
  scenario_label = "optimal_calories",
  metric_label   = "protein_inadequacy_EAR"
)

# 2) Global inadequacy (EAR), realistic calories (MM)
summary_real_EAR <- make_summary_blocks(
  scen_central,
  value_col      = "prot_inad_EAR_MM",
  scenario_label = "realistic_calories",
  metric_label   = "protein_inadequacy_EAR"
)

# 3) Global sub-optimal (OPT), optimal calories
summary_opt_OPT <- make_summary_blocks(
  prot2_central,
  value_col      = "prot_inad_OPT_exact",
  scenario_label = "optimal_calories",
  metric_label   = "protein_inadequacy_OPT"
)

# 4) Global sub-optimal (OPT), realistic calories (MM)
summary_real_OPT <- make_summary_blocks(
  scen_central,
  value_col      = "prot_inad_OPT_MM",
  scenario_label = "realistic_calories",
  metric_label   = "protein_inadequacy_OPT"
)

# 5) Global inadequacy, quality-adjusted EAR, optimal calories
summary_opt_EAR_QA <- make_summary_blocks(
  prot2_central,
  value_col      = "prot_inad_EAR_exact_QA",
  scenario_label = "optimal_calories",
  metric_label   = "protein_inadequacy_EAR_quality_adjusted"
)

# 6) Global inadequacy, quality-adjusted EAR, realistic calories (MM)
summary_real_EAR_QA <- make_summary_blocks(
  scen_central,
  value_col      = "prot_inad_EAR_MM_QA",
  scenario_label = "realistic_calories",
  metric_label   = "protein_inadequacy_EAR_quality_adjusted"
)

# ==============================================================
# 5b. Mean-only adequacy summaries (percent of requirement met)
# ==============================================================

mean_opt_EAR <- make_mean_blocks(
  prot2_central,
  value_col      = "prot_mean_pct_EAR_exact",
  scenario_label = "optimal_calories",
  metric_label   = "mean_protein_pct_EAR"
)

mean_real_EAR <- make_mean_blocks(
  scen_central,
  value_col      = "prot_mean_pct_EAR_MM",
  scenario_label = "realistic_calories",
  metric_label   = "mean_protein_pct_EAR"
)

mean_opt_OPT <- make_mean_blocks(
  prot2_central,
  value_col      = "prot_mean_pct_OPT_exact",
  scenario_label = "optimal_calories",
  metric_label   = "mean_protein_pct_OPT"
)

mean_real_OPT <- make_mean_blocks(
  scen_central,
  value_col      = "prot_mean_pct_OPT_MM",
  scenario_label = "realistic_calories",
  metric_label   = "mean_protein_pct_OPT"
)

mean_opt_EAR_QA <- make_mean_blocks(
  prot2_central,
  value_col      = "prot_mean_pct_EAR_exact_QA",
  scenario_label = "optimal_calories",
  metric_label   = "mean_protein_pct_EAR_quality_adjusted"
)

mean_real_EAR_QA <- make_mean_blocks(
  scen_central,
  value_col      = "prot_mean_pct_EAR_MM_QA",
  scenario_label = "realistic_calories",
  metric_label   = "mean_protein_pct_EAR_quality_adjusted"
)

# ==============================================================
# 5c. Prevalence summaries for LOW/MEDIUM/HIGH calorie scenarios
# ==============================================================

summary_cal_inad_3scen <- bind_rows(
  make_summary_blocks(filter(scen3_long, calorie_scenario == "low"),
                      "cal_inad", "realistic_calories_low", "calorie_inadequacy_MDER"),
  make_summary_blocks(filter(scen3_long, calorie_scenario == "medium"),
                      "cal_inad", "realistic_calories_medium", "calorie_inadequacy_MDER"),
  make_summary_blocks(filter(scen3_long, calorie_scenario == "high"),
                      "cal_inad", "realistic_calories_high", "calorie_inadequacy_MDER")
)

summary_prot_EAR_3scen <- bind_rows(
  make_summary_blocks(filter(scen3_long, calorie_scenario == "low"),
                      "prot_inad_EAR", "realistic_calories_low", "protein_inadequacy_EAR"),
  make_summary_blocks(filter(scen3_long, calorie_scenario == "medium"),
                      "prot_inad_EAR", "realistic_calories_medium", "protein_inadequacy_EAR"),
  make_summary_blocks(filter(scen3_long, calorie_scenario == "high"),
                      "prot_inad_EAR", "realistic_calories_high", "protein_inadequacy_EAR")
)

summary_prot_OPT_3scen <- bind_rows(
  make_summary_blocks(filter(scen3_long, calorie_scenario == "low"),
                      "prot_inad_OPT", "realistic_calories_low", "protein_inadequacy_OPT"),
  make_summary_blocks(filter(scen3_long, calorie_scenario == "medium"),
                      "prot_inad_OPT", "realistic_calories_medium", "protein_inadequacy_OPT"),
  make_summary_blocks(filter(scen3_long, calorie_scenario == "high"),
                      "prot_inad_OPT", "realistic_calories_high", "protein_inadequacy_OPT")
)

summary_prot_EAR_QA_3scen <- bind_rows(
  make_summary_blocks(filter(scen3_long, calorie_scenario == "low"),
                      "prot_inad_EAR_QA", "realistic_calories_low", "protein_inadequacy_EAR_quality_adjusted"),
  make_summary_blocks(filter(scen3_long, calorie_scenario == "medium"),
                      "prot_inad_EAR_QA", "realistic_calories_medium", "protein_inadequacy_EAR_quality_adjusted"),
  make_summary_blocks(filter(scen3_long, calorie_scenario == "high"),
                      "prot_inad_EAR_QA", "realistic_calories_high", "protein_inadequacy_EAR_quality_adjusted")
)


# ==============================================================
# 6. Build TWO separate tidy tables + export
# ==============================================================

summary_prev_WB <- bind_rows(
  summary_opt_EAR,
  summary_real_EAR,
  summary_opt_OPT,
  summary_real_OPT,
  summary_opt_EAR_QA,
  summary_real_EAR_QA,
  
  # --- NEW: low/medium/high calorie scenarios ---
  summary_cal_inad_3scen,
  summary_prot_EAR_3scen,
  summary_prot_OPT_3scen,
  summary_prot_EAR_QA_3scen
) %>%
  mutate(prevalence_pct = prevalence * 100) %>%
  select(group_level, region_WB, sex, metric, scenario, prevalence, prevalence_pct) %>%
  arrange(metric, scenario, group_level, region_WB, sex)


cat("\n\n=== Prevalence summaries (WB regions) ===\n")
print(summary_prev_WB, n = nrow(summary_prev_WB))

write_xlsx(
  summary_prev_WB,
  "./output/summary_WB_regions_quality_prevalences.xlsx"
)

# ---- B) Mean-only % of requirement met ----
summary_mean_WB <- bind_rows(
  mean_opt_EAR,
  mean_real_EAR,
  mean_opt_OPT,
  mean_real_OPT,
  mean_opt_EAR_QA,
  mean_real_EAR_QA
) %>%
  # rename "value" to something explicit
  rename(mean_pct_requirement_met = value) %>%
  select(group_level, region_WB, sex, metric, scenario, mean_pct_requirement_met) %>%
  arrange(metric, scenario, group_level, region_WB, sex)

cat("\n\n=== Mean-only summaries (WB regions) ===\n")
print(summary_mean_WB, n = nrow(summary_mean_WB))

write_xlsx(
  summary_mean_WB,
  "./output/summary_WB_regions_quality_mean_only.xlsx"
)

cat("\n✅ Script 7 completed: wrote TWO outputs (prevalences + mean-only).\n")


