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


# ==============================================================
# 4b. Helper: summaries by AGE GROUP
# ==============================================================

make_age_blocks <- function(df, value_col, scenario_label, metric_label) {
  
  by_age <- df %>%
    group_by(age_group) %>%
    summarise(
      prevalence = weighted.mean(.data[[value_col]], w = population, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    mutate(
      group_level = "age",
      sex         = NA_character_,
      region_WB   = NA_character_,
      scenario    = scenario_label,
      metric      = metric_label
    )
  
  by_age
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

# ==============================================================
# 5d. AGE summaries for LOW/MEDIUM/HIGH calorie scenarios (protein EAR)
# ==============================================================

summary_prot_EAR_age_3scen <- bind_rows(
  make_age_blocks(
    filter(scen3_long, calorie_scenario == "low"),
    value_col      = "prot_inad_EAR",
    scenario_label = "realistic_calories_low",
    metric_label   = "protein_inadequacy_EAR"
  ),
  make_age_blocks(
    filter(scen3_long, calorie_scenario == "medium"),
    value_col      = "prot_inad_EAR",
    scenario_label = "realistic_calories_medium",
    metric_label   = "protein_inadequacy_EAR"
  ),
  make_age_blocks(
    filter(scen3_long, calorie_scenario == "high"),
    value_col      = "prot_inad_EAR",
    scenario_label = "realistic_calories_high",
    metric_label   = "protein_inadequacy_EAR"
  )
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
  summary_prot_EAR_QA_3scen,
  
  # --- NEW: age block for Fig2 panel C ---
  summary_prot_EAR_age_3scen
) %>%
  mutate(prevalence_pct = prevalence * 100) %>%
  # keep the original columns, and also keep age_group if present
  select(group_level, region_WB, sex, age_group, metric, scenario, prevalence, prevalence_pct) %>%
  arrange(metric, scenario, group_level, region_WB, sex, age_group)





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



library(dplyr)
library(ggplot2)
library(forcats)

# ---- Region × sex data ----
fig1a_region <- summary_mean_WB %>%
  filter(
    scenario == "optimal_calories",
    metric == "mean_protein_pct_EAR",
    group_level == "region_sex"
  )

# ---- Global × sex data ----
fig1a_global <- summary_mean_WB %>%
  filter(
    scenario == "optimal_calories",
    metric == "mean_protein_pct_EAR",
    group_level == "sex"
  ) %>%
  mutate(region_WB = "Global")

# ---- Combine ----
fig1a_data <- bind_rows(fig1a_global, fig1a_region) %>%
  mutate(
    sex = case_when(
      sex %in% c("Female", "female", "F") ~ "Female",
      sex %in% c("Male", "male", "M")     ~ "Male",
      TRUE ~ sex
    )
  )

# ---- Order regions: Global first, others by mean ----
region_order <- fig1a_data %>%
  filter(region_WB != "Global") %>%
  group_by(region_WB) %>%
  summarise(order_val = mean(mean_pct_requirement_met, na.rm = TRUE), .groups = "drop") %>%
  arrange(order_val)

region_levels <- c("Global", region_order$region_WB)

fig1a_data <- fig1a_data %>%
  mutate(region_WB = factor(region_WB, levels = region_levels))

# ---- Plot ----
p_fig1a <- ggplot(
  fig1a_data,
  aes(x = region_WB, y = mean_pct_requirement_met, fill = sex)
) +
  geom_col(
    position = position_dodge(width = 0.55),
    width = 0.45
  ) +
  # EAR reference line ON TOP
  geom_hline(
    yintercept = 100,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.7
  ) +
  labs(
    x = NULL,
    y = "Mean protein intake (% of EAR)",
    fill = NULL
  ) +
  coord_cartesian(
    ylim = c(0, max(fig1a_data$mean_pct_requirement_met, na.rm = TRUE) * 1.05)
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "top"
  )

# ---- Save PNG only ----
ggsave(
  "./output/fig1_panel_a_wb_region_sex_bar_global.png",
  p_fig1a,
  width = 7.6,
  height = 5,
  dpi = 400
)

p_fig1a


# ============================================================
# Figure 1b — Within-region variability (IQR summary; no outliers)
#   Shows IQR of COUNTRY-level mean %EAR within each WB region
#   (Optional: whiskers set to p10–p90 instead of min–max)
#   Requires: prot2_central (iso3, region_WB, population, prot_mean_pct_EAR_exact)
# ============================================================

library(dplyr)
library(ggplot2)
library(forcats)

# ---- 1) Country-level mean %EAR (within country, weighted over strata) ----
country_means <- prot2_central %>%
  group_by(iso3, region_WB) %>%
  summarise(
    country_mean_pct_EAR = weighted.mean(
      prot_mean_pct_EAR_exact,
      w = population,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  filter(is.finite(country_mean_pct_EAR))

# ---- 2) Region boxes: p10–p90 whiskers + IQR box ----
region_boxes <- country_means %>%
  filter(!is.na(region_WB)) %>%
  group_by(region_WB) %>%
  summarise(
    ymin   = quantile(country_mean_pct_EAR, 0.10, na.rm = TRUE),  # p10
    lower  = quantile(country_mean_pct_EAR, 0.25, na.rm = TRUE),  # Q1
    middle = quantile(country_mean_pct_EAR, 0.50, na.rm = TRUE),  # median
    upper  = quantile(country_mean_pct_EAR, 0.75, na.rm = TRUE),  # Q3
    ymax   = quantile(country_mean_pct_EAR, 0.90, na.rm = TRUE),  # p90
    .groups = "drop"
  )

# ---- 3) Global box computed across ALL countries (no country points shown) ----
global_box <- country_means %>%
  summarise(
    ymin   = quantile(country_mean_pct_EAR, 0.10, na.rm = TRUE),
    lower  = quantile(country_mean_pct_EAR, 0.25, na.rm = TRUE),
    middle = quantile(country_mean_pct_EAR, 0.50, na.rm = TRUE),
    upper  = quantile(country_mean_pct_EAR, 0.75, na.rm = TRUE),
    ymax   = quantile(country_mean_pct_EAR, 0.90, na.rm = TRUE)
  ) %>%
  mutate(region_WB = "Global")

# ---- 4) Combine + order (Global first; rest match Panel a if available) ----
region_boxes2 <- bind_rows(global_box, region_boxes)

if (exists("region_levels")) {
  # Ensure "Global" exists and is first
  region_levels2 <- c("Global", setdiff(region_levels, "Global"))
  region_boxes2 <- region_boxes2 %>%
    mutate(region_WB = factor(region_WB, levels = region_levels2))
} else {
  # Global first, rest ordered by median
  region_levels2 <- c(
    "Global",
    region_boxes %>% arrange(middle) %>% pull(region_WB)
  )
  region_boxes2 <- region_boxes2 %>%
    mutate(region_WB = factor(region_WB, levels = region_levels2))
}

# ---- 5) Plot: thin boxes with p10–p90 whiskers ----
p_fig1b_boxes_global <- ggplot(region_boxes2, aes(x = region_WB)) +
  geom_boxplot(
    aes(
      ymin = ymin,
      lower = lower,
      middle = middle,
      upper = upper,
      ymax = ymax
    ),
    stat = "identity",
    width = 0.32,          # thinner
    linewidth = 0.6,
    fill = "grey85",
    color = "grey30"
  ) +
  geom_hline(
    yintercept = 100,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.7
  ) +
  labs(
    x = NULL,
    y = "Country mean protein intake (% of EAR)"
  ) +
  coord_cartesian(
    ylim = c(0, max(region_boxes2$ymax, na.rm = TRUE) * 1.05)
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

# ---- 6) Save PNG ----
ggsave(
  "./output/fig1b_region_boxes_p10_p90_with_global.png",
  p_fig1b_boxes_global,
  width = 7.6,
  height = 2.8,
  dpi = 400
)

p_fig1b_boxes_global




# ============================================================
# Figure 1 (a–b): side-by-side + shared y-axis limits
# ============================================================

library(patchwork)

# ---- Shared y-axis upper limit (use the plotted data sources) ----
ymax_common <- max(
  max(fig1a_data$mean_pct_requirement_met, na.rm = TRUE),
  max(region_boxes2$ymax, na.rm = TRUE)
) * 1.05

# ---- Rebuild plots with identical y-limits ----
p_fig1a_shared <- p_fig1a + coord_cartesian(ylim = c(0, ymax_common))
p_fig1b_shared <- p_fig1b_boxes_global + coord_cartesian(ylim = c(0, ymax_common))

# ---- Combine side-by-side, balanced widths, add panel tags ----
fig1_ab <- (p_fig1a_shared | p_fig1b_shared) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "a")

# ---- Save PNG (two equal half-page panels) ----
ggsave(
  "./output/Figure1_ab_side_by_side_shared_axes.png",
  fig1_ab,
  width = 15.2,   # ~2 × 7.6
  height = 5.0,
  dpi = 400
)

fig1_ab




# ============================================================
# Figure 3 — STANDALONE MAPS (NO STITCHING)
#   Produces THREE separate PNGs (so you can stitch later):
#     (1) EAR inadequacy map (realistic calories, medium)
#     (2) Proportion map: Protein inadequacy (EAR) relative to calories (binary 3-level)
#     (3) "Quality" map loaded as in Script 6 (ASF share fallback)
#
# NOTES
# - No density map.
# - No multi-panel composition in this block.
# - Uses Script 6 palettes + robust ISO parenting join logic.
# - Assumes Script 7 already created: dat, scen3_long, prot2_central (for optional later),
#   and your normalize_iso3 logic is defined here (self-contained).
# ============================================================

library(dplyr)
library(ggplot2)
library(scales)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)

# ---- Okabe–Ito palette (same as Script 6) ----
oi <- list(
  black   = "#000000",
  orange  = "#E69F00",
  sky     = "#56B4E9",
  green   = "#009E73",
  yellow  = "#F0E442",
  blue    = "#0072B2",
  vermil  = "#D55E00",
  purple  = "#CC79A7",
  greyNA  = "#E5E5E5"
)

na_gray <- "#8c8c8c"   # match Panel a NA color across all panels


# ============================================================
# 1) World geometry + robust parent-ISO3 join (from Script 6)
# ============================================================

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  mutate(
    iso3_iso_a3     = if ("iso_a3"     %in% names(.)) iso_a3     else NA_character_,
    iso3_iso_a3_eh  = if ("iso_a3_eh"  %in% names(.)) iso_a3_eh  else NA_character_,
    iso3_adm0_a3    = if ("adm0_a3"    %in% names(.)) adm0_a3    else NA_character_,
    iso3_wb_a3      = if ("wb_a3"      %in% names(.)) wb_a3      else NA_character_
  )

normalize_iso3 <- function(code) {
  dplyr::case_when(
    code %in% c("XKX","KOS") ~ "KOS",
    code == "SDS" ~ "SSD",
    code == "ROM" ~ "ROU",
    code == "ZAR" ~ "COD",
    code == "TMP" ~ "TLS",
    code == "WBG" ~ "PSE",
    TRUE ~ code
  )
}

territory_to_parent <- c(
  # France
  "GUF"="FRA","GLP"="FRA","MTQ"="FRA","REU"="FRA","MYT"="FRA","NCL"="FRA",
  "PYF"="FRA","WLF"="FRA","SPM"="FRA","BLM"="FRA","MAF"="FRA","ATF"="FRA","MFO"="FRA",
  # Norway
  "SJM"="NOR","BVT"="NOR",
  # Denmark
  "GRL"="DNK","FRO"="DNK",
  # Netherlands
  "ABW"="NLD","CUW"="NLD","BES"="NLD","SXM"="NLD",
  # UK
  "GGY"="GBR","JEY"="GBR","IMN"="GBR",
  # USA
  "PRI"="USA","VIR"="USA","GUM"="USA","MNP"="USA","ASM"="USA",
  # China SAR/Macao
  "HKG"="CHN","MAC"="CHN"
)

world <- world %>%
  mutate(
    ne_iso_raw  = coalesce(iso3_adm0_a3, iso3_iso_a3_eh, iso3_iso_a3, iso3_wb_a3),
    ne_iso_norm = normalize_iso3(ne_iso_raw),
    parent_iso3 = dplyr::case_when(
      ne_iso_norm %in% names(territory_to_parent) ~ territory_to_parent[ne_iso_norm],
      TRUE ~ ne_iso_norm
    )
  )

# ---- Shared tight viewport (match Script 6) ----
lat_lim <- c(-58, 85)
lon_lim <- c(-180, 180)

common_coord_sf <- function() {
  coord_sf(
    xlim   = lon_lim,
    ylim   = lat_lim,
    expand = FALSE,
    clip   = "on",
    datum  = NA
  )
}

tight_map_theme <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.text   = element_blank(),
      axis.title  = element_blank(),
      panel.grid  = element_blank(),
      plot.margin = margin(3, 4, 3, 4)
    )
}

# ---- Continuous red ramp (same as Script 6 continuous protein maps) ----
map_fill_cont_red <- function(cap, label_fun = percent_format(accuracy = 1)) {
  scale_fill_gradientn(
    colours  = c("#f0f0f0", "#fdbb84", "#e34a33"),
    limits   = c(0, cap),
    oob      = scales::squish,
    labels   = label_fun,
    na.value = na_gray
  )
}


# ---- Map saver to ./output (PNG) ----
save_map_output <- function(plot, filename, width = 12, height = 6.2, dpi = 400) {
  ggsave(
    filename = file.path("./output", paste0(filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
  message("✅ Saved: ", file.path("./output", paste0(filename, ".png")))
}

# ============================================================
# 2) Build country-level inputs from Script 7 objects
#    Requires: scen3_long created in your Script 7
# ============================================================

stopifnot(exists("scen3_long"))

# ---- EAR inadequacy (realistic calories, medium) ----
country_EAR_real <- scen3_long %>%
  filter(calorie_scenario == "medium") %>%
  group_by(iso3) %>%
  summarise(
    prot_EAR = weighted.mean(prot_inad_EAR, w = population, na.rm = TRUE),
    cal_MD   = weighted.mean(cal_inad,      w = population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(parent_iso3 = normalize_iso3(iso3))

map_EAR_df <- world %>% left_join(country_EAR_real, by = "parent_iso3")

# ============================================================
# 3) PANEL 1 — EAR inadequacy map (realistic, medium)
# ============================================================

p_map_EAR <- ggplot(map_EAR_df) +
  geom_sf(aes(fill = prot_EAR), color = NA) +
  map_fill_cont_red(cap = 0.30, label_fun = percent_format(accuracy = 1)) +
  common_coord_sf() +
  labs(title = "Protein inadequacy (EAR)") +
  tight_map_theme()

save_map_output(p_map_EAR, "Fig3_panel_EAR_inadequacy")

# ============================================================
# 4) PANEL 2 — Proportion map: Protein (EAR) relative to calories
#    2 bins only: ratio > 1 vs ratio < 1
#    (ties at exactly 1 treated as NA OR assigned to one side; choose below)
# ============================================================

# --- Two-class palette (keep the same “muted” feel as before) ---
ratio2_palette <- c(
  "<1" = "#2A9D8F",   # muted teal
  ">1" = "#4A4A4A"    # charcoal
)



map_ratio2_df <- map_EAR_df %>%
  mutate(
    ratio = prot_EAR / cal_MD,
    ratio_cat = case_when(
      is.na(ratio) ~ NA_character_,
      ratio < 1    ~ "<1",
      ratio > 1    ~ ">1",
      TRUE         ~ NA_character_   # exactly 1 -> NA (rare). If you prefer, set to "<1" or ">1".
    ),
    ratio_cat = factor(ratio_cat, levels = c("<1", ">1"))
  )

p_map_ratio2 <- ggplot(map_ratio2_df) +
  geom_sf(aes(fill = ratio_cat), color = NA) +
  scale_fill_manual(values = ratio2_palette, na.value = na_gray, drop = FALSE) +
  common_coord_sf() +
  labs(
    title = "Protein inadequacy (EAR) relative to caloric inadequacy",
    fill  = "Ratio"
  ) +
  tight_map_theme() +
  theme(legend.position = "bottom")

save_map_output(
  p_map_ratio2,
  "Fig3_panel_proportion_ratio_EAR_to_calories_2class",
  width = 12,
  height = 6.6
)

p_map_ratio2

# ============================================================
# 5) PANEL 3 — "Quality" map (loaded as in Script 6)
#     - tries TRUE protein quality if you join it
#     - otherwise falls back to ASF share (explicit label)
# ============================================================

protein_props <- readRDS("./output/protein_asf_props.rds") %>%
  filter(year == 2018) %>%
  mutate(parent_iso3 = normalize_iso3(iso3))

# OPTIONAL: if you have a TRUE quality file, load + join here.
# Example:
# protein_quality_file <- readRDS("./output/protein_quality_country.rds") %>%
#   mutate(parent_iso3 = normalize_iso3(iso3))
# protein_props <- protein_props %>%
#   left_join(protein_quality_file %>% select(parent_iso3, protein_quality), by = "parent_iso3")

if ("protein_quality" %in% names(protein_props)) {
  country_quality <- protein_props %>%
    group_by(parent_iso3) %>%
    summarise(protein_quality = mean(protein_quality, na.rm = TRUE), .groups = "drop")
  quality_title <- "Protein quality"
  quality_fill_scale <- scale_fill_viridis_c(option = "mako", direction = -1, na.value = na_gray)
  quality_labels <- waiver()
} else {
  country_quality <- protein_props %>%
    group_by(parent_iso3) %>%
    summarise(protein_quality = mean(prop_asf, na.rm = TRUE), .groups = "drop")
  quality_title <- "Protein from animal-source foods (% of total protein)"
  quality_fill_scale <- scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    limits = c(0, 1),
    labels = percent_format(accuracy = 1),
    na.value = na_gray
  )
  quality_labels <- percent_format(accuracy = 1)
}

map_quality_df <- world %>% left_join(country_quality, by = "parent_iso3")

p_map_quality <- ggplot(map_quality_df) +
  geom_sf(aes(fill = protein_quality), color = NA) +
  quality_fill_scale +
  common_coord_sf() +
  labs(title = quality_title, fill = NULL) +
  tight_map_theme()

save_map_output(p_map_quality, "Fig3_panel_quality_or_ASFshare", width = 12, height = 6.2)

# ---- Print to device (optional) ----
print(p_map_EAR)
print(p_map_ratio2)
print(p_map_quality)


# ============================================================
# PANEL (candidate) — Sex gap map (Female − Male)
#   Protein inadequacy (EAR) under MM calories:
#   Uses scen3_long filtered to calorie_scenario == "medium"
#   (No dependency on scen_calc / Script 6)
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ---- 0) Sanity: require objects that should already exist ----
stopifnot(exists("scen3_long"))
stopifnot(exists("world"))
stopifnot(exists("normalize_iso3"))
stopifnot(exists("common_coord_sf"))
stopifnot(exists("tight_map_theme"))
stopifnot(exists("oi"))
stopifnot(exists("save_map_output"))

# ---- 1) Country×sex protein inadequacy (EAR) in "MM" calories = medium ----
country_sex_MM_EAR <- scen3_long %>%
  filter(calorie_scenario == "medium") %>%
  group_by(iso3, sex) %>%
  summarise(
    prot_EAR = weighted.mean(prot_inad_EAR, w = population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sex = case_when(
      sex %in% c("Female", "female", "Females", "F") ~ "Female",
      sex %in% c("Male", "male", "Males", "M")       ~ "Male",
      TRUE ~ as.character(sex)
    )
  ) %>%
  filter(sex %in% c("Female", "Male"))

# ---- 2) Wide + compute gap (Female − Male) ----
country_gap <- country_sex_MM_EAR %>%
  pivot_wider(names_from = sex, values_from = prot_EAR) %>%
  mutate(
    gap_FminusM = Female - Male,
    parent_iso3 = normalize_iso3(iso3)
  )

# ---- 3) Join to map geometry ----
map_gap_df <- world %>%
  left_join(country_gap, by = "parent_iso3")

# ============================================================
# Sex gap map — continuous, data-informed full scale
#   Female − Male protein inadequacy (EAR), MM scenario
# ============================================================

# Choose data-informed limits based on observed distribution
gap_limits <- c(-0.01, 0.1501)   # −1 pp to +15 pp

p_map_sexgap <- ggplot(map_gap_df) +
  geom_sf(aes(fill = gap_FminusM), color = NA) +
  scale_fill_gradient2(
    low = "#0072B2",      # blue: Male higher (negative)
    mid = "grey95",
    high = "#D55E00",     # vermilion: Female higher
    midpoint = 0,
    limits = gap_limits,
    oob = squish,
    labels = percent_format(accuracy = 1),
    na.value = na_gray
  ) +
  common_coord_sf() +
  labs(
    title = "Sex gap in protein inadequacy (EAR), medium-calorie scenario",
    fill  = "Female − Male\n(pp)"
  ) +
  tight_map_theme() +
  theme(legend.position = "bottom")

save_map_output(
  p_map_sexgap,
  "Fig3_panel_sexgap_FminusM_proteinEAR_real_scale",
  width = 12,
  height = 6.6
)

p_map_sexgap


# ---- Inspect full country-level table ----
# country_gap_inspect <- country_gap %>%
#   arrange(desc(abs(gap_FminusM)))
# 
# print(country_gap_inspect, n = nrow(country_gap_inspect))
# summary(country_gap$gap_FminusM)
# 
# quantile(country_gap$gap_FminusM, probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = TRUE)



# ============================================================
# Figure 3 — 2×2 GRID (FOUR MAPS) + SAVE
#   Uses the FOUR ggplot objects you already created above:
#     p_map_EAR
#     p_map_ratio2
#     p_map_quality
#     p_map_sexgap
#
# Output:
#   ./output/Figure3_2x2_maps.png
#
# Notes:
# - No re-computation: just stitches existing plots.
# - Uses patchwork; keeps your tight map theme.
# - Legend handling:
#     * default = keep per-panel legends (often looks fine)
#     * optional = collect legends (can get messy since scales differ)
# ============================================================

library(patchwork)

# ---- (Optional) add panel tags + make titles consistent size ----
pA <- p_map_EAR     + labs(tag = "a") + theme(plot.title = element_text(size = 13, face = "bold"))
pB <- p_map_ratio2  + labs(tag = "b") + theme(plot.title = element_text(size = 13, face = "bold"))
pC <- p_map_sexgap  + labs(tag = "c") + theme(plot.title = element_text(size = 13, face = "bold"))
pD <- p_map_quality + labs(tag = "d") + theme(plot.title = element_text(size = 13, face = "bold"))

# ---- 2×2 grid (balanced) ----
fig3_2x2 <- (pA | pB) / (pC | pD) +
  plot_layout(widths = c(1, 1), heights = c(1, 1)) +
  plot_annotation(
    theme = theme(
      plot.margin = margin(6, 6, 6, 6)
    )
  )

# ---- Save (tweak size if you want more breathing room) ----
ggsave(
  filename = "./output/Figure3_2x2_maps.png",
  plot = fig3_2x2,
  width = 14.5,
  height = 10.0,
  dpi = 450,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

fig3_2x2


# ============================================================
# Fig2 — WB region protein inadequacy (EAR)
#   Main estimate = MEDIUM scenario (big point)
#   Sensitivity   = LOW–HIGH scenarios (range line)
#
# Requires: summary_prev_WB from Script 7
# Output: ./output/Fig2_WB_regions_proteinEAR_medium_point_with_low_high_range.png
# ============================================================

library(dplyr)
library(ggplot2)
library(forcats)
library(scales)
library(tidyr)

# ---- 1) Pull protein EAR prevalences for the 3 calorie scenarios ----
df3 <- summary_prev_WB %>%
  filter(
    metric == "protein_inadequacy_EAR",
    scenario %in% c("realistic_calories_low", "realistic_calories_medium", "realistic_calories_high"),
    group_level %in% c("global", "region")
  ) %>%
  mutate(
    region_WB = if_else(group_level == "global", "Global", region_WB),
    scen = recode(
      scenario,
      realistic_calories_low    = "low",
      realistic_calories_medium = "medium",
      realistic_calories_high   = "high"
    )
  ) %>%
  filter(!is.na(region_WB)) %>%
  group_by(region_WB, scen) %>%                 # safety against duplicates
  summarise(prevalence_pct = mean(prevalence_pct, na.rm = TRUE), .groups = "drop")

# ---- 2) Wide: one row per region with low/medium/high ----
df_wide <- df3 %>%
  pivot_wider(names_from = scen, values_from = prevalence_pct) %>%
  mutate(
    low    = as.numeric(low),
    medium = as.numeric(medium),
    high   = as.numeric(high)
  )

# ---- 3) Order regions by MEDIUM (Global first) ----
region_levels <- df_wide %>%
  filter(region_WB != "Global") %>%
  arrange(medium) %>%
  pull(region_WB)

df_wide <- df_wide %>%
  mutate(region_WB = factor(region_WB, levels = c("Global", region_levels)))

# ---- 4) Plot: low–high linerange + emphasized medium point ----
p_fig2 <- ggplot(df_wide, aes(x = region_WB)) +
  geom_linerange(
    aes(ymin = low, ymax = high),
    linewidth = 0.9,
    color = "grey60"
  ) +
  # optional little end-caps (makes it feel more “CI-like”)
  geom_point(aes(y = low),  shape = 95, size = 6, color = "grey60") +
  geom_point(aes(y = high), shape = 95, size = 6, color = "grey60") +
  # main estimate (medium) — visually dominant
  geom_point(
    aes(y = medium),
    size = 3.6,
    color = "grey15"
  ) +
  labs(
    x = NULL,
    y = "Protein inadequacy (EAR), % of population",
    title = "Protein inadequacy by World Bank region (medium scenario; low–high range)"
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title  = element_text(face = "bold"),
    plot.margin = margin(6, 6, 6, 6)
  )

ggsave(
  filename = "./output/Fig2_WB_regions_proteinEAR_medium_point_with_low_high_range.png",
  plot = p_fig2,
  width = 9.2,
  height = 4.8,
  dpi = 450,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

p_fig2



# ============================================================
# Fig2 — SEX protein inadequacy (EAR)
#   Main estimate = MEDIUM scenario (big point)
#   Sensitivity   = LOW–HIGH scenarios (range line)
#
# Requires: summary_prev_WB from Script 7
# Output: ./output/Fig2_sex_proteinEAR_medium_point_with_low_high_range.png
# ============================================================

library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)

# ---- 1) Pull protein EAR prevalences for the 3 calorie scenarios (by sex) ----
df3_sex <- summary_prev_WB %>%
  filter(
    metric == "protein_inadequacy_EAR",
    group_level == "sex",
    scenario %in% c("realistic_calories_low", "realistic_calories_medium", "realistic_calories_high")
  ) %>%
  mutate(
    sex = case_when(
      sex %in% c("Females", "Female", "female", "F") ~ "Female",
      sex %in% c("Males", "Male", "male", "M")       ~ "Male",
      TRUE ~ as.character(sex)
    ),
    scen = recode(
      scenario,
      realistic_calories_low    = "low",
      realistic_calories_medium = "medium",
      realistic_calories_high   = "high"
    )
  ) %>%
  filter(!is.na(sex)) %>%
  group_by(sex, scen) %>%  # safety against duplicates
  summarise(prevalence_pct = mean(prevalence_pct, na.rm = TRUE), .groups = "drop")

# ---- 2) Wide: one row per sex with low/medium/high ----
df_wide_sex <- df3_sex %>%
  pivot_wider(names_from = scen, values_from = prevalence_pct) %>%
  mutate(
    low    = as.numeric(low),
    medium = as.numeric(medium),
    high   = as.numeric(high),
    sex    = factor(sex, levels = c("Female", "Male"))
  )

# ---- 3) Plot: low–high linerange + emphasized medium point ----
p_fig2_sex <- ggplot(df_wide_sex, aes(x = sex)) +
  geom_linerange(
    aes(ymin = low, ymax = high),
    linewidth = 0.9,
    color = "grey60"
  ) +
  # optional little end-caps (makes it feel more “CI-like”)
  geom_point(aes(y = low),  shape = 95, size = 6, color = "grey60") +
  geom_point(aes(y = high), shape = 95, size = 6, color = "grey60") +
  # main estimate (medium) — visually dominant
  geom_point(
    aes(y = medium),
    size = 3.6,
    color = "grey15"
  ) +
  labs(
    x = NULL,
    y = "Protein inadequacy (EAR), % of population",
    title = "Protein inadequacy by sex (medium scenario; low–high range)"
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title  = element_text(face = "bold"),
    plot.margin = margin(6, 6, 6, 6)
  )

ggsave(
  filename = "./output/Fig2_sex_proteinEAR_medium_point_with_low_high_range.png",
  plot = p_fig2_sex,
  width = 5.2,
  height = 4.8,
  dpi = 450,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

p_fig2_sex


# ============================================================
# Fig2 — AGE protein inadequacy (EAR)
#   Main estimate = MEDIUM scenario (big point)
#   Sensitivity   = LOW–HIGH scenarios (range line)
#
# Requires: summary_prev_WB from Script 7
# Output: ./output/Fig2_age_proteinEAR_medium_point_with_low_high_range.png
# ============================================================

library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)

# ---- 0) Quick sanity: do we have an age column? ----
# If this errors or prints NULL for age-related fields, you need to add age summaries in Script 7.
cat("\nColumns in summary_prev_WB:\n")
print(names(summary_prev_WB))
cat("\nGroup levels available:\n")
print(sort(unique(summary_prev_WB$group_level)))

# ---- 1) Pull protein EAR prevalences for the 3 calorie scenarios (by age) ----
# EXPECTS:
#   group_level == "age"
#   and an age column named "age_group"
df3_age <- summary_prev_WB %>%
  filter(
    metric == "protein_inadequacy_EAR",
    group_level == "age",
    scenario %in% c("realistic_calories_low", "realistic_calories_medium", "realistic_calories_high")
  ) %>%
  mutate(
    scen = recode(
      scenario,
      realistic_calories_low    = "low",
      realistic_calories_medium = "medium",
      realistic_calories_high   = "high"
    )
  ) %>%
  filter(!is.na(age_group)) %>%
  group_by(age_group, scen) %>%  # safety against duplicates
  summarise(prevalence_pct = mean(prevalence_pct, na.rm = TRUE), .groups = "drop")

# ---- 2) Wide: one row per age group with low/medium/high ----
df_wide_age <- df3_age %>%
  pivot_wider(names_from = scen, values_from = prevalence_pct) %>%
  mutate(
    low    = as.numeric(low),
    medium = as.numeric(medium),
    high   = as.numeric(high)
  )

# ---- 3) Order age groups (keep natural order if age_group is already a factor) ----
# If age_group is character like "1-4","5-9",..., we preserve that order as it appears
# in the data; tweak here if you want a custom order.
if (!is.factor(df_wide_age$age_group)) {
  age_levels <- df_wide_age$age_group
  df_wide_age <- df_wide_age %>%
    mutate(age_group = factor(age_group, levels = unique(age_levels)))
}

# ---- 4) Plot: low–high linerange + emphasized medium point ----
p_fig2_age <- ggplot(df_wide_age, aes(x = age_group)) +
  geom_linerange(
    aes(ymin = low, ymax = high),
    linewidth = 0.9,
    color = "grey60"
  ) +
  geom_point(aes(y = low),  shape = 95, size = 6, color = "grey60") +
  geom_point(aes(y = high), shape = 95, size = 6, color = "grey60") +
  geom_point(
    aes(y = medium),
    size = 3.6,
    color = "grey15"
  ) +
  labs(
    x = NULL,
    y = "Protein inadequacy (EAR), % of population",
    title = "Protein inadequacy by age (medium scenario; low–high range)"
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title  = element_text(face = "bold"),
    plot.margin = margin(6, 6, 6, 6)
  )

ggsave(
  filename = "./output/Fig2_age_proteinEAR_medium_point_with_low_high_range.png",
  plot = p_fig2_age,
  width = 9.2,
  height = 4.8,
  dpi = 450,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

p_fig2_age


# ============================================================
# Fig2 panel C — AGE protein inadequacy (EAR)
#   Main estimate = MEDIUM scenario (big point)
#   Sensitivity   = LOW–HIGH scenarios (range line)
#
# Requires: summary_prev_WB now includes group_level == "age" and age_group
# Output: ./output/Fig2_age_proteinEAR_medium_point_with_low_high_range.png
# ============================================================

library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)

df3_age <- summary_prev_WB %>%
  filter(
    metric == "protein_inadequacy_EAR",
    group_level == "age",
    scenario %in% c("realistic_calories_low", "realistic_calories_medium", "realistic_calories_high")
  ) %>%
  mutate(
    scen = recode(
      scenario,
      realistic_calories_low    = "low",
      realistic_calories_medium = "medium",
      realistic_calories_high   = "high"
    )
  ) %>%
  filter(!is.na(age_group)) %>%
  group_by(age_group, scen) %>%
  summarise(prevalence_pct = mean(prevalence_pct, na.rm = TRUE), .groups = "drop")

df_wide_age <- df3_age %>%
  pivot_wider(names_from = scen, values_from = prevalence_pct) %>%
  mutate(
    low    = as.numeric(low),
    medium = as.numeric(medium),
    high   = as.numeric(high)
  )

# Keep your existing order from the analysis (usually already ordered nicely).
# If you want to force chronological order, set factor levels explicitly here.
df_wide_age <- df_wide_age %>%
  mutate(age_group = factor(age_group, levels = unique(age_group)))

# ---- Force chronological order by numeric start age ----
age_levels <- df_wide_age %>%
  distinct(age_group) %>%
  mutate(
    age_start = suppressWarnings(as.numeric(sub("^([0-9]+).*", "\\1", age_group)))
  ) %>%
  arrange(age_start) %>%
  pull(age_group)

df_wide_age <- df_wide_age %>%
  mutate(age_group = factor(age_group, levels = age_levels))


# ---- 3) Order age groups using the canonical order from `dat` ----
age_levels <- dat %>%
  distinct(age_group) %>%
  pull(age_group)



p_fig2_age <- ggplot(df_wide_age, aes(x = age_group)) +
  geom_linerange(
    aes(ymin = low, ymax = high),
    linewidth = 0.9,
    color = "grey60"
  ) +
  geom_point(aes(y = low),  shape = 95, size = 6, color = "grey60") +
  geom_point(aes(y = high), shape = 95, size = 6, color = "grey60") +
  geom_point(aes(y = medium), size = 3.6, color = "grey15") +
  labs(
    x = NULL,
    y = "Protein inadequacy (EAR), % of population",
    title = "Protein inadequacy by age (medium scenario; low–high range)"
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title  = element_text(face = "bold"),
    plot.margin = margin(6, 6, 6, 6)
  )

ggsave(
  filename = "./output/Fig2_age_proteinEAR_medium_point_with_low_high_range.png",
  plot = p_fig2_age,
  width = 9.2,
  height = 4.8,
  dpi = 450,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

p_fig2_age


# ---- Shared y-axis upper limit for FIGURE 2 (prevalence %) ----
ymax_common_fig2 <- max(
  df_wide$high,
  df_wide_sex$high,
  df_wide_age$high,
  na.rm = TRUE
) * 1.5

# ============================================================
# Figure 2 — final assembly (cleaned)
#   - Titles simplified
#   - Y-axis label removed from panel c
# ============================================================

library(patchwork)





# ---- Shared y-axis upper limit (already computed earlier) ----
# Uses ymax_common from your previous step

pA <- p_fig2 +
  coord_cartesian(ylim = c(0, ymax_common_fig2)) +
  labs(
    title = "Protein inadequacy by World Bank region",
    tag = "a"
  ) +
  theme(
    plot.title = element_text(size = 13, face = "bold")
  )

pB <- p_fig2_sex +
  coord_cartesian(ylim = c(0, ymax_common_fig2)) +
  labs(
    title = "Protein inadequacy by sex",
    tag = "b"
  ) +
  theme(
    plot.title = element_text(size = 13, face = "bold")
  )

pC <- p_fig2_age +
  coord_cartesian(ylim = c(0, ymax_common_fig2)) +
  labs(
    title = "Protein inadequacy by age",
    tag = "c",
    y = NULL              # <-- remove y-axis label here
  ) +
  theme(
    plot.title = element_text(size = 13, face = "bold")
  )

# ---- Layout: a on top, b–c below ----
fig2_abc <- pA / (pB | pC) +
  plot_layout(heights = c(1.1, 1)) +
  plot_annotation(
    theme = theme(
      plot.margin = margin(6, 6, 6, 6)
    )
  )

ggsave(
  filename = "./output/Figure2_abc_panels.png",
  plot = fig2_abc,
  width = 14.5,
  height = 9.2,
  dpi = 450,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

fig2_abc


