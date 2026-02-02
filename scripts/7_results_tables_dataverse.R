# ============================================================
# Dataverse export — Protein FULL + REDUCED (NO "optimal" scenario)
#   FULL:  low / medium / high calorie scenarios only
#   REDUCED: medium-only stored EAR inadequacy (as per Dataverse spec)
# ============================================================

library(dplyr)
library(tidyr)
library(readr)
library(countrycode)
library(stringr)

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); setwd("..")
options(scipen = 999)

# -----------------------------
# Load main dataset
# -----------------------------
dat <- readRDS("./output/dat_quality.rds") %>%
  distinct(iso3, sex, age_group, .keep_all = TRUE) %>%
  filter(age_group != "0-0.99")

# -----------------------------
# Helpers
# -----------------------------
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

format_age_endash <- function(x) {
  str_replace_all(as.character(x), "-", "\u2013")
}

format_sex <- function(x) {
  case_when(
    x %in% c("Females", "Female") ~ "females",
    x %in% c("Males", "Male")     ~ "males",
    TRUE ~ tolower(as.character(x))
  )
}

# -----------------------------
# Attach country + WB region (as in your Script 7)
#   NOTE: countrycode destination for WB region is "region"
# -----------------------------
region_lookup <- dat %>%
  distinct(iso3) %>%
  mutate(
    country = countrycode(iso3, origin = "iso3c", destination = "country.name"),
    region  = countrycode(iso3, origin = "iso3c", destination = "region")
  )

dat2 <- dat %>%
  left_join(region_lookup, by = "iso3") %>%
  mutate(
    sex = format_sex(sex),
    age_group = format_age_endash(age_group),
    
    # ---- FIX population once, upstream ----
    population = as.integer(round(population))
  )


# ============================================================
# FULL dataset (Dataverse "Full"): realistic scenarios only
#   calorie_scenario: low / medium / high
# ============================================================
protein_full <- dat2 %>%
  transmute(
    units = "g/day",
    iso3, country, sex, age_group,
    population,
    region,
    
    calorie_requirement = mder_stratum_kcal,
    
    best_dist_calorie, cv_calorie,
    best_dist_protein, cv_protein,
    
    protein_share = protein_kcal_share_mean,
    prop_asf,
    quality_factor = quality_factor_EAR,
    
    protein_requirement      = ear_mean_g_day,
    protein_requirement_opt  = opt_mean_g_day,
    protein_requirement_adj  = ear_mean_g_day_adj,
    
    kcal_mean_low, kcal_mean_medium, kcal_mean_high
  ) %>%
  pivot_longer(
    cols = c(kcal_mean_low, kcal_mean_medium, kcal_mean_high),
    names_to = "calorie_scenario",
    values_to = "mean_calorie_intake"
  ) %>%
  mutate(
    calorie_scenario = recode(
      calorie_scenario,
      kcal_mean_low    = "low",
      kcal_mean_medium = "medium",
      kcal_mean_high   = "high"
    )
  ) %>%
  rowwise() %>%
  mutate(
    mean_protein_intake = (mean_calorie_intake * protein_share) / 4,
    
    sev_calorie = calculate_inadequacy(
      mean_intake       = mean_calorie_intake,
      cv_intake         = cv_calorie,
      distribution_type = best_dist_calorie,
      requirement       = calorie_requirement
    ),
    
    sev = calculate_inadequacy(
      mean_intake       = mean_protein_intake,
      cv_intake         = cv_protein,
      distribution_type = best_dist_protein,
      requirement       = protein_requirement
    ),
    
    sev_opt = calculate_inadequacy(
      mean_intake       = mean_protein_intake,
      cv_intake         = cv_protein,
      distribution_type = best_dist_protein,
      requirement       = protein_requirement_opt
    ),
    
    sev_adj = calculate_inadequacy(
      mean_intake       = mean_protein_intake,
      cv_intake         = cv_protein,
      distribution_type = best_dist_protein,
      requirement       = protein_requirement_adj
    ),
    
    ndeficient     = sev     * population,
    ndeficient_opt = sev_opt * population,
    ndeficient_adj = sev_adj * population
  ) %>%
  ungroup() %>%
  select(
    units,
    iso3, country, sex, age_group, population, region,
    calorie_scenario,
    mean_calorie_intake,
    calorie_requirement,
    sev_calorie,
    mean_protein_intake,
    protein_share,
    protein_requirement,
    protein_requirement_opt,
    protein_requirement_adj,
    best_dist_protein, cv_protein,
    best_dist_calorie, cv_calorie,
    sev, sev_opt, sev_adj,
    ndeficient, ndeficient_opt, ndeficient_adj,
    prop_asf,
    quality_factor
  ) %>%
  arrange(iso3, sex, age_group, factor(calorie_scenario, levels = c("low","medium","high")))




# ============================================================
# REDUCED dataset (Dataverse "Reduced"): medium-only, stored EAR
# ============================================================
protein_reduced <- dat2 %>%
  transmute(
    iso3,
    country,
    sex,
    age_group,
    population,
    mean_protein_intake = protein_grams_true,
    protein_requirement = ear_mean_g_day,
    sev = prevalence_protein_inadequate_EAR,
    ndeficient = prevalence_protein_inadequate_EAR * population
  ) %>%
  arrange(iso3, sex, age_group)

# -----------------------------
# Write outputs
# -----------------------------
out_dir <- "./output/dataverse"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

write_csv(protein_full,    file.path(out_dir, "protein_full.csv"),    na = "")
write_csv(protein_reduced, file.path(out_dir, "protein_reduced.csv"), na = "")

cat("\n✅ Dataverse exports written:\n")
cat(" - ", normalizePath(file.path(out_dir, "protein_full.csv")), "\n")
cat(" - ", normalizePath(file.path(out_dir, "protein_reduced.csv")), "\n")

cat("\nFULL scenario counts:\n")
print(table(protein_full$calorie_scenario, useNA = "ifany"))

cat("\nFULL dimensions:\n")
print(dim(protein_full))

cat("\nREDUCED dimensions:\n")
print(dim(protein_reduced))
