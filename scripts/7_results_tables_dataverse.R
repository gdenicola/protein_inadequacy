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
library(writexl)

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); setwd("..")
options(scipen = 999)

# -----------------------------
# Load main dataset
# -----------------------------
dat <- readRDS("./output/dat_quality.rds") %>%
  distinct(iso3, sex, age_group, .keep_all = TRUE) %>%
  filter(age_group != "0-0.99")


age_levels <- c(
  "1–4","5–9","10–14","15–19","20–24","25–29","30–34","35–39",
  "40–44","45–49","50–54","55–59","60–64","65–69","70–74",
  "75–79","80–84","85–89","90–94","95–99"
)

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
    
    age_group = factor(
      format_age_endash(age_group),
      levels = age_levels,
      ordered = TRUE
    ),
    
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

# 1. CSV
write_csv(protein_full,    file.path(out_dir, "protein_full.csv"),    na = "")
write_csv(protein_reduced, file.path(out_dir, "protein_reduced.csv"), na = "")

# 2. RDS (R Data Serialization - best for reloading into R later)
saveRDS(protein_full,    file.path(out_dir, "protein_full.rds"))
saveRDS(protein_reduced, file.path(out_dir, "protein_reduced.rds"))

# 3. XLSX (Excel)
write_xlsx(protein_full,    file.path(out_dir, "protein_full.xlsx"))
write_xlsx(protein_reduced, file.path(out_dir, "protein_reduced.xlsx"))

cat("\n✅ Dataverse exports written:\n")
cat(" - CSVs:  ", normalizePath(out_dir), "(*.csv)\n")
cat(" - RDS:   ", normalizePath(out_dir), "(*.rds)\n")
cat(" - Excel: ", normalizePath(out_dir), "(*.xlsx)\n")




# ============================================================
# EXTRA EXPORT (v2): Country + WB-region + World rows (MEDIUM scenario)
#   - Country rows: one per iso3
#   - WB region rows: one per WB region (population-weighted)
#   - World row: global aggregate (population-weighted)
#   Columns:
#     iso3, country_name, population_2018, wb_region, world,
#     mean_protein_pct_EAR, pct_inadequate_protein, pct_inadequate_protein_quality_adj
# ============================================================

# Safety: out_dir exists (defined above in your script, but keep robust)
if (!exists("out_dir")) out_dir <- "./output/dataverse"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Base: medium scenario micro rows ----
pf_med <- protein_full %>%
  filter(calorie_scenario == "medium")

# ---- Country rows (185-ish) ----
country_rows <- pf_med %>%
  group_by(iso3, country, region) %>%
  summarise(
    population_2018 = sum(population, na.rm = TRUE),
    
    mean_protein_pct_EAR = 100 * sum(population * (mean_protein_intake / protein_requirement), na.rm = TRUE) /
      sum(population, na.rm = TRUE),
    
    pct_inadequate_protein = 100 * sum(ndeficient, na.rm = TRUE) / sum(population, na.rm = TRUE),
    
    pct_inadequate_protein_quality_adj = 100 * sum(ndeficient_adj, na.rm = TRUE) / sum(population, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  transmute(
    iso3 = iso3,
    country_name = country,
    population_2018 =  (round(population_2018)),
    wb_region = region,
    world = "World",
    mean_protein_pct_EAR = round(mean_protein_pct_EAR, 2),
    pct_inadequate_protein = round(pct_inadequate_protein, 2),
    pct_inadequate_protein_quality_adj = round(pct_inadequate_protein_quality_adj, 2)
  )

# ---- WB region rows (one per region) ----
region_rows <- pf_med %>%
  group_by(region) %>%
  summarise(
    population_2018 = sum(population, na.rm = TRUE),
    
    mean_protein_pct_EAR = 100 * sum(population * (mean_protein_intake / protein_requirement), na.rm = TRUE) /
      sum(population, na.rm = TRUE),
    
    pct_inadequate_protein = 100 * sum(ndeficient, na.rm = TRUE) / sum(population, na.rm = TRUE),
    
    pct_inadequate_protein_quality_adj = 100 * sum(ndeficient_adj, na.rm = TRUE) / sum(population, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  transmute(
    iso3 = NA_character_,
    country_name = region,          # region name goes in country_name field
    population_2018 =  (round(population_2018)),
    wb_region = region,             # keep region also here for consistency
    world = "World",
    mean_protein_pct_EAR = round(mean_protein_pct_EAR, 2),
    pct_inadequate_protein = round(pct_inadequate_protein, 2),
    pct_inadequate_protein_quality_adj = round(pct_inadequate_protein_quality_adj, 2)
  )

# ---- World row (global) ----
world_row <- pf_med %>%
  summarise(
    population_2018 = sum(population, na.rm = TRUE),
    
    mean_protein_pct_EAR = 100 * sum(population * (mean_protein_intake / protein_requirement), na.rm = TRUE) /
      sum(population, na.rm = TRUE),
    
    pct_inadequate_protein = 100 * sum(ndeficient, na.rm = TRUE) / sum(population, na.rm = TRUE),
    
    pct_inadequate_protein_quality_adj = 100 * sum(ndeficient_adj, na.rm = TRUE) / sum(population, na.rm = TRUE)
  ) %>%
  transmute(
    iso3 = NA_character_,
    country_name = "World",
    population_2018 =  (round(population_2018)),
    wb_region = "World",
    world = "World",
    mean_protein_pct_EAR = round(mean_protein_pct_EAR, 2),
    pct_inadequate_protein = round(pct_inadequate_protein, 2),
    pct_inadequate_protein_quality_adj = round(pct_inadequate_protein_quality_adj, 2)
  )

# ---- Combine: countries + regions + world ----
country_region_world_table_medium <- bind_rows(
  country_rows,
  region_rows,
  world_row
) %>%
  # Nice ordering: countries first (iso3 not NA), then regions, then world
  mutate(
    row_type = case_when(
      !is.na(iso3) ~ "country",
      country_name == "World" ~ "world",
      TRUE ~ "region"
    )
  ) %>%
  arrange(
    factor(row_type, levels = c("country", "region", "world")),
    wb_region,
    iso3,
    country_name
  ) %>%
  select(-row_type)

# -----------------------------
# Write outputs
# -----------------------------
write_csv(
  country_region_world_table_medium,
  file.path(out_dir, "protein_country_region_world_table_medium.csv"),
  na = ""
)

write_xlsx(
  country_region_world_table_medium,
  file.path(out_dir, "protein_country_region_world_table_medium.xlsx")
)

saveRDS(
  country_region_world_table_medium,
  file.path(out_dir, "protein_country_region_world_table_medium.rds")
)

cat("\n✅ Country + WB-region + World MEDIUM table written:\n")
cat(" - CSV:   ", file.path(normalizePath(out_dir), "protein_country_region_world_table_medium.csv"), "\n")
cat(" - Excel: ", file.path(normalizePath(out_dir), "protein_country_region_world_table_medium.xlsx"), "\n")
cat(" - RDS:   ", file.path(normalizePath(out_dir), "protein_country_region_world_table_medium.rds"), "\n")
