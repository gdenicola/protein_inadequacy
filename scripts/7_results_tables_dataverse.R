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
# EXTRA EXPORT (v4): Country + Region + World rows (MEDIUM scenario)
#   - Country rows: one per iso3
#   - Region rows: one per WB region (population-weighted)
#   - World row: global aggregate (population-weighted)
#   Columns:
#     iso3, country_name, population_2018_thousands,
#     mean_protein_pct_EAR,
#     pct_inadequate_protein,
#     pct_inadequate_protein_quality_adj
# ============================================================

# Safety: out_dir exists
if (!exists("out_dir")) out_dir <- "./output/dataverse"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Base: medium scenario micro rows ----
pf_med <- protein_full %>%
  filter(calorie_scenario == "medium")

# ---- Country rows ----
country_rows <- pf_med %>%
  group_by(iso3, country, region) %>%
  summarise(
    population_2018 = sum(population, na.rm = TRUE),
    
    mean_protein_pct_EAR =
      100 * sum(population * (mean_protein_intake / protein_requirement), na.rm = TRUE) /
      sum(population, na.rm = TRUE),
    
    pct_inadequate_protein =
      100 * sum(ndeficient, na.rm = TRUE) / sum(population, na.rm = TRUE),
    
    pct_inadequate_protein_quality_adj =
      100 * sum(ndeficient_adj, na.rm = TRUE) / sum(population, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  transmute(
    iso3 = iso3,
    country_name = country,
    population_2018_thousands = round(population_2018 / 1000, 0),
    mean_protein_pct_EAR = round(mean_protein_pct_EAR, 2),
    pct_inadequate_protein = round(pct_inadequate_protein, 2),
    pct_inadequate_protein_quality_adj = round(pct_inadequate_protein_quality_adj, 2)
  )

# ---- Region rows ----
region_rows <- pf_med %>%
  group_by(region) %>%
  summarise(
    population_2018 = sum(population, na.rm = TRUE),
    
    mean_protein_pct_EAR =
      100 * sum(population * (mean_protein_intake / protein_requirement), na.rm = TRUE) /
      sum(population, na.rm = TRUE),
    
    pct_inadequate_protein =
      100 * sum(ndeficient, na.rm = TRUE) / sum(population, na.rm = TRUE),
    
    pct_inadequate_protein_quality_adj =
      100 * sum(ndeficient_adj, na.rm = TRUE) / sum(population, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  transmute(
    iso3 = NA_character_,
    country_name = region,
    population_2018_thousands = round(population_2018 / 1000, 0),
    mean_protein_pct_EAR = round(mean_protein_pct_EAR, 2),
    pct_inadequate_protein = round(pct_inadequate_protein, 2),
    pct_inadequate_protein_quality_adj = round(pct_inadequate_protein_quality_adj, 2)
  )

# ---- World row ----
world_row <- pf_med %>%
  summarise(
    population_2018 = sum(population, na.rm = TRUE),
    
    mean_protein_pct_EAR =
      100 * sum(population * (mean_protein_intake / protein_requirement), na.rm = TRUE) /
      sum(population, na.rm = TRUE),
    
    pct_inadequate_protein =
      100 * sum(ndeficient, na.rm = TRUE) / sum(population, na.rm = TRUE),
    
    pct_inadequate_protein_quality_adj =
      100 * sum(ndeficient_adj, na.rm = TRUE) / sum(population, na.rm = TRUE)
  ) %>%
  transmute(
    iso3 = NA_character_,
    country_name = "World",
    population_2018_thousands = round(population_2018 / 1000, 0),
    mean_protein_pct_EAR = round(mean_protein_pct_EAR, 2),
    pct_inadequate_protein = round(pct_inadequate_protein, 2),
    pct_inadequate_protein_quality_adj = round(pct_inadequate_protein_quality_adj, 2)
  )

# ---- Combine + order: World -> Regions -> Countries ----
country_region_world_table_medium <- bind_rows(
  world_row,
  region_rows,
  country_rows
) %>%
  mutate(
    row_type = case_when(
      country_name == "World" ~ "world",
      is.na(iso3) ~ "region",
      TRUE ~ "country"
    )
  ) %>%
  arrange(
    factor(row_type, levels = c("world", "region", "country")),
    if_else(row_type == "region", country_name, NA_character_),
    if_else(row_type == "country", country_name, NA_character_)
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

cat("\n✅ Country + Region + World MEDIUM table written (clean version):\n")
cat(" - CSV:   ", file.path(normalizePath(out_dir), "protein_country_region_world_table_medium.csv"), "\n")
cat(" - Excel: ", file.path(normalizePath(out_dir), "protein_country_region_world_table_medium.xlsx"), "\n")
cat(" - RDS:   ", file.path(normalizePath(out_dir), "protein_country_region_world_table_medium.rds"), "\n")




# ============================================================
# EXTRA EXPORT: Country + Region + World rows (MEDIUM scenario) — OPTIMAL PROTEIN
#   Parallel to EAR table but using OPT requirement thresholds
#   - Country rows: one per iso3
#   - Region rows: one per WB region (population-weighted)
#   - World row: global aggregate (population-weighted)
#   Columns:
#     iso3, country_name, population_2018_thousands,
#     mean_protein_pct_OPT,
#     pct_inadequate_protein_OPT
# ============================================================

# Safety: out_dir exists
if (!exists("out_dir")) out_dir <- "./output/dataverse"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Base: medium scenario micro rows ----
pf_med <- protein_full %>%
  filter(calorie_scenario == "medium")

# ---- Country rows ----
country_rows_opt <- pf_med %>%
  group_by(iso3, country, region) %>%
  summarise(
    population_2018 = sum(population, na.rm = TRUE),
    
    mean_protein_pct_OPT =
      100 * sum(population * (mean_protein_intake / protein_requirement_opt), na.rm = TRUE) /
      sum(population, na.rm = TRUE),
    
    pct_inadequate_protein_OPT =
      100 * sum(ndeficient_opt, na.rm = TRUE) / sum(population, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  transmute(
    iso3 = iso3,
    country_name = country,
    population_2018_thousands = round(population_2018 / 1000, 0),
    mean_protein_pct_OPT = round(mean_protein_pct_OPT, 2),
    pct_inadequate_protein_OPT = round(pct_inadequate_protein_OPT, 2)
  )

# ---- Region rows ----
region_rows_opt <- pf_med %>%
  group_by(region) %>%
  summarise(
    population_2018 = sum(population, na.rm = TRUE),
    
    mean_protein_pct_OPT =
      100 * sum(population * (mean_protein_intake / protein_requirement_opt), na.rm = TRUE) /
      sum(population, na.rm = TRUE),
    
    pct_inadequate_protein_OPT =
      100 * sum(ndeficient_opt, na.rm = TRUE) / sum(population, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  transmute(
    iso3 = NA_character_,
    country_name = region,
    population_2018_thousands = round(population_2018 / 1000, 0),
    mean_protein_pct_OPT = round(mean_protein_pct_OPT, 2),
    pct_inadequate_protein_OPT = round(pct_inadequate_protein_OPT, 2)
  )

# ---- World row ----
world_row_opt <- pf_med %>%
  summarise(
    population_2018 = sum(population, na.rm = TRUE),
    
    mean_protein_pct_OPT =
      100 * sum(population * (mean_protein_intake / protein_requirement_opt), na.rm = TRUE) /
      sum(population, na.rm = TRUE),
    
    pct_inadequate_protein_OPT =
      100 * sum(ndeficient_opt, na.rm = TRUE) / sum(population, na.rm = TRUE)
  ) %>%
  transmute(
    iso3 = NA_character_,
    country_name = "World",
    population_2018_thousands = round(population_2018 / 1000, 0),
    mean_protein_pct_OPT = round(mean_protein_pct_OPT, 2),
    pct_inadequate_protein_OPT = round(pct_inadequate_protein_OPT, 2)
  )

# ---- Combine + order: World -> Regions -> Countries ----
country_region_world_table_medium_opt <- bind_rows(
  world_row_opt,
  region_rows_opt,
  country_rows_opt
) %>%
  mutate(
    row_type = case_when(
      country_name == "World" ~ "world",
      is.na(iso3) ~ "region",
      TRUE ~ "country"
    )
  ) %>%
  arrange(
    factor(row_type, levels = c("world", "region", "country")),
    if_else(row_type == "region", country_name, NA_character_),
    if_else(row_type == "country", country_name, NA_character_)
  ) %>%
  select(-row_type)

# -----------------------------
# Write outputs
# -----------------------------
write_csv(
  country_region_world_table_medium_opt,
  file.path(out_dir, "protein_country_region_world_table_medium_opt.csv"),
  na = ""
)

write_xlsx(
  country_region_world_table_medium_opt,
  file.path(out_dir, "protein_country_region_world_table_medium_opt.xlsx")
)

saveRDS(
  country_region_world_table_medium_opt,
  file.path(out_dir, "protein_country_region_world_table_medium_opt.rds")
)

cat("\n✅ Country + Region + World MEDIUM table written (OPTIMAL version):\n")
cat(" - CSV:   ", file.path(normalizePath(out_dir), "protein_country_region_world_table_medium_opt.csv"), "\n")
cat(" - Excel: ", file.path(normalizePath(out_dir), "protein_country_region_world_table_medium_opt.xlsx"), "\n")
cat(" - RDS:   ", file.path(normalizePath(out_dir), "protein_country_region_world_table_medium_opt.rds"), "\n")


#####exports for supplementary tables:

ear_table_doc <- country_region_world_table_medium %>%
  select(
    Country = country_name,
    `Population (thousands, 2018)` = population_2018_thousands,
    `Mean protein intake (% EAR)` = mean_protein_pct_EAR,
    `Protein inadequacy (%)` = pct_inadequate_protein,
    `Protein inadequacy, quality-adjusted (%)` = pct_inadequate_protein_quality_adj
  )



opt_table_doc <- country_region_world_table_medium_opt %>%
  select(
    Country = country_name,
    `Population (thousands, 2018)` = population_2018_thousands,
    `Mean protein intake (% OPT)` = mean_protein_pct_OPT,
    `Below optimal intake (%)` = pct_inadequate_protein_OPT
  )


write_xlsx(ear_table_doc, file.path(out_dir, "Supplementary_Table_1_EAR.xlsx"))
write_xlsx(opt_table_doc, file.path(out_dir, "Supplementary_Table_2_OPT.xlsx"))



# ============================================================
# EXTRA EXPORT (v3): Three macro tables (LOW/MED/HIGH)
#   Table 1 — EAR (low/med/high): mean %EAR + % inadequate (EAR)
#   Table 2 — OPT (low/med/high): mean %OPT + % below OPT
#   Table 3 — EAR quality-adjusted (low/med/high): % inadequate (EAR adj)
#
# Notes:
# - Aggregation is population-weighted across sex x age_group within each iso3/region/world
# - Population is reported in THOUSANDS (rounded integer)
# - Ordering: World, then WB regions, then countries (alphabetic)
# ============================================================

# Safety: out_dir exists (defined above in your script, but keep robust)
if (!exists("out_dir")) out_dir <- "./output/dataverse"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Helper: compute ordered rows after binding (world -> regions -> countries) ----
order_world_region_country <- function(df) {
  df %>%
    mutate(
      row_type = case_when(
        name == "World" ~ "world",
        is.na(iso3)     ~ "region",
        TRUE            ~ "country"
      )
    ) %>%
    arrange(
      factor(row_type, levels = c("world", "region", "country")),
      if_else(row_type == "region", name, NA_character_),
      if_else(row_type == "country", name, NA_character_)
    ) %>%
    select(-row_type)
}

# ---- Base: micro rows for all scenarios ----
pf_all <- protein_full %>%
  filter(calorie_scenario %in% c("low", "medium", "high"))

# ============================================================
# TABLE 1: EAR (low/med/high) — mean %EAR + % inadequate (EAR)
# ============================================================

# ---- Country rows ----
t1_country <- pf_all %>%
  group_by(calorie_scenario, iso3, country, region) %>%
  summarise(
    pop = sum(population, na.rm = TRUE),
    
    mean_pct_EAR = 100 * sum(population * (mean_protein_intake / protein_requirement), na.rm = TRUE) /
      sum(population, na.rm = TRUE),
    
    pct_inadequate_EAR = 100 * sum(ndeficient, na.rm = TRUE) / sum(population, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  transmute(
    iso3 = iso3,
    name = country,
    pop_thousands = as.integer(round(pop / 1000)),
    group = "country",
    calorie_scenario = calorie_scenario,
    mean_pct = round(mean_pct_EAR, 2),
    pct_inadequate = round(pct_inadequate_EAR, 2)
  )

# ---- Region rows ----
t1_region <- pf_all %>%
  group_by(calorie_scenario, region) %>%
  summarise(
    pop = sum(population, na.rm = TRUE),
    
    mean_pct_EAR = 100 * sum(population * (mean_protein_intake / protein_requirement), na.rm = TRUE) /
      sum(population, na.rm = TRUE),
    
    pct_inadequate_EAR = 100 * sum(ndeficient, na.rm = TRUE) / sum(population, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  transmute(
    iso3 = NA_character_,
    name = region,
    pop_thousands = as.integer(round(pop / 1000)),
    group = "region",
    calorie_scenario = calorie_scenario,
    mean_pct = round(mean_pct_EAR, 2),
    pct_inadequate = round(pct_inadequate_EAR, 2)
  )

# ---- World row ----
t1_world <- pf_all %>%
  group_by(calorie_scenario) %>%
  summarise(
    pop = sum(population, na.rm = TRUE),
    
    mean_pct_EAR = 100 * sum(population * (mean_protein_intake / protein_requirement), na.rm = TRUE) /
      sum(population, na.rm = TRUE),
    
    pct_inadequate_EAR = 100 * sum(ndeficient, na.rm = TRUE) / sum(population, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  transmute(
    iso3 = NA_character_,
    name = "World",
    pop_thousands = as.integer(round(pop / 1000)),
    group = "world",
    calorie_scenario = calorie_scenario,
    mean_pct = round(mean_pct_EAR, 2),
    pct_inadequate = round(pct_inadequate_EAR, 2)
  )

t1_long <- bind_rows(t1_country, t1_region, t1_world) %>%
  mutate(
    calorie_scenario = factor(calorie_scenario, levels = c("low","medium","high"))
  )

t1_wide <- t1_long %>%
  select(iso3, name, pop_thousands, group, calorie_scenario, mean_pct, pct_inadequate) %>%
  pivot_wider(
    names_from  = calorie_scenario,
    values_from = c(mean_pct, pct_inadequate),
    names_glue  = "{.value}_{calorie_scenario}"
  ) %>%
  order_world_region_country()

# ============================================================
# TABLE 2: OPT (low/med/high) — mean %OPT + % below OPT
# ============================================================

t2_country <- pf_all %>%
  group_by(calorie_scenario, iso3, country, region) %>%
  summarise(
    pop = sum(population, na.rm = TRUE),
    
    mean_pct_OPT = 100 * sum(population * (mean_protein_intake / protein_requirement_opt), na.rm = TRUE) /
      sum(population, na.rm = TRUE),
    
    pct_below_OPT = 100 * sum(ndeficient_opt, na.rm = TRUE) / sum(population, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  transmute(
    iso3 = iso3,
    name = country,
    pop_thousands = as.integer(round(pop / 1000)),
    group = "country",
    calorie_scenario = calorie_scenario,
    mean_pct = round(mean_pct_OPT, 2),
    pct_inadequate = round(pct_below_OPT, 2)
  )

t2_region <- pf_all %>%
  group_by(calorie_scenario, region) %>%
  summarise(
    pop = sum(population, na.rm = TRUE),
    
    mean_pct_OPT = 100 * sum(population * (mean_protein_intake / protein_requirement_opt), na.rm = TRUE) /
      sum(population, na.rm = TRUE),
    
    pct_below_OPT = 100 * sum(ndeficient_opt, na.rm = TRUE) / sum(population, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  transmute(
    iso3 = NA_character_,
    name = region,
    pop_thousands = as.integer(round(pop / 1000)),
    group = "region",
    calorie_scenario = calorie_scenario,
    mean_pct = round(mean_pct_OPT, 2),
    pct_inadequate = round(pct_below_OPT, 2)
  )

t2_world <- pf_all %>%
  group_by(calorie_scenario) %>%
  summarise(
    pop = sum(population, na.rm = TRUE),
    
    mean_pct_OPT = 100 * sum(population * (mean_protein_intake / protein_requirement_opt), na.rm = TRUE) /
      sum(population, na.rm = TRUE),
    
    pct_below_OPT = 100 * sum(ndeficient_opt, na.rm = TRUE) / sum(population, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  transmute(
    iso3 = NA_character_,
    name = "World",
    pop_thousands = as.integer(round(pop / 1000)),
    group = "world",
    calorie_scenario = calorie_scenario,
    mean_pct = round(mean_pct_OPT, 2),
    pct_inadequate = round(pct_below_OPT, 2)
  )

t2_long <- bind_rows(t2_country, t2_region, t2_world) %>%
  mutate(
    calorie_scenario = factor(calorie_scenario, levels = c("low","medium","high"))
  )

t2_wide <- t2_long %>%
  select(iso3, name, pop_thousands, group, calorie_scenario, mean_pct, pct_inadequate) %>%
  pivot_wider(
    names_from  = calorie_scenario,
    values_from = c(mean_pct, pct_inadequate),
    names_glue  = "{.value}_{calorie_scenario}"
  ) %>%
  order_world_region_country()

# ============================================================
# TABLE 3: EAR quality-adjusted (low/med/high)
#   Here we export: % inadequate (EAR adj) by scenario
#   (mean %EAR QA is ambiguous; prevalence is the key result)
# ============================================================

t3_country <- pf_all %>%
  group_by(calorie_scenario, iso3, country, region) %>%
  summarise(
    pop = sum(population, na.rm = TRUE),
    
    pct_inadequate_EAR_QA = 100 * sum(ndeficient_adj, na.rm = TRUE) / sum(population, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  transmute(
    iso3 = iso3,
    name = country,
    pop_thousands = as.integer(round(pop / 1000)),
    group = "country",
    calorie_scenario = calorie_scenario,
    pct_inadequate = round(pct_inadequate_EAR_QA, 2)
  )

t3_region <- pf_all %>%
  group_by(calorie_scenario, region) %>%
  summarise(
    pop = sum(population, na.rm = TRUE),
    
    pct_inadequate_EAR_QA = 100 * sum(ndeficient_adj, na.rm = TRUE) / sum(population, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  transmute(
    iso3 = NA_character_,
    name = region,
    pop_thousands = as.integer(round(pop / 1000)),
    group = "region",
    calorie_scenario = calorie_scenario,
    pct_inadequate = round(pct_inadequate_EAR_QA, 2)
  )

t3_world <- pf_all %>%
  group_by(calorie_scenario) %>%
  summarise(
    pop = sum(population, na.rm = TRUE),
    
    pct_inadequate_EAR_QA = 100 * sum(ndeficient_adj, na.rm = TRUE) / sum(population, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  transmute(
    iso3 = NA_character_,
    name = "World",
    pop_thousands = as.integer(round(pop / 1000)),
    group = "world",
    calorie_scenario = calorie_scenario,
    pct_inadequate = round(pct_inadequate_EAR_QA, 2)
  )

t3_long <- bind_rows(t3_country, t3_region, t3_world) %>%
  mutate(
    calorie_scenario = factor(calorie_scenario, levels = c("low","medium","high"))
  )

t3_wide <- t3_long %>%
  select(iso3, name, pop_thousands, group, calorie_scenario, pct_inadequate) %>%
  pivot_wider(
    names_from  = calorie_scenario,
    values_from = pct_inadequate,
    names_glue  = "pct_inadequate_QA_{calorie_scenario}"
  ) %>%
  order_world_region_country()

# -----------------------------
# Write outputs (CSV + XLSX + RDS)
# -----------------------------
write_csv(t1_wide, file.path(out_dir, "table1_EAR_low_med_high.csv"), na = "")
write_csv(t2_wide, file.path(out_dir, "table2_OPT_low_med_high.csv"), na = "")
write_csv(t3_wide, file.path(out_dir, "table3_EAR_quality_adj_low_med_high.csv"), na = "")

write_xlsx(list(
  table1_EAR_low_med_high = t1_wide,
  table2_OPT_low_med_high = t2_wide,
  table3_EAR_quality_adj_low_med_high = t3_wide
), file.path(out_dir, "supp_tables_EAR_OPT_QA_low_med_high.xlsx"))

saveRDS(t1_wide, file.path(out_dir, "table1_EAR_low_med_high.rds"))
saveRDS(t2_wide, file.path(out_dir, "table2_OPT_low_med_high.rds"))
saveRDS(t3_wide, file.path(out_dir, "table3_EAR_quality_adj_low_med_high.rds"))

cat("\n✅ Supplementary tables written (LOW/MED/HIGH):\n")
cat(" - Table 1 (EAR): ", file.path(normalizePath(out_dir), "table1_EAR_low_med_high.csv"), "\n")
cat(" - Table 2 (OPT): ", file.path(normalizePath(out_dir), "table2_OPT_low_med_high.csv"), "\n")
cat(" - Table 3 (EAR QA): ", file.path(normalizePath(out_dir), "table3_EAR_quality_adj_low_med_high.csv"), "\n")
cat(" - Combined Excel: ", file.path(normalizePath(out_dir), "supp_tables_EAR_OPT_QA_low_med_high.xlsx"), "\n")



# ============================================================
# GOOGLE DOC READY VERSIONS
# ============================================================

# ---- Table 1: EAR ----
table1_doc <- t1_wide %>%
  select(
    `Country / Region` = name,
    `Population (thousands, 2018)` = pop_thousands,
    `Mean protein intake (% EAR, Low)` = mean_pct_low,
    `Mean protein intake (% EAR, Medium)` = mean_pct_medium,
    `Mean protein intake (% EAR, High)` = mean_pct_high,
    `Protein inadequacy (% , Low)` = pct_inadequate_low,
    `Protein inadequacy (% , Medium)` = pct_inadequate_medium,
    `Protein inadequacy (% , High)` = pct_inadequate_high
  )

# ---- Table 2: OPT ----
table2_doc <- t2_wide %>%
  select(
    `Country / Region` = name,
    `Population (thousands, 2018)` = pop_thousands,
    `Mean protein intake (% OPT, Low)` = mean_pct_low,
    `Mean protein intake (% OPT, Medium)` = mean_pct_medium,
    `Mean protein intake (% OPT, High)` = mean_pct_high,
    `Below optimal intake (% , Low)` = pct_inadequate_low,
    `Below optimal intake (% , Medium)` = pct_inadequate_medium,
    `Below optimal intake (% , High)` = pct_inadequate_high
  )

# ---- Table 3: EAR Quality-Adjusted ----
table3_doc <- t3_wide %>%
  select(
    `Country / Region` = name,
    `Population (thousands, 2018)` = pop_thousands,
    `Protein inadequacy (% , Low, QA)` = pct_inadequate_QA_low,
    `Protein inadequacy (% , Medium, QA)` = pct_inadequate_QA_medium,
    `Protein inadequacy (% , High, QA)` = pct_inadequate_QA_high
  )

# Export clean Excel file with all 3 sheets
write_xlsx(
  list(
    `Table 1 - EAR` = table1_doc,
    `Table 2 - OPT` = table2_doc,
    `Table 3 - EAR QA` = table3_doc
  ),
  file.path(out_dir, "SI_tables_GoogleDoc_ready.xlsx")
)

cat("\n✅ Google Doc–ready SI tables exported:\n")
cat(" - ", file.path(normalizePath(out_dir), "SI_tables_GoogleDoc_ready.xlsx"), "\n")


#get proportions of protein sources from sheet data
library(dplyr)
library(readr)

# 1) Load ASF shares by country (from Script 5)
protein_props <- readRDS("./output/protein_asf_props.rds") %>%
  filter(year == 2018) %>%
  select(iso3, prop_asf)

# 2) Get 2018 country populations from your quality-adjusted stratum dataset
#    (robust: sum across sex/age strata to iso3 totals)
dat_quality <- readRDS("./output/dat_quality.rds")

pop_iso3_2018 <- dat_quality %>%
  group_by(iso3) %>%
  summarise(pop_2018 = sum(population, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(iso3), pop_2018 > 0)

# 3) Join and compute global population-weighted ASF share
global_asf <- pop_iso3_2018 %>%
  inner_join(protein_props, by = "iso3") %>%
  filter(!is.na(prop_asf)) %>%
  summarise(
    global_prop_asf = sum(pop_2018 * prop_asf) / sum(pop_2018),
    n_countries = n(),
    pop_covered = sum(pop_2018),
    .groups = "drop"
  ) %>%
  mutate(
    global_prop_plant = 1 - global_prop_asf,
    pct_animal = 100 * global_prop_asf,
    pct_plant  = 100 * global_prop_plant
  )

global_asf

