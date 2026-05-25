# ==============================================================
# Clean supplementary table output script
# Protein inadequacy manuscript
#
# Produces Google-Doc-ready Supplementary Tables 1-5:
#   Supplementary Table 1 - EAR protein adequacy/inadequacy
#   Supplementary Table 2 - Optimal protein intake
#   Supplementary Table 3 - Quality-adjusted EAR inadequacy
#   Supplementary Table 4 - Total dietary protein and ASF share
#   Supplementary Table 5 - Energy-adequate counterfactual at EER
#
# Supplementary Table 6 is kept hardcoded in the manuscript.
#
# Inputs expected:
#   output/dat_quality.rds
#   output/protein_asf_props.rds
#
# Outputs written to:
#   output/supplementary_tables_clean/
# ============================================================== 

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(countrycode)
  library(readr)
  library(writexl)
})

rm(list = ls())
options(scipen = 999)

# ---- Project root and outputs ----
if (!file.exists(file.path("output", "dat_quality.rds")) &&
    requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  if (basename(script_dir) %in% c("scripts", "code", "R")) {
    setwd(dirname(script_dir))
  } else {
    setwd(script_dir)
  }
}

output_dir <- file.path("output", "supplementary_tables_clean")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================
# Helpers
# ============================================================== 

clean_sex <- function(x) {
  dplyr::case_when(
    x %in% c("Female", "female", "Females", "F") ~ "Female",
    x %in% c("Male", "male", "Males", "M") ~ "Male",
    TRUE ~ as.character(x)
  )
}

calculate_inadequacy <- function(mean_intake, cv_intake, distribution_type, requirement) {
  if (is.na(requirement) || requirement <= 0 ||
      is.na(mean_intake) || is.na(cv_intake) ||
      mean_intake <= 0 || cv_intake <= 0 ||
      is.na(distribution_type)) {
    return(NA_real_)
  }

  distribution_type <- as.character(distribution_type)

  if (distribution_type == "gamma") {
    shape_k <- 1 / (cv_intake^2)
    scale_theta <- mean_intake * (cv_intake^2)
    return(pgamma(requirement, shape = shape_k, scale = scale_theta))
  }

  if (distribution_type == "log-normal") {
    meanlog <- log(mean_intake) - 0.5 * log(1 + cv_intake^2)
    sdlog <- sqrt(log(1 + cv_intake^2))
    return(plnorm(requirement, meanlog = meanlog, sdlog = sdlog))
  }

  NA_real_
}

weighted_mean_safe <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  weighted.mean(x[ok], w = w[ok])
}

format_percent_string <- function(x, digits = 1) {
  paste0(round(100 * x, digits), "%")
}

order_world_region_country <- function(df, name_col = "Country / Region") {
  df %>%
    mutate(
      row_type_order = case_when(
        row_type == "world" ~ 1L,
        row_type == "region" ~ 2L,
        row_type == "country" ~ 3L,
        TRUE ~ 4L
      )
    ) %>%
    arrange(row_type_order, .data[[name_col]]) %>%
    select(-row_type_order, -row_type)
}

# ==============================================================
# Load analytical data
# ============================================================== 

dat <- readRDS(file.path("output", "dat_quality.rds")) %>%
  distinct(iso3, sex, age_group, .keep_all = TRUE) %>%
  filter(age_group != "0-0.99") %>%
  mutate(
    sex = clean_sex(sex),
    country = countrycode(iso3, origin = "iso3c", destination = "country.name"),
    country = if_else(is.na(country), iso3, country),
    region = countrycode(iso3, origin = "iso3c", destination = "region")
  )

# ==============================================================
# Reconstruct scenario-level protein intake and inadequacy
# ============================================================== 

protein_full <- dat %>%
  transmute(
    iso3,
    country,
    region,
    sex,
    age_group,
    population,

    protein_share = protein_kcal_share_mean,

    protein_requirement_EAR = ear_mean_g_day,
    protein_requirement_OPT = opt_mean_g_day,
    protein_requirement_EAR_QA = ear_mean_g_day_adj,

    best_dist_protein,
    cv_protein,

    kcal_mean_low,
    kcal_mean_medium,
    kcal_mean_high
  ) %>%
  pivot_longer(
    cols = c(kcal_mean_low, kcal_mean_medium, kcal_mean_high),
    names_to = "energy_scenario",
    values_to = "mean_energy_intake"
  ) %>%
  mutate(
    energy_scenario = recode(
      energy_scenario,
      kcal_mean_low = "Low",
      kcal_mean_medium = "Central",
      kcal_mean_high = "High"
    )
  ) %>%
  rowwise() %>%
  mutate(
    mean_protein_intake = (mean_energy_intake * protein_share) / 4,

    prev_EAR = calculate_inadequacy(
      mean_intake = mean_protein_intake,
      cv_intake = cv_protein,
      distribution_type = best_dist_protein,
      requirement = protein_requirement_EAR
    ),

    prev_OPT = calculate_inadequacy(
      mean_intake = mean_protein_intake,
      cv_intake = cv_protein,
      distribution_type = best_dist_protein,
      requirement = protein_requirement_OPT
    ),

    prev_EAR_QA = calculate_inadequacy(
      mean_intake = mean_protein_intake,
      cv_intake = cv_protein,
      distribution_type = best_dist_protein,
      requirement = protein_requirement_EAR_QA
    )
  ) %>%
  ungroup()

# ==============================================================
# Generic table builders
# ============================================================== 

build_requirement_table <- function(df, requirement_col, prevalence_col) {
  country_rows <- df %>%
    group_by(energy_scenario, iso3, country, region) %>%
    summarise(
      mean_pct = 100 * weighted_mean_safe(
        mean_protein_intake / .data[[requirement_col]],
        w = .data[["population"]]
      ),
      prevalence_pct = 100 * weighted_mean_safe(
        .data[[prevalence_col]],
        w = .data[["population"]]
      ),
      population_total = sum(.data[["population"]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    transmute(
      row_type = "country",
      `Country / Region` = country,
      energy_scenario,
      population_thousands = population_total / 1000,
      mean_pct,
      prevalence_pct
    )

  region_rows <- df %>%
    filter(!is.na(region)) %>%
    group_by(energy_scenario, region) %>%
    summarise(
      mean_pct = 100 * weighted_mean_safe(
        mean_protein_intake / .data[[requirement_col]],
        w = .data[["population"]]
      ),
      prevalence_pct = 100 * weighted_mean_safe(
        .data[[prevalence_col]],
        w = .data[["population"]]
      ),
      population_total = sum(.data[["population"]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    transmute(
      row_type = "region",
      `Country / Region` = region,
      energy_scenario,
      population_thousands = population_total / 1000,
      mean_pct,
      prevalence_pct
    )

  world_rows <- df %>%
    group_by(energy_scenario) %>%
    summarise(
      mean_pct = 100 * weighted_mean_safe(
        mean_protein_intake / .data[[requirement_col]],
        w = .data[["population"]]
      ),
      prevalence_pct = 100 * weighted_mean_safe(
        .data[[prevalence_col]],
        w = .data[["population"]]
      ),
      population_total = sum(.data[["population"]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    transmute(
      row_type = "world",
      `Country / Region` = "World",
      energy_scenario,
      population_thousands = population_total / 1000,
      mean_pct,
      prevalence_pct
    )

  bind_rows(world_rows, region_rows, country_rows) %>%
    mutate(energy_scenario = factor(energy_scenario, levels = c("Low", "Central", "High"))) %>%
    pivot_wider(
      id_cols = c(row_type, `Country / Region`),
      names_from = energy_scenario,
      values_from = c(population_thousands, mean_pct, prevalence_pct),
      names_glue = "{.value}_{energy_scenario}"
    ) %>%
    mutate(
      population_thousands = coalesce(
        population_thousands_Central,
        population_thousands_Low,
        population_thousands_High
      )
    )
}

build_prevalence_only_table <- function(df, prevalence_col) {
  country_rows <- df %>%
    group_by(energy_scenario, iso3, country, region) %>%
    summarise(
      prevalence_pct = 100 * weighted_mean_safe(
        .data[[prevalence_col]],
        w = .data[["population"]]
      ),
      population_total = sum(.data[["population"]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    transmute(
      row_type = "country",
      `Country / Region` = country,
      energy_scenario,
      population_thousands = population_total / 1000,
      prevalence_pct
    )

  region_rows <- df %>%
    filter(!is.na(region)) %>%
    group_by(energy_scenario, region) %>%
    summarise(
      prevalence_pct = 100 * weighted_mean_safe(
        .data[[prevalence_col]],
        w = .data[["population"]]
      ),
      population_total = sum(.data[["population"]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    transmute(
      row_type = "region",
      `Country / Region` = region,
      energy_scenario,
      population_thousands = population_total / 1000,
      prevalence_pct
    )

  world_rows <- df %>%
    group_by(energy_scenario) %>%
    summarise(
      prevalence_pct = 100 * weighted_mean_safe(
        .data[[prevalence_col]],
        w = .data[["population"]]
      ),
      population_total = sum(.data[["population"]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    transmute(
      row_type = "world",
      `Country / Region` = "World",
      energy_scenario,
      population_thousands = population_total / 1000,
      prevalence_pct
    )

  bind_rows(world_rows, region_rows, country_rows) %>%
    mutate(energy_scenario = factor(energy_scenario, levels = c("Low", "Central", "High"))) %>%
    pivot_wider(
      id_cols = c(row_type, `Country / Region`),
      names_from = energy_scenario,
      values_from = c(population_thousands, prevalence_pct),
      names_glue = "{.value}_{energy_scenario}"
    ) %>%
    mutate(
      population_thousands = coalesce(
        population_thousands_Central,
        population_thousands_Low,
        population_thousands_High
      )
    )
}

# ==============================================================
# Supplementary Table 1
# Protein adequacy under minimum EAR thresholds
# ============================================================== 

supp_table_1_EAR <- build_requirement_table(
  df = protein_full,
  requirement_col = "protein_requirement_EAR",
  prevalence_col = "prev_EAR"
) %>%
  transmute(
    row_type,
    `Country / Region`,
    `Population (thousands)` = round(population_thousands, 0),
    `Mean protein intake (% EAR, Low)` = round(mean_pct_Low, 2),
    `Mean protein intake (% EAR, Central)` = round(mean_pct_Central, 2),
    `Mean protein intake (% EAR, High)` = round(mean_pct_High, 2),
    `Protein inadequacy (% , Low)` = round(prevalence_pct_Low, 2),
    `Protein inadequacy (% , Central)` = round(prevalence_pct_Central, 2),
    `Protein inadequacy (% , High)` = round(prevalence_pct_High, 2)
  ) %>%
  order_world_region_country()

# ==============================================================
# Supplementary Table 2
# Protein intake relative to optimal requirement thresholds
# ============================================================== 

supp_table_2_OPT <- build_requirement_table(
  df = protein_full,
  requirement_col = "protein_requirement_OPT",
  prevalence_col = "prev_OPT"
) %>%
  transmute(
    row_type,
    `Country / Region`,
    `Population (thousands)` = round(population_thousands, 0),
    `Mean protein intake (% OPT, Low)` = round(mean_pct_Low, 2),
    `Mean protein intake (% OPT, Central)` = round(mean_pct_Central, 2),
    `Mean protein intake (% OPT, High)` = round(mean_pct_High, 2),
    `Below optimal intake (% , Low)` = round(prevalence_pct_Low, 2),
    `Below optimal intake (% , Central)` = round(prevalence_pct_Central, 2),
    `Below optimal intake (% , High)` = round(prevalence_pct_High, 2)
  ) %>%
  order_world_region_country()

# ==============================================================
# Supplementary Table 3
# Protein inadequacy under quality-adjusted EAR thresholds
# ============================================================== 

supp_table_3_QA <- build_prevalence_only_table(
  df = protein_full,
  prevalence_col = "prev_EAR_QA"
) %>%
  transmute(
    row_type,
    `Country / Region`,
    `Population (thousands)` = round(population_thousands, 0),
    `Protein inadequacy (% , Low, QA)` = round(prevalence_pct_Low, 2),
    `Protein inadequacy (% , Central, QA)` = round(prevalence_pct_Central, 2),
    `Protein inadequacy (% , High, QA)` = round(prevalence_pct_High, 2)
  ) %>%
  order_world_region_country()

# ==============================================================
# Supplementary Table 4
# Total dietary protein and share derived from animal-source foods
# ============================================================== 

protein_props <- readRDS(file.path("output", "protein_asf_props.rds")) %>%
  filter(year == 2018)

country_population <- dat %>%
  group_by(iso3) %>%
  summarise(population_2018 = sum(population, na.rm = TRUE), .groups = "drop")

protein_country <- protein_props %>%
  left_join(country_population, by = "iso3") %>%
  mutate(
    country = countrycode(iso3, origin = "iso3c", destination = "country.name"),
    country = if_else(is.na(country), iso3, country),
    region = countrycode(iso3, origin = "iso3c", destination = "region")
  ) %>%
  filter(!is.na(population_2018), population_2018 > 0)

asf_world <- protein_country %>%
  summarise(
    total_protein = weighted_mean_safe(total_protein, w = population_2018),
    prop_asf = weighted_mean_safe(prop_asf, w = population_2018),
    .groups = "drop"
  ) %>%
  mutate(
    row_type = "world",
    `Country / Region` = "World"
  )

asf_region <- protein_country %>%
  filter(!is.na(region)) %>%
  group_by(region) %>%
  summarise(
    total_protein = weighted_mean_safe(total_protein, w = population_2018),
    prop_asf = weighted_mean_safe(prop_asf, w = population_2018),
    .groups = "drop"
  ) %>%
  mutate(
    row_type = "region",
    `Country / Region` = region
  ) %>%
  select(-region)

asf_country <- protein_country %>%
  transmute(
    row_type = "country",
    `Country / Region` = country,
    total_protein,
    prop_asf
  )

supp_table_4_ASF <- bind_rows(asf_world, asf_region, asf_country) %>%
  transmute(
    row_type,
    `Country / Region`,
    `Total dietary protein (g/day)` = round(total_protein, 1),
    `Protein from animal-source foods` = format_percent_string(prop_asf, digits = 1)
  ) %>%
  order_world_region_country()

# ==============================================================
# Supplementary Table 5
# Energy-adequate counterfactual: protein intake at EER
# ============================================================== 

counterfactual_full <- dat %>%
  transmute(
    iso3,
    country,
    region,
    sex,
    age_group,
    population,
    protein_share = protein_kcal_share_mean,
    eer_kcal = eer_kcal_marco_mean,
    protein_requirement_RDA = rda_mean_g_day
  ) %>%
  mutate(
    mean_protein_intake_eer = (eer_kcal * protein_share) / 4,
    mean_pct_RDA = 100 * mean_protein_intake_eer / protein_requirement_RDA
  )

summarise_counterfactual <- function(df, group_vars, row_type_value, name_var = NULL, name_value = NULL) {
  out <- df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      population_total = sum(population, na.rm = TRUE),
      mean_protein_g_day = weighted_mean_safe(
        mean_protein_intake_eer,
        w = population
      ),
      mean_pct_RDA = weighted_mean_safe(
        mean_pct_RDA,
        w = population
      ),
      .groups = "drop"
    ) %>%
    mutate(row_type = row_type_value)

  if (!is.null(name_var)) {
    out <- out %>%
      mutate(`Country / Region` = .data[[name_var]])
  } else {
    out <- out %>%
      mutate(`Country / Region` = name_value)
  }

  out %>%
    select(
      row_type,
      `Country / Region`,
      sex,
      population_total,
      mean_protein_g_day,
      mean_pct_RDA
    )
}

counterfactual_with_total_sex <- counterfactual_full %>%
  bind_rows(
    counterfactual_full %>% mutate(sex = "Total")
  ) %>%
  mutate(sex = factor(sex, levels = c("Total", "Female", "Male")))

counterfactual_world <- summarise_counterfactual(
  df = counterfactual_with_total_sex,
  group_vars = c("sex"),
  row_type_value = "world",
  name_value = "World"
)

counterfactual_region <- counterfactual_with_total_sex %>%
  filter(!is.na(region)) %>%
  summarise_counterfactual(
    group_vars = c("region", "sex"),
    row_type_value = "region",
    name_var = "region"
  )

counterfactual_country <- summarise_counterfactual(
  df = counterfactual_with_total_sex,
  group_vars = c("iso3", "country", "region", "sex"),
  row_type_value = "country",
  name_var = "country"
)

supp_table_5_EER_counterfactual <- bind_rows(
  counterfactual_world,
  counterfactual_region,
  counterfactual_country
) %>%
  pivot_wider(
    id_cols = c(row_type, `Country / Region`),
    names_from = sex,
    values_from = c(
      population_total,
      mean_protein_g_day,
      mean_pct_RDA
    ),
    names_glue = "{.value}_{sex}"
  ) %>%
  mutate(
    population_thousands = population_total_Total / 1000
  ) %>%
  transmute(
    row_type,
    `Country / Region`,
    `Population (thousands)` = round(population_thousands, 0),

    `Mean protein intake (g/day, Total)` = round(mean_protein_g_day_Total, 1),
    `Mean protein intake (g/day, Female)` = round(mean_protein_g_day_Female, 1),
    `Mean protein intake (g/day, Male)` = round(mean_protein_g_day_Male, 1),

    `Mean protein intake (% RDA, Total)` = round(mean_pct_RDA_Total, 1),
    `Mean protein intake (% RDA, Female)` = round(mean_pct_RDA_Female, 1),
    `Mean protein intake (% RDA, Male)` = round(mean_pct_RDA_Male, 1)
  ) %>%
  order_world_region_country()

# ==============================================================
# Export Google-Doc-ready files
# ============================================================== 

write_xlsx(
  list(
    `Supplementary Table 1` = supp_table_1_EAR,
    `Supplementary Table 2` = supp_table_2_OPT,
    `Supplementary Table 3` = supp_table_3_QA,
    `Supplementary Table 4` = supp_table_4_ASF,
    `Supplementary Table 5` = supp_table_5_EER_counterfactual
  ),
  file.path(output_dir, "supplementary_tables_1_to_5_GoogleDoc_ready.xlsx")
)

write_csv(
  supp_table_1_EAR,
  file.path(output_dir, "Supplementary_Table_1_EAR.csv"),
  na = ""
)

write_csv(
  supp_table_2_OPT,
  file.path(output_dir, "Supplementary_Table_2_OPT.csv"),
  na = ""
)

write_csv(
  supp_table_3_QA,
  file.path(output_dir, "Supplementary_Table_3_quality_adjusted_EAR.csv"),
  na = ""
)

write_csv(
  supp_table_4_ASF,
  file.path(output_dir, "Supplementary_Table_4_ASF_share.csv"),
  na = ""
)

write_csv(
  supp_table_5_EER_counterfactual,
  file.path(output_dir, "Supplementary_Table_5_EER_counterfactual.csv"),
  na = ""
)

cat(
  "Supplementary Tables 1-5 written to ",
  normalizePath(output_dir, mustWork = FALSE),
  "\n",
  sep = ""
)

# Suggested caption for the manuscript:
# Supplementary Table 5. Protein intake under the energy-adequate
# counterfactual by country, World Bank region, and globally (2018).
# Mean protein intake is estimated assuming each country-sex-age stratum
# consumes energy equal to its Estimated Energy Requirement (EER). Values
# are population-weighted across age-sex strata.
