# ==============================================================
# Script 5 — Plots & Maps (color-blind friendly)
#   PART A: Protein ribbon plots with "exact calories" (from Script 2)
#   PART B: Full 3×3 scenarios (from Script 4): ribbons, boxplots, maps
#   NEW:   Common y-axis; Robust joins; Prot vs Cal difference maps
#          + Optimal-calorie-world maps (EAR & OPT) with same categories
# ==============================================================

# ---- Libraries ----
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# ---- Setup ----
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); setwd("..")
options(scipen = 999)

# ---- Okabe–Ito palette (color-blind friendly) ----
oi <- list(
  black   = "#000000",
  orange  = "#E69F00",
  sky     = "#56B4E9",
  green   = "#009E73",
  yellow  = "#F0E442",
  blue    = "#0072B2",
  vermil  = "#D55E00",
  purple  = "#CC79A7",
  greyNA  = "#E5E5E5"   # for true NA on maps
)

# ---- Shared helper for inadequacy ----
calculate_inadequacy <- function(mean_intake, cv_intake, distribution_type, requirement) {
  if (is.na(requirement) || requirement <= 0 || is.na(mean_intake) || is.na(cv_intake) ||
      mean_intake <= 0 || cv_intake <= 0) return(NA_real_)
  if (distribution_type == "gamma") {
    shape_k <- 1/(cv_intake^2); if (is.infinite(shape_k)) return(NA_real_)
    scale_theta <- mean_intake * (cv_intake^2)
    pgamma(requirement, shape = shape_k, scale = scale_theta)
  } else if (distribution_type == "log-normal") {
    meanlog <- log(mean_intake) - 0.5*log(1 + cv_intake^2)
    sdlog   <- sqrt(log(1 + cv_intake^2)); if (is.na(sdlog) || sdlog <= 0) return(NA_real_)
    plnorm(requirement, meanlog = meanlog, sdlog = sdlog)
  } else NA_real_
}

age_order <- c("1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39",
               "40-44","45-49","50-54","55-59","60-64","65-69","70-74",
               "75-79","80-84","85-89","90-94","95-99")

# ==============================================================
# PART A — Ribbon plots for Protein assuming "exact calories"
#           (load Script 2 BEFORE Script 4)
# ==============================================================

prot2 <- readRDS("./output/script2_final_results.rds")

prot2_slim <- prot2 %>%
  select(any_of(c(
    "iso3","sex","age_group",
    "best_dist","cv",
    "ear_mean_g_day","opt_mean_g_day",
    "protein_kcal_share_lower","protein_kcal_share_mean","protein_kcal_share_upper",
    "kcal_consumed_optimal",
    "eer_kcal_marco_mean",
    "population"
  ))) %>%
  mutate(age_group = as.character(age_group)) %>%
  filter(age_group != "0-0.99") %>%
  rename(best_dist_protein = best_dist,
         cv_protein        = cv)

need_pop_from_s4 <- !"population" %in% names(prot2_slim) || all(is.na(prot2_slim$population))

prot2_long <- prot2_slim %>%
  mutate(exact_kcal = ifelse(!is.na(kcal_consumed_optimal), kcal_consumed_optimal, eer_kcal_marco_mean)) %>%
  select(iso3, sex, age_group, population, exact_kcal,
         best_dist_protein, cv_protein, ear_mean_g_day, opt_mean_g_day,
         protein_kcal_share_lower, protein_kcal_share_mean, protein_kcal_share_upper) %>%
  pivot_longer(starts_with("protein_kcal_share_"),
               names_to = "prot_level", values_to = "protein_share") %>%
  mutate(prot_level = recode(prot_level,
                             protein_kcal_share_lower = "lower",
                             protein_kcal_share_mean  = "mean",
                             protein_kcal_share_upper = "upper"))

if (need_pop_from_s4 || any(is.na(prot2_long$population))) {
  dat_pop <- readRDS("./output/final_analysis_with_all_inadequacy.rds") %>%
    distinct(iso3, sex, age_group, .keep_all = TRUE) %>%
    filter(age_group != "0-0.99") %>%
    select(iso3, sex, age_group, population)
  prot2_long <- prot2_long %>%
    select(-population) %>%
    left_join(dat_pop, by = c("iso3","sex","age_group"))
}

prot2_calc <- prot2_long %>%
  rowwise() %>%
  mutate(
    protein_grams_exact = (exact_kcal * protein_share) / 4,
    prot_inad_EAR_exact = calculate_inadequacy(
      protein_grams_exact, cv_protein, best_dist_protein, ear_mean_g_day),
    prot_inad_OPT_exact = calculate_inadequacy(
      protein_grams_exact, cv_protein, best_dist_protein, opt_mean_g_day)
  ) %>%
  ungroup()

global_by_age_exact <- prot2_calc %>%
  group_by(age_group, prot_level) %>%
  summarise(
    prot_EAR = weighted.mean(prot_inad_EAR_exact, w = population, na.rm = TRUE),
    prot_OPT = weighted.mean(prot_inad_OPT_exact, w = population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(age_group = factor(age_group, levels = age_order),
         prot_level = factor(prot_level, levels = c("lower","mean","upper")))

exact_long <- global_by_age_exact %>%
  pivot_longer(cols = c(prot_EAR, prot_OPT),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         prot_EAR = "Protein inadequacy (EAR) — exact calories",
                         prot_OPT = "Protein inadequacy (OPT) — exact calories"))

exact_ribbon <- exact_long %>%
  group_by(age_group, metric) %>%
  summarise(ymin = min(value, na.rm = TRUE),
            ymax = max(value, na.rm = TRUE), .groups = "drop")

exact_main <- exact_long %>% filter(prot_level == "mean")

# ==============================================================
# PART B — Full 3×3 scenarios (load Script 4 now)
# ==============================================================

dat <- readRDS("./output/final_analysis_with_all_inadequacy.rds") %>%
  distinct(iso3, sex, age_group, .keep_all = TRUE) %>%
  filter(age_group != "0-0.99")

prot_long <- dat %>%
  select(iso3, sex, age_group, population,
         cv_protein, best_dist_protein, ear_mean_g_day, opt_mean_g_day,
         cv_calorie, best_dist_calorie, mder_stratum_kcal,
         protein_kcal_share_lower, protein_kcal_share_mean, protein_kcal_share_upper) %>%
  pivot_longer(starts_with("protein_kcal_share_"),
               names_to = "prot_level", values_to = "protein_share") %>%
  mutate(prot_level = recode(prot_level,
                             protein_kcal_share_lower = "lower",
                             protein_kcal_share_mean  = "mean",
                             protein_kcal_share_upper = "upper"))

cal_long <- dat %>%
  select(iso3, sex, age_group, kcal_mean_low, kcal_mean_medium, kcal_mean_high) %>%
  pivot_longer(
    cols = c(kcal_mean_low, kcal_mean_medium, kcal_mean_high),
    names_to = "cal_scenario",
    values_to = "kcal_mean"
  ) %>%
  mutate(cal_scenario = sub("^kcal_mean_", "", cal_scenario))

scen <- prot_long %>%
  inner_join(cal_long, by = c("iso3","sex","age_group"),
             relationship = "many-to-many")

scen_calc <- scen %>%
  rowwise() %>%
  mutate(
    protein_grams = (kcal_mean * protein_share) / 4,
    prot_inad_EAR = calculate_inadequacy(protein_grams, cv_protein, best_dist_protein, ear_mean_g_day),
    prot_inad_OPT = calculate_inadequacy(protein_grams, cv_protein, best_dist_protein, opt_mean_g_day),
    cal_inad      = calculate_inadequacy(kcal_mean,     cv_calorie,  best_dist_calorie,  mder_stratum_kcal)
  ) %>%
  ungroup()

global_by_age_scen <- scen_calc %>%
  group_by(age_group, prot_level, cal_scenario) %>%
  summarise(
    prot_inad_EAR = weighted.mean(prot_inad_EAR, w = population, na.rm = TRUE),
    prot_inad_OPT = weighted.mean(prot_inad_OPT, w = population, na.rm = TRUE),
    cal_inad      = weighted.mean(cal_inad,      w = population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    age_group = factor(age_group, levels = age_order),
    prot_level = factor(prot_level, levels = c("lower","mean","upper")),
    cal_scenario = factor(cal_scenario, levels = c("low","medium","high"))
  )

plot_long <- global_by_age_scen %>%
  pivot_longer(cols = c(prot_inad_EAR, prot_inad_OPT, cal_inad),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         prot_inad_EAR = "Protein inadequacy (EAR)",
                         prot_inad_OPT = "Protein inadequacy (OPT)",
                         cal_inad      = "Calorie inadequacy"))

ribbon_df <- plot_long %>%
  group_by(age_group, metric) %>%
  summarise(ymin = min(value, na.rm = TRUE),
            ymax = max(value, na.rm = TRUE),
            .groups = "drop")

main_line <- plot_long %>%
  filter(prot_level == "mean", cal_scenario == "medium")

# ==============================================================
# COMMON Y-AXIS LIMIT FOR *ALL* RIBBON PLOTS (starts at zero)
# ==============================================================

ymax_exact <- max(exact_ribbon$ymax, na.rm = TRUE)
ymax_scen  <- max(ribbon_df$ymax,    na.rm = TRUE)
ymax_all   <- max(ymax_exact, ymax_scen, 0, na.rm = TRUE) * 1.05  # 5% headroom

# Helpers that **use the common y limit** and show all scenarios
plot_exact_with_ribbon <- function(metric_name, title_txt) {
  df_main   <- exact_main %>% filter(metric == metric_name)
  df_rib    <- exact_ribbon %>% filter(metric == metric_name)
  df_guides <- exact_long %>% filter(metric == metric_name, prot_level != "mean")
  ggplot() +
    geom_ribbon(data = df_rib, aes(x = age_group, ymin = ymin, ymax = ymax),
                fill = "grey85", alpha = 0.55) +
    geom_line(data = df_guides, aes(x = age_group, y = value, group = prot_level),
              linewidth = 0.7, alpha = 0.55, color = "grey35", linetype = 2) +
    geom_line(data = df_main, aes(x = age_group, y = value, group = 1),
              linewidth = 1.4, color = oi$black) +
    geom_point(data = df_main, aes(x = age_group, y = value),
               size = 2.2, color = oi$black) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, ymax_all),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = title_txt,
         subtitle = "Bold = mean protein share (exact calories); ribbon spans lower–upper shares",
         x = "Age group", y = "Prevalence") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank())
}

plot_age_with_ribbon <- function(metric_name, title_txt, subtitle_txt) {
  df_main   <- main_line %>% filter(metric == metric_name)
  df_rib    <- ribbon_df %>% filter(metric == metric_name)
  df_guides <- plot_long  %>% filter(metric == metric_name,
                                     !(prot_level == "mean" & cal_scenario == "medium")) %>%
    group_by(age_group, prot_level, cal_scenario) %>%
    summarise(value = first(value), .groups = "drop")
  ggplot() +
    geom_ribbon(data = df_rib,
                aes(x = age_group, ymin = ymin, ymax = ymax),
                fill = "grey85", alpha = 0.55) +
    geom_line(data = df_guides,
              aes(x = age_group, y = value,
                  group = interaction(prot_level, cal_scenario)),
              linewidth = 0.7, alpha = 0.55, color = "grey35", linetype = 2) +
    geom_line(data = df_main, aes(x = age_group, y = value, group = 1),
              linewidth = 1.4, color = oi$black) +
    geom_point(data = df_main, aes(x = age_group, y = value),
               size = 2.2, color = oi$black) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, ymax_all),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = title_txt, subtitle = subtitle_txt, x = "Age group", y = "Prevalence") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank())
}

# --- Print ribbon plots (common y-axis) ---
plot_exact_EAR <- plot_exact_with_ribbon("Protein inadequacy (EAR) — exact calories",
                                         "Protein Inadequacy by Age — Exact Calories (EAR)")
plot_exact_OPT <- plot_exact_with_ribbon("Protein inadequacy (OPT) — exact calories",
                                         "Protein Inadequacy by Age — Exact Calories (OPT)")
print(plot_exact_EAR); print(plot_exact_OPT)

plot_prot_EAR <- plot_age_with_ribbon(
  "Protein inadequacy (EAR)",
  "Global Protein Inadequacy by Age (EAR threshold)",
  "Central scenario (M–M) in black; other 8 scenario lines dashed grey; ribbon = min–max across scenarios"
)
plot_prot_OPT <- plot_age_with_ribbon(
  "Protein inadequacy (OPT)",
  "Global Protein Inadequacy by Age (OPT threshold)",
  "Central scenario (M–M) in black; other 8 scenario lines dashed grey; ribbon = min–max across scenarios"
)
plot_cal <- plot_age_with_ribbon(
  "Calorie inadequacy",
  "Global Calorie Inadequacy by Age (MDER)",
  "Central scenario (M–M) in black; other 8 scenario lines dashed grey; ribbon = min–max across scenarios"
)
print(plot_prot_EAR); print(plot_prot_OPT); print(plot_cal)

# ---- Boxplots by sex (MM scenario), color-blind friendly ----
prot_map <- c(lower="L", mean="M", upper="U")
cal_map  <- c(low="L",   medium="M", high="H")

scen_calc <- scen_calc %>%
  mutate(scenario_label = paste0(prot_map[prot_level], cal_map[cal_scenario]))

country_sex_MM <- scen_calc %>%
  filter(scenario_label == "MM") %>%
  group_by(iso3, sex) %>%
  summarise(
    prot_EAR = weighted.mean(prot_inad_EAR, w = population, na.rm = TRUE),
    prot_OPT = weighted.mean(prot_inad_OPT, w = population, na.rm = TRUE),
    cal_MD   = weighted.mean(cal_inad,      w = population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(prot_EAR, prot_OPT, cal_MD),
               names_to = "metric", values_to = "prevalence") %>%
  mutate(metric = factor(metric,
                         levels = c("prot_EAR","prot_OPT","cal_MD"),
                         labels = c("Protein (EAR)","Protein (OPT)","Calories (MDER)")))

plot_box_by_sex <- ggplot(country_sex_MM, aes(x = sex, y = prevalence, fill = sex)) +
  geom_boxplot(alpha = 0.9, outlier.alpha = 0.85, outlier.size = 1.8) +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.02, 0.08))) +
  scale_fill_manual(values = c("Females" = oi$blue, "Males" = oi$vermil), guide = "none") +
  labs(title = "Country-level Prevalence of Inadequacy by Sex (MM scenario)",
       x = "Sex", y = "Prevalence of inadequacy") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())
print(plot_box_by_sex)

# ---- Country totals for maps (all sexes, MM) ----
country_total_MM <- scen_calc %>%
  filter(scenario_label == "MM") %>%
  group_by(iso3) %>%
  summarise(
    prot_EAR = weighted.mean(prot_inad_EAR, w = population, na.rm = TRUE),
    prot_OPT = weighted.mean(prot_inad_OPT, w = population, na.rm = TRUE),
    cal_MD   = weighted.mean(cal_inad,      w = population, na.rm = TRUE),
    .groups = "drop"
  )

# ---- World geometry + robust joins (FRA/NOR territories, SSD, etc.) ----
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  mutate(
    iso3_iso_a3     = if ("iso_a3"     %in% names(.)) iso_a3     else NA_character_,
    iso3_iso_a3_eh  = if ("iso_a3_eh"  %in% names(.)) iso_a3_eh  else NA_character_,
    iso3_adm0_a3    = if ("adm0_a3"    %in% names(.)) adm0_a3    else NA_character_,
    iso3_wb_a3      = if ("wb_a3"      %in% names(.)) wb_a3      else NA_character_
  )

normalize_iso3 <- function(code) {
  code <- dplyr::case_when(
    code %in% c("XKX","KOS") ~ "KOS",
    code == "SDS" ~ "SSD",
    code == "ROM" ~ "ROU",
    code == "ZAR" ~ "COD",
    code == "TMP" ~ "TLS",
    code == "WBG" ~ "PSE",
    TRUE ~ code
  )
  code
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

country_total_MM_parented <- country_total_MM %>%
  mutate(parent_iso3 = normalize_iso3(iso3))

map_df <- world %>% left_join(country_total_MM_parented, by = "parent_iso3")

# ---- Six-band categories for maps (0–5, 5–10, 10–15, 15–20, 20–25, 25%+) ----
cat6_cut <- function(x) {
  cut(x,
      breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, Inf),
      labels = c("0–5%", "5–10%", "10–15%", "15–20%", "20–25%", "25%+"),
      right = FALSE, include.lowest = TRUE)
}

map_df <- map_df %>%
  mutate(
    cat6_prot_EAR = cat6_cut(prot_EAR),
    cat6_prot_OPT = cat6_cut(prot_OPT),
    cat6_cal_MD   = cat6_cut(cal_MD)
  )

# Color-blind friendly sequential-ish palette (low→high), NA distinct grey
cat6_palette <- c(
  "0–5%"   = oi$green,
  "5–10%"  = "#7FC97F",
  "10–15%" = oi$yellow,
  "15–20%" = oi$orange,
  "20–25%" = oi$vermil,
  "25%+"   = "#A50F15"
)

plot_cat6_map <- function(df, col, title_txt) {
  ggplot(df) +
    geom_sf(aes(fill = .data[[col]]), color = NA) +
    scale_fill_manual(values = cat6_palette, na.value = oi$greyNA, drop = FALSE) +
    coord_sf(clip = "on") +
    labs(title = title_txt, fill = "Prevalence") +
    theme_minimal(base_size = 14) +
    theme(axis.text = element_blank(), panel.grid = element_blank())
}

# ==============================================================
# NEW: "Optimal-calorie world" maps (EAR & OPT protein) FIRST
# ==============================================================

country_exact_mean <- prot2_calc %>%
  filter(prot_level == "mean") %>%
  group_by(iso3) %>%
  summarise(
    prot_EAR_exact = weighted.mean(prot_inad_EAR_exact, w = population, na.rm = TRUE),
    prot_OPT_exact = weighted.mean(prot_inad_OPT_exact, w = population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(parent_iso3 = normalize_iso3(iso3))

map_df_exact <- world %>%
  left_join(country_exact_mean, by = "parent_iso3") %>%
  mutate(
    cat6_prot_EAR_exact = cat6_cut(prot_EAR_exact),
    cat6_prot_OPT_exact = cat6_cut(prot_OPT_exact)
  )

print(plot_cat6_map(map_df_exact, "cat6_prot_EAR_exact", "Protein Inadequacy (EAR) — Optimal-Calorie World"))
print(plot_cat6_map(map_df_exact, "cat6_prot_OPT_exact", "Protein Inadequacy (OPT) — Optimal-Calorie World"))

# --- Then the MM scenario maps (as before) ---
print(plot_cat6_map(map_df, "cat6_prot_EAR", "Protein Inadequacy (EAR) — MM scenario"))
print(plot_cat6_map(map_df, "cat6_prot_OPT", "Protein Inadequacy (OPT) — MM scenario"))
print(plot_cat6_map(map_df, "cat6_cal_MD",   "Calorie Inadequacy (MDER) — MM scenario"))

# ==============================================================
# Comparison maps: Protein vs Calories (clear labels)
#   - Similar: |Prot - Cal| ≤ 5 pp
#   - Prot. inad. > Cal. inad.: (protein higher)
#   - Prot. inad. < Cal. inad.: (calories higher)
#   Palette: three qualitative colors (no green/yellow/red; no hot/cold cue)
# ==============================================================

comp_levels <- c(
  "Prot. inad. ≈ Cal. inad. (±5 pp)",
  "Prot. inad. > Cal. inad. (≥5 pp)",
  "Prot. inad. < Cal. inad. (≥5 pp)"
)

# Harmonious purple–cyan–blue combo (color-blind friendly, non-hot/cold)
comp_palette <- c(
  "Prot. inad. ≈ Cal. inad. (±5 pp)" = "#D9D9D9",  # neutral light gray
  "Prot. inad. > Cal. inad. (≥5 pp)" = "#636363",  # dark slate
  "Prot. inad. < Cal. inad. (≥5 pp)" = "#BDB76B"   # muted olive
)





compare_map_df <- function(df, prot_col, title_txt) {
  df %>%
    mutate(
      diff_pp = (cal_MD - .data[[prot_col]]) * 100,
      comp_cat = case_when(
        is.na(diff_pp)                     ~ NA_character_,
        abs(diff_pp) <= 5                  ~ comp_levels[1],
        diff_pp < -5                       ~ comp_levels[2], # protein inadequacy > calorie inadequacy
        diff_pp >  5                       ~ comp_levels[3]  # protein inadequacy < calorie inadequacy
      ),
      comp_cat = factor(comp_cat, levels = comp_levels)
    ) %>%
    {
      ggplot(.) +
        geom_sf(aes(fill = comp_cat), color = NA) +
        scale_fill_manual(values = comp_palette, na.value = oi$greyNA, drop = FALSE) +
        coord_sf(clip = "on") +
        labs(title = title_txt, fill = NULL) +
        theme_minimal(base_size = 14) +
        theme(axis.text = element_blank(), panel.grid = element_blank())
    }
}

print(compare_map_df(map_df, "prot_EAR", "Protein (EAR) vs Calories — Difference Map (MM scenario)"))
print(compare_map_df(map_df, "prot_OPT", "Protein (OPT) vs Calories — Difference Map (MM scenario)"))

# ---- Optional global totals (sanity) ----
global_scenarios <- scen_calc %>%
  group_by(prot_level, cal_scenario) %>%
  summarise(
    prot_inad_EAR = weighted.mean(prot_inad_EAR, w = population, na.rm = TRUE),
    prot_inad_OPT = weighted.mean(prot_inad_OPT, w = population, na.rm = TRUE),
    cal_inad      = weighted.mean(cal_inad,      w = population, na.rm = TRUE),
    .groups = "drop"
  ) %>% arrange(prot_level, cal_scenario)

cat("\n=== Global population-weighted averages (all scenarios) ===\n")
print(
  global_scenarios %>%
    mutate(across(c(prot_inad_EAR, prot_inad_OPT, cal_inad),
                  ~ scales::percent(.x, accuracy = 0.1)))
)



# ---- Helper to save ggplots to ../plots ----
save_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  out_dir <- file.path(getwd(), "plots")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  filepath <- file.path(out_dir, paste0(filename, ".png"))
  ggsave(filepath, plot = plot, width = width, height = height, dpi = dpi)
  message("✅ Saved: ", filepath)
}


save_plot(plot_exact_EAR, "ribbon_exact_EAR")
save_plot(plot_exact_OPT, "ribbon_exact_OPT")

save_plot(plot_prot_EAR, "ribbon_prot_EAR")
save_plot(plot_prot_OPT, "ribbon_prot_OPT")
save_plot(plot_cal, "ribbon_calorie")

save_plot(plot_box_by_sex, "box_by_sex")

save_plot(plot_cat6_map(map_df_exact, "cat6_prot_EAR_exact",
                        "Protein Inadequacy (EAR) — Optimal-Calorie World"),
          "map_protein_EAR_optcal")
save_plot(plot_cat6_map(map_df_exact, "cat6_prot_OPT_exact",
                        "Protein Inadequacy (OPT) — Optimal-Calorie World"),
          "map_protein_OPT_optcal")

save_plot(plot_cat6_map(map_df, "cat6_prot_EAR",
                        "Protein Inadequacy (EAR) — MM scenario"),
          "map_protein_EAR_MM")
save_plot(plot_cat6_map(map_df, "cat6_prot_OPT",
                        "Protein Inadequacy (OPT) — MM scenario"),
          "map_protein_OPT_MM")
save_plot(plot_cat6_map(map_df, "cat6_cal_MD",
                        "Calorie Inadequacy (MDER) — MM scenario"),
          "map_calorie_MM")

save_plot(compare_map_df(map_df, "prot_EAR",
                         "Protein (EAR) vs Calories — Difference Map (MM scenario)"),
          "diffmap_EAR_vs_calories")
save_plot(compare_map_df(map_df, "prot_OPT",
                         "Protein (OPT) vs Calories — Difference Map (MM scenario)"),
          "diffmap_OPT_vs_calories")



cat("\n✅ Script 5 finished: optimal-calorie-world maps first; new neutral 3-color palette for difference maps.\n")
