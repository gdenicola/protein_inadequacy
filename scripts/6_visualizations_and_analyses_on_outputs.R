# ==============================================================
# Script 6 — Plots & Maps (color-blind friendly)
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
library(viridis)

# ---- Setup ----
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); setwd("..")
options(scipen = 999)


# Load main stratum-level dataset (with all inadequacy + extras)
# This is the output of script 5
dat <- readRDS("./output/dat_quality.rds") %>%
  distinct(iso3, sex, age_group, .keep_all = TRUE)


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
# Use dat (dat_quality) directly for exact-calorie ribbons
prot2_long <- dat %>%
  filter(age_group != "0-0.99") %>%
  transmute(
    iso3,
    sex,
    age_group = as.character(age_group),
    population,
    best_dist_protein,
    cv_protein,
    ear_mean_g_day,
    ear_mean_g_day_adj,
    opt_mean_g_day,
    protein_kcal_share_lower,
    protein_kcal_share_mean,
    protein_kcal_share_upper,
    exact_kcal = ifelse(!is.na(kcal_consumed_optimal),
                        kcal_consumed_optimal,
                        eer_kcal_marco_mean)
  ) %>%
  pivot_longer(
    starts_with("protein_kcal_share_"),
    names_to  = "prot_level",
    values_to = "protein_share"
  ) %>%
  mutate(
    prot_level = recode(
      prot_level,
      protein_kcal_share_lower = "lower",
      protein_kcal_share_mean  = "mean",
      protein_kcal_share_upper = "upper"
    )
  )


prot2_calc <- prot2_long %>%
  rowwise() %>%
  mutate(
    protein_grams_exact = (exact_kcal * protein_share) / 4,
    # Original EAR-based inadequacy
    prot_inad_EAR_exact = calculate_inadequacy(
      protein_grams_exact, cv_protein, best_dist_protein, ear_mean_g_day
    ),
    # OPT-based inadequacy
    prot_inad_OPT_exact = calculate_inadequacy(
      protein_grams_exact, cv_protein, best_dist_protein, opt_mean_g_day
    ),
    # NEW: quality-adjusted EAR-based inadequacy
    prot_inad_EAR_exact_Q = calculate_inadequacy(
      protein_grams_exact, cv_protein, best_dist_protein, ear_mean_g_day_adj
    )
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
                         prot_OPT = "Below optimal protein intake (OPT) — exact calories"))


exact_ribbon <- exact_long %>%
  group_by(age_group, metric) %>%
  summarise(ymin = min(value, na.rm = TRUE),
            ymax = max(value, na.rm = TRUE), .groups = "drop")

exact_main <- exact_long %>% filter(prot_level == "mean")

# ==============================================================
# PART B — Full 3×3 scenarios (load Script 4 now)
# ==============================================================

# ==============================================================
# PART B — Full 3×3 scenarios (from dat_quality)
# ==============================================================

dat_scen <- dat %>%
  distinct(iso3, sex, age_group, .keep_all = TRUE) %>%
  filter(age_group != "0-0.99")

prot_long <- dat_scen %>%
  select(iso3, sex, age_group, population,
         cv_protein, best_dist_protein,
         ear_mean_g_day, ear_mean_g_day_adj,  # <-- NEW
         opt_mean_g_day,
         cv_calorie, best_dist_calorie, mder_stratum_kcal,
         protein_kcal_share_lower, protein_kcal_share_mean, protein_kcal_share_upper) %>%
  pivot_longer(starts_with("protein_kcal_share_"),
               names_to = "prot_level", values_to = "protein_share") %>%
  mutate(prot_level = recode(prot_level,
                             protein_kcal_share_lower = "lower",
                             protein_kcal_share_mean  = "mean",
                             protein_kcal_share_upper = "upper"))


cal_long <- dat_scen %>%
  select(iso3, sex, age_group, kcal_mean_low, kcal_mean_medium, kcal_mean_high) %>%
  pivot_longer(
    cols = c(kcal_mean_low, kcal_mean_medium, kcal_mean_high),
    names_to = "cal_scenario",
    values_to = "kcal_mean"
  ) %>%
  mutate(cal_scenario = sub("^kcal_mean_", "", cal_scenario))

scen <- prot_long %>%
  inner_join(
    cal_long,
    by = c("iso3", "sex", "age_group"),
    relationship = "many-to-many"
  )


scen_calc <- scen %>%
  rowwise() %>%
  mutate(
    protein_grams = (kcal_mean * protein_share) / 4,
    # Original EAR inadequacy
    prot_inad_EAR = calculate_inadequacy(
      protein_grams, cv_protein, best_dist_protein, ear_mean_g_day
    ),
    # OPT inadequacy
    prot_inad_OPT = calculate_inadequacy(
      protein_grams, cv_protein, best_dist_protein, opt_mean_g_day
    ),
    # Calorie inadequacy
    cal_inad      = calculate_inadequacy(
      kcal_mean,     cv_calorie,  best_dist_calorie,  mder_stratum_kcal
    ),
    # NEW: quality-adjusted EAR inadequacy (scenario world)
    prot_inad_EAR_Q = calculate_inadequacy(
      protein_grams, cv_protein, best_dist_protein, ear_mean_g_day_adj
    )
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
                         prot_inad_OPT = "Below optimal protein intake (OPT)",
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
plot_exact_EAR <- plot_exact_with_ribbon(
  "Protein inadequacy (EAR) — exact calories",
  "Protein Inadequacy by Age — Exact Calories (EAR)"
)

plot_exact_OPT <- plot_exact_with_ribbon(
  "Below optimal protein intake (OPT) — exact calories",
  "Proportion Below Optimal Protein Intake by Age — Exact Calories (OPT)"
)

print(plot_exact_EAR); print(plot_exact_OPT)

plot_prot_EAR <- plot_age_with_ribbon(
  "Protein inadequacy (EAR)",
  "Global Protein Inadequacy by Age (EAR threshold)",
  "Central scenario (M–M) in black; other 8 scenario lines dashed grey; ribbon = min–max across scenarios"
)

plot_prot_OPT <- plot_age_with_ribbon(
  "Below optimal protein intake (OPT)",
  "Global Share Below Optimal Protein Intake by Age (OPT threshold)",
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
  mutate(metric = factor(
    metric,
    levels = c("prot_EAR","prot_OPT","cal_MD"),
    labels = c("Protein (EAR)",
               "Below optimal protein intake (OPT)",
               "Calories (MDER)")
  ))


plot_box_by_sex <- ggplot(country_sex_MM, aes(x = sex, y = prevalence, fill = sex)) +
  geom_boxplot(alpha = 0.9, outlier.alpha = 0.85, outlier.size = 1.8) +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  # Flipped palette: Females = vermilion, Males = blue
  scale_fill_manual(
    values = c("Females" = oi$vermil,
               "Males" = oi$blue),
    guide = "none"
  ) +
  labs(
    title = "Country-level Prevalence of Inadequacy by Sex (MM scenario)",
    x = "Sex",
    y = "Prevalence of inadequacy"
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())


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


# ---- Shared, tight world viewport (kills top/bottom whitespace) ----
lat_lim <- c(-58, 85)         # adjust if needed
lon_lim <- c(-180, 180)

common_coord_sf <- function() {
  coord_sf(
    xlim   = lon_lim,
    ylim   = lat_lim,
    expand = FALSE,           # critical: removes padding
    clip   = "on",
    datum  = NA
  )
}

tight_map_theme <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.text   = element_blank(),
      panel.grid  = element_blank(),
      plot.margin = margin(4, 6, 4, 6)  # small margins
    )
}

# Map-specific saver with wide aspect to avoid letterboxing
save_map <- function(plot, filename, width = 12, height = 6, dpi = 300) {
  out_dir <- file.path(getwd(), "plots")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  filepath <- file.path(out_dir, paste0(filename, ".png"))
  ggsave(filepath, plot = plot, width = width, height = height,
         dpi = dpi, limitsize = FALSE, units = "in", bg = "white")
  message("✅ Saved (map): ", filepath)
}

# ==============================================================
# NEW: Sex-gap map (Female − Male) in protein inadequacy (EAR), MM scenario
#      Uses SAME continuous palette as EAR inadequacy maps
# ==============================================================

# 1) Country-level protein inadequacy (EAR) by sex, MM scenario
country_sex_MM_EAR <- scen_calc %>%
  filter(scenario_label == "MM", prot_level == "mean") %>%   # medium calories + mean protein share
  group_by(iso3, sex) %>%
  summarise(
    prot_EAR = weighted.mean(prot_inad_EAR, w = population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(names_from = sex, values_from = prot_EAR) %>%
  mutate(
    sexgap_pp = (Females - Males) * 100,       # percentage points
    parent_iso3 = normalize_iso3(iso3)
  )

# 2) Join to world geometry (same parent/territory logic as other maps)
map_df_sexgap <- world %>%
  left_join(country_sex_MM_EAR, by = "parent_iso3")

# 3) Plot (same palette as EAR inadequacy continuous maps)
sexgap_map <- ggplot(map_df_sexgap) +
  geom_sf(aes(fill = sexgap_pp), color = NA) +
  scale_fill_gradientn(
    colours = c("#f0f0f0", "#fdbb84", "#e34a33"),
    limits  = c(0, 15),     # match the visual you showed (0–15 pp)
    oob     = scales::squish,
    breaks  = c(0, 5, 10, 15),
    labels  = function(x) paste0(x, "%"),
    #na.value = oi$greyNA
  ) +
  common_coord_sf() +
  labs(
    title = "Sex gap in protein inadequacy (EAR), medium-calorie scenario",
    fill  = "Female \u2212 Male\n(pp)"
  ) +
  tight_map_theme()

print(sexgap_map)

# Save (goes to /plots via your save_map helper)
save_map(sexgap_map, "map_sexgap_protein_EAR_MM", width = 12, height = 6, dpi = 300)

plot_cat6_map <- function(df, col, title_txt) {
  ggplot(df) +
    geom_sf(aes(fill = .data[[col]]), color = NA) +
    scale_fill_manual(values = cat6_palette, na.value = oi$greyNA, drop = FALSE) +
    common_coord_sf() +
    labs(title = title_txt, fill = "Prevalence") +
    tight_map_theme()
}


# ------------------------------------------------------------
# Continuous plasma maps (capped at a chosen max)
# ------------------------------------------------------------
# Continuous protein maps (no transform, non-white baseline)
# Continuous protein maps (no transform, non-white baseline)
plot_protein_map_cont <- function(df, value_col, title_txt, cap) {
  ggplot(df) +
    geom_sf(aes(fill = .data[[value_col]]), color = NA) +
    scale_fill_gradientn(
      colours = c("#f0f0f0", "#fdbb84", "#e34a33"),
      limits  = c(0, cap),
      oob     = squish,
      labels  = percent_format(accuracy = 1),
      na.value = "#8c8c8c"   # <-- darker NA grey
    ) +
    common_coord_sf() +
    labs(title = title_txt, fill = "Prevalence") +
    tight_map_theme()
}

# Continuous calorie maps (same palette, capped at 30%)
plot_calorie_map_cont <- function(df, value_col, title_txt, cap = 0.301) {
  ggplot(df) +
    geom_sf(aes(fill = .data[[value_col]]), color = NA) +
    scale_fill_gradientn(
      colours = c("#f0f0f0", "#fdbb84", "#e34a33"),
      limits  = c(0, cap),
      oob     = squish,
      labels  = percent_format(accuracy = 1),
      na.value = "#8c8c8c"   # <-- darker NA grey
    ) +
    common_coord_sf() +
    labs(title = title_txt, fill = "Prevalence") +
    tight_map_theme()
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

print(plot_protein_map_cont(
  map_df_exact, "prot_EAR_exact",
  "Protein Inadequacy (EAR) — Optimal-Calorie World",
  cap = 0.301
))

print(plot_protein_map_cont(
  map_df_exact, "prot_OPT_exact",
  "Proportion Below Optimal Protein Intake (OPT) — Optimal-Calorie World",
  cap = 0.701
))


print(plot_protein_map_cont(
  map_df, "prot_EAR",
  "Protein Inadequacy (EAR) — MM scenario",
  cap = 0.301
))

print(plot_protein_map_cont(
  map_df, "prot_OPT",
  "Proportion Below Optimal Protein Intake (OPT) — MM scenario",
  cap = 0.70
))

print(plot_calorie_map_cont(
  map_df, "cal_MD",
  "Calorie Inadequacy (MDER) — MM scenario",
  cap = 0.30
))



# ==============================================================
# Ratio map: Protein-to-Calorie Inadequacy (MM scenario)
# ==============================================================

# ---- Ratio map palette (soft diverging with yellow middle) ----
ratio_palette <- c(
  "Protein inadequacy lower"  = "#8CBFE6",  # lighter, less heavy blue
  "About the same"            = "#F1D77E",  # muted yellow
  "Protein inadequacy higher" = "#FCAE91"   # soft red (kept)
)

# ==============================================================
# Ratio map: Protein-to-Energy Inadequacy (MM scenario) — 3-class
#   Robust to near-zero energy inadequacy
#   Labels are interpretive; thresholds live in caption
# ==============================================================

ratio_map_df <- function(df, prot_col, threshold_label) {
  
  title_txt <- if (threshold_label == "OPT") {
    "Protein inadequacy (OPT) relative to energy inadequacy"
  } else {
    paste0("Protein inadequacy (", threshold_label, ") relative to energy inadequacy")
  }
  
  df2 <- df %>%
    mutate(
      # Avoid divide-by-zero / near-zero explosions:
      # if energy inadequacy < 0.5 percentage points, treat ratio as undefined
      energy_safe = ifelse(is.na(cal_MD) | cal_MD < 0.005, NA_real_, cal_MD),
      ratio = .data[[prot_col]] / energy_safe,
      
      ratio_cat = case_when(
        is.na(ratio)        ~ NA_character_,
        ratio < 0.95        ~ "Protein inadequacy lower",
        ratio > 1.05        ~ "Protein inadequacy higher",
        TRUE                ~ "About the same"
      ),
      ratio_cat = factor(
        ratio_cat,
        levels = c("Protein inadequacy lower", "About the same", "Protein inadequacy higher")
      )
    )
  
  ggplot(df2) +
    geom_sf(aes(fill = ratio_cat), color = NA) +
    scale_fill_manual(
      values   = ratio_palette,
      #na.value = oi$greyNA,
      drop     = FALSE
    ) +
    common_coord_sf() +
    labs(
      title = title_txt,
      fill  = NULL
    ) +
    tight_map_theme() +
    theme(
      legend.position   = "bottom",
      legend.direction  = "horizontal",
      legend.key.height = unit(0.45, "cm"),
      legend.key.width  = unit(1.6, "cm"),
      legend.text       = element_text(size = 11)
    )
}


# Generate and save both maps
ratio_map_EAR <- ratio_map_df(map_df, "prot_EAR", "EAR")
ratio_map_OPT <- ratio_map_df(map_df, "prot_OPT", "OPT")

print(ratio_map_EAR)
print(ratio_map_OPT)



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


# ---- Global population-weighted averages (all scenarios) — QUALITY-ADJUSTED EAR ----
global_scenarios_Q <- scen_calc %>%
  group_by(prot_level, cal_scenario) %>%
  summarise(
    prot_inad_EAR_Q = weighted.mean(prot_inad_EAR_Q, w = population, na.rm = TRUE),
    prot_inad_OPT   = weighted.mean(prot_inad_OPT,   w = population, na.rm = TRUE),
    cal_inad        = weighted.mean(cal_inad,        w = population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(prot_level, cal_scenario) %>%
  mutate(
    across(
      c(prot_inad_EAR_Q, prot_inad_OPT, cal_inad),
      ~ scales::percent(.x, accuracy = 0.1)
    )
  )

# ---- Global population-weighted averages (all scenarios) — QUALITY-ADJUSTED EAR ----
global_scenarios_Q <- scen_calc %>%
  group_by(prot_level, cal_scenario) %>%
  summarise(
    prot_inad_EAR_Q = weighted.mean(prot_inad_EAR_Q, w = population, na.rm = TRUE),
    prot_inad_OPT   = weighted.mean(prot_inad_OPT,   w = population, na.rm = TRUE),
    cal_inad        = weighted.mean(cal_inad,        w = population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(prot_level, cal_scenario) %>%
  mutate(
    across(
      c(prot_inad_EAR_Q, prot_inad_OPT, cal_inad),
      ~ scales::percent(.x, accuracy = 0.1)
    )
  )

cat("\n=== Global population-weighted averages (all scenarios) — QUALITY-ADJUSTED EAR ===\n")
print(global_scenarios_Q)

# ---- Helper to save ggplots to ../plots ----
save_plot <- function(plot, filename, width = 10, height = 6, dpi = 300) {
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

# ---- Continuous protein maps (EAR + OPT) for optimal-calorie world ----
save_map(
  plot_protein_map_cont(
    map_df_exact, "prot_EAR_exact",
    "Protein Inadequacy (EAR) — Optimal-Calorie World",
    cap = 0.30
  ),
  "map_protein_EAR_optcal"
)

save_map(
  plot_protein_map_cont(
    map_df_exact, "prot_OPT_exact",
    "Proportion Below Optimal Protein Intake (OPT) — Optimal-Calorie World",
    cap = 0.70
  ),
  "map_protein_OPT_optcal"
)

# ---- Continuous protein maps (EAR + OPT) for realistic MM ----
save_map(
  plot_protein_map_cont(
    map_df, "prot_EAR",
    "Protein Inadequacy (EAR) — MM scenario",
    cap = 0.30
  ),
  "map_protein_EAR_MM"
)

save_map(
  plot_protein_map_cont(
    map_df, "prot_OPT",
    "Proportion Below Optimal Protein Intake (OPT) — MM scenario",
    cap = 0.70
  ),
  "map_protein_OPT_MM"
)

# ---- Continuous calorie map ----
save_map(
  plot_calorie_map_cont(
    map_df, "cal_MD",
    "Calorie Inadequacy (MDER) — MM scenario",
    cap = 0.30
  ),
  "map_calorie_MM"
)



# ---- Save ratio maps (Protein-to-Calorie Inadequacy, MM scenario) ----
# Build maps once with short threshold labels so titles are clean.
ratio_map_EAR <- ratio_map_df(map_df, "prot_EAR", "EAR")
ratio_map_OPT <- ratio_map_df(map_df, "prot_OPT", "OPT")


# Save
save_map(ratio_map_EAR, "map_ratio_proteinEAR_to_calorie_MM")
save_map(ratio_map_OPT, "map_ratio_proteinOPT_to_calorie_MM")


# ==============================================================
# NEW: Maps of share of protein from ASF and from seafood
# ==============================================================

# Load ASF / seafood protein proportions (from Script 5)
protein_props <- readRDS("./output/protein_asf_props.rds")

# Keep 2018 only (to match the rest of Script 6) and normalize ISO codes
protein_props_parented <- protein_props %>%
  filter(year == 2018) %>%
  mutate(parent_iso3 = normalize_iso3(iso3))

# Join to world geometry (using the same parent_iso3 logic as other maps)
map_asf_sea <- world %>%
  left_join(protein_props_parented, by = "parent_iso3")

# Interpreted as "higher share" vs "lower share", not good vs bad.
asf_map <- ggplot(map_asf_sea) +
  geom_sf(aes(fill = prop_asf), color = NA) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,   # now low = bright, high = dark
    limits   = c(0, 1),
    labels   = scales::percent_format(accuracy = 1),
    na.value = oi$greyNA
  ) +
  common_coord_sf() +
  labs(
    title = "Share of Total Protein from Animal-Source Foods (2018)",
    fill  = "ASF protein\n(% of total)"
  ) +
  tight_map_theme()

#alternative single gradient scale (blue)
# scale_fill_gradientn(
#   colors = c("#EFF3FF", "#BDD7E7", "#6BAED6", "#2171B5", "#08306B"),
#   limits = c(0,1),
#   labels = scales::percent_format(accuracy = 1),
#   na.value = oi$greyNA
# )

#alt_red
# scale_fill_gradientn(
#   colors = c("#FEE5D9", "#FCBBA1", "#FC9272", "#CC1F1A", "#66000A"),
#   limits = c(0,1),
#   labels = scales::percent_format(accuracy = 1),
#   na.value = oi$greyNA
# )

seafood_map <- ggplot(map_asf_sea) +
  geom_sf(aes(fill = prop_seafood), color = NA) +
  scale_fill_viridis_c(
    option = "mako",
    direction = -1,   
    #limits   = c(0, 1),
    labels   = scales::percent_format(accuracy = 1),
    na.value = oi$greyNA
  ) +
  common_coord_sf() +
  labs(
    title = "Share of Total Protein from Seafood (2018)",
    fill  = "Seafood protein\n(% of total)"
  ) +
  tight_map_theme()

# Print to screen
print(asf_map)
print(seafood_map)

# Save to ../plots using the existing helper
save_map(asf_map,     "map_share_protein_ASF")
save_map(seafood_map, "map_share_protein_seafood")

# ==============================================================
# NEW: Quality-adjusted EAR maps (protein) — EAR only
# ==============================================================

# --- Optimal-calorie world (exact calories, quality-adjusted EAR) ---
country_exact_mean_Q <- prot2_calc %>%
  filter(prot_level == "mean") %>%
  group_by(iso3) %>%
  summarise(
    prot_EAR_exact_Q = weighted.mean(prot_inad_EAR_exact_Q,
                                     w = population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(parent_iso3 = normalize_iso3(iso3))

map_df_exact_Q <- world %>%
  left_join(country_exact_mean_Q, by = "parent_iso3")

# Print
print(
  plot_protein_map_cont(
    map_df_exact_Q, "prot_EAR_exact_Q",
    "Protein Inadequacy (EAR, quality-adjusted) — Optimal-Calorie World",
    cap = 0.30   # same cap as non-quality-adjusted EAR map
  )
)

# Save
save_map(
  plot_protein_map_cont(
    map_df_exact_Q, "prot_EAR_exact_Q",
    "Protein Inadequacy (EAR, quality-adjusted) — Optimal-Calorie World",
    cap = 0.30
  ),
  "map_protein_EAR_Q_optcal"
)

# --- Realistic MM world (scenario = MM, quality-adjusted EAR) ---
country_total_MM_Q <- scen_calc %>%
  filter(scenario_label == "MM") %>%
  group_by(iso3) %>%
  summarise(
    prot_EAR_Q = weighted.mean(prot_inad_EAR_Q,
                               w = population, na.rm = TRUE),
    .groups = "drop"
  )

country_total_MM_parented_Q <- country_total_MM_Q %>%
  mutate(parent_iso3 = normalize_iso3(iso3))

map_df_Q <- world %>%
  left_join(country_total_MM_parented_Q, by = "parent_iso3")

# Print
print(
  plot_protein_map_cont(
    map_df_Q, "prot_EAR_Q",
    "Protein Inadequacy (EAR, quality-adjusted) — MM scenario",
    cap = 0.30   # same cap as non-quality-adjusted EAR MM map
  )
)

# Save
save_map(
  plot_protein_map_cont(
    map_df_Q, "prot_EAR_Q",
    "Protein Inadequacy (EAR, quality-adjusted) — MM scenario",
    cap = 0.30
  ),
  "map_protein_EAR_Q_MM"
)

# ==============================================================
# Fig 3 — 2×2 patch with SMALL gutter (no trim, no crop)
# ==============================================================

library(magick)

# ---- file paths (edit only if your filenames differ) ----
fig_paths <- list(
  a = file.path(getwd(), "plots", "map_protein_EAR_MM.png"),
  b = file.path(getwd(), "plots", "map_ratio_proteinEAR_to_calorie_MM.png"),
  c = file.path(getwd(), "plots", "map_sexgap_protein_EAR_MM.png"),
  d = file.path(getwd(), "plots", "map_share_protein_ASF.png")
)

# ---- sanity check ----
missing <- names(fig_paths)[!file.exists(unlist(fig_paths))]
if (length(missing) > 0) {
  stop("Missing files: ", paste(missing, collapse = ", "),
       "\nFix fig_paths to match the PNGs in your plots/ folder.")
}

# ---- read images (NO trimming) ----
A <- image_read(fig_paths$a)
B <- image_read(fig_paths$b)
C <- image_read(fig_paths$c)
D <- image_read(fig_paths$d)

# ---- add small panel label INSIDE each map (top-right corner) ----
add_label <- function(img, lab) {
  image_annotate(
    img, lab,
    gravity  = "northeast",
    location = "+35+25",
    size     = 40,
    color    = "black"
  )
}

# A <- add_label(A, "a")
# B <- add_label(B, "b")
# C <- add_label(C, "c")
# D <- add_label(D, "d")

# ---- helpers to match sizes without cropping ----
match_height <- function(left, right) {
  hl <- image_info(left)$height
  hr <- image_info(right)$height
  h  <- min(hl, hr)
  list(
    image_scale(left,  paste0("x", h)),
    image_scale(right, paste0("x", h))
  )
}

match_width <- function(top, bottom) {
  wt <- image_info(top)$width
  wb <- image_info(bottom)$width
  w  <- min(wt, wb)
  list(
    image_scale(top,    paste0(w, "x")),
    image_scale(bottom, paste0(w, "x"))
  )
}

# ---- build a small white spacer (gutter) ----
# tweak gutter_px if you want: 20–60 is usually perfect
gutter_px <- 40
make_h_gutter <- function(height, width = gutter_px, color = "white") {
  image_blank(width = width, height = height, color = color)
}
make_v_gutter <- function(width, height = gutter_px, color = "white") {
  image_blank(width = width, height = height, color = color)
}

# ---- assemble top row with gutter ----
ab <- match_height(A, B); A2 <- ab[[1]]; B2 <- ab[[2]]
h_top <- image_info(A2)$height
gH_top <- make_h_gutter(height = h_top)

top_row <- image_append(c(A2, gH_top, B2), stack = FALSE)

# ---- assemble bottom row with gutter ----
cd <- match_height(C, D); C2 <- cd[[1]]; D2 <- cd[[2]]
h_bot <- image_info(C2)$height
gH_bot <- make_h_gutter(height = h_bot)

bot_row <- image_append(c(C2, gH_bot, D2), stack = FALSE)

# ---- match row widths, then stack with gutter ----
rows <- match_width(top_row, bot_row)
top_row2 <- rows[[1]]
bot_row2 <- rows[[2]]

w_all <- image_info(top_row2)$width
gV <- make_v_gutter(width = w_all)

fig3_2x2 <- image_append(c(top_row2, gV, bot_row2), stack = TRUE)

# ---- save ----
out_png <- file.path(getwd(), "plots", "Fig3_2x2_gutter.png")
image_write(fig3_2x2, path = out_png, format = "png")
message("✅ Saved: ", out_png)

cat("\n✅ Script 6 finished: optimal-calorie-world maps first; new neutral 3-color palette for difference maps.\n")

