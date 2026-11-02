# ==============================================================
# Script 7: Clean figure output script for protein inadequacy manuscript
# Produces the main-text figures only:
#   Figure 1
#   Figure 2
#   Figure 3
# ==============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(countrycode)
  library(ggplot2)
  library(scales)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(viridis)
  library(patchwork)
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

output_dir <- file.path("output", "main_figures_clean")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

DPI <- 450

FIG1_WIDTH        <- 8.4
FIG1_HEIGHT       <- 2.9
FIG2_WIDTH        <- 7.2
FIG2_HEIGHT       <- 5.0
FIG3_WIDTH        <- 7.2
FIG3_HEIGHT       <- 4.3

BASE_SIZE   <- 8.5
AXIS_SIZE   <- 7.2
TITLE_SIZE  <- 8.6
TAG_SIZE    <- 9.2
LEGEND_SIZE <- 7.0

sex_pal <- c("Female" = "#CC79A7", "Male" = "#0072B2")
na_gray <- "#8c8c8c"
central_point <- "black"

# ---- Helpers ----
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
      mean_intake <= 0 || cv_intake <= 0) {
    return(NA_real_)
  }
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

save_png <- function(plot, filename, width, height) {
  ggsave(
    filename = file.path(output_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = DPI,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
}

fig_theme <- function(center_title = FALSE) {
  theme_classic(base_size = BASE_SIZE) +
    theme(
      text = element_text(family = "sans"),
      axis.title = element_text(size = AXIS_SIZE),
      axis.text = element_text(size = AXIS_SIZE),
      plot.title = element_text(
        size = TITLE_SIZE,
        face = "bold",
        hjust = ifelse(center_title, 0.5, 0),
        margin = margin(b = 3)
      ),
      plot.tag = element_text(size = TAG_SIZE, face = "plain"),
      plot.tag.position = c(0.01, 0.985),
      legend.title = element_text(size = LEGEND_SIZE),
      legend.text = element_text(size = LEGEND_SIZE),
      legend.position = "top",
      plot.margin = margin(4, 4, 4, 4)
    )
}

fig2_theme <- function() {
  fig_theme(center_title = TRUE) +
    theme(
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25),
      panel.grid.minor.y = element_blank(),
      axis.line = element_line(color = "grey35", linewidth = 0.3),
      axis.ticks = element_line(color = "grey35", linewidth = 0.25),
      axis.text.y = element_text(size = 6.3),
      axis.title.y = element_text(size = 6.7),
      plot.tag.position = "topleft"
    )
}

map_theme <- function() {
  theme_minimal(base_size = BASE_SIZE) +
    theme(
      text = element_text(family = "sans"),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = TITLE_SIZE, face = "bold", hjust = 0.5),
      plot.tag = element_text(size = TAG_SIZE, face = "plain"),
      plot.tag.position = c(0.01, 0.985),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = LEGEND_SIZE),
      legend.text = element_text(size = LEGEND_SIZE),
      legend.key.height = grid::unit(0.30, "cm"),
      legend.key.width = grid::unit(0.92, "cm"),
      plot.margin = margin(2, 2, 2, 2)
    )
}

percent_axis <- function(x) paste0(x, "%")

normalize_iso3 <- function(code) {
  dplyr::case_when(
    code %in% c("XKX", "KOS") ~ "KOS",
    code == "SDS" ~ "SSD",
    code == "ROM" ~ "ROU",
    code == "ZAR" ~ "COD",
    code == "TMP" ~ "TLS",
    code == "WBG" ~ "PSE",
    TRUE ~ code
  )
}

territory_to_parent <- c(
  "GUF" = "FRA", "GLP" = "FRA", "MTQ" = "FRA", "REU" = "FRA", "MYT" = "FRA",
  "NCL" = "FRA", "PYF" = "FRA", "WLF" = "FRA", "SPM" = "FRA", "BLM" = "FRA",
  "MAF" = "FRA", "ATF" = "FRA", "MFO" = "FRA",
  "SJM" = "NOR", "BVT" = "NOR",
  "GRL" = "DNK", "FRO" = "DNK",
  "ABW" = "NLD", "CUW" = "NLD", "BES" = "NLD", "SXM" = "NLD",
  "GGY" = "GBR", "JEY" = "GBR", "IMN" = "GBR",
  "PRI" = "USA", "VIR" = "USA", "GUM" = "USA", "MNP" = "USA", "ASM" = "USA",
  "HKG" = "CHN", "MAC" = "CHN"
)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  mutate(
    iso3_iso_a3 = if ("iso_a3" %in% names(.)) iso_a3 else NA_character_,
    iso3_iso_a3_eh = if ("iso_a3_eh" %in% names(.)) iso_a3_eh else NA_character_,
    iso3_adm0_a3 = if ("adm0_a3" %in% names(.)) adm0_a3 else NA_character_,
    iso3_wb_a3 = if ("wb_a3" %in% names(.)) wb_a3 else NA_character_,
    ne_iso_raw = coalesce(iso3_adm0_a3, iso3_iso_a3_eh, iso3_iso_a3, iso3_wb_a3),
    ne_iso_norm = normalize_iso3(ne_iso_raw),
    parent_iso3 = case_when(
      ne_iso_norm %in% names(territory_to_parent) ~ territory_to_parent[ne_iso_norm],
      TRUE ~ ne_iso_norm
    )
  )

common_coord_sf <- function() {
  coord_sf(
    xlim = c(-180, 180),
    ylim = c(-58, 85),
    expand = FALSE,
    clip = "on",
    datum = NA
  )
}

map_fill_cont_red <- function(cap, label_fun = percent_format(accuracy = 1)) {
  scale_fill_gradientn(
    colours = c("#f0f0f0", "#fdbb84", "#e34a33"),
    limits = c(0, cap),
    oob = scales::squish,
    labels = label_fun,
    na.value = na_gray,
    guide = guide_colorbar(
      title.position = "left",
      title.vjust = 0.8,
      title.hjust = 1,
      barwidth = grid::unit(2.2, "cm"),
      barheight = grid::unit(0.25, "cm")
    )
  )
}

make_summary_blocks <- function(df, value_col, scenario_label, metric_label) {
  bind_rows(
    df %>%
      summarise(prevalence = weighted.mean(.data[[value_col]], w = population, na.rm = TRUE)) %>%
      mutate(group_level = "global", sex = NA_character_, region_WB = NA_character_),
    df %>%
      group_by(sex) %>%
      summarise(prevalence = weighted.mean(.data[[value_col]], w = population, na.rm = TRUE), .groups = "drop") %>%
      mutate(group_level = "sex", region_WB = NA_character_),
    df %>%
      group_by(region_WB) %>%
      summarise(prevalence = weighted.mean(.data[[value_col]], w = population, na.rm = TRUE), .groups = "drop") %>%
      mutate(group_level = "region", sex = NA_character_),
    df %>%
      group_by(region_WB, sex) %>%
      summarise(prevalence = weighted.mean(.data[[value_col]], w = population, na.rm = TRUE), .groups = "drop") %>%
      mutate(group_level = "region_sex")
  ) %>%
    mutate(scenario = scenario_label, metric = metric_label)
}

make_age_blocks <- function(df, value_col, scenario_label, metric_label) {
  df %>%
    group_by(age_group) %>%
    summarise(prevalence = weighted.mean(.data[[value_col]], w = population, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      group_level = "age",
      sex = NA_character_,
      region_WB = NA_character_,
      scenario = scenario_label,
      metric = metric_label
    )
}

# ---- Load data and derive analysis tables ----
dat_raw <- readRDS(file.path("output", "dat_quality.rds")) %>%
  distinct(iso3, sex, age_group, .keep_all = TRUE) %>%
  filter(age_group != "0-0.99") %>%
  mutate(sex = clean_sex(sex))

region_lookup <- dat_raw %>%
  distinct(iso3) %>%
  mutate(region_WB = countrycode(iso3, "iso3c", "region"))

dat <- dat_raw %>%
  select(-any_of("region_WB")) %>%
  left_join(region_lookup, by = "iso3")

scen3_long <- dat %>%
  transmute(
    iso3, sex, age_group, population, region_WB,
    best_dist_protein, cv_protein,
    ear_mean_g_day, ear_mean_g_day_adj, opt_mean_g_day,
    protein_kcal_share_mean,
    best_dist_calorie, cv_calorie, mder_stratum_kcal,
    kcal_mean_low, kcal_mean_medium, kcal_mean_high
  ) %>%
  pivot_longer(
    cols = c(kcal_mean_low, kcal_mean_medium, kcal_mean_high),
    names_to = "energy_scenario",
    values_to = "kcal_mean"
  ) %>%
  mutate(
    energy_scenario = recode(
      energy_scenario,
      kcal_mean_low = "low",
      kcal_mean_medium = "central",
      kcal_mean_high = "high"
    )
  ) %>%
  rowwise() %>%
  mutate(
    protein_grams = (kcal_mean * protein_kcal_share_mean) / 4,
    prot_inad_EAR = calculate_inadequacy(protein_grams, cv_protein, best_dist_protein, ear_mean_g_day),
    prot_inad_OPT = calculate_inadequacy(protein_grams, cv_protein, best_dist_protein, opt_mean_g_day),
    prot_inad_EAR_QA = calculate_inadequacy(protein_grams, cv_protein, best_dist_protein, ear_mean_g_day_adj),
    energy_inad = calculate_inadequacy(kcal_mean, cv_calorie, best_dist_calorie, mder_stratum_kcal)
  ) %>%
  ungroup()

prot2_central <- dat %>%
  transmute(
    iso3, sex, age_group, population, region_WB,
    best_dist_protein, cv_protein,
    ear_mean_g_day, ear_mean_g_day_adj, rda_mean_g_day, opt_mean_g_day,
    protein_kcal_share_mean,
    exact_kcal = ifelse(!is.na(kcal_consumed_optimal), kcal_consumed_optimal, eer_kcal_marco_mean)
  ) %>%
  rowwise() %>%
  mutate(
    protein_grams_exact = (exact_kcal * protein_kcal_share_mean) / 4,
    prot_mean_pct_EAR_exact = 100 * protein_grams_exact / ear_mean_g_day,
    prot_mean_pct_RDA_exact = 100 * protein_grams_exact / rda_mean_g_day,
    prot_mean_pct_OPT_exact = 100 * protein_grams_exact / opt_mean_g_day,
    prot_mean_pct_EAR_exact_QA = 100 * protein_grams_exact / ear_mean_g_day_adj,
    prot_inad_EAR_exact = calculate_inadequacy(
      protein_grams_exact, cv_protein, best_dist_protein, ear_mean_g_day
    )
  ) %>%
  ungroup()

summary_prev_WB <- bind_rows(
  make_summary_blocks(filter(scen3_long, energy_scenario == "low"),
                      "prot_inad_EAR", "realistic_energy_low", "protein_inadequacy_EAR"),
  make_summary_blocks(filter(scen3_long, energy_scenario == "central"),
                      "prot_inad_EAR", "realistic_energy_central", "protein_inadequacy_EAR"),
  make_summary_blocks(filter(scen3_long, energy_scenario == "high"),
                      "prot_inad_EAR", "realistic_energy_high", "protein_inadequacy_EAR"),
  make_age_blocks(filter(scen3_long, energy_scenario == "low"),
                  "prot_inad_EAR", "realistic_energy_low", "protein_inadequacy_EAR"),
  make_age_blocks(filter(scen3_long, energy_scenario == "central"),
                  "prot_inad_EAR", "realistic_energy_central", "protein_inadequacy_EAR"),
  make_age_blocks(filter(scen3_long, energy_scenario == "high"),
                  "prot_inad_EAR", "realistic_energy_high", "protein_inadequacy_EAR")
) %>%
  mutate(prevalence_pct = 100 * prevalence) %>%
  select(group_level, region_WB, sex, age_group, metric, scenario, prevalence, prevalence_pct)

# ==============================================================
# Figure 1
# ==============================================================

fig1a_country <- prot2_central %>%
  group_by(iso3) %>%
  summarise(
    mean_protein_g_day = weighted.mean(protein_grams_exact, w = population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(parent_iso3 = normalize_iso3(iso3))

fig1a_map_df <- world %>% left_join(fig1a_country, by = "parent_iso3")

p_fig1a <- ggplot(fig1a_map_df) +
  geom_sf(aes(fill = mean_protein_g_day), color = NA) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    limits = c(0, max(fig1a_map_df$mean_protein_g_day, na.rm = TRUE)),
    oob = squish,
    na.value = na_gray,
    labels = label_number(accuracy = 1),
    guide = guide_colorbar(
      title.position = "left",
      title.vjust = 0.8,
      title.hjust = 1,
      barwidth = grid::unit(2.2, "cm"),
      barheight = grid::unit(0.25, "cm")
    )  ) +
  common_coord_sf() +
  labs(
    tag = "(a)",
    title = "Mean protein intake under energy-adequate assumption",
    fill = "g/day"
  ) +
  map_theme() +
  theme(
    plot.title = element_text(
      size = TITLE_SIZE,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = -12)
    ),
    legend.position = c(0.50, -0.2),
    legend.justification = c(0.5, 0),
    legend.background = element_rect(fill = scales::alpha("white", 0.85), color = NA),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin = margin(-2, 2, 10, 2)
  )

# Revision (R1.8): panel (b) is now expressed relative to the EAR (%EAR), for
# consistency with the EAR cut-point estimates in the rest of the paper.
# The global %RDA value is still reported in the text.
fig1b_region_sex <- prot2_central %>%
  group_by(iso3, region_WB, sex) %>%
  summarise(
    country_mean_pct_EAR = weighted.mean(prot_mean_pct_EAR_exact, w = population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(region_WB), !is.na(sex)) %>%
  mutate(sex = factor(clean_sex(sex), levels = c("Female", "Male")))

fig1b_plot_df <- bind_rows(
  fig1b_region_sex %>% mutate(region_WB = "Global"),
  fig1b_region_sex
)

region_order_fig1 <- fig1b_region_sex %>%
  group_by(region_WB) %>%
  summarise(region_median = median(country_mean_pct_EAR, na.rm = TRUE), .groups = "drop") %>%
  arrange(region_median) %>%
  pull(region_WB)

fig1b_plot_df <- fig1b_plot_df %>%
  mutate(region_WB = factor(region_WB, levels = c("Global", region_order_fig1)))

p_fig1b <- ggplot(fig1b_plot_df, aes(x = region_WB, y = country_mean_pct_EAR, fill = sex)) +
  geom_boxplot(
    position = position_dodge(width = 0.72),
    width = 0.58,
    # Revision: show countries beyond the whiskers (1.5 x IQR) as small open circles,
    # so the extreme cases named in the text (e.g. Ecuador) appear in the figure.
    outlier.shape = 1,
    outlier.size = 0.9,
    outlier.stroke = 0.35,
    outlier.colour = "grey30",
    linewidth = 0.38
  ) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey45", linewidth = 0.42) +
  scale_fill_manual(values = sex_pal, drop = FALSE) +
  coord_cartesian(ylim = c(0, 430)) +
  labs(
    tag = "(b)",
    title = "Country-level mean protein intake by region and sex",
    x = NULL,
    y = "Country mean protein intake (% of EAR)",
    fill = NULL
  ) +
  fig_theme(center_title = TRUE) +
  theme(
    axis.text.x = element_text(angle = 32, hjust = 1, vjust = 1, size = 6.7),
    legend.position = c(0.62, 0.96),
    plot.tag.position = c(0.055, 0.985),
    legend.justification = c(0.5, 1),
    legend.direction = "horizontal",
    legend.background = element_rect(
      fill = scales::alpha("white", 0.85),
      color = NA
    ),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )

fig1_ab <- (p_fig1a | p_fig1b) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(theme = theme(plot.margin = margin(2, 2, 2, 2)))

save_png(fig1_ab, "Figure1_ab_panels.png", FIG1_WIDTH, FIG1_HEIGHT)

# ==============================================================
# Figure 2
# ==============================================================

prep_range_df <- function(df, group_cols) {
  df %>%
    filter(
      metric == "protein_inadequacy_EAR",
      scenario %in% c("realistic_energy_low", "realistic_energy_central", "realistic_energy_high")
    ) %>%
    mutate(
      scenario_short = recode(
        scenario,
        realistic_energy_low = "low",
        realistic_energy_central = "central",
        realistic_energy_high = "high"
      )
    ) %>%
    group_by(across(all_of(group_cols)), scenario_short) %>%
    summarise(prevalence_pct = mean(prevalence_pct, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = scenario_short, values_from = prevalence_pct) %>%
    mutate(
      range_min = pmin(low, high, na.rm = TRUE),
      range_max = pmax(low, high, na.rm = TRUE)
    )
}

# Revision (R1.4, R2.m1, R4.4): Monte Carlo 95% uncertainty intervals for the
# central scenario (from Script 6), drawn as a light band behind the central
# estimate, with the low/high calorie scenarios as separate open triangles.
mc_file <- file.path("output", "revision", "mc_summary.csv")
mc_central <- if (file.exists(mc_file)) {
  utils::read.csv(mc_file, stringsAsFactors = FALSE) %>% filter(scenario == "central", metric == "EAR")
} else {
  message("Note: ", mc_file, " not found - run Script 6 first to add Monte Carlo intervals to Figure 2.")
  NULL
}
add_mc <- function(df, levels, key) {
  if (is.null(mc_central)) return(df %>% mutate(mc_lo = NA_real_, mc_hi = NA_real_))
  m <- mc_central %>%
    filter(group_level %in% levels) %>%
    transmute(.key = group, mc_lo = 100 * q025, mc_hi = 100 * q975)
  df %>%
    mutate(.key = as.character(.data[[key]])) %>%
    left_join(m, by = ".key") %>%
    select(-.key)
}

fig2_region <- summary_prev_WB %>%
  filter(group_level %in% c("global", "region")) %>%
  mutate(region_WB = if_else(group_level == "global", "Global", region_WB)) %>%
  filter(!is.na(region_WB)) %>%
  prep_range_df(group_cols = "region_WB") %>%
  add_mc(levels = c("global", "region"), key = "region_WB")

region_levels_fig2 <- fig2_region %>%
  filter(region_WB != "Global") %>%
  arrange(central) %>%
  pull(region_WB)

fig2_region <- fig2_region %>%
  mutate(region_WB = factor(region_WB, levels = c("Global", region_levels_fig2)))

fig2_sex <- summary_prev_WB %>%
  filter(group_level == "sex", !is.na(sex)) %>%
  mutate(sex = factor(clean_sex(sex), levels = c("Female", "Male"))) %>%
  prep_range_df(group_cols = "sex") %>%
  add_mc(levels = "sex", key = "sex")

fig2_age <- summary_prev_WB %>%
  filter(group_level == "age", !is.na(age_group)) %>%
  prep_range_df(group_cols = "age_group") %>%
  add_mc(levels = "age", key = "age_group")

age_levels <- fig2_age %>%
  distinct(age_group) %>%
  mutate(age_start = suppressWarnings(as.numeric(sub("^([0-9]+).*", "\\1", age_group)))) %>%
  arrange(age_start) %>%
  pull(age_group)

fig2_age <- fig2_age %>%
  mutate(age_group = factor(age_group, levels = age_levels))

# Revision (R1.4, R2.m1, R4.4): central scenario as a point (a, b) or line (c); Monte Carlo
# 95% uncertainty interval for the central scenario as a light band; low- and high-calorie
# scenarios as open triangles (up = low-calorie, down = high-calorie), joined by dashed lines
# in (c). The y-axis is fixed at 0-20% in all panels for comparability.
fig2_col_central <- "#1F2A44"
fig2_col_ui      <- "#B9C6D8"
fig2_col_scen    <- "grey45"
fig2_dot         <- 1.9
fig2_lab_central <- "Central scenario"
fig2_lab_ui      <- "95% uncertainty interval (central scenario)"
fig2_lab_scen    <- c(low = "Low-calorie scenario", high = "High-calorie scenario")
fig2_shape_scen  <- c("Low-calorie scenario" = 24, "High-calorie scenario" = 25)
fig2_y_label     <- "Protein inadequacy (EAR), % of population"

fig2_y_scale <- list(
  coord_cartesian(ylim = c(0, 20)),
  scale_y_continuous(
    breaks = c(0, 5, 10, 15, 20),
    labels = percent_axis,
    expand = expansion(mult = c(0, 0.02))
  )
)

fig2_legend <- function(show) {
  pos <- if (!show) {
    theme(legend.position = "none")
  } else if (utils::packageVersion("ggplot2") >= "3.5.0") {
    theme(legend.position = "inside", legend.position.inside = c(0.01, 0.99))
  } else {
    theme(legend.position = c(0.01, 0.99))
  }
  pos + theme(
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.text = element_text(size = 6.6),
    legend.key.height = unit(8, "pt"),
    legend.key.width = unit(12, "pt"),
    legend.spacing.y = unit(0, "pt"),
    legend.margin = margin(0, 0, 0, 0),
    legend.background = element_blank()
  )
}

fig2_point_plot <- function(df, xvar, title, legend = FALSE, y_label = fig2_y_label) {
  scen <- df %>%
    select(all_of(xvar), low, high) %>%
    pivot_longer(c(low, high), names_to = "scen", values_to = "prev") %>%
    mutate(scen_lab = unname(fig2_lab_scen[scen]))
  ggplot(df, aes(x = .data[[xvar]])) +
    geom_linerange(aes(ymin = mc_lo, ymax = mc_hi, colour = fig2_lab_ui),
                   linewidth = 3.4, key_glyph = "path", na.rm = TRUE) +
    geom_point(data = scen, aes(y = prev, shape = scen_lab),
               size = 1.05, stroke = 0.55, colour = fig2_col_scen, fill = "white", na.rm = TRUE) +
    geom_point(aes(y = central, colour = fig2_lab_central), size = fig2_dot) +
    scale_colour_manual(
      NULL,
      values = setNames(c(fig2_col_central, fig2_col_ui), c(fig2_lab_central, fig2_lab_ui)),
      breaks = c(fig2_lab_central, fig2_lab_ui),
      guide = guide_legend(order = 1, override.aes = list(
        linetype = c(0, 1), shape = c(16, NA), linewidth = c(0, 3), size = c(fig2_dot, 0)
      ))
    ) +
    scale_shape_manual(
      NULL,
      values = fig2_shape_scen,
      breaks = names(fig2_shape_scen),
      guide = guide_legend(order = 2, override.aes = list(size = 1.5, colour = fig2_col_scen, fill = "white"))
    ) +
    fig2_y_scale +
    labs(x = NULL, y = y_label, title = title) +
    fig2_theme() +
    fig2_legend(legend)
}

fig2_region_labels <- c(
  "North America" = "North\nAmerica",
  "Europe & Central Asia" = "Europe &\nCentral Asia",
  "South Asia" = "South\nAsia",
  "Middle East & North Africa" = "Middle East &\nNorth Africa",
  "East Asia & Pacific" = "East Asia &\nPacific",
  "Latin America & Caribbean" = "Latin America\n& Caribbean",
  "Sub-Saharan Africa" = "Sub-Saharan\nAfrica"
)

p_fig2_region <- fig2_point_plot(fig2_region, "region_WB", "Protein inadequacy by region", legend = TRUE) +
  scale_x_discrete(labels = function(x) ifelse(x %in% names(fig2_region_labels), fig2_region_labels[x], x)) +
  labs(tag = "(a)") +
  theme(axis.text.x = element_text(lineheight = 0.9))

p_fig2_sex <- fig2_point_plot(fig2_sex, "sex", "Protein inadequacy by sex") +
  labs(tag = "(b)") +
  theme(plot.title.position = "plot")

fig2_age_scen <- fig2_age %>%
  select(age_group, low, high) %>%
  pivot_longer(c(low, high), names_to = "scen", values_to = "prev") %>%
  mutate(shape = unname(fig2_shape_scen[fig2_lab_scen[scen]]))

p_fig2_age <- ggplot(fig2_age, aes(x = age_group, group = 1)) +
  geom_ribbon(aes(ymin = mc_lo, ymax = mc_hi), fill = fig2_col_ui, na.rm = TRUE) +
  geom_line(data = fig2_age_scen, aes(y = prev, group = scen),
            colour = fig2_col_scen, linewidth = 0.3, linetype = "22") +
  geom_point(data = fig2_age_scen, aes(y = prev, group = scen, shape = I(shape)),
             size = 0.95, stroke = 0.45, colour = fig2_col_scen, fill = "white") +
  geom_line(aes(y = central), colour = fig2_col_central, linewidth = 0.6) +
  geom_point(aes(y = central), colour = fig2_col_central, size = 1.3) +
  fig2_y_scale +
  scale_x_discrete(expand = expansion(add = 0.6)) +
  labs(x = NULL, y = NULL, title = "Protein inadequacy by age", tag = "(c)") +
  fig2_theme() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 34, hjust = 1, vjust = 1, size = 6.6)
  )

fig2_bottom <- (p_fig2_sex | p_fig2_age) + plot_layout(widths = c(0.45, 1.55))
fig2_abc <- p_fig2_region / fig2_bottom +
  plot_layout(heights = c(1.05, 1)) +
  plot_annotation(theme = theme(plot.margin = margin(3, 3, 3, 3)))

save_png(fig2_abc, "Figure2_abc_panels.png", FIG2_WIDTH, FIG2_HEIGHT)

# ==============================================================
# Figure 3
# ==============================================================

country_EAR_central <- scen3_long %>%
  filter(energy_scenario == "central") %>%
  group_by(iso3) %>%
  summarise(
    prot_EAR = weighted.mean(prot_inad_EAR, w = population, na.rm = TRUE),
    energy_MD = weighted.mean(energy_inad, w = population, na.rm = TRUE),
    protein_g_day = weighted.mean(protein_grams, w = population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(parent_iso3 = normalize_iso3(iso3))

map_central_df <- world %>% left_join(country_EAR_central, by = "parent_iso3")

p_map_protein_g <- ggplot(map_central_df) +
  geom_sf(aes(fill = protein_g_day), color = NA) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    limits = c(0, max(map_central_df$protein_g_day, na.rm = TRUE)),
    oob = squish,
    na.value = na_gray,
    labels = label_number(accuracy = 1),
    guide = guide_colorbar(
      title.position = "left",
      title.vjust = 0.8,
      title.hjust = 1,
      barwidth = grid::unit(2.2, "cm"),
      barheight = grid::unit(0.25, "cm")
    )  ) +
  common_coord_sf() +
  labs(tag = "(a)", title = "Mean protein intake", fill = "g/day") +
  map_theme()

p_map_EAR <- ggplot(map_central_df) +
  geom_sf(aes(fill = prot_EAR), color = NA) +
  map_fill_cont_red(cap = 0.30, label_fun = percent_format(accuracy = 1)) +
  common_coord_sf() +
  labs(tag = "(b)", title = "Protein inadequacy (EAR)", fill = "Prevalence") +
  map_theme()

# Revision (R2.7): panel (c) compares two separately estimated (marginal)
# prevalences; labels now say so. Overlap within individuals is not identified.
ratio_palette <- c(
  "Protein prevalence lower" = "#8CBFE6",
  "Similar" = "#F1D77E",
  "Protein prevalence higher" = "#FCAE91"
)

map_ratio_df <- map_central_df %>%
  mutate(
    energy_safe = ifelse(is.na(energy_MD) | energy_MD < 0.005, NA_real_, energy_MD),
    ratio = prot_EAR / energy_safe,
    ratio_cat = case_when(
      is.na(ratio) ~ NA_character_,
      ratio < 0.95 ~ "Protein prevalence lower",
      ratio > 1.05 ~ "Protein prevalence higher",
      TRUE ~ "Similar"
    ),
    ratio_cat = factor(ratio_cat, levels = c("Protein prevalence lower", "Similar", "Protein prevalence higher"))
  )

p_map_ratio <- ggplot(map_ratio_df) +
  geom_sf(aes(fill = ratio_cat), color = NA) +
  scale_fill_manual(
    values = ratio_palette,
    na.value = na_gray,
    drop = FALSE,
    na.translate = FALSE,
    guide = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  common_coord_sf() +
  labs(tag = "(c)", title = "Prevalence of protein vs energy inadequacy", fill = NULL) +
  map_theme() +
  theme(
    legend.text = element_text(size = 5.8),
    legend.key.width = grid::unit(0.34, "cm"),
    legend.spacing.x = grid::unit(0.06, "cm")
  )

protein_props <- readRDS(file.path("output", "protein_asf_props.rds")) %>%
  filter(year == 2018) %>%
  mutate(parent_iso3 = normalize_iso3(iso3))

country_quality <- protein_props %>%
  group_by(parent_iso3) %>%
  summarise(asf_share = mean(prop_asf, na.rm = TRUE), .groups = "drop")

map_quality_df <- world %>% left_join(country_quality, by = "parent_iso3")

p_map_quality <- ggplot(map_quality_df) +
  geom_sf(aes(fill = asf_share), color = NA) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    limits = c(0, 1),
    labels = percent_format(accuracy = 1),
    na.value = na_gray,
    guide = guide_colorbar(
      title.position = "left",
      title.vjust = 0.8,
      title.hjust = 1,
      barwidth = grid::unit(2.2, "cm"),
      barheight = grid::unit(0.25, "cm")
    )  ) +
  common_coord_sf() +
  labs(tag = "(d)", title = "Protein from animal-source foods", fill = "ASF protein share") +
  map_theme()

fig3_2x2 <- (p_map_protein_g | p_map_EAR) / (p_map_ratio | p_map_quality) +
  plot_layout(widths = c(1, 1), heights = c(1, 1)) +
  plot_annotation(theme = theme(plot.margin = margin(0, 2, 0, 2))) &
  theme(
    plot.margin = margin(0, 1, 0, 1),
    legend.margin = margin(-2, 0, -2, 0),
    legend.box.margin = margin(-3, 0, -3, 0)
  )

save_png(fig3_2x2, "Figure3_2x2_maps.png", FIG3_WIDTH, FIG3_HEIGHT)

cat("Main figures written to ", normalizePath(output_dir, mustWork = FALSE), "\n", sep = "")
