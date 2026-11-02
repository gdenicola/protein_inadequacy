# ==============================================================
# Script 6: Uncertainty and sensitivity analyses (Nature Food revision)
#
# Run after Script 5 and BEFORE Scripts 7 and 8
# (Figure 2 and Supplementary Table 1 read the Monte Carlo output).
#
# Input:  output/dat_quality.rds (+ logs and donor-matching inputs written by
#         Scripts 1, 3, 4 to output/revision/)
# Output: output/revision/  (csv files + revision_tables.xlsx)
#
#  A. Monte Carlo uncertainty analysis ........................ R1.4, R1.11, R2.m1, R4.4
#  B. Full probability approach vs EAR cut-point .............. R1.3
#  C. Distribution family (all gamma / all log-normal) ........ R1.12
#  D. Protein-quality adjustment, sourced parameter grid ...... R1.6, R2.3
#  E. Higher intake levels proposed for adults (exploratory) .. R1.5, R1.15, R2.6
#  F. nutriR matching tiers and donor countries ............... R1.9, R2.4, R4.3
#  G. Donor bootstrap (5 nearest donor countries) ............. R4.3, R4.6, R4.4
#  H. Protein-to-energy (P/E) ratios .......................... R1.17
#  I. Countries where protein exceeds energy inadequacy ....... R2.m4, R2.7
#  J. Imputation log .......................................... R2.4
#  K. Calibration check ....................................... R1.16, R2.1
# ==============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(countrycode)
})

rm(list = ls())
options(scipen = 999)

# ---- Project root ----
if (!file.exists(file.path("output", "dat_quality.rds")) &&
    requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  if (basename(script_dir) %in% c("scripts", "code", "R")) setwd(dirname(script_dir)) else setwd(script_dir)
}

out_dir <- file.path("output", "revision")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Settings ----
N_DRAWS      <- 10000   # Monte Carlo draws
CHUNK        <- 500     # draws processed at a time (memory)
SEED         <- 20261020
CV_REQ       <- 0.12    # CV of the protein requirement distribution (Rand et al. 2003; IOM)
Z95          <- qnorm(0.975)
SCENARIOS    <- c(low = "kcal_mean_low", central = "kcal_mean_medium", high = "kcal_mean_high")

# Protein-quality grid (R1.6): see Script 5 for sources
QF_GRID <- tidyr::expand_grid(
  qf_max    = c(0.11, 0.15, 0.25),
  asf_upper = c(0.33, 0.40, 0.45),
  asf_lower = c(0.15, 0.20)
)

# ==============================================================
# Helpers
# ==============================================================

# P(intake < requirement) for gamma / log-normal intake with given mean and CV.
# Vectorised; req and mu can be vectors or matrices (strata in rows), cv and
# dist are per-stratum vectors and are recycled down the columns.
p_below <- function(req, mu, cv, dist) {
  dims <- dim(mu)
  req <- as.vector(req); mu <- as.vector(mu)
  n <- length(mu)
  cv <- rep_len(cv, n); dist <- rep_len(as.character(dist), n)
  out <- rep(NA_real_, n)
  ok <- !is.na(req) & !is.na(mu) & !is.na(cv) & req > 0 & mu > 0 & cv > 0
  g <- ok & dist == "gamma"
  l <- ok & dist == "log-normal"
  out[g] <- pgamma(req[g], shape = 1 / cv[g]^2, scale = mu[g] * cv[g]^2)
  out[l] <- plnorm(req[l], meanlog = log(mu[l]) - 0.5 * log(1 + cv[l]^2),
                   sdlog = sqrt(log(1 + cv[l]^2)))
  if (!is.null(dims)) dim(out) <- dims
  out
}

# Gauss-Hermite nodes/weights (Golub-Welsch), for E[g(Z)], Z ~ N(0, 1)
gauss_hermite <- function(n = 30) {
  i <- seq_len(n - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- sqrt(i / 2)
  J[cbind(i + 1, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE)
  list(z = sqrt(2) * e$values, w = e$vectors[1, ]^2)   # weights sum to 1
}
GH <- gauss_hermite(30)

# Probability approach: requirement ~ log-normal with median EAR and CV_REQ
p_below_probability <- function(ear, mu, cv, dist, cv_req = CV_REQ) {
  sd_req <- sqrt(log(1 + cv_req^2))
  Reduce(`+`, lapply(seq_along(GH$z), function(k) {
    GH$w[k] * p_below(ear * exp(sd_req * GH$z[k]), mu, cv, dist)
  }))
}

wmean <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

clean_sex <- function(x) dplyr::case_when(
  x %in% c("Female", "female", "Females", "F") ~ "Female",
  x %in% c("Male", "male", "Males", "M") ~ "Male",
  TRUE ~ as.character(x)
)

age_start <- function(a) suppressWarnings(as.numeric(sub("^([0-9.]+).*", "\\1", a)))

tier_group <- function(tier) dplyr::case_when(
  tier == "exact" ~ "Exact country-sex-age match",
  tier %in% c("nearest_age", "opposite_sex") ~ "Same country (nearest age or opposite sex)",
  !is.na(tier) ~ "Donor country",
  TRUE ~ NA_character_
)

# Revision (R1.9, R4.3): the donor step can pick a country as its own nearest
# donor (it has nutriR data for other age-sex groups). These strata use the
# country's own data, so the SI reports them as same country. Reporting only:
# the donor bootstrap (G) still treats them as donor-matched strata.
tier_group_own <- function(tier, iso3, donor) dplyr::if_else(
  tier_group(tier) %in% "Donor country" & !is.na(donor) & donor == iso3,
  "Same country (nearest age or opposite sex)",
  tier_group(tier)
)

# ==============================================================
# Load data
# ==============================================================
dat <- readRDS(file.path("output", "dat_quality.rds")) %>%
  distinct(iso3, sex, age_group, .keep_all = TRUE) %>%
  filter(age_group != "0-0.99") %>%
  mutate(
    sex_label = clean_sex(sex),
    region  = countrycode(iso3, origin = "iso3c", destination = "region"),
    country = countrycode(iso3, origin = "iso3c", destination = "country.name"),
    country = if_else(is.na(country), iso3, country),
    adult   = !is.na(opt_mean_g_day)
  ) %>%
  arrange(iso3, sex, age_start(age_group))

stopifnot(all(c("gdd_mean", "gdd_lower", "gdd_upper", "weight_mean_kg", "weight_low_kg",
                "weight_high_kg", "eer_kcal_marco_mean", "eer_kcal_marco_low",
                "eer_kcal_marco_high", "protein_kcal_share_mean", "ear_g_per_kg",
                "tgt16_g_per_kg", "iaao_ear_g_per_kg", "opt_g_per_kg") %in% names(dat)))

N <- nrow(dat)
pop <- dat$population
cv_p <- dat$cv_protein; dist_p <- dat$best_dist_protein

# ---- Aggregation matrices: rows = groups, columns = strata, entries = population
make_W <- function(groups_df, keep = rep(TRUE, N)) {
  # groups_df: one row per stratum with columns group_level, group (long format)
  lev <- groups_df %>% distinct(group_level, group)
  W <- matrix(0, nrow(lev), N)
  key <- paste(lev$group_level, lev$group)
  idx <- match(paste(groups_df$group_level, groups_df$group), key)
  W[cbind(idx, groups_df$stratum)] <- pop[groups_df$stratum] * keep[groups_df$stratum]
  list(W = W, groups = lev)
}

strata_groups <- function(levels = c("global", "region", "sex", "age", "country")) {
  s <- seq_len(N)
  bind_rows(
    if ("global"  %in% levels) tibble(stratum = s, group_level = "global",  group = "Global"),
    if ("region"  %in% levels) tibble(stratum = s, group_level = "region",  group = dat$region),
    if ("sex"     %in% levels) tibble(stratum = s, group_level = "sex",     group = dat$sex_label),
    if ("age"     %in% levels) tibble(stratum = s, group_level = "age",     group = dat$age_group),
    if ("country" %in% levels) tibble(stratum = s, group_level = "country", group = dat$iso3)
  ) %>% filter(!is.na(group))
}

aggregate_W <- function(Wl, P) {
  P[is.na(P)] <- 0
  G <- Wl$W %*% P
  den <- rowSums(Wl$W)
  G / ifelse(den > 0, den, NA_real_)
}

W_all   <- make_W(strata_groups())
W_small <- make_W(strata_groups(c("global", "region", "sex", "age")))
W_small_adult <- make_W(strata_groups(c("global", "region", "sex", "age")), keep = dat$adult)

summarise_strata <- function(p, keep = rep(TRUE, N), levels = c("global", "region", "sex", "age")) {
  Wl <- make_W(strata_groups(levels), keep = keep)
  tibble(Wl$groups, value = as.vector(aggregate_W(Wl, matrix(p, ncol = 1))))
}

# ==============================================================
# Central (deterministic) inputs per scenario
# ==============================================================
protein_mean <- function(kcal) kcal * dat$protein_kcal_share_mean / 4

central_p <- lapply(SCENARIOS, function(col) {
  mu <- protein_mean(dat[[col]])
  list(
    EAR  = p_below(dat$ear_mean_g_day,      mu, cv_p, dist_p),
    IAAO = p_below(dat$iaao_ear_mean_g_day, mu, cv_p, dist_p),
    T12  = p_below(dat$opt_mean_g_day,      mu, cv_p, dist_p),
    T16  = p_below(dat$tgt16_mean_g_day,    mu, cv_p, dist_p),
    energy = p_below(dat$mder_stratum_kcal, dat[[col]], dat$cv_calorie, dat$best_dist_calorie),
    mu = mu
  )
})

# ==============================================================
# A. Monte Carlo uncertainty analysis (R1.4, R1.11)
#
# Each draw samples, for every country, one standard-normal deviate for GDD
# protein intake and one for body weight. Within a country the same deviates
# apply to all sex-age strata (estimates for a country share the same surveys
# and model); countries are independent.
#   - Protein energy share: log-normal around the GDD median, sd on the log scale
#     from the GDD 95% uncertainty interval.
#   - Body weight: normal around the NCD-RisC mean, sd from its 95% interval.
#   - EER: normal around the mean, sd from its 95% interval, driven by the SAME
#     deviate as body weight, so a heavier draw raises both the requirement and
#     the EER used to allocate national energy intake (Stu's point).
# National energy intake per scenario stays fixed (it is the scenario); within
# a country it is reallocated across strata in proportion to the drawn EERs.
# Not included (assessed separately): donor choice for intake distributions,
# CVs, distribution family, waste function / calorie scenario.
# ==============================================================
set.seed(SEED)

iso_idx <- match(dat$iso3, unique(dat$iso3))
n_c <- max(iso_idx)

sdlog_gdd <- ifelse(dat$gdd_lower > 0 & dat$gdd_upper > 0,
                    (log(dat$gdd_upper) - log(dat$gdd_lower)) / (2 * Z95), 0)
sd_w   <- pmax(0, (dat$weight_high_kg - dat$weight_low_kg) / (2 * Z95))
sd_eer <- pmax(0, (dat$eer_kcal_marco_high - dat$eer_kcal_marco_low) / (2 * Z95))
sdlog_gdd[is.na(sdlog_gdd)] <- 0; sd_w[is.na(sd_w)] <- 0; sd_eer[is.na(sd_eer)] <- 0

# country-level population-weighted mean EER (central), excluding infants
Wc <- matrix(0, n_c, N); Wc[cbind(iso_idx, seq_len(N))] <- pop
avg_eer_c <- as.vector(Wc %*% dat$eer_kcal_marco_mean) / rowSums(Wc)

metrics <- c("EAR", "IAAO", "T12", "T16")
mc_store <- setNames(lapply(names(SCENARIOS), function(s) setNames(vector("list", length(metrics)), metrics)),
                     names(SCENARIOS))

n_chunks <- ceiling(N_DRAWS / CHUNK)
cat("Monte Carlo:", N_DRAWS, "draws in", n_chunks, "chunks\n")
for (ch in seq_len(n_chunks)) {
  B <- min(CHUNK, N_DRAWS - (ch - 1) * CHUNK)
  zg <- matrix(rnorm(n_c * B), n_c, B)[iso_idx, , drop = FALSE]
  zw <- matrix(rnorm(n_c * B), n_c, B)[iso_idx, , drop = FALSE]

  share_d <- dat$protein_kcal_share_mean * exp(sdlog_gdd * zg)
  w_d     <- pmax(dat$weight_mean_kg + sd_w * zw, 0.5 * dat$weight_mean_kg)
  eer_d   <- pmax(dat$eer_kcal_marco_mean + sd_eer * zw, 0.5 * dat$eer_kcal_marco_mean)
  avg_d   <- (Wc %*% eer_d) / rowSums(Wc)                      # n_c x B
  realloc <- (eer_d / dat$eer_kcal_marco_mean) * (avg_eer_c[iso_idx] / avg_d[iso_idx, , drop = FALSE])

  req <- list(
    EAR  = dat$ear_g_per_kg      * w_d,
    IAAO = dat$iaao_ear_g_per_kg * w_d,
    T12  = dat$opt_g_per_kg      * w_d,
    T16  = dat$tgt16_g_per_kg    * w_d
  )
  for (s in names(SCENARIOS)) {
    mu_d <- dat[[SCENARIOS[[s]]]] * realloc * share_d / 4
    for (m in metrics) {
      P <- p_below(req[[m]], mu_d, cv_p, dist_p)
      G <- if (m == "EAR") aggregate_W(W_all, P) else aggregate_W(W_small_adult, P)
      mc_store[[s]][[m]] <- cbind(mc_store[[s]][[m]], G)
    }
  }
  if (ch %% 5 == 0 || ch == n_chunks) cat("  chunk", ch, "of", n_chunks, "done\n")
}

point_estimates <- bind_rows(lapply(names(SCENARIOS), function(s) {
  bind_rows(
    tibble(W_all$groups, scenario = s, metric = "EAR",
           point = as.vector(aggregate_W(W_all, matrix(central_p[[s]]$EAR, ncol = 1)))),
    bind_rows(lapply(c("IAAO", "T12", "T16"), function(m) {
      tibble(W_small_adult$groups, scenario = s, metric = m,
             point = as.vector(aggregate_W(W_small_adult, matrix(central_p[[s]][[m]], ncol = 1))))
    }))
  )
}))

mc_summary_all <- bind_rows(lapply(names(SCENARIOS), function(s) {
  bind_rows(lapply(metrics, function(m) {
    G <- mc_store[[s]][[m]]
    groups <- if (m == "EAR") W_all$groups else W_small_adult$groups
    tibble(groups, scenario = s, metric = m,
           mc_mean = rowMeans(G, na.rm = TRUE),
           q025 = apply(G, 1, quantile, probs = 0.025, na.rm = TRUE),
           q975 = apply(G, 1, quantile, probs = 0.975, na.rm = TRUE))
  }))
})) %>%
  left_join(point_estimates, by = c("group_level", "group", "scenario", "metric")) %>%
  select(scenario, metric, group_level, group, point, mc_mean, q025, q975)

mc_summary <- mc_summary_all %>% filter(group_level != "country")
mc_country <- mc_summary_all %>%
  filter(group_level == "country") %>%
  rename(iso3 = group) %>%
  mutate(country = countrycode(iso3, "iso3c", "country.name"),
         country = if_else(is.na(country), iso3, country)) %>%
  select(scenario, metric, iso3, country, point, mc_mean, q025, q975)

write_csv(mc_summary, file.path(out_dir, "mc_summary.csv"))
write_csv(mc_country, file.path(out_dir, "mc_country.csv"))

cat("\nGlobal protein inadequacy (EAR cut-point) with Monte Carlo 95% UI:\n")
print(mc_summary %>% filter(group_level == "global") %>%
        mutate(across(c(point, mc_mean, q025, q975), ~ round(100 * .x, 2))))

# ==============================================================
# B. Full probability approach vs EAR cut-point (R1.3)
# ==============================================================
prob_approach <- bind_rows(lapply(names(SCENARIOS), function(s) {
  mu <- central_p[[s]]$mu
  p_prob <- p_below_probability(dat$ear_mean_g_day, mu, cv_p, dist_p)
  left_join(
    summarise_strata(central_p[[s]]$EAR) %>% rename(cut_point = value),
    summarise_strata(p_prob) %>% rename(probability_approach = value),
    by = c("group_level", "group")
  ) %>% mutate(scenario = s)
})) %>%
  mutate(difference_pp = 100 * (probability_approach - cut_point),
         across(c(cut_point, probability_approach), ~ round(100 * .x, 2)),
         difference_pp = round(difference_pp, 2)) %>%
  select(scenario, group_level, group, cut_point, probability_approach, difference_pp)
write_csv(prob_approach, file.path(out_dir, "si_probability_approach.csv"))

# ==============================================================
# C. Distribution family (R1.12)
# ==============================================================
dist_family <- bind_rows(lapply(names(SCENARIOS), function(s) {
  mu <- central_p[[s]]$mu
  reduce(list(
    summarise_strata(central_p[[s]]$EAR, levels = c("global", "region", "sex")) %>% rename(nutriR_assigned = value),
    summarise_strata(p_below(dat$ear_mean_g_day, mu, cv_p, "gamma"), levels = c("global", "region", "sex")) %>% rename(all_gamma = value),
    summarise_strata(p_below(dat$ear_mean_g_day, mu, cv_p, "log-normal"), levels = c("global", "region", "sex")) %>% rename(all_lognormal = value)
  ), left_join, by = c("group_level", "group")) %>% mutate(scenario = s)
})) %>%
  mutate(across(c(nutriR_assigned, all_gamma, all_lognormal), ~ round(100 * .x, 2))) %>%
  select(scenario, group_level, group, nutriR_assigned, all_gamma, all_lognormal)
write_csv(dist_family, file.path(out_dir, "si_distribution_family.csv"))

family_counts <- dat %>% count(best_dist_protein, wt = population, name = "population") %>%
  mutate(share_population = round(100 * population / sum(population), 1))
write_csv(family_counts, file.path(out_dir, "si_distribution_family_counts.csv"))

# ==============================================================
# D. Protein-quality adjustment grid (R1.6, R2.3)
# ==============================================================
quality_factor_fun <- function(prop_asf, qf_max, asf_upper, asf_lower) {
  f <- 1 + qf_max * pmin(1, pmax(0, (asf_upper - prop_asf) / (asf_upper - asf_lower)))
  ifelse(is.na(f), 1, f)
}
quality_grid <- bind_rows(lapply(seq_len(nrow(QF_GRID)), function(i) {
  g <- QF_GRID[i, ]
  qf <- quality_factor_fun(dat$prop_asf, g$qf_max, g$asf_upper, g$asf_lower)
  bind_rows(lapply(names(SCENARIOS), function(s) {
    summarise_strata(p_below(dat$ear_mean_g_day * qf, central_p[[s]]$mu, cv_p, dist_p),
                     levels = c("global", "region")) %>%
      mutate(scenario = s)
  })) %>% mutate(qf_max = g$qf_max, asf_upper = g$asf_upper, asf_lower = g$asf_lower)
})) %>%
  mutate(prevalence_pct = round(100 * value, 2)) %>%
  select(qf_max, asf_upper, asf_lower, scenario, group_level, group, prevalence_pct) %>%
  left_join(
    bind_rows(lapply(names(SCENARIOS), function(s)
      summarise_strata(central_p[[s]]$EAR, levels = c("global", "region")) %>% mutate(scenario = s))) %>%
      transmute(scenario, group_level, group, unadjusted_pct = round(100 * value, 2)),
    by = c("scenario", "group_level", "group")
  )
write_csv(quality_grid, file.path(out_dir, "si_quality_grid.csv"))

# ==============================================================
# E. Higher intake levels proposed for adults (R1.5, R1.15, R2.6)
# ==============================================================
target_summary <- mc_summary %>%
  filter(metric %in% c("IAAO", "T12", "T16")) %>%
  mutate(metric = recode(metric,
                         IAAO = "Below 0.93 g/kg/d (IAAO EAR-type), adults",
                         T12  = "Below 1.2 g/kg/d (proposed target intake), adults",
                         T16  = "Below 1.6 g/kg/d (proposed target intake), adults"),
         across(c(point, mc_mean, q025, q975), ~ round(100 * .x, 2)))
write_csv(target_summary, file.path(out_dir, "si_target_intakes_summary.csv"))

adult_pop <- sum(pop[dat$adult])
target_headcount <- point_estimates %>%
  filter(group_level == "global", metric %in% c("IAAO", "T12", "T16")) %>%
  mutate(adult_population_billion = adult_pop / 1e9,
         adults_below_billion = point * adult_pop / 1e9)
write_csv(target_headcount, file.path(out_dir, "si_target_intakes_headcount.csv"))

# ==============================================================
# F. nutriR matching tiers and donors (R1.9, R2.4, R4.3)
# ==============================================================
tiers_long <- bind_rows(
  dat %>% transmute(nutrient = "Protein", iso3, country, region, population, tier = match_tier_protein, donor = donor_iso3_protein),
  dat %>% transmute(nutrient = "Energy",  iso3, country, region, population, tier = match_tier_energy,  donor = donor_iso3_energy)
) %>% mutate(tier_group = tier_group_own(tier, iso3, donor))

tier_summary <- bind_rows(
  tiers_long %>% mutate(region = "World"),
  tiers_long %>% filter(!is.na(region))
) %>%
  group_by(nutrient, region, tier_group, tier) %>%
  summarise(strata = n(), population = sum(population), .groups = "drop_last") %>%
  group_by(nutrient, region) %>%
  mutate(pct_strata = round(100 * strata / sum(strata), 1),
         pct_population = round(100 * population / sum(population), 1)) %>%
  ungroup()
write_csv(tier_summary, file.path(out_dir, "si_matching_tiers.csv"))

matching_by_country <- dat %>%
  mutate(grp_protein = tier_group_own(match_tier_protein, iso3, donor_iso3_protein),
         grp_energy  = tier_group_own(match_tier_energy,  iso3, donor_iso3_energy)) %>%
  group_by(iso3, country, region) %>%
  summarise(
    exact_protein_pct = round(100 * sum(population[match_tier_protein == "exact"]) / sum(population), 1),
    donor_country_protein_pct = round(100 * sum(population[grp_protein == "Donor country"]) / sum(population), 1),
    protein_donor = paste(sort(unique(donor_iso3_protein[grp_protein == "Donor country"])), collapse = ", "),
    exact_energy_pct = round(100 * sum(population[match_tier_energy == "exact"]) / sum(population), 1),
    donor_country_energy_pct = round(100 * sum(population[grp_energy == "Donor country"]) / sum(population), 1),
    energy_donor = paste(sort(unique(donor_iso3_energy[grp_energy == "Donor country"])), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(region, country)
write_csv(matching_by_country, file.path(out_dir, "si_matching_by_country.csv"))

# Additions to Supplementary Table 1 (read by Script 8)
exact_share <- function(df) tibble(
  `Exact nutriR match, protein (% of population)` = round(100 * sum(df$population[df$match_tier_protein == "exact"]) / sum(df$population), 1),
  `Exact nutriR match, energy (% of population)`  = round(100 * sum(df$population[df$match_tier_energy  == "exact"]) / sum(df$population), 1)
)
ui_lookup <- mc_summary_all %>%
  filter(scenario == "central", metric == "EAR", group_level %in% c("global", "region", "country")) %>%
  mutate(`Country / Region` = case_when(
    group_level == "global" ~ "World",
    group_level == "region" ~ group,
    TRUE ~ if_else(is.na(countrycode(group, "iso3c", "country.name")), group,
                   countrycode(group, "iso3c", "country.name"))
  ),
  `Protein inadequacy (% , Central, 95% UI)` = sprintf("%.2f-%.2f", 100 * q025, 100 * q975)) %>%
  select(`Country / Region`, `Protein inadequacy (% , Central, 95% UI)`)
si_table1_additions <- bind_rows(
  exact_share(dat) %>% mutate(`Country / Region` = "World"),
  dat %>% filter(!is.na(region)) %>% group_by(`Country / Region` = region) %>% group_modify(~ exact_share(.x)) %>% ungroup(),
  dat %>% group_by(`Country / Region` = country) %>% group_modify(~ exact_share(.x)) %>% ungroup()
) %>%
  left_join(ui_lookup, by = "Country / Region") %>%
  select(`Country / Region`, `Protein inadequacy (% , Central, 95% UI)`, everything())
write_csv(si_table1_additions, file.path(out_dir, "si_table1_additions.csv"))

# ==============================================================
# G. Donor bootstrap (R4.3, R4.6, R4.4)
# Strata without their own nutriR distribution take the shape of the single
# nearest donor country (Scripts 1 and 3). Each bootstrap draw instead gives
# every country one of its DONOR_K = 5 nearest donor countries, drawn with equal
# probabilities (the same donor for all of a country's donor-matched strata).
# The waste functions are scenario assumptions fixed once in Script 3 and are
# not changed here. Central scenario.
#   1. Protein donors drawn: protein inadequacy (EAR).
#   2. Energy donors drawn: energy inadequacy around the 8% benchmark. Protein
#      inadequacy does not change, because mean energy intake (supply and
#      waste) is unchanged and protein intake depends on energy only through it.
# ==============================================================
N_BOOT <- 1000

donor_boot <- tryCatch({
  # inputs saved by Scripts 1 and 3
  source(file.path("scripts", "helpers_donor_matching.R"))
  inp_p <- readRDS(file.path(out_dir, "donor_inputs_protein.rds"))
  inp_e <- readRDS(file.path(out_dir, "donor_inputs_energy.rds"))
  dat_key <- paste(dat$iso3, dat$sex, dat$age_group)

  # N x DONOR_K matrices of CV and family from the r-th nearest donor;
  # column 1 is the main analysis
  rank_matrices <- function(inp, cv_main, dist_main, tier, donor_main, protein) {
    don <- tier_group(tier) %in% "Donor country"
    first <- tapply(donor_main[don], dat$iso3[don], function(x) x[1])
    rs <- donor_rank_shapes(inp, k = DONOR_K, first = first)
    if (protein) {
      # Script 1 builds the 1-4 group from the 2-4 shape
      rs <- rs %>%
        filter(!age_group %in% c("0-0.99", "1-1.99")) %>%
        mutate(age_group = if_else(age_group == "2-4", "1-4", age_group))
    }
    rs <- rs %>% mutate(row = match(paste(iso3, sex, age_group), dat_key)) %>%
      filter(!is.na(row), don[row])
    CV <- matrix(cv_main, N, DONOR_K)
    DS <- matrix(as.character(dist_main), N, DONOR_K)
    for (r in 2:DONOR_K) {
      x <- rs[rs$donor_rank == r, ]
      CV[x$row, r] <- x$cv; DS[x$row, r] <- x$best_dist
    }
    list(CV = CV, DS = DS, n_donor_strata = sum(don))
  }
  RP <- rank_matrices(inp_p, cv_p, dist_p, dat$match_tier_protein, dat$donor_iso3_protein, protein = TRUE)
  RE <- rank_matrices(inp_e, dat$cv_calorie, dat$best_dist_calorie, dat$match_tier_energy, dat$donor_iso3_energy, protein = FALSE)
  cat(sprintf("\nDonor bootstrap: %d protein and %d energy donor-matched strata\n",
              RP$n_donor_strata, RE$n_donor_strata))

  mu_c <- central_p$central$mu
  Wg <- make_W(strata_groups(c("global", "region")))
  agg <- function(p) as.vector(aggregate_W(Wg, matrix(p, ncol = 1)))
  pick <- function(M, rk) M[cbind(seq_len(N), rk)]

  set.seed(SEED + 9)
  Gp <- matrix(NA_real_, nrow(Wg$groups), N_BOOT)
  Ge <- matrix(NA_real_, nrow(Wg$groups), N_BOOT)
  for (b in seq_len(N_BOOT)) {
    rp <- sample.int(DONOR_K, n_c, replace = TRUE)[iso_idx]
    re <- sample.int(DONOR_K, n_c, replace = TRUE)[iso_idx]
    Gp[, b] <- agg(p_below(dat$ear_mean_g_day, mu_c, pick(RP$CV, rp), pick(RP$DS, rp)))
    Ge[, b] <- agg(p_below(dat$mder_stratum_kcal, dat$kcal_mean_medium, pick(RE$CV, re), pick(RE$DS, re)))
  }

  summarise_boot <- function(G, main, variant_label, outcome_label) {
    below <- if (outcome_label == "Protein inadequacy (EAR)") {
      ifelse(Wg$groups$group_level == "global", round(100 * rowMeans(G < inp_e$target), 1), NA_real_)
    } else NA_real_
    tibble(Wg$groups,
           variant = variant_label, outcome = outcome_label,
           main_pct   = round(100 * main, 2),
           median_pct = round(100 * apply(G, 1, median), 2),
           q025_pct   = round(100 * apply(G, 1, quantile, probs = 0.025), 2),
           q975_pct   = round(100 * apply(G, 1, quantile, probs = 0.975), 2),
           pct_draws_below_energy_benchmark = below)
  }
  bind_rows(
    summarise_boot(Gp, agg(central_p$central$EAR), "Protein donors drawn", "Protein inadequacy (EAR)"),
    summarise_boot(Ge, agg(central_p$central$energy), "Energy donors drawn", "Energy inadequacy")
  ) %>%
    select(variant, outcome, group_level, group, everything())
}, error = function(e) {
  message("Donor bootstrap skipped: ", conditionMessage(e))
  tibble(note = paste("Donor bootstrap skipped:", conditionMessage(e)))
})
write_csv(donor_boot, file.path(out_dir, "si_donor_bootstrap.csv"))

if ("variant" %in% names(donor_boot)) {
  cat("\nDonor bootstrap, central scenario (global, %):\n")
  print(as.data.frame(donor_boot %>% filter(group_level == "global")))
}

# ==============================================================
# H. Protein-to-energy ratios (R1.17)
#   Diet P/E     = protein energy share reconstructed from GDD (by stratum)
#   Reference P/E = 4 x requirement (g/d) / EER (kcal/d): the P/E at which mean
#   intake equals the requirement when energy intake equals the EER. With the
#   RDA this is the "safe" P/E in the sense of Millward & Jackson (2004).
#   Mean intake / RDA in the energy-adequate counterfactual = diet P/E / safe P/E.
# ==============================================================
pe <- dat %>%
  mutate(
    pe_diet = protein_kcal_share_mean,
    pe_ref_EAR = 4 * ear_mean_g_day / eer_kcal_marco_mean,
    pe_ref_RDA = 4 * rda_mean_g_day / eer_kcal_marco_mean
  )
pe_summary <- pe %>%
  group_by(age_group, sex = sex_label) %>%
  summarise(
    pe_diet_mean   = wmean(pe_diet, population),
    pe_diet_p10    = quantile(pe_diet, 0.10, na.rm = TRUE),
    pe_diet_p90    = quantile(pe_diet, 0.90, na.rm = TRUE),
    pe_ref_EAR_mean = wmean(pe_ref_EAR, population),
    pe_ref_RDA_mean = wmean(pe_ref_RDA, population),
    pct_pop_diet_below_ref_RDA = 100 * sum(population[pe_diet < pe_ref_RDA], na.rm = TRUE) / sum(population),
    .groups = "drop"
  ) %>%
  arrange(sex, age_start(age_group)) %>%
  mutate(across(starts_with("pe_"), ~ round(.x, 3)), pct_pop_diet_below_ref_RDA = round(pct_pop_diet_below_ref_RDA, 1))
write_csv(pe_summary, file.path(out_dir, "si_pe_ratios.csv"))

# ==============================================================
# I. Countries where protein inadequacy prevalence exceeds energy
#    inadequacy prevalence, central scenario (R2.m4, R2.7)
#    Same rule as Figure 3c (ratio > 1.05; energy prevalence >= 0.5%)
# ==============================================================
country_central <- dat %>%
  mutate(p_prot = central_p$central$EAR, p_energy = central_p$central$energy) %>%
  group_by(iso3, country, region) %>%
  summarise(protein_inadequacy = wmean(p_prot, population),
            energy_inadequacy = wmean(p_energy, population),
            asf_share = first(prop_asf), .groups = "drop") %>%
  mutate(ratio = protein_inadequacy / ifelse(energy_inadequacy < 0.005, NA, energy_inadequacy))
protein_exceeds_energy <- country_central %>%
  filter(!is.na(ratio), ratio > 1.05) %>%
  left_join(mc_country %>% filter(scenario == "central", metric == "EAR") %>% select(iso3, q025, q975), by = "iso3") %>%
  transmute(country, region,
            `Protein inadequacy (%)` = round(100 * protein_inadequacy, 1),
            `Protein inadequacy 95% UI` = sprintf("%.1f-%.1f", 100 * q025, 100 * q975),
            `Energy inadequacy (%)` = round(100 * energy_inadequacy, 1),
            `Ratio protein / energy` = round(ratio, 2),
            `Protein from animal-source foods (%)` = round(100 * asf_share, 1)) %>%
  arrange(desc(`Ratio protein / energy`))
write_csv(protein_exceeds_energy, file.path(out_dir, "si_protein_exceeds_energy.csv"))

# ==============================================================
# J. Imputation log (R2.4)
# ==============================================================
read_if <- function(f) if (file.exists(f)) read_csv(f, show_col_types = FALSE, col_types = cols(.default = "c")) else NULL
imputation_log <- bind_rows(
  read_if(file.path(out_dir, "imputation_log_weights.csv")),
  dat %>% filter(!is.na(eer_imputed_source)) %>%
    transmute(input = "Estimated energy requirement (Springmann 2025)", iso3, sex, age_group,
              imputation = eer_imputed_source),
  read_if(file.path(out_dir, "imputation_log_country_inputs_script3.csv")),
  read_if(file.path(out_dir, "imputation_log_mder.csv"))
) %>%
  mutate(country = countrycode(iso3, "iso3c", "country.name")) %>%
  select(input, iso3, country, any_of(c("sex", "age_group")), imputation) %>%
  arrange(input, iso3)
write_csv(imputation_log, file.path(out_dir, "si_imputation_log.csv"))

# ==============================================================
# K. Calibration check (R1.16, R2.1)
# ==============================================================
calibration <- tibble(
  scenario = names(SCENARIOS),
  global_energy_inadequacy_pct = sapply(names(SCENARIOS), function(s) round(100 * wmean(central_p[[s]]$energy, pop), 2)),
  global_protein_inadequacy_pct = sapply(names(SCENARIOS), function(s) round(100 * wmean(central_p[[s]]$EAR, pop), 2))
)
calib_file <- file.path(out_dir, "calibration_parameters.csv")
if (file.exists(calib_file)) {
  cat("\nWaste-function parameters (from Script 3):\n")
  print(read_csv(calib_file, show_col_types = FALSE))
}
write_csv(calibration, file.path(out_dir, "si_calibration_check.csv"))
cat("\nRealised global energy and protein inadequacy by scenario:\n")
print(calibration)

# ==============================================================
# One workbook with everything
# ==============================================================
if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(list(
    mc_summary = mc_summary,
    mc_country = mc_country,
    probability_approach = prob_approach,
    distribution_family = dist_family,
    quality_grid = quality_grid,
    target_intakes = target_summary,
    matching_tiers = tier_summary,
    matching_by_country = matching_by_country,
    donor_bootstrap = donor_boot,
    pe_ratios = pe_summary,
    protein_exceeds_energy = protein_exceeds_energy,
    imputation_log = imputation_log,
    calibration_check = calibration
  ), file.path(out_dir, "revision_tables.xlsx"))
}

cat("\nScript 6 outputs written to ", normalizePath(out_dir, mustWork = FALSE), "\n", sep = "")
