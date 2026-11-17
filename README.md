# Model-based estimates of inadequate protein intake worldwide

This repository contains the analysis code for the manuscript:

**Model-based estimates of inadequate protein intake worldwide** *(under review)*  
Giacomo De Nicola, Ty Beal, Stella Nordhagen, Stuart M. Phillips, Simone Passarelli, Christopher D. Golden

The analysis estimates the prevalence of inadequate protein intake across 184 countries in 2018, disaggregated by sex and age group. It combines modeled dietary intake data, protein requirement thresholds, caloric intake scenarios, intake distribution parameters, food supply data, body-weight data, and population estimates to reconstruct country-sex-age-specific protein and energy intake distributions. All estimates are model-based: they rely on reconstructed rather than observed intake distributions.

## Overview

For each country-sex-age stratum, the analysis reconstructs mean protein intake, specifies an intake distribution, applies a protein requirement, and estimates inadequacy as the probability that usual intake falls below the requirement (EAR cut-point method).

The analysis has two main components.

First, we estimate mean protein intake under a deterministic **energy-adequate counterfactual**, where each stratum is assumed to consume energy equal to its Estimated Energy Requirement (EER). This isolates the role of dietary protein density from total energy insufficiency.

Second, we estimate protein and energy inadequacy under three **realistic caloric intake scenarios**, where country-level energy intake is derived from FAO Food Balance Sheet energy supply adjusted for modeled downstream food waste. The central scenario is calibrated to 8% global energy inadequacy (FAO prevalence of undernourishment for 2018), the low-calorie scenario to 16%, and the high-calorie scenario assumes no downstream waste.

Additional analyses: a Monte Carlo uncertainty analysis, the full probability approach, alternative distribution families, a donor-country bootstrap, an exploratory analysis of higher intake levels proposed for adults (0.93, 1.2 and 1.6 g/kg/d), an exploratory low-ASF scenario for protein quality, and protein-to-energy ratios.

## Quick start

1. Install R (version 4.x) and RStudio, and the packages listed under *Software requirements*.
2. Obtain the input data (see *Input data*) and place it in `data/` with the folder structure shown below.
3. Open the project in RStudio and run the scripts in this order, each from a clean session:

```r
source("scripts/1_main_protein.R")
source("scripts/2_protein_inadequacy_under_optimal_caloric_intake.R")
source("scripts/3_caloric_intake_estimation.R")
source("scripts/4_caloric_and_protein_inadequacy_calculation_and_overlap.R")
source("scripts/5_protein_composition.R")
source("scripts/6_revision_uncertainty_and_sensitivity.R")
source("scripts/7_output_figures.R")
source("scripts/8_create_supplementary_tables.R")
source("scripts/9_dataverse_tables.R")
```

Scripts 1-5 set the working directory with `rstudioapi`, so run them from RStudio (or set the working directory to the repository root and remove the `setwd()` lines). Script 3 downloads GDP data from the World Bank (`WDI` package) on its first run, which needs an internet connection, and saves it as `data/WDI_gdp_per_capita_ppp_2018.csv`; later runs read that file. Script 6 must run before Scripts 7 and 8, which read its Monte Carlo output.

## Input data

The `data/` folder is not tracked in Git because of file size. All inputs are publicly available from the sources below; the assembled folder is also available from the corresponding author on request.

```text
data/
  GDD_FinalEstimates_01102022/Country-level estimates/
      v23_cnty.csv                         Global Dietary Database 2018, total protein (Script 1)
      v01_cnty.csv ... v18_cnty.csv, v57_cnty.csv
                                           GDD 2018 food-group intakes (Script 5)
  data_extract_040725.xlsx                 sheet 1: Estimated Energy Requirements (Springmann 2025)
                                           sheet 2: body weight by country, sex and age (NCD-RisC)
  OWID_daily-per-capita-caloric-supply.csv FAO Food Balance Sheets dietary energy supply (via Our World in Data)
  OWID_minimum-requirement-calories.csv    FAO Minimum Dietary Energy Requirement (via Our World in Data)
  USDA Food Composition Data (FNDDS)/
      food.csv, nutrient.csv, food_nutrient.csv, wweia_food_category.csv
                                           USDA Food and Nutrient Database for Dietary Studies 2021-2023
```

Sources:

- Global Dietary Database, 2018 estimates: https://www.globaldietarydatabase.org
- Energy requirements: Springmann M. Supplementary Datafile for "Estimates of energy intake, requirements, and imbalances based on anthropometric measurements at global, regional, and national levels and for sociodemographic groups: a modelling study". Zenodo, 2025. doi:10.5281/zenodo.16760810
- Body weight: NCD Risk Factor Collaboration (NCD-RisC), https://www.ncdrisc.org
- FAO Food Balance Sheets and Minimum Dietary Energy Requirement: https://www.fao.org/faostat (accessed via https://ourworldindata.org)
- USDA FNDDS: https://www.ars.usda.gov/northeast-area/beltsville-md-bhnrc/beltsville-human-nutrition-research-center/food-surveys-research-group/
- Downloaded at run time through R packages: intake distribution shapes and CVs (`nutriR`), population (`wpp2024`, UN World Population Prospects 2024), GDP per capita PPP (`WDI`, World Bank).

All external datasets remain subject to their own licenses and terms of use.

## What each script does and writes

| Script | Purpose | Main outputs (in `output/`) |
|---|---|---|
| `1_main_protein.R` | GDD protein intake, nutriR protein distributions (with matching tier and donor country), body weight, protein requirements (EAR, RDA, and the adult-only 0.93 / 1.2 / 1.6 g/kg/d levels), population | `gdd_prepped_data_for_caloric_adjustment.rds`, `wpp2024_population_aggregated_2018.rds`, `revision/imputation_log_weights.csv`, `revision/donor_inputs_protein.rds` |
| `2_protein_inadequacy_under_optimal_caloric_intake.R` | Protein energy shares; energy-adequate counterfactual (energy = EER) | `script2_final_results.rds` |
| `3_caloric_intake_estimation.R` | Energy supply, GDP-linked waste functions calibrated to 8% and 16% global energy inadequacy, disaggregation to strata, nutriR energy distributions (donor countries ranked on energy supply before waste; the same shapes in all scenarios) | `final_calorie_distributions_wide.rds`, `revision/calibration_parameters.csv`, `revision/imputation_log_country_inputs_script3.csv`, `revision/donor_inputs_energy.rds` |
| `4_caloric_and_protein_inadequacy_calculation_and_overlap.R` | MDER, protein and energy inadequacy (central scenario), calibration checks | `final_analysis_with_all_inadequacy.rds`, `revision/script4_revision_checks.csv` |
| `5_protein_composition.R` | Share of protein from animal-source foods; adjusted EAR for the exploratory low-ASF scenario (Supplementary Note 1) | `dat_quality.rds`, `protein_asf_props.rds` |
| `6_revision_uncertainty_and_sensitivity.R` | Monte Carlo uncertainty intervals, probability approach, distribution family, quality grid, higher intake levels, matching tiers, donor bootstrap, P/E ratios, imputation log | `revision/*.csv`, `revision/revision_tables.xlsx` |
| `7_output_figures.R` | Main figures | `main_figures_clean/` |
| `8_create_supplementary_tables.R` | Supplementary tables | `supplementary_tables_clean/` |
| `9_dataverse_tables.R` | Dataverse export files | `dataverse/` |
| `helpers_donor_matching.R` | Helpers for the donor bootstrap (shapes from each country's 5 nearest donor countries); sourced by Script 6 | none |
| `999_extra_visualizations.R` | Exploratory plots (not needed for the paper) | `extra_plots/` |

The main analytical dataset is `output/dat_quality.rds`: one row per country-sex-age stratum, with protein and energy intake, requirement thresholds, distribution parameters, matching tiers, inadequacy estimates, and the animal-source-food protein share.

## Reproducing the figures and tables without the raw data

The repository includes the processed stratum-level dataset (`output/dat_quality.rds`) and the protein-source shares (`output/protein_asf_props.rds`), so Scripts 6, 7 and 8 can be run directly after cloning to reproduce the figures and tables (a few minutes on a standard laptop). Two further sets of files will be archived with the stratum-level estimates in the Harvard Dataverse (doi:10.7910/DVN/EDULN9) at publication:

- `donor_inputs_protein.rds` and `donor_inputs_energy.rds`, written by Scripts 1 and 3 and needed for the donor bootstrap (place in `output/revision/`)
- `imputation_log_country_inputs_script3.csv`, `imputation_log_mder.csv` and `imputation_log_weights.csv`, needed for the table of imputed inputs (place in `output/revision/`)

Without them, Script 6 skips the donor bootstrap and the table of imputed inputs is incomplete.

## Software requirements

R packages: `dplyr`, `tidyr`, `readr`, `readxl`, `purrr`, `tibble`, `stringr`, `rlang`, `ggplot2`, `scales`, `patchwork`, `viridis`, `sf`, `rnaturalearth`, `rnaturalearthdata`, `countrycode`, `WDI`, `writexl`, `rstudioapi`, plus `nutriR` (https://github.com/cfree14/nutriR) and `wpp2024` (https://github.com/PPgp/wpp2024), which are installed from GitHub. The results use nutriR commit 49977ef; later versions add countries and change the donor pool, so install that commit: `remotes::install_github("cfree14/nutriR@49977ef35e08c5e846c52fcc847487a9cc8f35e8")`.

The package versions used for the published results are listed in `session_info.txt`. The analysis was run on macOS 26.6 (Apple silicon) with R 4.4.3; no non-standard hardware is needed. Installing the packages takes a few minutes, and running Scripts 1-9 takes about 5 minutes on a standard laptop.

## Archiving

The code is version-controlled on GitHub; the version used for the published results will be archived with a DOI at publication.

## Citation

If using this repository, please cite the associated manuscript:

De Nicola G, Beal T, Nordhagen S, Phillips SM, Passarelli S, Golden CD. *Model-based estimates of inadequate protein intake worldwide.* Under review.

Please also cite the relevant primary data sources listed above.

## Contact

Giacomo De Nicola, Harvard T.H. Chan School of Public Health, gdenicola@hsph.harvard.edu
