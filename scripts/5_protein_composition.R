# Script 5: Proportion of protein from different food types
# ---------------------------------------------------------
# Goal:
# - Derive representative protein density (g / 100 g, as eaten)
#   for each GDD food group using USDA FNDDS + WWEIA categories.

# Load packages -------------------------------------------------------------
library(dplyr)
library(readr)
library(nutriR)
library(ggplot2)
library(tidyr)
library(purrr)
library(countrycode)

# Setup ---------------------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
options(scipen = 999)
rm(list = ls())

# Read preprocessed GDD data (for later use in this script) -----------------
dat <- readRDS("./output/final_analysis_with_all_inadequacy.rds") %>%
  distinct(iso3, sex, age_group, .keep_all = TRUE) %>%
  filter(age_group != "0-0.99")

# Load FNDDS / WWEIA food composition data ---------------------------------
food_crosswalk     <- read_csv("data/USDA Food Composition Data (FNDDS)/food.csv")
nutrient_crosswalk <- read_csv("data/USDA Food Composition Data (FNDDS)/nutrient.csv")
nutrient_amounts   <- read_csv("data/USDA Food Composition Data (FNDDS)/food_nutrient.csv")
food_categories    <- read_csv("data/USDA Food Composition Data (FNDDS)/wweia_food_category.csv")

# Derive protein (g / 100 g) per FNDDS food --------------------------------
# Protein = nutrient_nbr 203. In food_nutrient, this corresponds to nutrient_id == 203.
# Units are g per 100 g edible portion (as consumed).

# sanity check (not printed): this will be a single row for Protein
protein_def <- nutrient_crosswalk %>%
  filter(nutrient_nbr == 203)

protein_amounts <- nutrient_amounts %>%
  filter(nutrient_id == 203) %>%
  select(fdc_id, protein_g_100g = amount)

# Attach protein to each FNDDS food ----------------------------------------
food_with_protein <- food_crosswalk %>%
  left_join(protein_amounts, by = "fdc_id")
# food_with_protein: one row per FNDDS food with protein_g_100g (may be NA for some)

# Attach WWEIA food category labels ----------------------------------------
foods_with_wweia <- food_with_protein %>%
  left_join(
    food_categories,
    by = c("food_category_id" = "wweia_food_category")
  )
# Now we have: description, protein_g_100g, and wweia_food_category_description per food.


# ---------------------------------------------------------------------------
# Map WWEIA food categories to GDD food groups
# ---------------------------------------------------------------------------
# Initial mapping: hand-crafted dictionary from WWEIA category descriptions
# to GDD categories. This is the conceptual backbone of the whole exercise.

wweia_to_gdd <- tribble(
  ~wweia_food_category_description,                  ~GDD_category,
  
  # -------- FRUITS --------
  "Bananas",                                         "Fruits",
  "Apples",                                          "Fruits",
  "Citrus fruits",                                   "Fruits",
  "Blueberries and other berries",                   "Fruits",
  "Strawberries",                                    "Fruits",
  "Melons",                                          "Fruits",
  "Grapes",                                          "Fruits",
  "Pears",                                           "Fruits",
  "Pineapple",                                       "Fruits",
  "Mango and papaya",                                "Fruits",
  "Peaches and nectarines",                          "Fruits",
  "Other fruits and fruit salads",                   "Fruits",
  "Dried fruits",                                    "Fruits",
  "Baby food: fruit",                                "Fruits",
  
  # -------- NON-STARCHY VEGETABLES --------
  "Broccoli",                                        "Non-starchy vegetables",
  "Spinach",                                         "Non-starchy vegetables",
  "Lettuce and lettuce salads",                      "Non-starchy vegetables",
  "Carrots",                                         "Non-starchy vegetables",
  "Tomatoes",                                        "Non-starchy vegetables",
  "Onions",                                          "Non-starchy vegetables",
  "String beans",                                    "Non-starchy vegetables",
  "Cabbage",                                         "Non-starchy vegetables",
  "Other red and orange vegetables",                 "Non-starchy vegetables",
  "Other dark green vegetables",                     "Non-starchy vegetables",
  "Other vegetables and combinations",               "Non-starchy vegetables",
  
  # -------- POTATOES --------
  "White potatoes, baked or boiled",                 "Potatoes",
  "Mashed potatoes and white potato mixtures",       "Potatoes",
  
  # -------- OTHER STARCHY VEGETABLES --------
  "Corn",                                            "Other starchy vegetables",
  "Other starchy vegetables",                        "Other starchy vegetables",
  "Plantains",                                       "Other starchy vegetables",
  "Yams",                                            "Other starchy vegetables",
  "Cassava",                                         "Other starchy vegetables",
  "Fried vegetables",                                "Other starchy vegetables", # will refine below
  
  # -------- LEGUMES --------
  "Beans, peas, legumes",                            "Beans and legumes",
  "Bean, pea, legume dishes",                        "Beans and legumes",
  
  # -------- NUTS & SEEDS --------
  "Nuts and seeds",                                  "Nuts and seeds",
  
  # -------- REFINED GRAINS (initial) --------
  "Yeast breads",                                    "Refined grains",
  "Rolls and buns",                                  "Refined grains",
  "Doughnuts, sweet rolls, pastries",                "Refined grains",
  "Bagels and English muffins",                      "Refined grains",
  "Crackers, excludes saltines",                     "Refined grains",
  "Saltine crackers",                                "Refined grains",
  "Pretzels/snack mix",                              "Refined grains",
  "Pancakes, waffles, French toast",                 "Refined grains",
  "Tortillas",                                       "Refined grains",
  
  # -------- WHOLE GRAINS (initial) --------
  "Oatmeal",                                         "Whole grains",
  "Grits and other cooked cereals",                  "Whole grains",
  "Ready-to-eat cereal, lower sugar (=<21.2g/100g)", "Whole grains",
  "Ready-to-eat cereal, higher sugar (>21.2g/100g)", "Whole grains",
  "Rice",                                            "Whole grains",       # refined vs whole handled in refinement
  "Pasta, noodles, cooked grains",                   "Whole grains",
  
  # -------- PROCESSED MEAT --------
  "Cold cuts and cured meats",                       "Total processed meats",
  "Frankfurters",                                    "Total processed meats",
  "Sausages",                                        "Total processed meats",
  "Bacon",                                           "Total processed meats",
  "Deli and cured meat sandwiches",                  "Total processed meats",
  
  # -------- UNPROCESSED RED MEAT --------
  "Beef, excludes ground",                           "Unprocessed red meats",
  "Ground beef",                                     "Unprocessed red meats",
  "Pork",                                            "Unprocessed red meats",
  "Lamb, goat, game",                                "Unprocessed red meats",
  
  # -------- SEAFOOD (initial) --------
  "Fish",                                            "Total seafood",
  "Shellfish",                                       "Total seafood",
  "Seafood mixed dishes",                            "Total seafood",
  "Seafood sandwiches",                              "Total seafood",
  
  # -------- EGGS --------
  "Eggs and omelets",                                "Eggs",
  
  # -------- CHEESE --------
  "Cheese",                                          "Cheese",
  "Cottage/ricotta cheese",                          "Cheese",
  
  # -------- YOGURT / FERMENTED MILK --------
  "Yogurt, regular",                                 "Yogurt including fermented milk",
  "Yogurt, Greek",                                   "Yogurt including fermented milk",
  "Plant-based yogurt",                              "Yogurt including fermented milk",
  
  # -------- MILK --------
  "Milk, whole",                                     "Total milk",
  "Milk, reduced fat",                               "Total milk",
  "Milk, lowfat",                                    "Total milk",
  "Milk, nonfat",                                    "Total milk",
  "Flavored milk, whole",                            "Total milk",
  "Flavored milk, reduced fat",                      "Total milk",
  "Flavored milk, lowfat",                           "Total milk",
  "Flavored milk, nonfat",                           "Total milk",
  
  # -------- JUICE --------
  "Citrus juice",                                    "Fruit juices",
  "Apple juice",                                     "Fruit juices",
  "Other fruit juice",                               "Fruit juices",
  
  # -------- SSB --------
  "Soft drinks",                                     "Sugar-sweetened beverages",
  "Fruit drinks",                                    "Sugar-sweetened beverages",
  "Sport and energy drinks",                         "Sugar-sweetened beverages",
  
  # -------- COFFEE & TEA --------
  "Coffee",                                          "Coffee",
  "Tea",                                             "Tea"
)

# Refine mapping for specific categories ------------------------------------
# - Other starchy vegetables: drop "Fried vegetables"
# - Refined grains: cooked refined grains + tortillas + yeast breads
# - Total seafood: fish + shellfish only (no sandwiches/mixed dishes)
# - Coffee: keep as consumed (no change needed)

wweia_to_gdd_refined <- wweia_to_gdd %>%
  mutate(
    GDD_category = case_when(
      # ---- OTHER STARCHY VEGETABLES ----
      wweia_food_category_description %in% c(
        "Corn",
        "Other starchy vegetables"
      ) ~ "Other starchy vegetables",
      wweia_food_category_description == "Fried vegetables" ~ NA_character_,
      
      # ---- REFINED GRAINS (A = 2) ----
      # cooked refined grains + tortillas + yeast breads
      wweia_food_category_description %in% c(
        "Pasta, noodles, cooked grains",
        "Rice",
        "Grits and other cooked cereals",
        "Tortillas",
        "Yeast breads"
      ) ~ "Refined grains",
      
      # remove clearly non-staple refined grain products from this group
      wweia_food_category_description %in% c(
        "Bagels and English muffins",
        "Crackers, excludes saltines",
        "Doughnuts, sweet rolls, pastries",
        "Pancakes, waffles, French toast",
        "Pretzels/snack mix",
        "Rolls and buns",
        "Saltine crackers"
      ) ~ NA_character_,
      
      # ---- TOTAL SEAFOOD ----
      # keep pure fish & shellfish, drop sandwiches & mixed dishes
      wweia_food_category_description %in% c("Fish", "Shellfish") ~ "Total seafood",
      wweia_food_category_description %in% c("Seafood mixed dishes",
                                             "Seafood sandwiches") ~ NA_character_,
      
      # ---- EVERYTHING ELSE: keep previous mapping ----
      TRUE ~ GDD_category
    )
  )

# Attach refined GDD categories to foods ------------------------------------
foods_with_gdd <- foods_with_wweia %>%
  left_join(
    wweia_to_gdd_refined,
    by = "wweia_food_category_description"
  )

# Compute median protein per 100 g for each GDD category --------------------
protein_per_gdd <- foods_with_gdd %>%
  filter(!is.na(GDD_category),
         !is.na(protein_g_100g)) %>%
  group_by(GDD_category) %>%
  summarise(
    n_foods             = n(),
    median_protein_100g = median(protein_g_100g, na.rm = TRUE),
    mean_protein_100g   = mean(protein_g_100g,   na.rm = TRUE),
    .groups             = "drop"
  ) %>%
  arrange(GDD_category)

# protein_per_gdd is the main output:
# one row per GDD_category with median and mean protein g/100 g as eaten.


# OPTIONAL: diagnostic object with foods closest to the median --------------
# (not printed; useful for sanity checks or supplementary tables)

closest_to_median <- foods_with_gdd %>%
  filter(!is.na(GDD_category),
         !is.na(protein_g_100g)) %>%
  group_by(GDD_category) %>%
  mutate(
    median_value = median(protein_g_100g, na.rm = TRUE),
    distance     = abs(protein_g_100g - median_value)
  ) %>%
  arrange(GDD_category, distance) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  select(
    GDD_category,
    description,
    wweia_food_category_description,
    protein_g_100g,
    median_value
  )
# You can inspect closest_to_median interactively if needed.




## now to load the GDD data for each category
# v01  Fruits  
# v02  Non-starchy vegetables  
# v03  Potatoes  
# v04  Other starchy vegetables  
# v05  Beans and legumes  
# v06  Nuts and seeds  
# v07  Refined grains  
# v08  Whole grains  
# v09  Total processed meats  
# v10  Unprocessed red meats  
# v11  Total seafood  
# v12  Eggs  
# v13  Cheese  
# v14  Yogurt (including fermented milk)  
# v15  Sugar-sweetened beverages  
# v16  Fruit juices  
# v17  Coffee  
# v18  Tea  
# v57  Total milk


gdd_path <- "./data/GDD_FinalEstimates_01102022/Country-level estimates/"

# which variables we want
target_vars <- c(sprintf("v%02d", 1:18), "v57")

# build full paths (e.g. v01_cnty.csv, v02_cnty.csv, ..., v57_cnty.csv)
gdd_files <- tibble(
  var_code = target_vars,
  file     = file.path(gdd_path, paste0(var_code, "_cnty.csv"))
)

# load them ALL as-is, no filtering yet
gdd_raw <- gdd_files %>%
  mutate(
    data = map(file, ~ read_csv(.x))  # keep full structure for now
  ) %>%
  { set_names(.$data, .$var_code) }

# now you have:
# gdd_raw[["v01"]]
# gdd_raw[["v02"]]
# ...
# gdd_raw[["v18"]]
# gdd_raw[["v57"]]


## 1) Map v-codes to GDD category names --------------------------------------
var_labels <- tribble(
  ~var_code, ~GDD_category,
  "v01", "Fruits",
  "v02", "Non-starchy vegetables",
  "v03", "Potatoes",
  "v04", "Other starchy vegetables",
  "v05", "Beans and legumes",
  "v06", "Nuts and seeds",
  "v07", "Refined grains",
  "v08", "Whole grains",
  "v09", "Total processed meats",
  "v10", "Unprocessed red meats",
  "v11", "Total seafood",
  "v12", "Eggs",
  "v13", "Cheese",
  "v14", "Yogurt including fermented milk",
  "v15", "Sugar-sweetened beverages",
  "v16", "Fruit juices",
  "v17", "Coffee",
  "v18", "Tea",
  "v57", "Total milk"
)

## 2) Helper to extract national, all-age, all-sex, all-edu/urban rows ------
extract_nat_all_2018 <- function(df) {
  df %>%
    filter(
      year   == 2018,
      urban  == 999,
      edu    == 999,
      female == 999,
      age    == 999
    ) %>%
    select(iso3, year, median)
}

## 3) Stack all v's long, then pivot wide by GDD category --------------------
gdd_long_2018 <- imap_dfr(
  gdd_raw,
  ~ extract_nat_all_2018(.x) %>%
    mutate(var_code = .y)  # .y is e.g. "v01"
)

gdd_intake_2018_wide <- gdd_long_2018 %>%
  left_join(var_labels, by = "var_code") %>%
  select(iso3, year, GDD_category, median) %>%
  pivot_wider(
    names_from  = GDD_category,
    values_from = median
  )

# Result:
# gdd_intake_2018_wide = one row per country (iso3),
#                        columns = grams/day (median) for each GDD food group.

## --------------------------------------------------------------------
## Compute country-level protein from ASFs and seafood (proportions)
## Requires:
##   - protein_per_gdd: GDD_category, median_protein_100g
##   - gdd_intake_2018_wide: iso3, year, GDD_category columns (grams/day)
## --------------------------------------------------------------------


# 1) Define ASF and seafood groups ------------------------------------
asf_groups <- c(
  "Total processed meats",
  "Unprocessed red meats",
  "Total seafood",
  "Eggs",
  "Cheese",
  "Total milk",
  "Yogurt including fermented milk"
)

seafood_group <- "Total seafood"

# 2) Protein density lookup (g protein per g food) --------------------
protein_density_lookup <- protein_per_gdd %>%
  mutate(protein_per_g = median_protein_100g / 100) %>%
  select(GDD_category, protein_per_g)

# 3) Merge intake (g/day) with protein density ------------------------
gdd_protein_long <- gdd_intake_2018_wide %>%
  pivot_longer(
    cols      = -c(iso3, year),
    names_to  = "GDD_category",
    values_to = "intake_g_day"
  ) %>%
  left_join(protein_density_lookup, by = "GDD_category") %>%
  mutate(
    protein_g_day = intake_g_day * protein_per_g
  )

# 4) Total protein per country ----------------------------------------
protein_totals <- gdd_protein_long %>%
  group_by(iso3, year) %>%
  summarise(
    total_protein = sum(protein_g_day, na.rm = TRUE),
    .groups = "drop"
  )

# 5) Protein from animal source foods (ASFs) --------------------------
protein_asf <- gdd_protein_long %>%
  filter(GDD_category %in% asf_groups) %>%
  group_by(iso3, year) %>%
  summarise(
    asf_protein = sum(protein_g_day, na.rm = TRUE),
    .groups = "drop"
  )

# 6) Protein from seafood only ----------------------------------------
protein_seafood <- gdd_protein_long %>%
  filter(GDD_category == seafood_group) %>%
  group_by(iso3, year) %>%
  summarise(
    seafood_protein = sum(protein_g_day, na.rm = TRUE),
    .groups = "drop"
  )

# 7) Combine and compute proportions ----------------------------------
protein_props <- protein_totals %>%
  left_join(protein_asf,     by = c("iso3", "year")) %>%
  left_join(protein_seafood, by = c("iso3", "year")) %>%
  mutate(
    prop_asf     = asf_protein     / total_protein,
    prop_seafood = seafood_protein / total_protein
  )

# protein_props:
#   iso3, year, total_protein (g/d),
#   asf_protein (g/d), seafood_protein (g/d),
#   prop_asf (share of protein from ASFs),
#   prop_seafood (share of protein from seafood).




# Save ASF / seafood protein proportions for use in Script 6
saveRDS(protein_props, "./output/protein_asf_props.rds")


protein_props_named <- protein_props %>%
  mutate(
    country = countrycode(iso3, origin = "iso3c", destination = "country.name")
  ) %>%
  select(country, everything()) %>%
  arrange(desc(prop_asf))   # <-- sort by ASF share, highest first


# library(writexl)
# 
# write_xlsx(
#   protein_props_named,
#   "./output/protein_asf_props.xlsx"
# )
