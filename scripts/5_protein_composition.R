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
