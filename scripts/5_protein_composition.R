#Script 5: Proportion of protein from different food types

# Load required packages
library(dplyr)
library(readr)
library(nutriR)
library(ggplot2)
library(tidyr)
library(purrr)
library(countrycode)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
options(scipen=999)
rm(list = ls())

library(dplyr)
library(readr)


#read in preprocessed GDD data
dat <- readRDS("./output/final_analysis_with_all_inadequacy.rds") %>%
  distinct(iso3, sex, age_group, .keep_all = TRUE) %>%
  filter(age_group != "0-0.99")


#load food composition data
food_crosswalk      <- read_csv("data/USDA Food Composition Data (FNDDS)/food.csv")
nutrient_crosswalk  <- read_csv("data/USDA Food Composition Data (FNDDS)/nutrient.csv")
nutrient_amounts    <- read_csv("data/USDA Food Composition Data (FNDDS)/food_nutrient.csv")
food_categories <- read_csv("data/USDA Food Composition Data (FNDDS)/wweia_food_category.csv")

# sanity check: protein definition
nutrient_crosswalk %>% 
  filter(nutrient_nbr == 203)
# id = 1003, name = "Protein", unit_name = "G"

# 1) protein rows in food_nutrient: nutrient_id == 203
protein_amounts <- nutrient_amounts %>%
  filter(nutrient_id == 203) %>%          # <-- use 203, not 1003
  select(fdc_id, protein_g_100g = amount)

# 2) join to foods
food_with_protein <- food_crosswalk %>%
  left_join(protein_amounts, by = "fdc_id")

##food_with_protein now contains all foods
##with the corresponding amount of protein per 100g (as eaten)


#we need to select top foods for each of the gdd cateogories, to then create a 
#value of protein per 100g for each category.


## 1) Attach WWEIA food category labels to each food -------------------------

foods_with_wweia <- food_with_protein %>%
  left_join(
    food_categories,
    by = c("food_category_id" = "wweia_food_category")
  )



# Create an empty tibble of all WWEIA categories
wweia_list <- foods_with_wweia %>%
  distinct(wweia_food_category_description) %>%
  arrange(wweia_food_category_description)

# Define an initial mapping (partial, editable)
# Fill in only the ones you know for sure; others remain NA
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
  # Note: fried potatoes may be excluded depending on your rule
  
  # -------- OTHER STARCHY VEGETABLES --------
  "Corn",                                            "Other starchy vegetables",
  "Other starchy vegetables",                        "Other starchy vegetables",
  "Plantains",                                       "Other starchy vegetables",
  "Yams",                                            "Other starchy vegetables",
  "Cassava",                                         "Other starchy vegetables",
  "Fried vegetables",                                "Other starchy vegetables", # optional
  
  # -------- LEGUMES --------
  "Beans, peas, legumes",                             "Beans and legumes",
  "Bean, pea, legume dishes",                         "Beans and legumes",
  
  # -------- NUTS & SEEDS --------
  "Nuts and seeds",                                  "Nuts and seeds",
  
  # -------- REFINED GRAINS --------
  "Yeast breads",                                    "Refined grains",
  "Rolls and buns",                                  "Refined grains",
  "Doughnuts, sweet rolls, pastries",                "Refined grains",
  "Bagels and English muffins",                      "Refined grains",
  "Crackers, excludes saltines",                     "Refined grains",
  "Saltine crackers",                                "Refined grains",
  "Pretzels/snack mix",                              "Refined grains",
  "Pancakes, waffles, French toast",                 "Refined grains",
  "Tortillas",                                       "Refined grains",
  
  # -------- WHOLE GRAINS --------
  "Oatmeal",                                         "Whole grains",
  "Grits and other cooked cereals",                  "Whole grains",
  "Ready-to-eat cereal, lower sugar (=<21.2g/100g)", "Whole grains",
  "Ready-to-eat cereal, higher sugar (>21.2g/100g)", "Whole grains",
  "Rice",                                            "Whole grains",       # verify white vs brown
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
  
  # -------- SEAFOOD --------
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
  
  # -------- COFFEE --------
  "Coffee",                                          "Coffee",
  
  # -------- TEA --------
  "Tea",                                             "Tea"
)


foods_with_gdd <- foods_with_wweia %>%
  left_join(wweia_to_gdd, by = "wweia_food_category_description")


# 4) Compute the MEDIAN protein per 100 g for each GDD category
protein_per_gdd <- foods_with_gdd %>%
  filter(!is.na(GDD_category),
         !is.na(protein_g_100g)) %>%
  group_by(GDD_category) %>%
  summarise(
    n_foods = n(),
    median_protein_100g = median(protein_g_100g, na.rm = TRUE),
    mean_protein_100g   = mean(protein_g_100g, na.rm = TRUE),  # optional
    .groups = "drop"
  ) %>%
  arrange(GDD_category)

# inspect
protein_per_gdd



# refine mapping: adjust 3 groups (other starchy veg, refined grains, seafood)
wweia_to_gdd_refined <- wweia_to_gdd %>%
  mutate(
    GDD_category = case_when(
      # ---- OTHER STARCHY VEGETABLES ----
      wweia_food_category_description %in% c(
        "Corn",
        "Other starchy vegetables"
      ) ~ "Other starchy vegetables",
      wweia_food_category_description == "Fried vegetables" ~ NA_character_,  # drop fried veg from this group
      
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
      wweia_food_category_description %in% c("Seafood mixed dishes", "Seafood sandwiches") ~ NA_character_,
      
      # ---- EVERYTHING ELSE: keep previous mapping ----
      TRUE ~ GDD_category
    )
  )


# function: get foods closest to median per GDD category
closest_to_median <- foods_with_gdd %>%
  filter(!is.na(GDD_category), !is.na(protein_g_100g)) %>%
  group_by(GDD_category) %>%
  mutate(
    median_value = median(protein_g_100g, na.rm = TRUE),
    distance = abs(protein_g_100g - median_value)
  ) %>%
  arrange(GDD_category, distance) %>%
  # keep the top 5 closest foods to the median for inspection
  slice_head(n = 5) %>%
  ungroup() %>%
  select(
    GDD_category,
    description,
    wweia_food_category_description,
    protein_g_100g,
    median_value
  )


# rebuild foods_with_gdd using refined mapping
foods_with_gdd <- foods_with_wweia %>%
  left_join(
    wweia_to_gdd_refined,
    by = "wweia_food_category_description"
  )

# recompute median & mean protein per GDD category
protein_per_gdd <- foods_with_gdd %>%
  filter(!is.na(GDD_category),
         !is.na(protein_g_100g)) %>%
  group_by(GDD_category) %>%
  summarise(
    n_foods = n(),
    median_protein_100g = median(protein_g_100g, na.rm = TRUE),
    mean_protein_100g   = mean(protein_g_100g,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(GDD_category)

protein_per_gdd


closest_to_median <- foods_with_gdd %>%
  filter(!is.na(GDD_category),
         !is.na(protein_g_100g)) %>%
  group_by(GDD_category) %>%
  mutate(
    median_value = median(protein_g_100g, na.rm = TRUE),
    distance = abs(protein_g_100g - median_value)
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



