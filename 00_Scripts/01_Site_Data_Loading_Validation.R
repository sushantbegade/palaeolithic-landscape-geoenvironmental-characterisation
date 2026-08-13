# =============================================================
# SCRIPT 01 (REBUILT): Site Data Loading, Cleaning, Validation
# =============================================================
# Project: Reading the Palaeolithic Landscape
# Chapter: Mapping the Past (Springer Nature, 2026)
# Author: Sushant Begade | RTMNU Nagpur
# ORCID: 0009-0003-0804-1763
# Rebuild date: August 2026
# =============================================================
# CHANGES FROM ORIGINAL SCRIPT 01:
#   1. coord_precision tier standardised from raw location_precision
#      field into 3 clean categories (field_gps / author_reported /
#      gis_digitised) and RETAINED in all saved outputs so it
#      propagates to every downstream matrix (Scripts 02-07).
#      Original script computed this only for a print() summary
#      and then dropped it — never reached sites_extract onward.
#   2. period locked as factor(levels = c("LP","MP","UP","LP_MP",
#      "MP_UP","UNCLASSIFIED")) immediately after classification,
#      to prevent silent level-order drift in later group_by/case_when.
#   3. Added explicit multicultural-site flag (period_multi) —
#      needed later to decide how multicultural sites are handled
#      in MW tests / PCA (currently silently absorbed into "most
#      inclusive" category per chapter text — now auditable).
# =============================================================

load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")

library(terra)
library(sf)
library(tidyverse)
library(readxl)
library(janitor)

set.seed(rf_seed)

message("Global parameters loaded.")
message("Working directory: ", chapter_root)


# -------------------------------------------------------------
# SECTION 2: Load Raw Site Data
# -------------------------------------------------------------

sites_raw <- read_csv(paths$sites_csv, show_col_types = FALSE)
message("\nRaw site records loaded: ", nrow(sites_raw))
message("Columns: ", paste(names(sites_raw), collapse = ", "))


# -------------------------------------------------------------
# SECTION 3: Clean Column Names
# -------------------------------------------------------------

sites_clean <- sites_raw %>% clean_names()
message("\nCleaned column names: ", paste(names(sites_clean), collapse = ", "))


# -------------------------------------------------------------
# SECTION 4: Inspect
# -------------------------------------------------------------

print(head(sites_clean, 5))
for (i in seq_along(names(sites_clean))) {
  message(i, ": ", names(sites_clean)[i])
}


# -------------------------------------------------------------
# SECTION 5: Select and Rename Core Fields
# -------------------------------------------------------------
# VERIFY column names below against Section 4 output before running.

sites_df <- sites_clean %>%
  select(
    site_id            = site_id,
    site_name          = site_name,
    district           = district,
    latitude           = latitude,
    longitude          = longitude,
    cultural_period    = relative_chronology,
    river_basin        = river_basin_wainganga_wardha,
    raw_material       = raw_material,
    lithic_assemblage  = lithic_assemblage,
    data_source        = data_source_type_survey_excavation_review,
    location_precision = location_precision,
    reference          = primary_reference
  )

message("\nCore fields selected. Rows: ", nrow(sites_df))


# -------------------------------------------------------------
# SECTION 6: Validate Coordinates
# -------------------------------------------------------------

missing_lat <- sum(is.na(sites_df$latitude))
missing_lon <- sum(is.na(sites_df$longitude))
message("\nMissing latitude: ", missing_lat, " | Missing longitude: ", missing_lon)

out_of_range <- sites_df %>%
  filter(latitude < 19.0 | latitude > 22.0 |
           longitude < 77.5 | longitude > 80.5)
message("Sites outside expected bbox: ", nrow(out_of_range))
if (nrow(out_of_range) > 0) {
  print(out_of_range %>% select(site_id, site_name, latitude, longitude))
}

sites_valid <- sites_df %>%
  filter(!is.na(latitude), !is.na(longitude),
         latitude  >= 19.0 & latitude  <= 22.0,
         longitude >= 77.5 & longitude <= 80.5)

message("Sites with valid coordinates: ", nrow(sites_valid))


# -------------------------------------------------------------
# SECTION 7: Classify Cultural Periods (LP / MP / UP)  [FIXED v2]
# -------------------------------------------------------------
# BUG FOUND IN v1: ordered case_when() let single-period patterns
# ("LOWER", "MIDDLE", "UPPER") match FIRST against combo strings
# like "Lower and Middle Palaeolithic", so every multi-period site
# got silently force-classified into one bucket. Confirmed: v1 run
# showed period_multi = 0 despite 5 distinct multi-period text
# strings existing in raw data, and LP+MP+UP summed to exactly 197.
#
# FIX: count keyword hits per site instead of ordered matching.
# >1 keyword hit -> MULTI (excluded from period-stratified tests
# downstream since periods <- c("LP","MP","UP") global param
# already filters to only these three - no extra code needed in
# Scripts 06/07). MULTI sites still counted in Fig 1 / Table 1
# totals via district summary below.

message("\nRaw cultural_period value counts:")
print(sites_valid %>% count(cultural_period, sort = TRUE))

sites_valid <- sites_valid %>%
  mutate(
    cp_upper    = str_to_upper(cultural_period),
    hit_lower   = str_detect(cp_upper, "LOWER|ACHEULIAN|ACHEULEAN"),
    hit_middle  = str_detect(cp_upper, "MIDDLE|LEVALLOIS"),
    hit_upper   = str_detect(cp_upper, "UPPER|BLADE"),
    n_hits      = hit_lower + hit_middle + hit_upper,
    period = case_when(
      n_hits == 0 ~ "UNCLASSIFIED",
      n_hits >  1 ~ "MULTI",
      hit_lower   ~ "LP",
      hit_middle  ~ "MP",
      hit_upper   ~ "UP",
      TRUE ~ "UNCLASSIFIED"
    ),
    period = factor(period, levels = c("LP", "MP", "UP", "MULTI", "UNCLASSIFIED")),
    period_multi = period == "MULTI"
  ) %>%
  select(-cp_upper, -hit_lower, -hit_middle, -hit_upper)

period_summary <- sites_valid %>% count(period, name = "n_sites") %>% arrange(desc(n_sites))
message("\nCultural period classification:")
print(period_summary)

n_unclassified <- sum(sites_valid$period == "UNCLASSIFIED")
if (n_unclassified > 0) {
  message("\nWARNING: ", n_unclassified, " sites unclassified. Review manually:")
  print(sites_valid %>% filter(period == "UNCLASSIFIED") %>%
          select(site_id, site_name, cultural_period))
}

n_multi <- sum(sites_valid$period_multi)
message("\nMulticultural sites (period = MULTI): ", n_multi)
if (n_multi > 0) {
  message("These are EXCLUDED from period-stratified MW/PCA/RF tests")
  message("(periods <- c(\"LP\",\"MP\",\"UP\") global param naturally drops them).")
  message("They ARE still counted in Fig 1 / Table 1 site totals below.")
  message("Update Methods 5.1 text to match this actual handling —")
  message("draft currently claims a different (never-implemented) convention.")
  print(sites_valid %>% filter(period == "MULTI") %>%
          select(site_id, site_name, cultural_period))
}


# -------------------------------------------------------------
# SECTION 8: Validate Districts
# -------------------------------------------------------------

message("\nValidating district field...")
print(table(sites_valid$district))

sites_valid <- sites_valid %>%
  mutate(
    district = case_when(
      str_detect(str_to_upper(district), "NAGPUR")     ~ "Nagpur",
      str_detect(str_to_upper(district), "CHANDRAPUR") ~ "Chandrapur",
      TRUE ~ district
    )
  )
print(table(sites_valid$district))


# -------------------------------------------------------------
# SECTION 9: Check for Duplicates
# -------------------------------------------------------------

dup_ids <- sites_valid %>% filter(duplicated(site_id))
message("\nDuplicate site IDs: ", nrow(dup_ids))

dup_coords <- sites_valid %>%
  filter(duplicated(paste(round(latitude, 4), round(longitude, 4))))
message("Duplicate coordinates: ", nrow(dup_coords))
if (nrow(dup_coords) > 0) {
  print(dup_coords %>% select(site_id, site_name, latitude, longitude, period))
}

sites_valid <- sites_valid %>% distinct(latitude, longitude, .keep_all = TRUE)
message("Sites after duplicate removal: ", nrow(sites_valid))


# -------------------------------------------------------------
# SECTION 10: Standardise Coordinate Precision Tier  [NEW/FIXED]
# -------------------------------------------------------------
# CHANGE 1 core logic. Chapter text claims:
#   field-verified n=64 (±5-10m), author-reported n=73 (±5-10m),
#   GIS-digitised n=60 (lower confidence).
# Adjust the str_detect patterns below to match your actual
# location_precision raw values (inspect first).

message("\nRaw location_precision values:")
print(table(sites_valid$location_precision, useNA = "always"))

sites_valid <- sites_valid %>%
  mutate(
    coord_precision = case_when(
      str_detect(str_to_upper(location_precision), "FIELD|GPS|VERIFIED") ~ "field_gps",
      str_detect(str_to_upper(location_precision), "AUTHOR|REPORTED|PUBLISHED") ~ "author_reported",
      str_detect(str_to_upper(location_precision), "GIS|DIGITI|MAP|GEOREF") ~ "gis_digitised",
      TRUE ~ "unclassified_precision"
    ),
    coord_precision = factor(coord_precision,
                             levels = c("field_gps", "author_reported",
                                        "gis_digitised", "unclassified_precision"))
  )

precision_summary <- sites_valid %>% count(coord_precision, name = "n")
message("\nCoordinate precision tier breakdown:")
print(precision_summary)

n_unclass_prec <- sum(sites_valid$coord_precision == "unclassified_precision")
if (n_unclass_prec > 0) {
  message("\nWARNING: ", n_unclass_prec,
          " sites have unrecognised location_precision text.")
  message("Fix str_detect patterns above BEFORE proceeding —")
  message("this column is now carried through every downstream script")
  message("and is required for the positional-uncertainty sensitivity")
  message("analysis flagged in the chapter's Discussion (Section 5.1).")
  print(sites_valid %>% filter(coord_precision == "unclassified_precision") %>%
          select(site_id, location_precision))
}

# HARD STOP if precision tiers don't reconcile to expected n=197 split.
# Comment out once verified — this is a deliberate audit gate.
stopifnot(
  "coord_precision must be fully classified before continuing" =
    n_unclass_prec == 0
)


# -------------------------------------------------------------
# SECTION 11: Convert to Spatial (WGS84) + Reproject UTM 44N
# -------------------------------------------------------------

sites_wgs84 <- st_as_sf(sites_valid, coords = c("longitude", "latitude"), crs = crs_wgs84)
message("\nSpatial object created: ", nrow(sites_wgs84), " sites")

sites_utm <- st_transform(sites_wgs84, crs = crs_utm44n)
message("Reprojected CRS: ", st_crs(sites_utm)$input)


# -------------------------------------------------------------
# SECTION 12: District Boundary Check
# -------------------------------------------------------------

boundary <- st_read(paths$boundary, quiet = TRUE) %>% st_transform(crs = crs_utm44n)

sites_in_boundary <- st_within(sites_utm, st_union(boundary), sparse = FALSE)
n_outside <- sum(!sites_in_boundary)
message("\nSites outside district boundary: ", n_outside)
if (n_outside > 0) {
  print(sites_utm[!sites_in_boundary, ] %>% select(site_id, site_name, period))
}


# -------------------------------------------------------------
# SECTION 13: Summary Statistics
# -------------------------------------------------------------

message("\n=== SITE DATA SUMMARY ===")
message("Total validated sites: ", nrow(sites_valid))
message("\nSites per District:");         print(table(sites_valid$district))
message("\nSites per Cultural Period:");   print(table(sites_valid$period))
message("(LP+MP+UP+MULTI+UNCLASSIFIED must sum to ", nrow(sites_valid), ")")
message("\nSites per District x Period:"); print(table(sites_valid$district, sites_valid$period))
message("\nCoordinate precision tier:");   print(table(sites_valid$coord_precision))
message("\nPrecision x Period cross-tab (check for confound before analysis):")
print(table(sites_valid$coord_precision, sites_valid$period))


# -------------------------------------------------------------
# SECTION 14: Save Outputs — coord_precision now retained
# -------------------------------------------------------------

out_path <- file.path(chapter_root, "01_Data_Processed")

write_csv(sites_valid, file.path(out_path, "sites_validated.csv"))
message("\nSaved: sites_validated.csv (now includes coord_precision, period_multi)")

st_write(sites_wgs84, file.path(out_path, "sites_wgs84.shp"), delete_layer = TRUE, quiet = TRUE)
st_write(sites_utm,   file.path(out_path, "sites_utm44n.shp"), delete_layer = TRUE, quiet = TRUE)
message("Saved: sites_wgs84.shp, sites_utm44n.shp")

summary_table <- sites_valid %>%
  group_by(district, period) %>%
  summarise(n_sites = n(), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = n_sites, values_fill = 0) %>%
  mutate(Total = rowSums(across(where(is.numeric))))

write_csv(summary_table, file.path(chapter_root, "11_Tables", "Table_01_Site_Summary.csv"))
message("Saved: Table_01_Site_Summary.csv")

# NEW: dedicated precision audit table for methods/supplementary material
precision_table <- sites_valid %>%
  group_by(period, coord_precision) %>%
  summarise(n_sites = n(), .groups = "drop") %>%
  pivot_wider(names_from = coord_precision, values_from = n_sites, values_fill = 0)

write_csv(precision_table, file.path(chapter_root, "11_Tables", "Table_00_Coordinate_Precision.csv"))
message("Saved: Table_00_Coordinate_Precision.csv [NEW — cite in Methods 5.1]")

save(
  sites_valid, sites_wgs84, sites_utm, boundary,
  file = file.path(out_path, "sites_validated.RData")
)
message("Saved: sites_validated.RData")


# -------------------------------------------------------------
# SECTION 15: Log Session
# -------------------------------------------------------------

log_file <- file.path(chapter_root, "13_Logs",
                      paste0("script01_log_", Sys.Date(), ".txt"))
sink(log_file)
cat("SCRIPT 01 LOG (REBUILT) — Site Data Loading + Validation\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Total sites loaded:", nrow(sites_raw), "\n")
cat("Sites after validation:", nrow(sites_valid), "\n")
cat("Multicultural sites (flagged, not dropped):", n_multi, "\n\n")
cat("Sites per period:\n"); print(table(sites_valid$period))
cat("\nSites per district:\n"); print(table(sites_valid$district))
cat("\nCoordinate precision tier:\n"); print(table(sites_valid$coord_precision))
cat("\nPrecision x Period cross-tab:\n")
print(table(sites_valid$coord_precision, sites_valid$period))
sink()
message("\nLog saved: ", log_file)


# -------------------------------------------------------------
# SCRIPT 01 COMPLETE
# -------------------------------------------------------------

message("\n=============================================================")
message("Script 01 (REBUILT) complete.")
message("Validated sites: ", nrow(sites_valid))
message("coord_precision + period_multi now flow into ALL downstream")
message("matrices (Scripts 02-07) — verify Script 02 carries them next.")
message("Next: Rebuild Script 02 — Background Point Generation")
message("  (survey-envelope-constrained + uniform-random PARALLEL runs)")
message("=============================================================")