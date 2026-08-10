# =============================================================
# SCRIPT 01: Site Data Loading, Cleaning, and Validation
# =============================================================
# Project: Reading the Palaeolithic Landscape
# Chapter: Mapping the Past (Springer Nature, 2026)
# Author: Sushant Begade | RTMNU Nagpur
# ORCID: 0009-0003-0804-1763
# Date: August 2026
# =============================================================
# PURPOSE: Load 197 Palaeolithic site records, clean and
# validate all fields, classify cultural periods (LP/MP/UP),
# reproject to UTM 44N, generate summary statistics,
# and save validated spatial + tabular outputs.
# =============================================================


# -------------------------------------------------------------
# SECTION 1: Load Global Parameters
# -------------------------------------------------------------

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

message("\nLoading site data...")

# Load from CSV (primary)
sites_raw <- read_csv(paths$sites_csv, show_col_types = FALSE)

message("Raw site records loaded: ", nrow(sites_raw))
message("Columns: ", paste(names(sites_raw), collapse = ", "))


# -------------------------------------------------------------
# SECTION 3: Clean Column Names
# -------------------------------------------------------------

sites_clean <- sites_raw %>%
  clean_names()  # janitor: lowercase, underscores, no special chars

message("\nCleaned column names: ", paste(names(sites_clean), collapse = ", "))


# -------------------------------------------------------------
# SECTION 4: Identify Key Columns
# -------------------------------------------------------------

# Print first 5 rows to inspect structure
print(head(sites_clean, 5))

# Print column names with indices for reference
for (i in seq_along(names(sites_clean))) {
  message(i, ": ", names(sites_clean)[i])
}


# -------------------------------------------------------------
# SECTION 5: Select and Rename Core Fields
# -------------------------------------------------------------

# IMPORTANT: Adjust column names below based on Section 4 output
# These are expected column names after clean_names()
# Verify against your actual output and edit if needed

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

message("\nValidating coordinates...")

# Check for missing
missing_lat <- sum(is.na(sites_df$latitude))
missing_lon <- sum(is.na(sites_df$longitude))
message("Missing latitude: ", missing_lat)
message("Missing longitude: ", missing_lon)

# Check coordinate ranges (Nagpur-Chandrapur approx bbox)
# Lat: 19.0 - 22.0 | Lon: 77.5 - 80.5
out_of_range <- sites_df %>%
  filter(
    latitude  < 19.0 | latitude  > 22.0 |
      longitude < 77.5 | longitude > 80.5
  )

message("Sites outside expected bbox: ", nrow(out_of_range))
if (nrow(out_of_range) > 0) {
  print(out_of_range %>% select(site_id, site_name, latitude, longitude))
}

# Remove sites with missing or out-of-range coordinates
sites_valid <- sites_df %>%
  filter(
    !is.na(latitude),
    !is.na(longitude),
    latitude  >= 19.0 & latitude  <= 22.0,
    longitude >= 77.5 & longitude <= 80.5
  )

message("Sites with valid coordinates: ", nrow(sites_valid))


# -------------------------------------------------------------
# SECTION 7: Classify Cultural Periods (LP / MP / UP)
# -------------------------------------------------------------

message("\nClassifying cultural periods...")

# Inspect unique values in cultural_period column
message("Unique cultural period values:")
print(unique(sites_valid$cultural_period))

# Standardise to LP / MP / UP
# Adjust patterns based on actual values in your data
sites_valid <- sites_valid %>%
  mutate(
    period = case_when(
      str_detect(str_to_upper(cultural_period),
                 "LOWER|ACHEULIAN|ACHEULEAN|LP") ~ "LP",
      str_detect(str_to_upper(cultural_period),
                 "MIDDLE|LEVALLOIS|MP") ~ "MP",
      str_detect(str_to_upper(cultural_period),
                 "UPPER|BLADE|UP") ~ "UP",
      str_detect(str_to_upper(cultural_period),
                 "MIDDLE.*UPPER|MP.*UP|UP.*MP") ~ "MP_UP",
      str_detect(str_to_upper(cultural_period),
                 "LOWER.*MIDDLE|LP.*MP|MP.*LP") ~ "LP_MP",
      TRUE ~ "UNCLASSIFIED"
    )
  )

# Period summary
period_summary <- sites_valid %>%
  count(period, name = "n_sites") %>%
  arrange(desc(n_sites))

message("\nCultural period classification:")
print(period_summary)

# Flag unclassified
n_unclassified <- sum(sites_valid$period == "UNCLASSIFIED")
if (n_unclassified > 0) {
  message("\nWARNING: ", n_unclassified, " sites unclassified.")
  message("Review these sites manually:")
  print(sites_valid %>%
          filter(period == "UNCLASSIFIED") %>%
          select(site_id, site_name, cultural_period))
}


# -------------------------------------------------------------
# SECTION 8: Validate Districts
# -------------------------------------------------------------

message("\nValidating district field...")
print(table(sites_valid$district))

# Standardise district names
sites_valid <- sites_valid %>%
  mutate(
    district = case_when(
      str_detect(str_to_upper(district), "NAGPUR")     ~ "Nagpur",
      str_detect(str_to_upper(district), "CHANDRAPUR") ~ "Chandrapur",
      TRUE ~ district
    )
  )

message("District breakdown after standardisation:")
print(table(sites_valid$district))


# -------------------------------------------------------------
# SECTION 9: Check for Duplicates
# -------------------------------------------------------------

message("\nChecking for duplicate sites...")

# Duplicate site IDs
dup_ids <- sites_valid %>%
  filter(duplicated(site_id))
message("Duplicate site IDs: ", nrow(dup_ids))

# Spatially duplicate (same lat/lon)
dup_coords <- sites_valid %>%
  filter(duplicated(paste(round(latitude, 4), round(longitude, 4))))
message("Duplicate coordinates: ", nrow(dup_coords))

if (nrow(dup_coords) > 0) {
  message("Review duplicate coordinate sites:")
  print(dup_coords %>% select(site_id, site_name, latitude, longitude, period))
}

# Remove exact coordinate duplicates (keep first)
sites_valid <- sites_valid %>%
  distinct(latitude, longitude, .keep_all = TRUE)

message("Sites after duplicate removal: ", nrow(sites_valid))


# -------------------------------------------------------------
# SECTION 10: Convert to Spatial Object (WGS84)
# -------------------------------------------------------------

message("\nConverting to spatial object (WGS84)...")

sites_wgs84 <- st_as_sf(
  sites_valid,
  coords = c("longitude", "latitude"),
  crs    = crs_wgs84
)

message("Spatial object created: ", nrow(sites_wgs84), " sites")
message("CRS: ", st_crs(sites_wgs84)$input)


# -------------------------------------------------------------
# SECTION 11: Reproject to UTM Zone 44N
# -------------------------------------------------------------

message("\nReprojecting to UTM Zone 44N...")

sites_utm <- st_transform(sites_wgs84, crs = crs_utm44n)

message("Reprojected CRS: ", st_crs(sites_utm)$input)


# -------------------------------------------------------------
# SECTION 12: Load and Validate District Boundary
# -------------------------------------------------------------

message("\nLoading district boundary...")

boundary <- st_read(paths$boundary, quiet = TRUE) %>%
  st_transform(crs = crs_utm44n)

message("Boundary CRS: ", st_crs(boundary)$input)

# Spatial join — confirm sites fall within boundary
sites_in_boundary <- st_within(sites_utm, st_union(boundary),
                               sparse = FALSE)
n_outside <- sum(!sites_in_boundary)
message("Sites outside district boundary: ", n_outside)

if (n_outside > 0) {
  message("Review sites outside boundary:")
  print(sites_utm[!sites_in_boundary, ] %>%
          select(site_id, site_name, period))
}


# -------------------------------------------------------------
# SECTION 13: Generate Summary Statistics
# -------------------------------------------------------------

message("\n=== SITE DATA SUMMARY ===")

message("\nTotal validated sites: ", nrow(sites_valid))

message("\nSites per District:")
print(table(sites_valid$district))

message("\nSites per Cultural Period:")
print(table(sites_valid$period))

message("\nSites per District × Period:")
print(table(sites_valid$district, sites_valid$period))

message("\nLocation Precision breakdown:")
print(table(sites_valid$location_precision))

message("\nRaw Material data availability:")
n_rm_reported <- sum(!str_detect(str_to_upper(sites_valid$raw_material),
                                 "NOT REPORTED|NA|UNKNOWN"),
                     na.rm = TRUE)
message("Sites with raw material data: ", n_rm_reported, " / ", nrow(sites_valid))

message("\nCoordinate ranges:")
message("Latitude:  ", round(min(sites_valid$latitude), 4),
        " to ", round(max(sites_valid$latitude), 4))
message("Longitude: ", round(min(sites_valid$longitude), 4),
        " to ", round(max(sites_valid$longitude), 4))


# -------------------------------------------------------------
# SECTION 14: Save Outputs
# -------------------------------------------------------------

message("\nSaving outputs...")

out_path <- file.path(chapter_root, "01_Data_Processed")

# 1. Cleaned CSV
write_csv(
  sites_valid,
  file.path(out_path, "sites_validated.csv")
)
message("Saved: sites_validated.csv")

# 2. Spatial shapefile — WGS84
st_write(
  sites_wgs84,
  file.path(out_path, "sites_wgs84.shp"),
  delete_layer = TRUE,
  quiet = TRUE
)
message("Saved: sites_wgs84.shp")

# 3. Spatial shapefile — UTM 44N
st_write(
  sites_utm,
  file.path(out_path, "sites_utm44n.shp"),
  delete_layer = TRUE,
  quiet = TRUE
)
message("Saved: sites_utm44n.shp")

# 4. Summary table
summary_table <- sites_valid %>%
  group_by(district, period) %>%
  summarise(n_sites = n(), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = n_sites, values_fill = 0) %>%
  mutate(Total = rowSums(across(where(is.numeric))))

write_csv(
  summary_table,
  file.path(chapter_root, "11_Tables", "Table_01_Site_Summary.csv")
)
message("Saved: Table_01_Site_Summary.csv")

# 5. Save as RData for use in subsequent scripts
save(
  sites_valid,
  sites_wgs84,
  sites_utm,
  boundary,
  file = file.path(out_path, "sites_validated.RData")
)
message("Saved: sites_validated.RData")


# -------------------------------------------------------------
# SECTION 15: Log Session
# -------------------------------------------------------------

log_file <- file.path(chapter_root, "13_Logs",
                      paste0("script01_log_", Sys.Date(), ".txt"))
sink(log_file)
cat("SCRIPT 01 LOG — Site Data Loading + Validation\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Total sites loaded:", nrow(sites_raw), "\n")
cat("Sites after validation:", nrow(sites_valid), "\n")
cat("Sites per period:\n")
print(table(sites_valid$period))
cat("\nSites per district:\n")
print(table(sites_valid$district))
sink()

message("\nLog saved: ", log_file)


# -------------------------------------------------------------
# SCRIPT 01 COMPLETE
# -------------------------------------------------------------

message("\n=============================================================")
message("Script 01 complete.")
message("Validated sites: ", nrow(sites_valid))
message("Outputs saved to: ", out_path)
message("Next: Run Script 02 — Background Point Generation")
message("=============================================================")