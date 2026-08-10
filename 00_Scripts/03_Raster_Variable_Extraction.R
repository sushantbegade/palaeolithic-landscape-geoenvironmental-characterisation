# =============================================================
# SCRIPT 03: Raster Variable Extraction
# =============================================================
# Project: Reading the Palaeolithic Landscape
# Chapter: Mapping the Past (Springer Nature, 2026)
# Author: Sushant Begade | RTMNU Nagpur
# ORCID: 0009-0003-0804-1763
# Date: August 2026
# =============================================================
# PURPOSE: Extract all raster variable values at 197 site
# points and 1000 background points. Covers:
# - Sentinel-2A bands (B01-B12, B8A)
# - Sentinel-2A derived indices (NDVI, NDWI, MNDWI, NDBI,
#   BSI, SAVI, MSAVI)
# - Soil properties (Depth, Erosion, Productivity,
#   Slope, Texture)
# - DEM
# - CHELSA modern climate (bio01, bio12, bio15)
# - CHELSA TraCE21k palaeoclimate (bio01, bio12, bio15)
# - WorldClim (bio01, bio12, bio15)
# Output: master extraction matrix for all analyses.
# =============================================================


# -------------------------------------------------------------
# SECTION 1: Load Global Parameters + Data
# -------------------------------------------------------------

load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")

load(file.path(chapter_root, "01_Data_Processed/sites_validated.RData"))
load(file.path(chapter_root, "02_Background_Points/background_pts.RData"))

library(terra)
library(sf)
library(tidyverse)

set.seed(rf_seed)

message("Parameters and data loaded.")
message("Sites: ", nrow(sites_utm))
message("Background points: ", nrow(background_pts))


# -------------------------------------------------------------
# SECTION 2: Prepare Extraction Points
# -------------------------------------------------------------

message("\nPreparing extraction points...")

# Inspect what columns are available
message("Sites columns: ", paste(names(sites_utm), collapse = ", "))
message("Background columns: ", paste(names(background_pts), collapse = ", "))

# Build site extraction dataframe from available columns
site_coords <- st_coordinates(sites_utm)

sites_extract <- data.frame(
  point_id   = paste0("SITE_", sprintf("%03d", sites_utm$site_id)),
  point_type = "site",
  period     = sites_utm$period,
  easting    = site_coords[, 1],
  northing   = site_coords[, 2]
)

# Build background extraction dataframe
bg_coords <- st_coordinates(background_pts)

bg_extract <- data.frame(
  point_id   = paste0("BG_", sprintf("%04d", 1:nrow(background_pts))),
  point_type = "background",
  period     = "background",
  easting    = bg_coords[, 1],
  northing   = bg_coords[, 2]
)

# Combined
all_pts_sf <- bind_rows(sites_extract, bg_extract)

message("Sites prepared: ",     nrow(sites_extract))
message("Background prepared: ", nrow(bg_extract))
message("Total points: ",        nrow(all_pts_sf))

# Convert to SpatVector for terra extraction
all_vect <- vect(
  all_pts_sf,
  geom = c("easting", "northing"),
  crs  = crs_utm44n
)

# Also keep site and background vects separately
sites_vect     <- vect(sites_extract,
                       geom = c("easting", "northing"),
                       crs  = crs_utm44n)
background_vect <- vect(bg_extract,
                        geom = c("easting", "northing"),
                        crs  = crs_utm44n)

message("SpatVectors created.")


# -------------------------------------------------------------
# SECTION 3: Helper Function — Safe Raster Extract
# -------------------------------------------------------------

safe_extract <- function(raster_path, points_vect,
                         var_names, reproject = TRUE) {
  tryCatch({
    r <- rast(raster_path)
    if (reproject) {
      r <- terra::project(r, crs_utm44n)
    }
    vals <- terra::extract(r, points_vect, ID = FALSE)
    if (ncol(vals) == length(var_names)) {
      names(vals) <- var_names
    } else {
      # More bands than expected — name sequentially
      names(vals) <- paste0(var_names[1], "_", seq_len(ncol(vals)))
    }
    message("  Extracted: ", paste(var_names, collapse = ", "))
    return(vals)
  }, error = function(e) {
    message("  ERROR extracting ", var_names[1], ": ", e$message)
    # Return NA columns
    df <- as.data.frame(matrix(NA, nrow = nrow(all_pts_sf),
                               ncol = length(var_names)))
    names(df) <- var_names
    return(df)
  })
}


# -------------------------------------------------------------
# SECTION 4: Extract DEM
# -------------------------------------------------------------

message("\nExtracting DEM...")

ext_dem <- safe_extract(
  paths$dem, all_vect, "elevation"
)


# -------------------------------------------------------------
# SECTION 5: Extract Sentinel-2A Bands
# -------------------------------------------------------------

message("\nExtracting Sentinel-2A bands...")

s2_bands <- list(
  list(path = paths$B01, name = "S2_B01"),
  list(path = paths$B02, name = "S2_B02"),
  list(path = paths$B03, name = "S2_B03"),
  list(path = paths$B04, name = "S2_B04"),
  list(path = paths$B05, name = "S2_B05"),
  list(path = paths$B06, name = "S2_B06"),
  list(path = paths$B07, name = "S2_B07"),
  list(path = paths$B08, name = "S2_B08"),
  list(path = paths$B09, name = "S2_B09"),
  list(path = paths$B11, name = "S2_B11"),
  list(path = paths$B12, name = "S2_B12"),
  list(path = paths$B8A, name = "S2_B8A")
)

ext_s2_bands <- lapply(s2_bands, function(b) {
  safe_extract(b$path, all_vect, b$name)
})
ext_s2_bands <- bind_cols(ext_s2_bands)


# -------------------------------------------------------------
# SECTION 6: Extract Sentinel-2A Derived Indices
# -------------------------------------------------------------

message("\nExtracting Sentinel-2A derived indices...")

s2_indices <- list(
  list(path = paths$NDVI,  name = "NDVI"),
  list(path = paths$NDWI,  name = "NDWI"),
  list(path = paths$MNDWI, name = "MNDWI"),
  list(path = paths$NDBI,  name = "NDBI"),
  list(path = paths$BSI,   name = "BSI"),
  list(path = paths$SAVI,  name = "SAVI"),
  list(path = paths$MSAVI, name = "MSAVI")
)

ext_indices <- lapply(s2_indices, function(idx) {
  safe_extract(idx$path, all_vect, idx$name)
})
ext_indices <- bind_cols(ext_indices)


# -------------------------------------------------------------
# SECTION 7: Extract Soil Properties
# -------------------------------------------------------------

message("\nExtracting soil properties...")

soil_layers <- list(
  list(path = paths$soil_depth,       name = "soil_depth"),
  list(path = paths$soil_erosion,     name = "soil_erosion"),
  list(path = paths$soil_productivity,name = "soil_productivity"),
  list(path = paths$soil_slope,       name = "soil_slope"),
  list(path = paths$soil_texture,     name = "soil_texture")
)

ext_soil <- lapply(soil_layers, function(s) {
  safe_extract(s$path, all_vect, s$name)
})
ext_soil <- bind_cols(ext_soil)


# -------------------------------------------------------------
# SECTION 8: Extract CHELSA Modern Climate
# -------------------------------------------------------------

message("\nExtracting CHELSA modern climate...")

ext_chelsa_modern <- bind_cols(
  safe_extract(paths$chelsa_bio01_modern, all_vect, "chelsa_bio01_modern"),
  safe_extract(paths$chelsa_bio12_modern, all_vect, "chelsa_bio12_modern"),
  safe_extract(paths$chelsa_bio15_modern, all_vect, "chelsa_bio15_modern")
)


# -------------------------------------------------------------
# SECTION 9: Extract CHELSA TraCE21k Palaeoclimate (~20ka LGM)
# -------------------------------------------------------------

message("\nExtracting CHELSA TraCE21k palaeoclimate (LGM ~20ka)...")

ext_chelsa_lgm <- bind_cols(
  safe_extract(paths$chelsa_bio01_lgm, all_vect, "chelsa_bio01_lgm"),
  safe_extract(paths$chelsa_bio12_lgm, all_vect, "chelsa_bio12_lgm"),
  safe_extract(paths$chelsa_bio15_lgm, all_vect, "chelsa_bio15_lgm")
)


# -------------------------------------------------------------
# SECTION 10: Extract WorldClim 2.1
# -------------------------------------------------------------

message("\nExtracting WorldClim 2.1...")

ext_worldclim <- bind_cols(
  safe_extract(paths$wc_bio01, all_vect, "wc_bio01"),
  safe_extract(paths$wc_bio12, all_vect, "wc_bio12"),
  safe_extract(paths$wc_bio15, all_vect, "wc_bio15")
)


# -------------------------------------------------------------
# SECTION 11: Assemble Master Extraction Matrix
# -------------------------------------------------------------

message("\nAssembling master extraction matrix...")

master_matrix <- bind_cols(
  all_pts_sf,         # point_id, point_type, period, easting, northing
  ext_dem,            # elevation
  ext_s2_bands,       # S2_B01 to S2_B8A (12 bands)
  ext_indices,        # NDVI, NDWI, MNDWI, NDBI, BSI, SAVI, MSAVI
  ext_soil,           # soil_depth, erosion, productivity, slope, texture
  ext_chelsa_modern,  # chelsa_bio01/12/15_modern
  ext_chelsa_lgm,     # chelsa_bio01/12/15_lgm
  ext_worldclim       # wc_bio01/12/15
)

message("Master matrix dimensions: ",
        nrow(master_matrix), " rows × ", ncol(master_matrix), " columns")
message("Columns: ", paste(names(master_matrix), collapse = ", "))


# -------------------------------------------------------------
# SECTION 12: Check for Missing Values
# -------------------------------------------------------------

message("\nChecking for missing values...")

# Count NAs per column
na_summary <- master_matrix %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(),
               names_to  = "variable",
               values_to = "n_missing") %>%
  filter(n_missing > 0) %>%
  arrange(desc(n_missing))

if (nrow(na_summary) == 0) {
  message("No missing values detected.")
} else {
  message("Variables with missing values:")
  print(na_summary)
  message("\nTotal rows with any NA: ",
          sum(!complete.cases(master_matrix)))
}

# Split into site and background for inspection
sites_matrix    <- master_matrix %>% filter(point_type == "site")
bg_matrix       <- master_matrix %>% filter(point_type == "background")

message("\nSite rows: ",       nrow(sites_matrix))
message("Background rows: ",  nrow(bg_matrix))
message("Complete site rows: ", sum(complete.cases(sites_matrix)))


# -------------------------------------------------------------
# SECTION 13: Basic Descriptive Statistics
# -------------------------------------------------------------

message("\nDescriptive statistics for key variables at site locations:")

site_stats <- sites_matrix %>%
  select(period, elevation, NDVI, NDWI, MNDWI, BSI,
         soil_depth, soil_productivity,
         chelsa_bio01_lgm, chelsa_bio12_lgm) %>%
  group_by(period) %>%
  summarise(across(where(is.numeric),
                   list(mean = ~round(mean(., na.rm = TRUE), 3),
                        sd   = ~round(sd(., na.rm = TRUE), 3)),
                   .names = "{.col}_{.fn}"),
            .groups = "drop")

print(site_stats)


# -------------------------------------------------------------
# SECTION 14: Save Outputs
# -------------------------------------------------------------

message("\nSaving outputs...")

out_path <- file.path(chapter_root, "03_Extracted_Values")

# 1. Full master matrix
write_csv(
  master_matrix,
  file.path(out_path, "master_extraction_matrix.csv")
)
message("Saved: master_extraction_matrix.csv")

# 2. Sites only
write_csv(
  sites_matrix,
  file.path(out_path, "sites_extraction_matrix.csv")
)
message("Saved: sites_extraction_matrix.csv")

# 3. Background only
write_csv(
  bg_matrix,
  file.path(out_path, "background_extraction_matrix.csv")
)
message("Saved: background_extraction_matrix.csv")

# 4. RData — primary input for all downstream scripts
save(
  master_matrix,
  sites_matrix,
  bg_matrix,
  file = file.path(out_path, "master_extraction_matrix.RData")
)
message("Saved: master_extraction_matrix.RData")

# 5. NA summary table
write_csv(
  na_summary,
  file.path(chapter_root, "11_Tables",
            "Table_QC_missing_values.csv")
)
message("Saved: Table_QC_missing_values.csv")


# -------------------------------------------------------------
# SECTION 15: Log Session
# -------------------------------------------------------------

log_file <- file.path(chapter_root, "13_Logs",
                      paste0("script03_log_", Sys.Date(), ".txt"))
sink(log_file)
cat("SCRIPT 03 LOG — Raster Variable Extraction\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Total points: ", nrow(master_matrix), "\n")
cat("Sites: ",        nrow(sites_matrix), "\n")
cat("Background: ",   nrow(bg_matrix), "\n")
cat("Variables extracted: ", ncol(master_matrix), "\n\n")
cat("Variables with missing values:\n")
print(na_summary)
sink()

message("Log saved: ", log_file)


# -------------------------------------------------------------
# SCRIPT 03 COMPLETE
# -------------------------------------------------------------

message("\n=============================================================")
message("Script 03 complete.")
message("Master matrix: ", nrow(master_matrix), " rows × ",
        ncol(master_matrix), " columns")
message("Outputs saved to: ", out_path)
message("Next: Run Script 04 — Vector Variable Extraction")
message("  (Structural geology proximity + geochemistry joins)")
message("=============================================================")