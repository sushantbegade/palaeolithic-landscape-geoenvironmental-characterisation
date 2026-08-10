# =============================================================
# SCRIPT 05: Taphonomic Bias Correction
# =============================================================
# Project: Reading the Palaeolithic Landscape
# Chapter: Mapping the Past (Springer Nature, 2026)
# Author: Sushant Begade | RTMNU Nagpur
# ORCID: 0009-0003-0804-1763
# Date: August 2026
# =============================================================
# PURPOSE: Model sediment burial depth across the study area
# using borehole drilling data and IDW interpolation.
# Output: continuous burial depth surface classified into
# low / medium / high burial risk tiers.
# Overlay with site distribution to identify zones of
# survey-induced site absence vs genuine hominin avoidance.
# This is Enhancement 2 — preempts primary methodological
# criticism of surface survey bias.
# =============================================================


# -------------------------------------------------------------
# SECTION 1: Load Global Parameters + Data
# -------------------------------------------------------------

load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")

load(file.path(chapter_root, "01_Data_Processed/sites_validated.RData"))
load(file.path(chapter_root, "02_Background_Points/background_pts.RData"))
load(file.path(chapter_root, "03_Extracted_Values/master_extraction_matrix.RData"))

library(terra)
library(sf)
library(tidyverse)
library(gstat)
library(readxl)

set.seed(rf_seed)

message("All data loaded.")


# -------------------------------------------------------------
# SECTION 2: Load Borehole Data
# -------------------------------------------------------------

message("\nLoading borehole data...")

# Load shapefile
borehole_sf <- st_read(paths$borehole_shp, quiet = TRUE) %>%
  st_transform(crs = crs_utm44n)

message("Borehole shapefile loaded: ", nrow(borehole_sf), " points")
message("Borehole columns: ",
        paste(names(borehole_sf), collapse = ", "))

# Load Excel for additional attributes
borehole_xl <- read_xlsx(paths$borehole_xlsx)
message("Borehole Excel loaded: ", nrow(borehole_xl),
        " rows × ", ncol(borehole_xl), " columns")
message("Excel columns: ",
        paste(names(borehole_xl), collapse = ", "))

# Print first few rows to inspect depth column
print(head(borehole_sf %>% st_drop_geometry(), 5))


# -------------------------------------------------------------
# SECTION 3: Identify Depth Column
# -------------------------------------------------------------

message("\nDepth column = length_m (borehole total length in metres)")
message("Stored as TEXT in shapefile — will coerce to numeric in Section 4.")

# Print all column names + sample values
borehole_df <- borehole_sf %>% st_drop_geometry()
for (col in names(borehole_df)) {
  vals <- borehole_df[[col]]
  if (is.numeric(vals)) {
    message("  NUMERIC — ", col, ": range ",
            round(min(vals, na.rm = TRUE), 2), " to ",
            round(max(vals, na.rm = TRUE), 2))
  } else {
    message("  TEXT — ", col, ": sample = ",
            paste(head(vals, 2), collapse = ", "))
  }
}


# -------------------------------------------------------------
# SECTION 4: Clean Borehole Data
# -------------------------------------------------------------

message("\nCleaning borehole data...")

# Extract coordinates
bh_coords <- st_coordinates(borehole_sf)

borehole_clean <- borehole_sf %>%
  st_drop_geometry() %>%
  mutate(
    easting  = bh_coords[, 1],
    northing = bh_coords[, 2]
  )

# MANUAL: length_m is the borehole depth column
# Stored as TEXT — coerce to numeric
borehole_clean$depth_m <- suppressWarnings(
  as.numeric(borehole_clean$length_m)
)

depth_col <- "length_m"
message("Depth column manually set: ", depth_col)

# Remove missing or zero depth
borehole_clean <- borehole_clean %>%
  filter(!is.na(depth_m), depth_m > 0)

message("Valid borehole records: ", nrow(borehole_clean))
message("Depth range: ",
        round(min(borehole_clean$depth_m), 2), " to ",
        round(max(borehole_clean$depth_m), 2), " m")

# -------------------------------------------------------------
# SECTION 5: Convert to Spatial for IDW
# -------------------------------------------------------------

message("\nPreparing spatial borehole data for IDW...")

# Convert to sf
borehole_sf_clean <- st_as_sf(
  borehole_clean,
  coords = c("easting", "northing"),
  crs    = crs_utm44n
)

# Convert to sp for gstat IDW
borehole_sp <- as(borehole_sf_clean, "Spatial")

message("Borehole spatial object ready: ",
        nrow(borehole_sp), " points")


# -------------------------------------------------------------
# SECTION 6: Create Interpolation Grid
# -------------------------------------------------------------

message("\nCreating interpolation grid...")

# Use DEM as grid template — 30m resolution
dem <- rast(paths$dem)
dem <- terra::project(dem, "EPSG:32644")

# Load study area
boundary <- st_read(paths$boundary, quiet = TRUE) %>%
  st_transform(crs = crs_utm44n)
study_area <- st_union(boundary)
study_area_vect <- vect(study_area)

# Crop DEM to study area
dem_masked <- dem %>%
  crop(study_area_vect) %>%
  mask(study_area_vect)

# Convert DEM to interpolation grid (sp format)
dem_r <- raster::raster(dem_masked)
grid_sp <- as(dem_r, "SpatialPixelsDataFrame")

message("Interpolation grid created: ",
        nrow(grid_sp), " pixels at 30m resolution")


# -------------------------------------------------------------
# SECTION 7: IDW Interpolation
# -------------------------------------------------------------

message("\nRunning IDW interpolation (power = 2)...")

# IDW formula — depth_m as target variable
idw_formula <- as.formula(paste("depth_m ~ 1"))

# Run IDW
idw_result <- gstat::idw(
  formula  = idw_formula,
  locations = borehole_sp,
  newdata  = grid_sp,
  idp      = 2      # inverse distance power = 2 (standard)
)

message("IDW interpolation complete.")

# Convert result to terra raster
burial_depth_r <- rast(raster::raster(idw_result["var1.pred"]))
terra::crs(burial_depth_r) <- "EPSG:32644"

# Mask to study area
burial_depth_r <- mask(burial_depth_r, study_area_vect)

message("Burial depth surface created.")
message("Depth range in surface: ",
        round(minmax(burial_depth_r)[1], 2), " to ",
        round(minmax(burial_depth_r)[2], 2), " m")


# -------------------------------------------------------------
# SECTION 8: Classify Burial Risk Tiers
# -------------------------------------------------------------

message("\nClassifying burial risk tiers...")

# Calculate quantile breaks
depth_vals <- values(burial_depth_r, na.rm = TRUE)
q33 <- quantile(depth_vals, 0.33, na.rm = TRUE)
q66 <- quantile(depth_vals, 0.66, na.rm = TRUE)

message("Burial depth quantile breaks:")
message("  33rd percentile (Low/Medium boundary): ",
        round(q33, 2), " m")
message("  66th percentile (Medium/High boundary): ",
        round(q66, 2), " m")

# Classify: 1 = Low, 2 = Medium, 3 = High burial risk
burial_risk_r <- classify(
  burial_depth_r,
  rcl = matrix(
    c(-Inf, q33,  1,
      q33, q66,  2,
      q66,  Inf, 3),
    ncol = 3, byrow = TRUE
  )
)

message("Burial risk classes:")
message("  1 = Low  (< ", round(q33, 2), " m)")
message("  2 = Medium (", round(q33, 2), "–",
        round(q66, 2), " m)")
message("  3 = High (> ", round(q66, 2), " m)")


# -------------------------------------------------------------
# SECTION 9: Extract Burial Depth at Site + Background Points
# -------------------------------------------------------------

message("\nExtracting burial depth at all points...")

# Rebuild all points SpatVector
all_pts_df <- master_matrix %>%
  select(point_id, point_type, period, easting, northing)

all_vect <- vect(
  all_pts_df,
  geom = c("easting", "northing"),
  crs  = crs_utm44n
)

# Extract burial depth
burial_vals <- terra::extract(burial_depth_r, all_vect,
                              ID = FALSE)
names(burial_vals) <- "burial_depth_m"

# Extract burial risk class
risk_vals <- terra::extract(burial_risk_r, all_vect,
                            ID = FALSE)
names(risk_vals) <- "burial_risk_class"

# Assemble taphonomy matrix
taphonomy_matrix <- bind_cols(
  all_pts_df,
  burial_vals,
  risk_vals
) %>%
  mutate(
    burial_risk_label = case_when(
      burial_risk_class == 1 ~ "Low",
      burial_risk_class == 2 ~ "Medium",
      burial_risk_class == 3 ~ "High",
      TRUE ~ NA_character_
    )
  )

message("Taphonomy matrix: ", nrow(taphonomy_matrix),
        " rows × ", ncol(taphonomy_matrix), " columns")


# -------------------------------------------------------------
# SECTION 10: Taphonomic Bias Analysis
# -------------------------------------------------------------

message("\n=== TAPHONOMIC BIAS ANALYSIS ===")

# Site burial risk distribution
message("\nSite distribution by burial risk tier:")
site_risk <- taphonomy_matrix %>%
  filter(point_type == "site") %>%
  count(period, burial_risk_label) %>%
  pivot_wider(names_from  = burial_risk_label,
              values_from = n,
              values_fill = 0)
print(site_risk)

# Background burial risk distribution
message("\nBackground distribution by burial risk tier:")
bg_risk <- taphonomy_matrix %>%
  filter(point_type == "background") %>%
  count(burial_risk_label) %>%
  mutate(pct = round(100 * n / sum(n), 1))
print(bg_risk)

# Key question: are sites underrepresented in high burial zones?
message("\nSite density by burial risk tier:")
site_pct <- taphonomy_matrix %>%
  filter(point_type == "site") %>%
  count(burial_risk_label) %>%
  mutate(pct_sites = round(100 * n / sum(n), 1))

bg_pct <- taphonomy_matrix %>%
  filter(point_type == "background") %>%
  count(burial_risk_label) %>%
  mutate(pct_bg = round(100 * n / sum(n), 1))

comparison <- left_join(site_pct, bg_pct,
                        by = "burial_risk_label") %>%
  mutate(
    ratio = round(pct_sites / pct_bg, 2),
    interpretation = case_when(
      ratio < 0.7  ~ "UNDERREPRESENTED — possible burial bias",
      ratio > 1.3  ~ "OVERREPRESENTED — genuine preference",
      TRUE         ~ "Proportional — no strong bias"
    )
  )

message("\nSite vs Background burial risk comparison:")
print(comparison)

# Mean burial depth per period
message("\nMean burial depth at site locations by period:")
taphonomy_matrix %>%
  filter(point_type == "site") %>%
  group_by(period) %>%
  summarise(
    mean_depth_m = round(mean(burial_depth_m, na.rm = TRUE), 2),
    sd_depth_m   = round(sd(burial_depth_m,   na.rm = TRUE), 2),
    n_high_risk  = sum(burial_risk_label == "High", na.rm = TRUE),
    .groups      = "drop"
  ) %>%
  print()


# -------------------------------------------------------------
# SECTION 11: Save Raster Outputs
# -------------------------------------------------------------

message("\nSaving raster outputs...")

out_path <- file.path(chapter_root, "06_Taphonomy")

# Burial depth surface
writeRaster(
  burial_depth_r,
  file.path(out_path, "burial_depth_surface.tif"),
  overwrite = TRUE
)
message("Saved: burial_depth_surface.tif")

# Burial risk classification
writeRaster(
  burial_risk_r,
  file.path(out_path, "burial_risk_classification.tif"),
  overwrite = TRUE
)
message("Saved: burial_risk_classification.tif")

# Taphonomy matrix CSV
write_csv(
  taphonomy_matrix,
  file.path(out_path, "taphonomy_matrix.csv")
)
message("Saved: taphonomy_matrix.csv")

# Bias comparison table
write_csv(
  comparison,
  file.path(chapter_root, "11_Tables",
            "Table_taphonomic_bias_comparison.csv")
)
message("Saved: Table_taphonomic_bias_comparison.csv")

# RData
save(
  burial_depth_r,
  burial_risk_r,
  taphonomy_matrix,
  comparison,
  q33, q66,
  file = file.path(out_path, "taphonomy_data.RData")
)
message("Saved: taphonomy_data.RData")


# -------------------------------------------------------------
# SECTION 12: Quick Validation Plot
# -------------------------------------------------------------

message("\nGenerating taphonomy validation plot...")

# Reload boundary for plotting
study_area_sf <- st_union(boundary) %>% st_as_sf()
sites_sf <- st_as_sf(
  taphonomy_matrix %>% filter(point_type == "site"),
  coords = c("easting", "northing"),
  crs    = crs_utm44n
)

png(
  file.path(chapter_root, "10_Figures",
            "QC_burial_depth_surface.png"),
  width = fig_width, height = fig_height,
  units = "mm", res = fig_dpi
)

plot(burial_depth_r,
     main = "Sediment Burial Depth Surface (IDW)",
     sub  = "Borehole data interpolated | UTM Zone 44N",
     col  = rev(terrain.colors(100)))
plot(st_geometry(study_area_sf), add = TRUE,
     border = "black", lwd = 1.5)
plot(st_geometry(sites_sf), add = TRUE,
     pch = 16, cex = 0.5, col = "red")

dev.off()
message("QC plot saved.")


# -------------------------------------------------------------
# SECTION 13: Log Session
# -------------------------------------------------------------

log_file <- file.path(chapter_root, "13_Logs",
                      paste0("script05_log_", Sys.Date(), ".txt"))
sink(log_file)
cat("SCRIPT 05 LOG — Taphonomic Bias Correction\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("IDW power:", 2, "\n")
cat("Depth column used:", depth_col, "\n")
cat("Valid borehole points:", nrow(borehole_clean), "\n")
cat("Depth range:", round(min(borehole_clean$depth_m), 2),
    "to", round(max(borehole_clean$depth_m), 2), "m\n")
cat("Low/Medium boundary (q33):", round(q33, 2), "m\n")
cat("Medium/High boundary (q66):", round(q66, 2), "m\n\n")
cat("Bias comparison table:\n")
print(comparison)
sink()

message("Log saved: ", log_file)


# -------------------------------------------------------------
# SCRIPT 05 COMPLETE
# -------------------------------------------------------------

message("\n=============================================================")
message("Script 05 complete.")
message("Burial depth surface: IDW interpolated from borehole data")
message("Burial risk tiers: Low / Medium / High")
message("Taphonomy matrix: ", nrow(taphonomy_matrix), " points")
message("Outputs saved to: ", out_path)
message("Next: Run Script 06 — Mann-Whitney U Tests")
message("  (Period-specific significance matrix)")
message("=============================================================")