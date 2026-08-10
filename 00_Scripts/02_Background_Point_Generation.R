# =============================================================
# SCRIPT 02: Background Point Generation
# =============================================================
# Project: Reading the Palaeolithic Landscape
# Chapter: Mapping the Past (Springer Nature, 2026)
# Author: Sushant Begade | RTMNU Nagpur
# ORCID: 0009-0003-0804-1763
# Date: August 2026
# =============================================================
# PURPOSE: Generate 1000 spatially filtered pseudo-absence
# background points across the study area extent.
# Spatial filtering ensures minimum 2km separation between
# background points to avoid spatial autocorrelation.
# Points constrained within district boundary.
# Save outputs for use in all extraction scripts.
# =============================================================


# -------------------------------------------------------------
# SECTION 1: Load Global Parameters + Validated Sites
# -------------------------------------------------------------

load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")

load(file.path(chapter_root, "01_Data_Processed/sites_validated.RData"))

library(terra)
library(sf)
library(tidyverse)
library(spatstat)

set.seed(rf_seed)

message("Parameters and site data loaded.")
message("Sites loaded: ", nrow(sites_utm))


# -------------------------------------------------------------
# SECTION 2: Prepare Study Area Mask
# -------------------------------------------------------------

message("\nPreparing study area mask...")

# Union district boundary — single polygon
study_area <- st_union(boundary) %>%
  st_as_sf()

# Confirm CRS
message("Study area CRS: ", st_crs(study_area)$input)
message("Study area approx. bbox:")
print(st_bbox(study_area))

# Convert to terra SpatVector for raster masking
study_area_vect <- vect(study_area)


# -------------------------------------------------------------
# SECTION 3: Load Reference Raster for Spatial Grid
# -------------------------------------------------------------

message("\nLoading reference raster (DEM)...")

dem <- rast(paths$dem)

# Force terra::project to avoid masking conflict
dem <- terra::project(dem, "EPSG:32644")

# Crop + mask to study area
dem_masked <- dem %>%
  crop(study_area_vect) %>%
  mask(study_area_vect)

message("DEM loaded and masked.")
message("Resolution: ", res(dem_masked)[1], " m")
message("Extent: ", paste(round(ext(dem_masked)[1:4], 0), collapse = ", "))


# -------------------------------------------------------------
# SECTION 4: Generate Candidate Random Points
# -------------------------------------------------------------

message("\nGenerating candidate background points...")

# Generate large pool — will filter down to n_background
# Oversample by 10x to ensure enough survive spatial filtering
n_candidate <- n_background * 10

# Sample random points within study area boundary
candidates <- st_sample(
  study_area,
  size = n_candidate,
  type = "random"
) %>%
  st_as_sf() %>%
  st_set_crs(crs_utm44n)

message("Candidate points generated: ", nrow(candidates))


# -------------------------------------------------------------
# SECTION 5: Spatial Filtering — Minimum Distance Between Points
# -------------------------------------------------------------

message("\nApplying spatial filter (min ", min_dist_bg, "m between points)...")

# Convert to matrix for distance filtering
coords_mat <- st_coordinates(candidates)

# Sequential spatial thinning
selected_idx <- c(1)  # Always keep first point

for (i in 2:nrow(candidates)) {
  # Calculate distances from candidate i to all already-selected points
  dists <- sqrt(
    (coords_mat[i, 1] - coords_mat[selected_idx, 1])^2 +
      (coords_mat[i, 2] - coords_mat[selected_idx, 2])^2
  )
  # Keep point only if minimum distance exceeded
  if (min(dists) >= min_dist_bg) {
    selected_idx <- c(selected_idx, i)
  }
  # Stop once we have enough
  if (length(selected_idx) >= n_background) break
}

background_pts <- candidates[selected_idx, ]

message("Background points after spatial filtering: ", nrow(background_pts))

# Check if we have enough
if (nrow(background_pts) < n_background) {
  warning("Only ", nrow(background_pts), " points generated after filtering.",
          " Consider reducing min_dist_bg or increasing n_candidate.")
} else {
  message("Target of ", n_background, " background points achieved.")
  # Trim to exactly n_background
  background_pts <- background_pts[1:n_background, ]
}


# -------------------------------------------------------------
# SECTION 6: Confirm All Points Within Boundary
# -------------------------------------------------------------

message("\nValidating background points within boundary...")

within_check <- st_within(background_pts, study_area, sparse = FALSE)
n_outside <- sum(!within_check)
message("Points outside boundary: ", n_outside)

# Remove any outside (should be zero)
background_pts <- background_pts[as.vector(within_check), ]
message("Final background points: ", nrow(background_pts))


# -------------------------------------------------------------
# SECTION 7: Add Metadata Columns
# -------------------------------------------------------------

background_pts <- background_pts %>%
  mutate(
    point_id   = paste0("BG_", sprintf("%04d", 1:n())),
    point_type = "background",
    period     = "background"
  )

# Extract coordinates as columns
bg_coords <- st_coordinates(background_pts)
background_pts <- background_pts %>%
  mutate(
    easting  = bg_coords[, 1],
    northing = bg_coords[, 2]
  )

message("\nBackground point metadata added.")


# -------------------------------------------------------------
# SECTION 8: Add Site Point Metadata for Combined Dataset
# -------------------------------------------------------------

message("\nPreparing combined site + background dataset...")

# Add coordinates to site points
site_coords <- st_coordinates(sites_utm)

# Build minimal site dataframe — avoid geometry conflict
sites_minimal <- data.frame(
  point_id   = paste0("SITE_", sprintf("%03d", sites_utm$site_id)),
  point_type = "site",
  period     = sites_utm$period,
  easting    = site_coords[, 1],
  northing   = site_coords[, 2]
)

# Build minimal background dataframe
background_minimal <- data.frame(
  point_id   = background_pts$point_id,
  point_type = "background",
  period     = "background",
  easting    = background_pts$easting,
  northing   = background_pts$northing
)

# Combine — plain dataframe, no geometry conflict
combined_pts <- bind_rows(sites_minimal, background_minimal)

message("Combined dataset: ", nrow(combined_pts), " points")
message("  Sites: ", sum(combined_pts$point_type == "site"))
message("  Background: ", sum(combined_pts$point_type == "background"))


# -------------------------------------------------------------
# SECTION 9: Minimum Distance Check — Sites vs Background
# -------------------------------------------------------------

message("\nChecking minimum distance between sites and background points...")

site_coords_mat <- site_coords        # already computed in Section 8
bg_coords_mat   <- as.matrix(
  background_pts %>% st_drop_geometry() %>% select(easting, northing)
)

min_dists_to_sites <- apply(bg_coords_mat, 1, function(bg) {
  min(sqrt(
    (bg[1] - site_coords_mat[, 1])^2 +
      (bg[2] - site_coords_mat[, 2])^2
  ))
})

message("Background points < 1km from any site: ",
        sum(min_dists_to_sites < 1000))
message("Min distance (bg to nearest site): ",
        round(min(min_dists_to_sites), 0), " m")
message("Mean distance (bg to nearest site): ",
        round(mean(min_dists_to_sites), 0), " m")


# -------------------------------------------------------------
# SECTION 10: Summary Statistics
# -------------------------------------------------------------

message("\n=== BACKGROUND POINT SUMMARY ===")
message("Total background points: ", nrow(background_pts))
message("Spatial filter applied: ", min_dist_bg, " m minimum separation")
message("Seed used: ", rf_seed)
message("Points within boundary: ", nrow(background_pts))
message("\nSpatial extent of background points:")
print(st_bbox(background_pts))


# -------------------------------------------------------------
# SECTION 11: Save Outputs
# -------------------------------------------------------------

message("\nSaving outputs...")

out_path <- file.path(chapter_root, "02_Background_Points")

# 1. Background points shapefile — UTM 44N
st_write(
  background_pts,
  file.path(out_path, "background_pts_utm44n.shp"),
  delete_layer = TRUE,
  quiet = TRUE
)
message("Saved: background_pts_utm44n.shp")

# 2. Background points — WGS84
background_wgs84 <- st_transform(background_pts, crs = crs_wgs84)
st_write(
  background_wgs84,
  file.path(out_path, "background_pts_wgs84.shp"),
  delete_layer = TRUE,
  quiet = TRUE
)
message("Saved: background_pts_wgs84.shp")

# 3. Combined site + background CSV
write_csv(
  combined_pts,
  file.path(out_path, "combined_sites_background.csv")
)
message("Saved: combined_sites_background.csv")

# RData
save(
  background_pts,
  background_wgs84,
  combined_pts,
  study_area,
  study_area_vect,
  dem_masked,
  file = file.path(out_path, "background_pts.RData")
)

# 4. RData for subsequent scripts
save(
  background_pts,
  background_wgs84,
  combined_pts,
  study_area,
  study_area_vect,
  dem_masked,
  file = file.path(out_path, "background_pts.RData")
)
message("Saved: background_pts.RData")


# -------------------------------------------------------------
# SECTION 12: Quick Visual Check
# -------------------------------------------------------------

message("\nGenerating quick validation plot...")

fig_path <- file.path(chapter_root, "10_Figures")

png(
  file.path(fig_path, "QC_background_points_distribution.png"),
  width  = fig_width,
  height = fig_height,
  units  = "mm",
  res    = fig_dpi
)

plot(st_geometry(study_area), col = "grey95", border = "grey40",
     main = "Background Points + Site Locations",
     sub  = paste0("Background n=", nrow(background_pts),
                   " | Sites n=", nrow(sites_utm)))
plot(st_geometry(background_pts), add = TRUE,
     pch = 16, cex = 0.3, col = "steelblue")
plot(st_geometry(sites_utm), add = TRUE,
     pch = 17, cex = 0.6, col = "red")
legend("bottomright",
       legend = c("Background", "Palaeolithic Sites"),
       pch    = c(16, 17),
       col    = c("steelblue", "red"),
       cex    = 0.7)

dev.off()
message("QC plot saved.")


# -------------------------------------------------------------
# SECTION 13: Log Session
# -------------------------------------------------------------

log_file <- file.path(chapter_root, "13_Logs",
                      paste0("script02_log_", Sys.Date(), ".txt"))
sink(log_file)
cat("SCRIPT 02 LOG — Background Point Generation\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Seed:", rf_seed, "\n")
cat("Target background points:", n_background, "\n")
cat("Min separation distance:", min_dist_bg, "m\n")
cat("Final background points:", nrow(background_pts), "\n")
cat("Total combined dataset:", nrow(combined_pts), "\n")
sink()

message("Log saved: ", log_file)


# -------------------------------------------------------------
# SCRIPT 02 COMPLETE
# -------------------------------------------------------------

message("\n=============================================================")
message("Script 02 complete.")
message("Background points generated: ", nrow(background_pts))
message("Combined dataset: ", nrow(combined_pts), " points")
message("Outputs saved to: ", out_path)
message("Next: Run Script 03 — Raster Variable Extraction")
message("=============================================================")