# =============================================================
# SCRIPT 02 (REBUILT): Background Point Generation
# =============================================================
# Project: Reading the Palaeolithic Landscape
# Chapter: Mapping the Past (Springer Nature, 2026)
# Author: Sushant Begade | RTMNU Nagpur
# ORCID: 0009-0003-0804-1763
# Rebuild date: August 2026
# =============================================================
# CHANGES FROM ORIGINAL SCRIPT 02:
#   1. TWO parallel background samples generated instead of one:
#      (a) background_uniform  — original method, random across
#          full district boundary (~22,000 km2). Represents
#          "landscape availability."
#      (b) background_envelope — NEW. Constrained to a buffered
#          convex hull around known site locations. Represents
#          "surveyed landscape only." Directly tests whether RF/MW
#          results in Script 06/07 depend on background extent
#          (i.e. survey-effort bias, reviewer critical problem #4).
#   2. Both use min_dist_bg = 2000m spatial thinning, same seed.
#   3. period now pulled from rebuilt Script 01 output — LP/MP/UP
#      counts are 48/62/27 (n=137 single-period), MULTI (n=60)
#      carried as separate category, NOT used for MW/PCA/RF period
#      stratification but included in combined_pts for Fig 1 /
#      descriptive stats and for the new Persistent-Places check.
#   4. coord_precision carried through into sites_minimal so it's
#      available in every downstream matrix without re-joining.
# =============================================================

load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")

load(file.path(chapter_root, "01_Data_Processed/sites_validated.RData"))

library(terra)
library(sf)
library(tidyverse)
library(spatstat)

set.seed(rf_seed)

message("Parameters and site data loaded.")
message("Sites loaded: ", nrow(sites_utm))
message("Period breakdown: LP=48 MP=62 UP=27 MULTI=60 (see Script 01 log)")


# -------------------------------------------------------------
# SECTION 2: Prepare Study Area Mask (full boundary)
# -------------------------------------------------------------

study_area <- st_union(boundary) %>% st_as_sf()
study_area_vect <- vect(study_area)

message("\nFull study area (uniform bg universe):")
print(st_bbox(study_area))


# -------------------------------------------------------------
# SECTION 2b: Survey-Envelope Polygon  [NEW, FIXED v2]
# -------------------------------------------------------------
# Rationale: uniform-random background across the FULL 22,000 km2
# boundary compares sites against landscape that was never
# meaningfully surveyed for Palaeolithic remains. That confounds
# "environmental preference" with "where archaeologists looked."
#
# v1 BUG: convex hull of all 197 sites + 5km buffer covered 84% of
# the full boundary — a hull bridges straight over unsurveyed
# interior/ridge terrain lying between site clusters (sites cluster
# along valleys, not across the whole district), so it barely
# restricts anything and the sensitivity check is close to useless.
#
# FIX: union of small buffers around EACH individual site instead
# of one hull. No bridging — only actual site vicinities count as
# "surveyed." envelope_buffer_m set well under the mean nearest-
# neighbour site spacing (~8.6km per Script 02 uniform-bg check)
# so buffers mostly stay disjoint rather than re-merging into one
# blob. ADJUST if you have real survey transect/block polygons —
# use those instead of this proxy if they exist. Flag as
# approximation in Methods regardless.

envelope_buffer_m <- 3000  # 3 km buffer around EACH site, then dissolved

survey_envelope <- sites_utm %>%
  st_buffer(dist = envelope_buffer_m) %>%
  st_union() %>%
  st_intersection(study_area) %>%
  st_as_sf()

envelope_area_km2   <- as.numeric(st_area(survey_envelope)) / 1e6
full_area_km2        <- as.numeric(st_area(study_area)) / 1e6

message("\nSurvey-envelope polygon created (", envelope_buffer_m/1000,
        "km buffer around each site, dissolved).")
message("Envelope area: ", round(envelope_area_km2, 0), " km2")
message("Full boundary area: ", round(full_area_km2, 0), " km2")
message("Envelope = ", round(100 * envelope_area_km2 / full_area_km2, 1),
        "% of full study area.")
if ((envelope_area_km2 / full_area_km2) > 0.5) {
  message("WARNING: envelope still >50% of full area — site buffers are")
  message("merging extensively. Consider reducing envelope_buffer_m further")
  message("or check for tight site clustering inflating union size.")
}

envelope_vect <- vect(survey_envelope)


# -------------------------------------------------------------
# SECTION 3: Reference Raster for Spatial Grid
# -------------------------------------------------------------

dem <- rast(paths$dem)
dem <- terra::project(dem, "EPSG:32644")
dem_masked <- dem %>% crop(study_area_vect) %>% mask(study_area_vect)
message("\nDEM loaded and masked. Resolution: ", res(dem_masked)[1], " m")


# -------------------------------------------------------------
# SECTION 4: Helper — Generate + Spatially Thin Background Set
# -------------------------------------------------------------

generate_background <- function(sampling_polygon, n_target, min_dist,
                                label, seed_val) {
  
  set.seed(seed_val)
  n_candidate <- n_target * 10
  
  message("\n--- Generating background: ", label, " ---")
  candidates <- st_sample(sampling_polygon, size = n_candidate, type = "random") %>%
    st_as_sf() %>% st_set_crs(crs_utm44n)
  message("  Candidates generated: ", nrow(candidates))
  
  coords_mat <- st_coordinates(candidates)
  selected_idx <- c(1)
  for (i in 2:nrow(candidates)) {
    dists <- sqrt((coords_mat[i, 1] - coords_mat[selected_idx, 1])^2 +
                    (coords_mat[i, 2] - coords_mat[selected_idx, 2])^2)
    if (min(dists) >= min_dist) selected_idx <- c(selected_idx, i)
    if (length(selected_idx) >= n_target) break
  }
  
  pts <- candidates[selected_idx, ]
  message("  After ", min_dist, "m spatial thinning: ", nrow(pts))
  
  if (nrow(pts) < n_target) {
    warning("  Only ", nrow(pts), " points survived thinning for '", label,
            "'. Consider larger n_candidate or smaller min_dist.")
  } else {
    pts <- pts[1:n_target, ]
  }
  
  within_check <- st_within(pts, study_area, sparse = FALSE)
  n_outside <- sum(!within_check)
  if (n_outside > 0) {
    message("  Removing ", n_outside, " points outside district boundary...")
    pts <- pts[as.vector(within_check), ]
  }
  
  pts <- pts %>%
    mutate(
      point_id     = paste0("BG_", toupper(substr(label,1,3)), "_", sprintf("%04d", 1:n())),
      point_type   = "background",
      period       = "background",
      bg_source    = label
    )
  bg_coords <- st_coordinates(pts)
  pts <- pts %>% mutate(easting = bg_coords[,1], northing = bg_coords[,2])
  
  message("  Final ", label, " background points: ", nrow(pts))
  return(pts)
}


# -------------------------------------------------------------
# SECTION 5: Generate BOTH Background Sets
# -------------------------------------------------------------

background_uniform <- generate_background(
  sampling_polygon = study_area,
  n_target  = n_background,
  min_dist  = min_dist_bg,
  label     = "uniform",
  seed_val  = rf_seed
)

background_envelope <- generate_background(
  sampling_polygon = survey_envelope,
  n_target  = n_background,
  min_dist  = min_dist_bg,
  label     = "envelope",
  seed_val  = rf_seed + 1   # different seed, same method, avoid identical draw
)


# -------------------------------------------------------------
# SECTION 6: Site Metadata (carry coord_precision + period_multi)
# -------------------------------------------------------------

site_coords <- st_coordinates(sites_utm)

sites_minimal <- data.frame(
  point_id        = paste0("SITE_", sprintf("%03d", sites_utm$site_id)),
  point_type      = "site",
  period          = as.character(sites_utm$period),      # LP/MP/UP/MULTI
  coord_precision = as.character(sites_utm$coord_precision),
  easting         = site_coords[, 1],
  northing        = site_coords[, 2]
)

message("\nSite records prepared: ", nrow(sites_minimal))
message("  By period: "); print(table(sites_minimal$period))


# -------------------------------------------------------------
# SECTION 7: Combined Datasets (one per background variant)
# -------------------------------------------------------------

bg_uniform_minimal <- background_uniform %>%
  st_drop_geometry() %>%
  transmute(point_id, point_type, period,
            coord_precision = NA_character_, easting, northing, bg_source = "uniform")

bg_envelope_minimal <- background_envelope %>%
  st_drop_geometry() %>%
  transmute(point_id, point_type, period,
            coord_precision = NA_character_, easting, northing, bg_source = "envelope")

sites_minimal <- sites_minimal %>% mutate(bg_source = NA_character_)

combined_uniform  <- bind_rows(sites_minimal, bg_uniform_minimal)
combined_envelope <- bind_rows(sites_minimal, bg_envelope_minimal)

message("\nCombined UNIFORM dataset: ", nrow(combined_uniform), " points")
message("Combined ENVELOPE dataset: ", nrow(combined_envelope), " points")


# -------------------------------------------------------------
# SECTION 8: Minimum Distance Checks (both variants)
# -------------------------------------------------------------

check_min_dist <- function(bg_sf, site_coords_mat, label) {
  bg_mat <- st_coordinates(bg_sf)
  min_d <- apply(bg_mat, 1, function(bg) {
    min(sqrt((bg[1] - site_coords_mat[,1])^2 + (bg[2] - site_coords_mat[,2])^2))
  })
  message("\n[", label, "] bg points <1km from any site: ", sum(min_d < 1000))
  message("[", label, "] min dist to nearest site: ", round(min(min_d), 0), " m")
  message("[", label, "] mean dist to nearest site: ", round(mean(min_d), 0), " m")
}

check_min_dist(background_uniform,  site_coords, "uniform")
check_min_dist(background_envelope, site_coords, "envelope")


# -------------------------------------------------------------
# SECTION 9: Save Outputs
# -------------------------------------------------------------

out_path <- file.path(chapter_root, "02_Background_Points")

st_write(background_uniform,
         file.path(out_path, "background_uniform_utm44n.shp"),
         delete_layer = TRUE, quiet = TRUE)
st_write(background_envelope,
         file.path(out_path, "background_envelope_utm44n.shp"),
         delete_layer = TRUE, quiet = TRUE)
message("\nSaved: background_uniform_utm44n.shp, background_envelope_utm44n.shp")

st_write(survey_envelope,
         file.path(out_path, "survey_envelope_polygon.shp"),
         delete_layer = TRUE, quiet = TRUE)
message("Saved: survey_envelope_polygon.shp  [cite in Methods / Fig 1 inset]")

write_csv(combined_uniform,  file.path(out_path, "combined_sites_background_uniform.csv"))
write_csv(combined_envelope, file.path(out_path, "combined_sites_background_envelope.csv"))
message("Saved: combined_sites_background_uniform.csv, combined_sites_background_envelope.csv")

save(
  background_uniform, background_envelope,
  combined_uniform, combined_envelope,
  study_area, study_area_vect,
  survey_envelope, envelope_vect,
  envelope_area_km2, full_area_km2,
  dem_masked,
  file = file.path(out_path, "background_pts.RData")
)
message("Saved: background_pts.RData")


# -------------------------------------------------------------
# SECTION 10: QC Plot — Both Variants Side by Side
# -------------------------------------------------------------

fig_path <- file.path(chapter_root, "10_Figures")

png(file.path(fig_path, "QC_background_uniform_vs_envelope.png"),
    width = fig_width * 2, height = fig_height, units = "mm", res = fig_dpi)
par(mfrow = c(1, 2))

plot(st_geometry(study_area), col = "grey95", border = "grey40",
     main = "Uniform background", sub = paste0("n=", nrow(background_uniform)))
plot(st_geometry(background_uniform), add = TRUE, pch = 16, cex = 0.3, col = "steelblue")
plot(st_geometry(sites_utm), add = TRUE, pch = 17, cex = 0.5, col = "red")

plot(st_geometry(study_area), col = "grey95", border = "grey40",
     main = "Envelope background", sub = paste0("n=", nrow(background_envelope)))
plot(st_geometry(survey_envelope), add = TRUE, border = "darkorange", lwd = 1.5, col = NA)
plot(st_geometry(background_envelope), add = TRUE, pch = 16, cex = 0.3, col = "steelblue")
plot(st_geometry(sites_utm), add = TRUE, pch = 17, cex = 0.5, col = "red")

dev.off()
message("\nQC comparison plot saved.")


# -------------------------------------------------------------
# SECTION 11: Log Session
# -------------------------------------------------------------

log_file <- file.path(chapter_root, "13_Logs",
                      paste0("script02_log_", Sys.Date(), ".txt"))
sink(log_file)
cat("SCRIPT 02 LOG (REBUILT) — Background Point Generation\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Uniform background: n=", nrow(background_uniform), "\n")
cat("Envelope background: n=", nrow(background_envelope), "\n")
cat("Survey envelope area:", round(envelope_area_km2,0), "km2 (",
    round(100*envelope_area_km2/full_area_km2,1), "% of full boundary)\n")
cat("Envelope buffer distance:", envelope_buffer_m, "m around each site (dissolved)\n")
sink()
message("Log saved: ", log_file)


# -------------------------------------------------------------
# SCRIPT 02 COMPLETE
# -------------------------------------------------------------

message("\n=============================================================")
message("Script 02 (REBUILT) complete.")
message("TWO background sets now exist: uniform (n=", nrow(background_uniform),
        ") + envelope (n=", nrow(background_envelope), ")")
message("Envelope covers ", round(100*envelope_area_km2/full_area_km2,1),
        "% of full study area.")
message("Scripts 03-07 must now run TWICE (once per bg variant) OR run once")
message("on uniform (primary) + envelope as a Discussion sensitivity check —")
message("decide before Script 03. Recommend: uniform = primary throughout,")
message("envelope = rerun of Script 06/07 ONLY, reported as robustness check.")
message("Next: Script 03 — Raster Variable Extraction")
message("  (fix background matrix numeric type bug + CHELSA-LGM handling)")
message("=============================================================")