# =============================================================
# SCRIPT 04: Vector Variable Extraction
# =============================================================
# Project: Reading the Palaeolithic Landscape
# Chapter: Mapping the Past (Springer Nature, 2026)
# Author: Sushant Begade | RTMNU Nagpur
# ORCID: 0009-0003-0804-1763
# Date: August 2026
# =============================================================
# PURPOSE: Extract vector-based variables at all 1197 points:
# PART A — Structural Geological Proximity:
#   Distance from each point to nearest:
#   fault, dyke, lineament, shear zone, mineral deposit
# PART B — Geochemistry Buffer Joins:
#   Mean elemental concentrations within 1km, 2km, 5km
#   buffers around each point from:
#   stream sediment, soil horizon, soil regolith data
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

set.seed(rf_seed)

message("All data loaded.")
message("Master matrix: ", nrow(master_matrix), " rows × ",
        ncol(master_matrix), " columns")


# -------------------------------------------------------------
# SECTION 2: Rebuild Spatial Points from Master Matrix
# -------------------------------------------------------------

message("\nRebuilding spatial points from master matrix...")

all_pts_sf <- st_as_sf(
  master_matrix,
  coords = c("easting", "northing"),
  crs    = crs_utm44n
)

message("Spatial points rebuilt: ", nrow(all_pts_sf))


# -------------------------------------------------------------
# SECTION 3: Helper Function — Safe Vector Distance
# -------------------------------------------------------------

safe_distance <- function(vector_path, points_sf,
                          var_name, geom_type = NULL) {
  tryCatch({
    vec <- st_read(vector_path, quiet = TRUE) %>%
      st_transform(crs = crs_utm44n)
    
    # Handle geometry type if needed
    if (!is.null(geom_type)) {
      vec <- vec %>%
        st_cast(geom_type, warn = FALSE)
    }
    
    # Remove empty geometries
    vec <- vec[!st_is_empty(vec), ]
    
    # Calculate nearest distance
    dists <- st_distance(points_sf, vec) %>%
      apply(1, min) %>%
      as.numeric()
    
    df <- data.frame(dists)
    names(df) <- var_name
    message("  Extracted distance: ", var_name,
            " | Range: ", round(min(dists), 0),
            "–", round(max(dists), 0), " m")
    return(df)
    
  }, error = function(e) {
    message("  ERROR: ", var_name, " — ", e$message)
    df <- data.frame(rep(NA, nrow(points_sf)))
    names(df) <- var_name
    return(df)
  })
}


# =============================================================
# PART A: STRUCTURAL GEOLOGICAL PROXIMITY
# =============================================================

message("\n--- PART A: Structural Geological Proximity ---")


# -------------------------------------------------------------
# SECTION 4: Distance to Fault
# -------------------------------------------------------------

message("\nCalculating distance to nearest fault...")
dist_fault <- safe_distance(paths$fault, all_pts_sf, "dist_fault")


# -------------------------------------------------------------
# SECTION 5: Distance to Dyke
# -------------------------------------------------------------

message("\nCalculating distance to nearest dyke...")
dist_dyke <- safe_distance(paths$dyke, all_pts_sf, "dist_dyke")


# -------------------------------------------------------------
# SECTION 6: Distance to Lineament
# -------------------------------------------------------------

message("\nCalculating distance to nearest lineament...")
dist_lineament <- safe_distance(paths$lineament, all_pts_sf,
                                "dist_lineament")


# -------------------------------------------------------------
# SECTION 7: Distance to Shear Zone
# -------------------------------------------------------------

message("\nCalculating distance to nearest shear zone...")
dist_shear <- safe_distance(paths$shear, all_pts_sf, "dist_shear")


# -------------------------------------------------------------
# SECTION 8: Distance to Mineral Deposit
# -------------------------------------------------------------

message("\nCalculating distance to nearest mineral deposit...")
dist_mineral <- safe_distance(paths$mineral, all_pts_sf, "dist_mineral")


# -------------------------------------------------------------
# SECTION 9: Assemble Structural Geology Matrix
# -------------------------------------------------------------

message("\nAssembling structural geology matrix...")

structural_matrix <- bind_cols(
  master_matrix %>% select(point_id, point_type, period),
  dist_fault,
  dist_dyke,
  dist_lineament,
  dist_shear,
  dist_mineral
)

message("Structural matrix: ", nrow(structural_matrix),
        " rows × ", ncol(structural_matrix), " columns")

# Quick summary at site locations
message("\nStructural proximity — site means by period:")
structural_matrix %>%
  filter(point_type == "site") %>%
  group_by(period) %>%
  summarise(
    fault_km     = round(mean(dist_fault,     na.rm = TRUE) / 1000, 2),
    dyke_km      = round(mean(dist_dyke,      na.rm = TRUE) / 1000, 2),
    lineament_km = round(mean(dist_lineament, na.rm = TRUE) / 1000, 2),
    shear_km     = round(mean(dist_shear,     na.rm = TRUE) / 1000, 2),
    mineral_km   = round(mean(dist_mineral,   na.rm = TRUE) / 1000, 2),
    .groups = "drop"
  ) %>%
  print()


# =============================================================
# PART B: GEOCHEMISTRY BUFFER JOINS
# =============================================================

message("\n--- PART B: Geochemistry Buffer Joins ---")


# -------------------------------------------------------------
# SECTION 10: Load Geochemistry Shapefiles
# -------------------------------------------------------------

message("\nLoading geochemistry data...")

# Stream Sediment
geoch_stream <- st_read(paths$geoch_stream, quiet = TRUE) %>%
  st_transform(crs = crs_utm44n)
message("Stream sediment points: ", nrow(geoch_stream))
message("Stream sediment columns: ",
        paste(names(geoch_stream), collapse = ", "))

# Soil Horizon
geoch_horizon <- st_read(paths$geoch_horizon, quiet = TRUE) %>%
  st_transform(crs = crs_utm44n)
message("Soil horizon points: ", nrow(geoch_horizon))
message("Soil horizon columns: ",
        paste(names(geoch_horizon), collapse = ", "))

# Soil Regolith
geoch_regolith <- st_read(paths$geoch_regolith, quiet = TRUE) %>%
  st_transform(crs = crs_utm44n)
message("Soil regolith points: ", nrow(geoch_regolith))
message("Soil regolith columns: ",
        paste(names(geoch_regolith), collapse = ", "))


# -------------------------------------------------------------
# SECTION 11: Helper Function — Buffer Geochemistry Join
# -------------------------------------------------------------

buffer_geoch_join <- function(points_sf, geoch_sf,
                              buffer_m, prefix) {
  tryCatch({
    # Get numeric geochemistry columns only
    geoch_numeric_cols <- geoch_sf %>%
      st_drop_geometry() %>%
      select(where(is.numeric)) %>%
      names()
    
    if (length(geoch_numeric_cols) == 0) {
      message("  No numeric columns found in geochemistry layer.")
      return(NULL)
    }
    
    message("  Numeric geochemical elements found: ",
            length(geoch_numeric_cols))
    
    # Buffer each point
    pts_buffered <- st_buffer(points_sf, dist = buffer_m)
    
    # Spatial join — geoch points within each buffer
    joined <- st_join(
      pts_buffered %>% select(point_id),
      geoch_sf %>% select(all_of(geoch_numeric_cols)),
      join = st_intersects
    )
    
    # Mean per point per element
    result <- joined %>%
      st_drop_geometry() %>%
      group_by(point_id) %>%
      summarise(
        across(all_of(geoch_numeric_cols),
               ~mean(., na.rm = TRUE)),
        n_geoch_pts = n(),
        .groups = "drop"
      )
    
    # Rename with prefix + buffer distance
    buf_label <- paste0("_", prefix, "_", buffer_m / 1000, "km")
    result <- result %>%
      rename_with(
        ~paste0(., buf_label),
        .cols = all_of(geoch_numeric_cols)
      ) %>%
      rename_with(
        ~paste0("n_pts_", prefix, "_", buffer_m / 1000, "km"),
        .cols = "n_geoch_pts"
      )
    
    # Left join back to all points
    output <- points_sf %>%
      st_drop_geometry() %>%
      select(point_id) %>%
      left_join(result, by = "point_id")
    
    message("  Buffer ", buffer_m / 1000, "km join complete. ",
            "Points with geoch data: ",
            sum(!is.na(output[[2]])), "/", nrow(output))
    
    return(output)
    
  }, error = function(e) {
    message("  ERROR in buffer join: ", e$message)
    return(NULL)
  })
}


# -------------------------------------------------------------
# SECTION 12: Run Buffer Joins — All Distances + All Layers
# -------------------------------------------------------------

message("\nRunning geochemistry buffer joins...")
message("Buffer distances: 1km, 2km, 5km")
message("Layers: stream sediment, soil horizon, soil regolith")

geoch_results <- list()

for (buf in buffer_distances) {
  buf_km <- buf / 1000
  message("\n--- Buffer: ", buf_km, "km ---")
  
  # Stream sediment
  message("  Stream sediment...")
  res_stream <- buffer_geoch_join(
    all_pts_sf, geoch_stream, buf, "stream"
  )
  if (!is.null(res_stream)) {
    geoch_results[[paste0("stream_", buf_km, "km")]] <- res_stream
  }
  
  # Soil horizon
  message("  Soil horizon...")
  res_horizon <- buffer_geoch_join(
    all_pts_sf, geoch_horizon, buf, "horizon"
  )
  if (!is.null(res_horizon)) {
    geoch_results[[paste0("horizon_", buf_km, "km")]] <- res_horizon
  }
  
  # Soil regolith
  message("  Soil regolith...")
  res_regolith <- buffer_geoch_join(
    all_pts_sf, geoch_regolith, buf, "regolith"
  )
  if (!is.null(res_regolith)) {
    geoch_results[[paste0("regolith_", buf_km, "km")]] <- res_regolith
  }
}

message("\nTotal geochemistry result tables: ", length(geoch_results))


# -------------------------------------------------------------
# SECTION 13: Assemble Geochemistry Matrix
# -------------------------------------------------------------

message("\nAssembling geochemistry matrix...")

# Start with point IDs
geoch_matrix <- all_pts_sf %>%
  st_drop_geometry() %>%
  select(point_id, point_type, period)

# Join all buffer results
for (name in names(geoch_results)) {
  res <- geoch_results[[name]]
  if (!is.null(res)) {
    geoch_matrix <- geoch_matrix %>%
      left_join(res, by = "point_id")
  }
}

message("Geochemistry matrix: ", nrow(geoch_matrix),
        " rows × ", ncol(geoch_matrix), " columns")

# Check coverage at site locations
sites_geoch <- geoch_matrix %>% filter(point_type == "site")
message("Site rows with geochemistry data (5km buffer): ",
        sum(!is.na(sites_geoch[, 4])), "/", nrow(sites_geoch))


# -------------------------------------------------------------
# SECTION 14: Assemble Full Vector Matrix
# -------------------------------------------------------------

message("\nAssembling full vector extraction matrix...")

vector_matrix <- structural_matrix %>%
  left_join(
    geoch_matrix %>% select(-point_type, -period),
    by = "point_id"
  )

message("Full vector matrix: ", nrow(vector_matrix),
        " rows × ", ncol(vector_matrix), " columns")


# -------------------------------------------------------------
# SECTION 15: Save Outputs
# -------------------------------------------------------------

message("\nSaving outputs...")

# Structural geology
write_csv(
  structural_matrix,
  file.path(chapter_root, "05_Structural_Geology",
            "structural_proximity_matrix.csv")
)
message("Saved: structural_proximity_matrix.csv")

# Geochemistry
write_csv(
  geoch_matrix,
  file.path(chapter_root, "04_Geochemistry",
            "geochemistry_buffer_matrix.csv")
)
message("Saved: geochemistry_buffer_matrix.csv")

# Full vector matrix
write_csv(
  vector_matrix,
  file.path(chapter_root, "05_Structural_Geology",
            "vector_extraction_matrix.csv")
)
message("Saved: vector_extraction_matrix.csv")

# RData
save(
  structural_matrix,
  geoch_matrix,
  vector_matrix,
  file = file.path(chapter_root, "05_Structural_Geology",
                   "vector_extraction_matrix.RData")
)
message("Saved: vector_extraction_matrix.RData")


# -------------------------------------------------------------
# SECTION 16: Log Session
# -------------------------------------------------------------

log_file <- file.path(chapter_root, "13_Logs",
                      paste0("script04_log_", Sys.Date(), ".txt"))
sink(log_file)
cat("SCRIPT 04 LOG — Vector Variable Extraction\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Structural variables: dist_fault, dist_dyke,",
    "dist_lineament, dist_shear, dist_mineral\n")
cat("Geochemistry buffers: 1km, 2km, 5km\n")
cat("Geochemistry layers: stream sediment,",
    "soil horizon, soil regolith\n\n")
cat("Structural matrix dimensions:",
    nrow(structural_matrix), "x", ncol(structural_matrix), "\n")
cat("Geochemistry matrix dimensions:",
    nrow(geoch_matrix), "x", ncol(geoch_matrix), "\n")
cat("Full vector matrix dimensions:",
    nrow(vector_matrix), "x", ncol(vector_matrix), "\n")
sink()

message("Log saved: ", log_file)


# -------------------------------------------------------------
# SCRIPT 04 COMPLETE
# -------------------------------------------------------------

message("\n=============================================================")
message("Script 04 complete.")
message("Structural proximity: 5 distance variables")
message("Geochemistry: ", length(geoch_results),
        " buffer-layer combinations")
message("Full vector matrix: ", nrow(vector_matrix),
        " rows × ", ncol(vector_matrix), " columns")
message("Next: Run Script 05 — Taphonomic Bias Correction")
message("  (Borehole IDW interpolation + burial depth surface)")
message("=============================================================")