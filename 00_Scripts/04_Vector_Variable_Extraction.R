# =============================================================
# SCRIPT 04 (REBUILT): Vector Variable Extraction
# =============================================================
# Project: Reading the Palaeolithic Landscape
# Chapter: Mapping the Past (Springer Nature, 2026)
# Author: Sushant Begade | RTMNU Nagpur
# ORCID: 0009-0003-0804-1763
# Rebuild date: August 2026
# =============================================================
# CHANGES FROM ORIGINAL SCRIPT 04:
#   1. Only loads master_extraction_matrix.RData — original loaded
#      sites_validated.RData + background_pts.RData too but never
#      used their contents. Dropping those avoids reloading stale
#      terra SpatVector/SpatRaster objects across sessions (the
#      "external pointer is not valid" crash reported after Script
#      03 traces to exactly this pattern).
#   2. PART A (structural geology proximity) — UNCHANGED, was already
#      correct/verified against mann_whitney_results_full.csv.
#   3. PART B (geochemistry) — MAJOR CHANGE. Original extracted
#      stream/horizon/regolith buffer means at 1/2/5km (633 columns)
#      and saved them, but NEVER joined them into the analysis matrix
#      used by Script 06 (MW tests) or Script 07 (PCA/RF). Abstract
#      claims geochemistry is integrated; it wasn't. FIX:
#        a. Keep buffer-join extraction (unchanged, still useful raw
#           output / supplementary material).
#        b. NEW: CLR (centred log-ratio) transform each layer's
#           elemental composition at the 5km buffer (primary scale,
#           best coverage per original text) — raw ppm/pct
#           concentrations are compositional/closed data, PCA on
#           raw values is not valid (Aitchison 1986, already in refs).
#        c. NEW: PCA on CLR-transformed data per layer (stream,
#           horizon, regolith), extract first 2 PCs each = 6 new
#           variables: geochem_stream_PC1/2, geochem_horizon_PC1/2,
#           geochem_regolith_PC1/2.
#        d. NEW: loadings tables saved per layer so Discussion can
#           name which elements actually drive each geochem PC
#           (required for interpretability, not just a black box).
#        e. These 6 variables get joined into vector_matrix and flow
#           into Script 06/07 analysis_vars — geochemistry finally
#           enters the actual statistical results.
# =============================================================

load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")
load(file.path(chapter_root, "03_Extracted_Values/master_extraction_matrix.RData"))

library(terra)
library(sf)
library(tidyverse)

set.seed(rf_seed)

message("Master matrix loaded: ", nrow(master_matrix), " rows x ", ncol(master_matrix), " cols")


# -------------------------------------------------------------
# SECTION 2: Rebuild Spatial Points from Master Matrix
# -------------------------------------------------------------

all_pts_sf <- st_as_sf(master_matrix, coords = c("easting", "northing"), crs = crs_utm44n)
message("Spatial points rebuilt: ", nrow(all_pts_sf))


# -------------------------------------------------------------
# SECTION 3: Helper — Safe Vector Distance (unchanged)
# -------------------------------------------------------------

safe_distance <- function(vector_path, points_sf, var_name, geom_type = NULL) {
  tryCatch({
    vec <- st_read(vector_path, quiet = TRUE) %>% st_transform(crs = crs_utm44n)
    if (!is.null(geom_type)) vec <- vec %>% st_cast(geom_type, warn = FALSE)
    vec <- vec[!st_is_empty(vec), ]
    dists <- st_distance(points_sf, vec) %>% apply(1, min) %>% as.numeric()
    df <- data.frame(dists); names(df) <- var_name
    message("  Extracted distance: ", var_name, " | Range: ",
            round(min(dists), 0), "-", round(max(dists), 0), " m")
    return(df)
  }, error = function(e) {
    message("  ERROR: ", var_name, " - ", e$message)
    df <- data.frame(rep(NA, nrow(points_sf))); names(df) <- var_name
    return(df)
  })
}


# =============================================================
# PART A: STRUCTURAL GEOLOGICAL PROXIMITY  (unchanged, verified)
# =============================================================

message("\n--- PART A: Structural Geological Proximity ---")

dist_fault     <- safe_distance(paths$fault,     all_pts_sf, "dist_fault")
dist_dyke      <- safe_distance(paths$dyke,      all_pts_sf, "dist_dyke")
dist_lineament <- safe_distance(paths$lineament, all_pts_sf, "dist_lineament")
dist_shear     <- safe_distance(paths$shear,     all_pts_sf, "dist_shear")
dist_mineral   <- safe_distance(paths$mineral,   all_pts_sf, "dist_mineral")

structural_matrix <- bind_cols(
  master_matrix %>% select(point_id, point_type, period, coord_precision, bg_source),
  dist_fault, dist_dyke, dist_lineament, dist_shear, dist_mineral
)

message("\nStructural matrix: ", nrow(structural_matrix), " rows x ", ncol(structural_matrix), " cols")


# =============================================================
# PART B: GEOCHEMISTRY — buffer extraction (unchanged) + CLR-PCA (NEW)
# =============================================================

message("\n--- PART B: Geochemistry ---")

geoch_stream   <- st_read(paths$geoch_stream,   quiet = TRUE) %>% st_transform(crs = crs_utm44n)
geoch_horizon  <- st_read(paths$geoch_horizon,  quiet = TRUE) %>% st_transform(crs = crs_utm44n)
geoch_regolith <- st_read(paths$geoch_regolith, quiet = TRUE) %>% st_transform(crs = crs_utm44n)
message("Stream sediment points: ", nrow(geoch_stream))
message("Soil horizon points: ",    nrow(geoch_horizon))
message("Soil regolith points: ",   nrow(geoch_regolith))


# -------------------------------------------------------------
# SECTION: Helper — Buffer Geochemistry Join (unchanged logic)
# -------------------------------------------------------------

buffer_geoch_join <- function(points_sf, geoch_sf, buffer_m, prefix) {
  tryCatch({
    # FIX (Bug 1): whitelist real elemental columns only. Raw shapefile
    # numeric columns include objectid/gid/longitude/latitude/n_pts —
    # none of these are chemical species and must never enter a CLR
    # compositional transform. Known element/oxide symbol set from
    # your own data dictionary (major oxides + trace/REE suite):
    element_pattern <- paste0(
      "^(sio2|al2o3|fe2o3|tio2|cao|mgo|mno|na2o|k2o|p2o5|loi|",
      "ba|ga|sc|v|th|pb|ni|co|rb|sr|y|zr|nb|cr|cu|zn|au|li|cs|",
      "as_{0,2}|sb|bi|se|ag|be|ge|mo|sn|la|ce|pr|nd|sm|eu|tb|gd|",
      "dy|ho|er|tm|yb|lu|hf|ta|w|u|pt|pd|f|in_{0,2}|tl|te|cd|hg)$"
    )
    
    geoch_numeric_cols <- geoch_sf %>% st_drop_geometry() %>%
      select(where(is.numeric)) %>% names()
    geoch_numeric_cols <- geoch_numeric_cols[
      str_detect(str_to_lower(geoch_numeric_cols), element_pattern)
    ]
    
    if (length(geoch_numeric_cols) == 0) {
      message("  WARNING: no columns matched element whitelist for ", prefix,
              " — check element_pattern against actual column names: ",
              paste(names(geoch_sf %>% st_drop_geometry() %>% select(where(is.numeric))),
                    collapse=", "))
      return(NULL)
    }
    
    pts_buffered <- st_buffer(points_sf, dist = buffer_m)
    joined <- st_join(pts_buffered %>% select(point_id),
                      geoch_sf %>% select(all_of(geoch_numeric_cols)),
                      join = st_intersects)
    
    result <- joined %>% st_drop_geometry() %>% group_by(point_id) %>%
      summarise(across(all_of(geoch_numeric_cols), ~mean(., na.rm = TRUE)),
                n_geoch_pts = n(), .groups = "drop") %>%
      # mean(all-NA, na.rm=TRUE) silently returns NaN not NA — fix so
      # complete.cases() downstream actually catches these rows
      mutate(across(all_of(geoch_numeric_cols), ~if_else(is.nan(.), NA_real_, .)))
    
    buf_label <- paste0("_", prefix, "_", buffer_m / 1000, "km")
    result <- result %>%
      rename_with(~paste0(., buf_label), .cols = all_of(geoch_numeric_cols)) %>%
      rename_with(~paste0("n_pts_", prefix, "_", buffer_m / 1000, "km"), .cols = "n_geoch_pts")
    
    output <- points_sf %>% st_drop_geometry() %>% select(point_id) %>%
      left_join(result, by = "point_id")
    
    message("  Buffer ", buffer_m / 1000, "km (", prefix, "): coverage ",
            sum(!is.na(output[[2]])), "/", nrow(output), " (",
            round(100 * sum(!is.na(output[[2]])) / nrow(output), 1), "%) | ",
            length(geoch_numeric_cols), " element cols")
    return(output)
  }, error = function(e) {
    message("  ERROR in buffer join: ", e$message)
    return(NULL)
  })
}


# -------------------------------------------------------------
# SECTION: Run Buffer Joins — All Distances + All Layers (unchanged)
# -------------------------------------------------------------

message("\nRunning geochemistry buffer joins (1/2/5km x stream/horizon/regolith)...")

geoch_results <- list()
for (buf in buffer_distances) {
  buf_km <- buf / 1000
  message("\n--- Buffer: ", buf_km, "km ---")
  for (layer_name in c("stream", "horizon", "regolith")) {
    geoch_sf <- switch(layer_name, stream = geoch_stream,
                       horizon = geoch_horizon, regolith = geoch_regolith)
    res <- buffer_geoch_join(all_pts_sf, geoch_sf, buf, layer_name)
    if (!is.null(res)) geoch_results[[paste0(layer_name, "_", buf_km, "km")]] <- res
  }
}

geoch_matrix <- all_pts_sf %>% st_drop_geometry() %>% select(point_id, point_type, period)
for (name in names(geoch_results)) {
  res <- geoch_results[[name]]
  if (!is.null(res)) geoch_matrix <- geoch_matrix %>% left_join(res, by = "point_id")
}
message("\nGeochemistry buffer matrix (raw, all scales): ", nrow(geoch_matrix),
        " rows x ", ncol(geoch_matrix), " cols")


# =============================================================
# NEW SECTION: CLR Transform + PCA per Layer at 5km (primary scale)
# =============================================================
# Raw elemental concentrations are compositional (closed, sum to a
# constant, non-negative) — Euclidean PCA on raw ppm/pct values is
# statistically invalid (Aitchison 1986; Filzmoser et al 2009, both
# already in your reference list). CLR transform first, then PCA.

message("\n=== GEOCHEMISTRY CLR + PCA (5km buffer, primary scale) ===")

clr_transform <- function(mat) {
  # mat: numeric matrix, columns = elements, rows = samples. Values > 0 assumed
  # after zero-replacement (below). Returns CLR-transformed matrix.
  log_mat <- log(mat)
  gm <- rowMeans(log_mat, na.rm = TRUE)
  sweep(log_mat, 1, gm, "-")
}

winsorize_cols <- function(mat, probs = c(0.01, 0.99)) {
  # Cap extreme concentrations per column before CLR. A CLR ratio blows
  # up when a trace element is near-zero except a few spikes — one
  # anomalous sample can otherwise dominate an entire PC axis, which is
  # exactly the signature seen in the first run (PC1 = 73-78% driven by
  # ultra-trace elements Cs/Mo/Se/W/Au and Te/In/Tl/Bi/Sb — not a
  # geologically typical grain-size/major-oxide PC1). Winsorizing caps
  # this without discarding the samples entirely.
  for (j in seq_len(ncol(mat))) {
    col <- mat[, j]
    qs <- quantile(col, probs = probs, na.rm = TRUE)
    col[col < qs[1]] <- qs[1]
    col[col > qs[2]] <- qs[2]
    mat[, j] <- col
  }
  mat
}

replace_zeros <- function(mat, min_positive_frac = 0.5) {
  # Zero-replacement (multiplicative, simple alternative to
  # zCompositions::multRepl). FIX (Bug 2): if a column's fraction of
  # positive values is below min_positive_frac (likely below-detection-
  # dominant trace element, e.g. Au/Pt/Pd/Ag/In/Te/Cd/Hg/Tl/Bi), DROP
  # the column entirely rather than silently letting pos_min=Inf force
  # log(0)=-Inf into the SVD. Dropped columns are reported by name.
  dropped <- character(0)
  keep <- rep(TRUE, ncol(mat))
  for (j in seq_len(ncol(mat))) {
    col <- mat[, j]
    frac_pos <- mean(col > 0, na.rm = TRUE)
    if (frac_pos < min_positive_frac) {
      keep[j] <- FALSE
      dropped <- c(dropped, colnames(mat)[j])
      next
    }
    pos_min <- min(col[col > 0], na.rm = TRUE)
    mat[!is.na(col) & col <= 0, j] <- pos_min / 2
  }
  if (length(dropped) > 0) {
    message("    Dropped (below-detection-dominant, <",
            round(min_positive_frac*100), "% positive values): ",
            paste(dropped, collapse = ", "))
  }
  mat[, keep, drop = FALSE]
}

MAJOR_OXIDES <- c("sio2","al2o3","fe2o3","tio2","cao","mgo","mno","na2o","k2o","p2o5")

run_geochem_clr_pca <- function(layer_name, buffer_km = 5, element_group = "major") {
  
  buf_label <- paste0("_", layer_name, "_", buffer_km, "km")
  layer_cols <- names(geoch_matrix)[str_ends(names(geoch_matrix), buf_label)]
  layer_cols <- setdiff(layer_cols, paste0("n_pts", buf_label))
  
  # FIX: split major oxides (well-behaved, standard lithology/weathering
  # signal, archaeologically interpretable — raw-material/landscape
  # geochemical character) from trace elements (near-detection-limit
  # log-ratios naturally dominate variance in mixed CLR-PCA; Reimann et
  # al 2008, Filzmoser et al 2009 — both already in refs). Confirmed
  # empirically: winsorizing barely changed PC1 (74.2% vs 73.5% stream,
  # 78.3% vs 77.6% horizon) — this is real compositional-data behaviour,
  # not an outlier artifact. element_group = "major" is PRIMARY (feeds
  # Script 06/07 analysis_vars). "trace" is SUPPLEMENTARY EXPLORATORY
  # ONLY — saved but not joined into vector_matrix.
  element_base <- str_remove(layer_cols, buf_label)
  if (element_group == "major") {
    layer_cols <- layer_cols[element_base %in% MAJOR_OXIDES]
  } else {
    layer_cols <- layer_cols[!(element_base %in% MAJOR_OXIDES)]
  }
  
  message("\n--- ", layer_name, " @ ", buffer_km, "km [", element_group, "]: ",
          length(layer_cols), " elemental columns ---")
  
  if (length(layer_cols) == 0) {
    message("  No columns found for this layer/buffer — skipping.")
    return(NULL)
  }
  
  sub <- geoch_matrix %>% select(point_id, all_of(layer_cols))
  complete_idx <- complete.cases(sub %>% select(-point_id))
  n_complete <- sum(complete_idx)
  message("  Complete-case coverage: ", n_complete, "/", nrow(sub),
          " (", round(100 * n_complete / nrow(sub), 1), "%)")
  
  if (n_complete < 30) {
    message("  WARNING: fewer than 30 complete cases — PCA unstable, skipping.")
    return(NULL)
  }
  
  mat <- as.matrix(sub[complete_idx, layer_cols])
  mat <- winsorize_cols(mat)  # cap outlier concentrations BEFORE zero-replace/CLR
  mat <- replace_zeros(mat)   # may drop below-detection-dominant columns
  if (ncol(mat) < 3) {
    message("  WARNING: fewer than 3 usable elements after dropping ",
            "below-detection columns — skipping this layer.")
    return(NULL)
  }
  clr_mat <- clr_transform(mat)
  stopifnot("CLR matrix must be finite before PCA" = all(is.finite(clr_mat)))
  
  pca_fit <- prcomp(clr_mat, center = TRUE, scale. = FALSE)
  var_pct <- round(100 * summary(pca_fit)$importance[2, 1:2], 1)
  message("  PC1 = ", var_pct[1], "% | PC2 = ", var_pct[2], "% variance explained")
  if (element_group == "major" && var_pct[1] > 70) {
    message("  Note: PC1 >70% with only ~10 major oxides is plausible")
    message("  (fewer variables = naturally higher PC1 share) — check")
    message("  loadings are geologically sensible rather than flagging by")
    message("  variance threshold alone.")
  } else if (element_group == "trace" && var_pct[1] > 55) {
    message("  Expected: trace-element PC1 often dominated by near-")
    message("  detection-limit log-ratio variance (Reimann et al 2008).")
    message("  Reported as supplementary/exploratory only, not in main model.")
  }
  
  # Loadings — top 5 elements per PC for Discussion/interpretation
  loadings <- as.data.frame(pca_fit$rotation[, 1:2]) %>%
    rownames_to_column("element") %>%
    mutate(element = str_remove(element, buf_label))
  top_pc1 <- loadings %>% arrange(desc(abs(PC1))) %>% head(5)
  top_pc2 <- loadings %>% arrange(desc(abs(PC2))) %>% head(5)
  message("  Top 5 elements on PC1: ", paste(top_pc1$element, collapse = ", "))
  message("  Top 5 elements on PC2: ", paste(top_pc2$element, collapse = ", "))
  
  scores <- as.data.frame(pca_fit$x[, 1:2])
  names(scores) <- paste0("geochem_", layer_name, "_", element_group, "_PC", 1:2)
  scores$point_id <- sub$point_id[complete_idx]
  
  list(
    scores   = scores,
    loadings = loadings,
    var_pct  = var_pct,
    n_complete = n_complete,
    element_group = element_group
  )
}

geochem_pca_results <- list(
  stream_major   = run_geochem_clr_pca("stream",   5, "major"),
  horizon_major  = run_geochem_clr_pca("horizon",  5, "major"),
  regolith_major = run_geochem_clr_pca("regolith", 5, "major"),
  stream_trace   = run_geochem_clr_pca("stream",   5, "trace"),
  horizon_trace  = run_geochem_clr_pca("horizon",  5, "trace"),
  regolith_trace = run_geochem_clr_pca("regolith", 5, "trace")
)

# PRIMARY variables (feed Script 06/07) = major-oxide PCs only.
# Trace-element PCs saved to CSV for supplementary/exploratory use,
# NOT joined into vector_matrix / analysis_vars.
primary_geochem_names <- c("stream_major", "horizon_major", "regolith_major")

geochem_pc_matrix <- all_pts_sf %>% st_drop_geometry() %>% select(point_id)
for (nm in primary_geochem_names) {
  res <- geochem_pca_results[[nm]]
  if (!is.null(res)) {
    geochem_pc_matrix <- geochem_pc_matrix %>% left_join(res$scores, by = "point_id")
  }
}

geochem_pc_vars <- setdiff(names(geochem_pc_matrix), "point_id")
message("\nPRIMARY geochemistry PC variables (major oxides, feed Script 06/07): ",
        paste(geochem_pc_vars, collapse = ", "))
message("Coverage per PC variable:")
print(geochem_pc_matrix %>% summarise(across(all_of(geochem_pc_vars), ~sum(!is.na(.)))))

# Supplementary trace-element PCs — saved separately, exploratory only
trace_geochem_matrix <- all_pts_sf %>% st_drop_geometry() %>% select(point_id)
for (nm in c("stream_trace", "horizon_trace", "regolith_trace")) {
  res <- geochem_pca_results[[nm]]
  if (!is.null(res)) {
    trace_geochem_matrix <- trace_geochem_matrix %>% left_join(res$scores, by = "point_id")
  }
}
write_csv(trace_geochem_matrix,
          file.path(chapter_root, "04_Geochemistry", "geochem_trace_PCs_SUPPLEMENTARY.csv"))
message("Saved: geochem_trace_PCs_SUPPLEMENTARY.csv [exploratory only, NOT in main model]")


# -------------------------------------------------------------
# SECTION: Assemble Full Vector Matrix — structural + geochem PCs
# -------------------------------------------------------------

message("\nAssembling full vector extraction matrix (structural + geochem PCs)...")

vector_matrix <- structural_matrix %>%
  left_join(geochem_pc_matrix, by = "point_id")

message("Full vector matrix: ", nrow(vector_matrix), " rows x ", ncol(vector_matrix), " cols")
message("This now includes ", length(geochem_pc_vars),
        " geochemistry PC variables ready for Script 06/07 analysis_vars —")
message("geochemistry is no longer extracted-and-abandoned.")


# -------------------------------------------------------------
# SECTION: Save Outputs
# -------------------------------------------------------------

write_csv(structural_matrix,
          file.path(chapter_root, "05_Structural_Geology", "structural_proximity_matrix.csv"))
write_csv(geoch_matrix,
          file.path(chapter_root, "04_Geochemistry", "geochemistry_buffer_matrix.csv"))
write_csv(vector_matrix,
          file.path(chapter_root, "05_Structural_Geology", "vector_extraction_matrix.csv"))
message("\nSaved: structural_proximity_matrix.csv, geochemistry_buffer_matrix.csv, vector_extraction_matrix.csv")

# NEW: geochem PCA diagnostics for Methods/Supplementary
geochem_variance_table <- bind_rows(lapply(names(geochem_pca_results), function(nm) {
  res <- geochem_pca_results[[nm]]
  if (is.null(res)) return(NULL)
  tibble(layer_group = nm, PC1_var_pct = res$var_pct[1], PC2_var_pct = res$var_pct[2],
         n_complete = res$n_complete, element_group = res$element_group,
         is_primary = nm %in% primary_geochem_names)
}))
write_csv(geochem_variance_table,
          file.path(chapter_root, "04_Geochemistry", "geochem_PCA_variance_explained.csv"))
message("Saved: geochem_PCA_variance_explained.csv [cite in Methods 5.3 — is_primary",
        " column marks major-oxide PCs used in main analysis]")

geochem_loadings_table <- bind_rows(lapply(names(geochem_pca_results), function(nm) {
  res <- geochem_pca_results[[nm]]
  if (is.null(res)) return(NULL)
  res$loadings %>% mutate(layer_group = nm, .before = 1)
}))
write_csv(geochem_loadings_table,
          file.path(chapter_root, "04_Geochemistry", "geochem_PCA_loadings.csv"))
message("Saved: geochem_PCA_loadings.csv [use for Discussion — name the elements]")

save(structural_matrix, geoch_matrix, vector_matrix,
     geochem_pca_results, geochem_pc_matrix, geochem_pc_vars,
     geochem_variance_table, geochem_loadings_table,
     file = file.path(chapter_root, "05_Structural_Geology", "vector_extraction_matrix.RData"))
message("Saved: vector_extraction_matrix.RData")


# -------------------------------------------------------------
# SECTION: Log Session
# -------------------------------------------------------------

log_file <- file.path(chapter_root, "13_Logs", paste0("script04_log_", Sys.Date(), ".txt"))
sink(log_file)
cat("SCRIPT 04 LOG (REBUILT) — Vector Variable Extraction\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Structural variables: dist_fault, dist_dyke, dist_lineament, dist_shear, dist_mineral\n\n")
cat("Geochemistry CLR-PCA (5km buffer) — major oxides = PRIMARY, trace = supplementary:\n")
print(geochem_variance_table)
cat("\nTop loadings per layer/group:\n")
print(geochem_loadings_table %>% group_by(layer_group) %>% slice_max(abs(PC1), n = 3))
sink()
message("Log saved: ", log_file)


# -------------------------------------------------------------
# SCRIPT 04 COMPLETE
# -------------------------------------------------------------

message("\n=============================================================")
message("Script 04 (REBUILT) complete.")
message("Structural proximity: 5 distance variables (unchanged)")
message("Geochemistry: ", length(geochem_pc_vars),
        " PRIMARY major-oxide CLR-PCA variables NOW LIVE")
message("  (+ 6 supplementary trace-element PCs saved separately, exploratory only)")
message("  (previously extracted but never analysed — abstract/results")
message("  mismatch from earlier assessment is now resolved)")
message("Full vector matrix: ", nrow(vector_matrix), " rows x ", ncol(vector_matrix), " cols")
message("Next: Script 05 — Taphonomic Bias Correction")
message("  (add LOOCV validation for IDW burial-depth surface,")
message("  verify length_m really means burial depth before trusting it)")
message("=============================================================")