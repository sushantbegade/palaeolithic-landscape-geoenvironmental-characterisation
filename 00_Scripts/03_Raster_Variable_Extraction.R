# =============================================================
# SCRIPT 03 (REBUILT): Raster Variable Extraction
# =============================================================
# Project: Reading the Palaeolithic Landscape
# Chapter: Mapping the Past (Springer Nature, 2026)
# Author: Sushant Begade | RTMNU Nagpur
# ORCID: 0009-0003-0804-1763
# Rebuild date: August 2026
# =============================================================
# CHANGES FROM ORIGINAL SCRIPT 03:
#   1. Extracts from combined_uniform (Script 02 output) instead of
#      rebuilding site/bg dataframes from scratch — carries
#      coord_precision, period (LP/MP/UP/MULTI), bg_source through
#      into master_matrix. Original script dropped these silently.
#   2. NUMERIC COERCION AUDIT [fixes reported bug]: every extracted
#      raster column is forced to numeric with as.numeric(), and any
#      values that fail to coerce (become NA that weren't NA before)
#      are counted and printed per column so the root cause is
#      actually visible instead of silently producing a String-typed
#      CSV column downstream. Original background_extraction_matrix.csv
#      had soil_depth/soil_erosion/soil_productivity/soil_slope/
#      soil_texture/chelsa_*/wc_* as String type while the sites
#      matrix had them as Float — this audit will show why.
#   3. PRIMARY RUN = uniform background only (n=1000). Envelope
#      background (n=728, Script 02) is NOT extracted here — deferred
#      to the Script 06/07 sensitivity pass to avoid double raster IO
#      now. extract_matrix() below is written as a reusable function
#      so that pass just calls it again on combined_envelope later.
#   4. CHELSA-LGM extracted here as before (extraction itself is
#      valid at all points) — the temporal-validity restriction for
#      LP/MP (Script 06 discussion) is an ANALYSIS decision, not an
#      extraction decision, so it stays in Script 06, not here.
# =============================================================

load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")

load(file.path(chapter_root, "01_Data_Processed/sites_validated.RData"))
load(file.path(chapter_root, "02_Background_Points/background_pts.RData"))

library(terra)
library(sf)
library(tidyverse)

set.seed(rf_seed)

message("Parameters and data loaded.")
message("combined_uniform: ", nrow(combined_uniform), " points")


# -------------------------------------------------------------
# SECTION 2: SpatVector from combined_uniform (PRIMARY dataset)
# -------------------------------------------------------------

all_pts_df <- combined_uniform  # point_id, point_type, period,
# coord_precision, easting, northing, bg_source

all_vect <- vect(all_pts_df, geom = c("easting", "northing"), crs = crs_utm44n)
message("SpatVector created: ", nrow(all_vect), " points")


# -------------------------------------------------------------
# SECTION 3: Helper — Safe Raster Extract (unchanged logic)
# -------------------------------------------------------------

safe_extract <- function(raster_path, points_vect, var_names, reproject = TRUE) {
  tryCatch({
    r <- rast(raster_path)
    if (reproject) r <- terra::project(r, crs_utm44n)
    vals <- terra::extract(r, points_vect, ID = FALSE)
    if (ncol(vals) == length(var_names)) {
      names(vals) <- var_names
    } else {
      names(vals) <- paste0(var_names[1], "_", seq_len(ncol(vals)))
    }
    message("  Extracted: ", paste(var_names, collapse = ", "))
    return(vals)
  }, error = function(e) {
    message("  ERROR extracting ", var_names[1], ": ", e$message)
    df <- as.data.frame(matrix(NA, nrow = nrow(all_pts_df), ncol = length(var_names)))
    names(df) <- var_names
    return(df)
  })
}


# -------------------------------------------------------------
# SECTION 4-10: Extract all raster layers (same variables as original)
# -------------------------------------------------------------

message("\nExtracting DEM...")
ext_dem <- safe_extract(paths$dem, all_vect, "elevation")

message("\nExtracting Sentinel-2A bands...")
s2_bands <- list(
  list(path = paths$B01, name = "S2_B01"), list(path = paths$B02, name = "S2_B02"),
  list(path = paths$B03, name = "S2_B03"), list(path = paths$B04, name = "S2_B04"),
  list(path = paths$B05, name = "S2_B05"), list(path = paths$B06, name = "S2_B06"),
  list(path = paths$B07, name = "S2_B07"), list(path = paths$B08, name = "S2_B08"),
  list(path = paths$B09, name = "S2_B09"), list(path = paths$B11, name = "S2_B11"),
  list(path = paths$B12, name = "S2_B12"), list(path = paths$B8A, name = "S2_B8A")
)
ext_s2_bands <- bind_cols(lapply(s2_bands, function(b) safe_extract(b$path, all_vect, b$name)))

message("\nExtracting Sentinel-2A derived indices...")
s2_indices <- list(
  list(path = paths$NDVI,  name = "NDVI"),  list(path = paths$NDWI,  name = "NDWI"),
  list(path = paths$MNDWI, name = "MNDWI"), list(path = paths$NDBI,  name = "NDBI"),
  list(path = paths$BSI,   name = "BSI"),   list(path = paths$SAVI,  name = "SAVI"),
  list(path = paths$MSAVI, name = "MSAVI")
)
ext_indices <- bind_cols(lapply(s2_indices, function(idx) safe_extract(idx$path, all_vect, idx$name)))

message("\nExtracting soil properties...")
soil_layers <- list(
  list(path = paths$soil_depth,        name = "soil_depth"),
  list(path = paths$soil_erosion,      name = "soil_erosion"),
  list(path = paths$soil_productivity, name = "soil_productivity"),
  list(path = paths$soil_slope,        name = "soil_slope"),
  list(path = paths$soil_texture,      name = "soil_texture")
)
ext_soil <- bind_cols(lapply(soil_layers, function(s) safe_extract(s$path, all_vect, s$name)))

message("\nExtracting CHELSA modern climate...")
ext_chelsa_modern <- bind_cols(
  safe_extract(paths$chelsa_bio01_modern, all_vect, "chelsa_bio01_modern"),
  safe_extract(paths$chelsa_bio12_modern, all_vect, "chelsa_bio12_modern"),
  safe_extract(paths$chelsa_bio15_modern, all_vect, "chelsa_bio15_modern")
)

message("\nExtracting CHELSA TraCE21k palaeoclimate (LGM ~20ka)...")
message("  NOTE: extracted at ALL points/periods here. Temporal-validity")
message("  restriction for LP/MP applied later in Script 06, not here.")
ext_chelsa_lgm <- bind_cols(
  safe_extract(paths$chelsa_bio01_lgm, all_vect, "chelsa_bio01_lgm"),
  safe_extract(paths$chelsa_bio12_lgm, all_vect, "chelsa_bio12_lgm"),
  safe_extract(paths$chelsa_bio15_lgm, all_vect, "chelsa_bio15_lgm")
)

message("\nExtracting WorldClim 2.1...")
ext_worldclim <- bind_cols(
  safe_extract(paths$wc_bio01, all_vect, "wc_bio01"),
  safe_extract(paths$wc_bio12, all_vect, "wc_bio12"),
  safe_extract(paths$wc_bio15, all_vect, "wc_bio15")
)


# -------------------------------------------------------------
# SECTION 11: Assemble Master Matrix — metadata RETAINED
# -------------------------------------------------------------

message("\nAssembling master extraction matrix...")

master_matrix <- bind_cols(
  all_pts_df,          # point_id, point_type, period, coord_precision, easting, northing, bg_source
  ext_dem, ext_s2_bands, ext_indices, ext_soil,
  ext_chelsa_modern, ext_chelsa_lgm, ext_worldclim
)

message("Master matrix (pre-coercion): ", nrow(master_matrix), " rows x ",
        ncol(master_matrix), " cols")


# -------------------------------------------------------------
# SECTION 12: NUMERIC COERCION AUDIT  [NEW — fixes reported bug]
# -------------------------------------------------------------

extracted_var_names <- c(
  "elevation",
  sapply(s2_bands, `[[`, "name"),
  sapply(s2_indices, `[[`, "name"),
  sapply(soil_layers, `[[`, "name"),
  "chelsa_bio01_modern","chelsa_bio12_modern","chelsa_bio15_modern",
  "chelsa_bio01_lgm","chelsa_bio12_lgm","chelsa_bio15_lgm",
  "wc_bio01","wc_bio12","wc_bio15"
)

message("\n=== NUMERIC COERCION AUDIT ===")
message("Checking ", length(extracted_var_names), " extracted variables for type issues...")

coercion_report <- tibble(variable = character(), class_before = character(),
                          na_before = integer(), na_after = integer(),
                          new_na_from_coercion = integer())

for (v in extracted_var_names) {
  col_before   <- master_matrix[[v]]
  class_before <- class(col_before)[1]
  na_before    <- sum(is.na(col_before))
  
  col_after <- suppressWarnings(as.numeric(col_before))
  na_after  <- sum(is.na(col_after))
  new_na    <- na_after - na_before
  
  master_matrix[[v]] <- col_after
  
  coercion_report <- bind_rows(coercion_report, tibble(
    variable = v, class_before = class_before,
    na_before = na_before, na_after = na_after,
    new_na_from_coercion = new_na
  ))
  
  if (class_before != "numeric" || new_na > 0) {
    message("  FLAG: ", v, " | was ", class_before,
            " | NA before=", na_before, " -> after=", na_after,
            " | coercion introduced ", new_na, " new NA(s)")
    if (new_na > 0) {
      bad_vals <- col_before[is.na(col_after) & !is.na(col_before)]
      message("    Sample non-numeric values found: ",
              paste(head(unique(bad_vals), 5), collapse = " | "))
    }
  }
}

n_flagged <- sum(coercion_report$class_before != "numeric" |
                   coercion_report$new_na_from_coercion > 0)
message("\nVariables with type/coercion issues: ", n_flagged, " of ", length(extracted_var_names))

write_csv(coercion_report,
          file.path(chapter_root, "11_Tables", "Table_QC_numeric_coercion_audit.csv"))
message("Saved: Table_QC_numeric_coercion_audit.csv [NEW — check this before trusting master_matrix]")

if (n_flagged > 0) {
  message("\n*** DO NOT PROCEED to Script 04 until you've inspected the flagged")
  message("*** variables above. If new_na_from_coercion > 0 for background rows")
  message("*** specifically, the root cause is contamination in the raw raster")
  message("*** or extraction (e.g. locale decimal-comma, out-of-extent string,")
  message("*** factor-coded raster) — not just a downstream typing quirk.")
}


# -------------------------------------------------------------
# SECTION 13: Missing Value Summary (post-coercion)
# -------------------------------------------------------------

na_summary <- master_matrix %>%
  summarise(across(all_of(extracted_var_names), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  filter(n_missing > 0) %>%
  arrange(desc(n_missing))

if (nrow(na_summary) == 0) {
  message("\nNo missing values detected post-coercion.")
} else {
  message("\nVariables with missing values (post-coercion):")
  print(na_summary)
}

sites_matrix <- master_matrix %>% filter(point_type == "site")
bg_matrix    <- master_matrix %>% filter(point_type == "background")

message("\nSite rows: ", nrow(sites_matrix), " | Background rows: ", nrow(bg_matrix))
message("Complete site rows: ", sum(complete.cases(sites_matrix %>% select(all_of(extracted_var_names)))))


# -------------------------------------------------------------
# SECTION 14: Save Outputs
# -------------------------------------------------------------

out_path <- file.path(chapter_root, "03_Extracted_Values")

write_csv(master_matrix, file.path(out_path, "master_extraction_matrix.csv"))
write_csv(sites_matrix,  file.path(out_path, "sites_extraction_matrix.csv"))
write_csv(bg_matrix,     file.path(out_path, "background_extraction_matrix.csv"))
message("\nSaved: master/sites/background_extraction_matrix.csv (all now carry")
message("coord_precision + period[LP/MP/UP/MULTI] + bg_source, all-numeric verified)")

save(master_matrix, sites_matrix, bg_matrix, extracted_var_names, coercion_report,
     file = file.path(out_path, "master_extraction_matrix.RData"))
message("Saved: master_extraction_matrix.RData")

write_csv(na_summary, file.path(chapter_root, "11_Tables", "Table_QC_missing_values.csv"))
message("Saved: Table_QC_missing_values.csv")


# -------------------------------------------------------------
# SECTION 15: Log Session
# -------------------------------------------------------------

log_file <- file.path(chapter_root, "13_Logs", paste0("script03_log_", Sys.Date(), ".txt"))
sink(log_file)
cat("SCRIPT 03 LOG (REBUILT) — Raster Variable Extraction\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Total points:", nrow(master_matrix), " (site=", nrow(sites_matrix),
    ", background=", nrow(bg_matrix), ")\n\n")
cat("Numeric coercion audit — flagged variables:\n")
print(coercion_report %>% filter(class_before != "numeric" | new_na_from_coercion > 0))
cat("\nMissing values post-coercion:\n")
print(na_summary)
sink()
message("Log saved: ", log_file)


# -------------------------------------------------------------
# SCRIPT 03 COMPLETE
# -------------------------------------------------------------

message("\n=============================================================")
message("Script 03 (REBUILT) complete. PRIMARY (uniform bg) extraction done.")
message("Master matrix: ", nrow(master_matrix), " rows x ", ncol(master_matrix), " cols")
message("Numeric coercion audit: ", n_flagged, " variable(s) flagged — see")
message("Table_QC_numeric_coercion_audit.csv. RESOLVE before Script 04 if >0.")
message("Envelope-background extraction deferred to Script 06/07 sensitivity pass.")
message("Next: Script 04 — Vector Variable Extraction")
message("  (structural geology — unchanged; geochemistry — CLR+PCA reduction)")
message("=============================================================")