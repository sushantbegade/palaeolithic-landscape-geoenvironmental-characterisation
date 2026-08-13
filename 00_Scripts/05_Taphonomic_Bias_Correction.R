# =============================================================
# SCRIPT 05 (REBUILT): Taphonomic Bias Correction
# =============================================================
# Project: Reading the Palaeolithic Landscape
# Chapter: Mapping the Past (Springer Nature, 2026)
# Author: Sushant Begade | RTMNU Nagpur
# ORCID: 0009-0003-0804-1763
# Rebuild date: August 2026
# =============================================================
# CHANGES FROM ORIGINAL SCRIPT 05:
#   1. HARD VERIFICATION CHECKPOINT (Section 3) — original script had
#      a comment asserting "length_m is the borehole depth column"
#      with no actual verification. This is a load-bearing assumption
#      for the chapter's headline Section 6.4 finding (LP mean burial
#      depth 137m vs MP 44.7m). Script now prints ALL borehole
#      attribute columns + summary stats and STOPS with an explicit
#      manual-confirmation gate before proceeding — you must check
#      this against actual CGWB borehole log metadata (a report/
#      logsheet describing what length_m records: total drilled
#      depth? depth to bedrock? depth to water strike? sediment
#      thickness specifically? These are NOT interchangeable).
#   2. LOOCV VALIDATION [NEW] — original used IDW power=2 with zero
#      validation, explicitly flagged as a problem in Discussion 5.6
#      text ("uncertainty mapping would be required... future work")
#      but never done. Now: leave-one-out CV via gstat::krige.cv() for
#      IDW, PLUS ordinary kriging with a fitted variogram as an
#      alternative, LOOCV compared, lower-RMSE method used for final
#      surface. This is the "future work" finally executed, not
#      deferred again.
#   3. Kriging (if selected) additionally provides a genuine per-pixel
#      prediction VARIANCE surface — something IDW structurally cannot
#      provide. This directly enables the uncertainty surface the
#      Discussion says is missing.
# =============================================================

load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")
load(file.path(chapter_root, "03_Extracted_Values/master_extraction_matrix.RData"))

library(terra)
library(sf)
library(tidyverse)
library(gstat)
library(readxl)

set.seed(rf_seed)

message("Master matrix loaded: ", nrow(master_matrix), " rows")


# -------------------------------------------------------------
# SECTION 2: Load Borehole Data
# -------------------------------------------------------------

borehole_sf <- st_read(paths$borehole_shp, quiet = TRUE) %>% st_transform(crs = crs_utm44n)
message("\nBorehole shapefile loaded: ", nrow(borehole_sf), " points")

borehole_xl <- read_xlsx(paths$borehole_xlsx)
message("Borehole Excel loaded: ", nrow(borehole_xl), " rows x ", ncol(borehole_xl), " cols")
message("Excel columns: ", paste(names(borehole_xl), collapse = ", "))


# -------------------------------------------------------------
# SECTION 3: SOURCE-DATA VERIFICATION — RESOLVED
# -------------------------------------------------------------
# FINDING (from column audit, August 2026): this dataset is NOT CGWB
# groundwater borehole data as originally assumed in the chapter text.
# Evidence: commodity_ = "Basemetal, Copper"; prospect_n/block_name =
# named mineral exploration prospects (TAMBEGADI PATHARI, MINJHARI
# AREA); geologist_ = named GSI/exploration geologists; stage_of_i =
# G2/G3 (standard Indian mineral-exploration reporting stage codes);
# source_of_ = GSI-style project codes (M2A_MEP/NC/CR/...); cl_inclina
# (45/50 degrees) + bearing = inclined diamond-drillhole survey
# parameters; water_leve is a column but entirely NA (never used for
# groundwater logging).
#
# CONCLUSION: this is GSI mineral-exploration diamond drillhole (DDH)
# data for Cu/base-metal prospecting, not CGWB water-borehole data.
# length_m is total DOWN-HOLE drilled length along an INCLINED
# trajectory — not sediment burial depth, not vertical depth, and not
# collected for any stratigraphic/taphonomic purpose.
#
# CORRECTIVE ACTIONS TAKEN:
#   1. Corrected to true VERTICAL depth using length_m * sin(inclination)
#      — still not "sediment burial depth," but geometrically honest.
#   2. Renamed throughout from "burial_depth" to "borehole_depth" /
#      "penetration_depth" so the reframing can't silently slide back
#      into overclaiming taphonomic validity downstream.
#   3. Chapter Table 1 (Data Sources) and Acknowledgements MUST be
#      corrected: source agency is GSI (mineral exploration drilling
#      records), not CGWB. Flag this for the prose-rewrite pass.
#   4. Discussion/Conclusion language must be softened substantially —
#      this is now an exploratory subsurface-penetration proxy, NOT a
#      validated taphonomic bias correction. The "motivates deep-
#      testing investment" claim from the original draft is NOT
#      defensible on this data and must be removed or heavily hedged.

message("\n=== SOURCE-DATA FINDING (RESOLVED) ===")
message("This is GSI mineral-exploration diamond drillhole data (Cu/base-metal")
message("prospecting), not CGWB groundwater boreholes — confirmed via")
message("commodity_/prospect_n/geologist_/stage_of_i/cl_inclina fields.")
message("length_m = down-hole drilled length along an INCLINED trajectory.")
message("Proceeding with vertical-depth correction + reframed, heavily-")
message("caveated construct: 'borehole penetration depth' NOT 'burial depth'.")

depth_col <- "length_m"
incl_col  <- "cl_inclina"


# -------------------------------------------------------------
# SECTION 4: Clean Borehole Data
# -------------------------------------------------------------

bh_coords <- st_coordinates(borehole_sf)
borehole_clean <- borehole_sf %>% st_drop_geometry() %>%
  mutate(easting = bh_coords[, 1], northing = bh_coords[, 2])

# FIX (was missed in prior version): apply the promised vertical-depth
# correction. length_m is down-hole drilled length along an INCLINED
# trajectory (cl_inclina = angle from horizontal, degrees). True
# vertical penetration = length_m * sin(inclination). Records with
# missing/non-numeric inclination CANNOT be corrected and are excluded
# rather than silently treated as vertical (would understate correction).

borehole_clean$length_m_raw <- suppressWarnings(as.numeric(borehole_clean[[depth_col]]))
borehole_clean$incl_deg     <- suppressWarnings(as.numeric(borehole_clean[[incl_col]]))

n_missing_incl <- sum(is.na(borehole_clean$incl_deg) & !is.na(borehole_clean$length_m_raw))
message("\nRecords with valid length_m but missing/non-numeric inclination: ",
        n_missing_incl, " — these CANNOT be vertical-corrected and are excluded")
message("(treating as vertical would understate their true penetration depth).")

borehole_clean <- borehole_clean %>%
  mutate(depth_m = length_m_raw * sin(incl_deg * pi / 180)) %>%
  filter(!is.na(depth_m), depth_m > 0)

message("Valid borehole records after vertical-depth correction: ", nrow(borehole_clean))
message("Vertical depth range: ", round(min(borehole_clean$depth_m), 2), " to ",
        round(max(borehole_clean$depth_m), 2), " m")
message("(compare to raw length_m range — correction should meaningfully reduce values")
message("since inclination angles were 45-50 degrees in sampled records, sin(45-50deg)",
        " ~ 0.71-0.77, i.e. ~25-30% reduction from raw drilled length)")

if (nrow(borehole_clean) < 30) {
  warning("Fewer than 30 valid borehole records (n=", nrow(borehole_clean),
          "). LOOCV and any interpolation will be unstable at this density.")
}

borehole_sf_clean <- st_as_sf(borehole_clean, coords = c("easting", "northing"), crs = crs_utm44n)
borehole_sp <- as(borehole_sf_clean, "Spatial")
message("Borehole spatial object ready: ", nrow(borehole_sp), " points")


# -------------------------------------------------------------
# SECTION 5: Interpolation Grid
# -------------------------------------------------------------

dem <- rast(paths$dem) %>% terra::project("EPSG:32644")
boundary <- st_read(paths$boundary, quiet = TRUE) %>% st_transform(crs = crs_utm44n)
study_area <- st_union(boundary)
study_area_vect <- vect(study_area)

dem_masked <- dem %>% crop(study_area_vect) %>% mask(study_area_vect)
dem_r <- raster::raster(dem_masked)
grid_sp <- as(dem_r, "SpatialPixelsDataFrame")
message("\nInterpolation grid: ", nrow(grid_sp), " pixels at ", res(dem_masked)[1], "m")


# =============================================================
# SECTION 6: LOOCV — IDW vs Ordinary Kriging  [NEW]
# =============================================================

message("\n=== LEAVE-ONE-OUT CROSS-VALIDATION: IDW vs Ordinary Kriging ===")

depth_formula <- as.formula("depth_m ~ 1")
n_bh <- nrow(borehole_sp)

# --- IDW LOOCV (power = 2, as original) ---
message("\nRunning LOOCV for IDW (power=2)...")
idw_cv <- krige.cv(depth_formula, locations = borehole_sp, set = list(idp = 2), nfold = n_bh)
idw_rmse <- sqrt(mean(idw_cv$residual^2, na.rm = TRUE))
idw_mae  <- mean(abs(idw_cv$residual), na.rm = TRUE)
message("IDW (power=2) LOOCV: RMSE = ", round(idw_rmse, 2),
        " m | MAE = ", round(idw_mae, 2), " m")

# --- Empirical variogram + Ordinary Kriging ---
message("\nFitting variogram for Ordinary Kriging...")
emp_vgm <- variogram(depth_formula, borehole_sp)

# Try a few common models, pick best fit by SSErr
vgm_candidates <- c("Sph", "Exp", "Gau")
fitted_vgms <- lapply(vgm_candidates, function(model_type) {
  tryCatch(
    fit.variogram(emp_vgm, model = vgm(model_type, nugget = TRUE)),
    error = function(e) NULL
  )
})
names(fitted_vgms) <- vgm_candidates
fitted_vgms <- fitted_vgms[!sapply(fitted_vgms, is.null)]

if (length(fitted_vgms) == 0) {
  message("WARNING: no variogram model converged — Ordinary Kriging unavailable,")
  message("falling back to IDW as final method regardless of LOOCV comparison.")
  ok_available <- FALSE
} else {
  sserr <- sapply(fitted_vgms, function(v) attr(v, "SSErr"))
  best_vgm_name <- names(fitted_vgms)[which.min(sserr)]
  best_vgm <- fitted_vgms[[best_vgm_name]]
  message("Best-fit variogram model: ", best_vgm_name, " (SSErr = ",
          round(min(sserr), 2), ")")
  print(best_vgm)
  
  message("\nRunning LOOCV for Ordinary Kriging (", best_vgm_name, " model)...")
  ok_cv <- tryCatch(
    krige.cv(depth_formula, locations = borehole_sp, model = best_vgm, nfold = n_bh),
    error = function(e) { message("OK LOOCV failed: ", e$message); NULL }
  )
  
  if (!is.null(ok_cv)) {
    ok_rmse <- sqrt(mean(ok_cv$residual^2, na.rm = TRUE))
    ok_mae  <- mean(abs(ok_cv$residual), na.rm = TRUE)
    message("Ordinary Kriging LOOCV: RMSE = ", round(ok_rmse, 2),
            " m | MAE = ", round(ok_mae, 2), " m")
    ok_available <- TRUE
  } else {
    ok_available <- FALSE
  }
}

# --- Model Selection ---
cv_comparison <- tibble(
  method = c("IDW (power=2)", if (ok_available) paste0("Ordinary Kriging (", best_vgm_name, ")")),
  RMSE_m = c(idw_rmse, if (ok_available) ok_rmse),
  MAE_m  = c(idw_mae,  if (ok_available) ok_mae)
)
message("\n=== LOOCV COMPARISON ===")
print(cv_comparison)

use_kriging <- ok_available && (ok_rmse < idw_rmse)
message("\nSelected method for final surface: ",
        if (use_kriging) paste0("Ordinary Kriging (", best_vgm_name, ") — lower LOOCV RMSE")
        else "IDW (power=2) — lower or comparable LOOCV RMSE, or kriging unavailable")

write_csv(cv_comparison, file.path(chapter_root, "11_Tables", "Table_burial_depth_LOOCV_comparison.csv"))
message("Saved: Table_burial_depth_LOOCV_comparison.csv [REQUIRED for Methods 5.6 —")
message("report this RMSE, do not present the surface without it]")


# -------------------------------------------------------------
# SECTION 7: Generate Final Surface (+ uncertainty if kriging)
# -------------------------------------------------------------

if (use_kriging) {
  message("\nGenerating Ordinary Kriging surface + prediction variance...")
  message("Applying maxdist limit = variogram range (", round(best_vgm$range[2],0),
          "m) so pixels beyond reliable data support return NA instead of a")
  message("fake-precise reversion to the global mean. Given 93 boreholes over")
  message("21,262 km2, MOST of the study area will be NA — this is the honest")
  message("reflection of sparse coverage, not a bug.")
  krige_maxdist <- best_vgm$range[2]
  ok_result <- krige(depth_formula, locations = borehole_sp, newdata = grid_sp,
                     model = best_vgm, maxdist = krige_maxdist)
  burial_depth_r <- rast(raster::raster(ok_result["var1.pred"]))
  burial_var_r   <- rast(raster::raster(ok_result["var1.var"]))
  terra::crs(burial_depth_r) <- "EPSG:32644"
  terra::crs(burial_var_r)   <- "EPSG:32644"
  burial_depth_r <- mask(burial_depth_r, study_area_vect)
  burial_var_r   <- mask(burial_var_r,   study_area_vect)
  
  pct_covered <- round(100 * sum(!is.na(values(burial_depth_r))) /
                         sum(!is.na(values(dem_masked))), 1)
  message("Surface coverage within maxdist of a borehole: ", pct_covered, "% of study area")
  
  writeRaster(burial_var_r, file.path(chapter_root, "06_Taphonomy", "borehole_depth_variance.tif"),
              overwrite = TRUE)
  message("Saved: borehole_depth_variance.tif")
} else {
  message("\nGenerating IDW surface (power=2), maxdist-limited...")
  krige_maxdist <- 12000  # fallback if kriging unavailable — conservative
  idw_result <- gstat::idw(depth_formula, locations = borehole_sp, newdata = grid_sp,
                           idp = 2, maxdist = krige_maxdist)
  burial_depth_r <- rast(raster::raster(idw_result["var1.pred"]))
  terra::crs(burial_depth_r) <- "EPSG:32644"
  burial_depth_r <- mask(burial_depth_r, study_area_vect)
  burial_var_r <- NULL
  pct_covered <- round(100 * sum(!is.na(values(burial_depth_r))) /
                         sum(!is.na(values(dem_masked))), 1)
  message("Surface coverage within maxdist of a borehole: ", pct_covered, "% of study area")
  message("NOTE: IDW has no native prediction-variance surface. Report LOOCV")
  message("RMSE as the global uncertainty estimate in Methods/Discussion.")
}

if (pct_covered < 30) {
  message("\n*** MAJOR CAVEAT: fewer than 30% of the study area falls within")
  message("*** reliable data support. Most sites will have NO borehole-depth")
  message("*** estimate. Section 6.4 in the chapter MUST be rewritten as a")
  message("*** small-sample exploratory note, NOT a taphonomic bias")
  message("*** correction applied across the full site record. Do not report")
  message("*** period means (LP/MP/UP) unless each period retains a")
  message("*** reasonable n of covered sites — check Section 10 output below.")
}

message("Surface range: ", round(minmax(burial_depth_r)[1], 2), " to ",
        round(minmax(burial_depth_r)[2], 2), " m")


# -------------------------------------------------------------
# SECTION 8: Classify Burial Risk Tiers (unchanged logic)
# -------------------------------------------------------------

depth_vals <- values(burial_depth_r, na.rm = TRUE)
message("\nPixels with a valid (maxdist-covered) estimate: ", length(depth_vals))
q33 <- quantile(depth_vals, 0.33, na.rm = TRUE)
q66 <- quantile(depth_vals, 0.66, na.rm = TRUE)
message("Burial depth quantile breaks: Low/Med=", round(q33,2),
        "m | Med/High=", round(q66,2), "m")

if (isTRUE(all.equal(as.numeric(q33), as.numeric(q66)))) {
  stop("Tertile breaks are degenerate (q33 == q66) — the covered-pixel ",
       "distribution is still spiked at a single value. Do not proceed to ",
       "risk-tier classification; investigate further before continuing.")
}

burial_risk_r <- classify(burial_depth_r, rcl = matrix(
  c(-Inf, q33, 1,  q33, q66, 2,  q66, Inf, 3), ncol = 3, byrow = TRUE))


# -------------------------------------------------------------
# SECTION 9: Extract at All Points (unchanged logic)
# -------------------------------------------------------------

all_pts_df <- master_matrix %>% select(point_id, point_type, period, coord_precision, easting, northing)
all_vect <- vect(all_pts_df, geom = c("easting", "northing"), crs = crs_utm44n)

burial_vals <- terra::extract(burial_depth_r, all_vect, ID = FALSE)
names(burial_vals) <- "burial_depth_m"
risk_vals <- terra::extract(burial_risk_r, all_vect, ID = FALSE)
names(risk_vals) <- "burial_risk_class"

taphonomy_matrix <- bind_cols(all_pts_df, burial_vals, risk_vals) %>%
  mutate(burial_risk_label = case_when(
    burial_risk_class == 1 ~ "Low", burial_risk_class == 2 ~ "Medium",
    burial_risk_class == 3 ~ "High", TRUE ~ NA_character_))

if (use_kriging) {
  var_vals <- terra::extract(burial_var_r, all_vect, ID = FALSE)
  names(var_vals) <- "burial_depth_variance"
  taphonomy_matrix <- bind_cols(taphonomy_matrix, var_vals)
}

message("\nTaphonomy matrix: ", nrow(taphonomy_matrix), " rows x ", ncol(taphonomy_matrix), " cols")


# -------------------------------------------------------------
# SECTION 10: Taphonomic Bias Analysis by Period (unchanged logic)
# -------------------------------------------------------------

message("\n=== TAPHONOMIC BIAS ANALYSIS ===")
message("Note: period here includes LP/MP/UP/MULTI — MULTI sites shown")
message("descriptively but excluded from period-stratified interpretation")
message("per the Script 01 rebuild decision.")

site_risk <- taphonomy_matrix %>% filter(point_type == "site") %>%
  count(period, burial_risk_label) %>%
  pivot_wider(names_from = burial_risk_label, values_from = n, values_fill = 0)
print(site_risk)

# Coverage check — how many sites even have a valid (non-NA) estimate?
site_coverage <- taphonomy_matrix %>% filter(point_type == "site") %>%
  group_by(period) %>%
  summarise(n_total = n(), n_covered = sum(!is.na(burial_depth_m)),
            pct_covered = round(100 * n_covered / n_total, 1), .groups = "drop")
message("\nSite coverage by borehole-depth surface (maxdist-limited):")
print(site_coverage)
message("\nIf pct_covered is low for any period (especially small-n UP),")
message("period-mean comparisons below are based on very few sites and")
message("must be reported with that n explicitly, not as a full-sample stat.")

period_burial_summary <- taphonomy_matrix %>%
  filter(point_type == "site", period %in% c("LP","MP","UP")) %>%
  group_by(period) %>%
  summarise(
    mean_depth_m = round(mean(burial_depth_m, na.rm = TRUE), 2),
    sd_depth_m   = round(sd(burial_depth_m,   na.rm = TRUE), 2),
    pct_high     = round(100 * mean(burial_risk_label == "High", na.rm = TRUE), 1),
    .groups = "drop"
  )
message("\nBurial depth by period (LP/MP/UP only):")
print(period_burial_summary)


# -------------------------------------------------------------
# SECTION 11: Save Outputs
# -------------------------------------------------------------

out_path <- file.path(chapter_root, "06_Taphonomy")

writeRaster(burial_depth_r, file.path(out_path, "burial_depth_surface.tif"), overwrite = TRUE)
writeRaster(burial_risk_r,  file.path(out_path, "burial_risk_classification.tif"), overwrite = TRUE)
write_csv(taphonomy_matrix, file.path(out_path, "taphonomy_matrix.csv"))
write_csv(period_burial_summary, file.path(chapter_root, "11_Tables", "Table_burial_by_period.csv"))
message("\nSaved: burial_depth_surface.tif, burial_risk_classification.tif, taphonomy_matrix.csv")

save(burial_depth_r, burial_risk_r, taphonomy_matrix, q33, q66,
     cv_comparison, use_kriging, period_burial_summary, site_coverage, pct_covered,
     file = file.path(out_path, "taphonomy_data.RData"))
message("Saved: taphonomy_data.RData")


# -------------------------------------------------------------
# SECTION 12: Log Session
# -------------------------------------------------------------

log_file <- file.path(chapter_root, "13_Logs", paste0("script05_log_", Sys.Date(), ".txt"))
sink(log_file)
cat("SCRIPT 05 LOG (REBUILT) — Taphonomic Bias Correction\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("depth_col used:", depth_col, "(manually confirmed by user before running)\n")
cat("Valid borehole points:", nrow(borehole_clean), "\n\n")
cat("LOOCV comparison:\n"); print(cv_comparison)
cat("\nMethod selected:", if (use_kriging) paste0("Ordinary Kriging (", best_vgm_name, ")") else "IDW (power=2)", "\n")
cat("\nBurial depth by period:\n"); print(period_burial_summary)
sink()
message("Log saved: ", log_file)


# -------------------------------------------------------------
# SCRIPT 05 COMPLETE
# -------------------------------------------------------------

message("\n=============================================================")
message("Script 05 (REBUILT) complete.")
message("Method used: ", if (use_kriging) "Ordinary Kriging" else "IDW", " (selected by LOOCV)")
message("LOOCV RMSE now reported — see Table_burial_depth_LOOCV_comparison.csv")
message("This is a REQUIRED citation in Methods 5.6 going forward.")
message("Next: Script 06 — Mann-Whitney U Tests")
message("  (join geochem PCs into analysis_vars, add spatial-block")
message("  permutation test for structural geology + burial depth,")
message("  add Persistent-Places MULTI-vs-background comparison)")
message("=============================================================")