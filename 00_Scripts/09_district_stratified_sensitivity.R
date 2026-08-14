# =============================================================
# SCRIPT 09 (NEW): District-Stratified Sensitivity Analysis
# =============================================================
# Project: Reading the Palaeolithic Landscape
# Chapter: Mapping the Past (Springer Nature, 2026)
# Author: Sushant Begade | RTMNU Nagpur
# ORCID: 0009-0003-0804-1763
# Created: August 2026
# =============================================================
# PURPOSE: Priority 1 outstanding limitation (manuscript Section
# 7.7). LP is 93.8% concentrated in Chandrapur District, MP is
# 88.7% concentrated in Nagpur District (Script 01). District may
# correlate with underlying geology, structural density, soils,
# drainage, and survey history independent of cultural period.
# This tests whether the headline LP-MP structural reversal
# (Script 06, Section 6.2 of manuscript) survives when period and
# district are disentangled by within-district stratification,
# comparing:
#   Analysis A: LP sites in Chandrapur vs background in Chandrapur
#   Analysis B: MP sites in Nagpur    vs background in Nagpur
# against the pooled (cross-district) result already reported.
#
# CRITICAL: background points (Script 02) were generated across the
# FULL study area boundary and were never assigned a district field.
# Section 1 of this script performs that assignment for the first
# time via spatial join against the district boundary polygon.
# =============================================================

load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")
load(file.path(chapter_root, "07_Statistics/mann_whitney_results.RData"))

library(sf)
library(tidyverse)

set.seed(rf_seed)

message("Loaded full_matrix_coords: ", nrow(full_matrix_coords), " rows")
message("Analysis variables: ", length(analysis_vars))


# -------------------------------------------------------------
# SECTION 1: Assign District to Every Point via Spatial Join
# -------------------------------------------------------------
# Sites already have a district field from Script 01 (administrative
# source), but background points do not. Both are reassigned here
# from the SAME boundary polygon so site and background district
# labels are on an identical, verifiable basis (not mixing an
# administrative source-record label for sites with a spatial-join
# label for background, which would itself be a subtle confound).

boundary_raw <- st_read(paths$boundary, quiet = TRUE) %>% st_transform(crs_utm44n)

message("\nBoundary shapefile columns: ", paste(names(boundary_raw), collapse = ", "))

# Auto-detect the district field name. EDIT district_field_name below
# manually if none of these guesses match your actual column — the
# printed column list above tells you the real name.
district_field_candidates <- c("district", "District", "DISTRICT", "District_N",
                               "DIST_NAME", "dist_name", "NAME_2", "name")
district_field_name <- intersect(district_field_candidates, names(boundary_raw))[1]

if (is.na(district_field_name)) {
  stop("Could not auto-detect district field name. Boundary columns are: ",
       paste(names(boundary_raw), collapse = ", "),
       " — set district_field_name manually to the correct column and rerun.")
}
message("Using district field: '", district_field_name, "'")

boundary_raw <- boundary_raw %>% rename(district_label = all_of(district_field_name))
message("District values found: ", paste(unique(boundary_raw$district_label), collapse = ", "))

all_pts_sf <- full_matrix_coords %>%
  st_as_sf(coords = c("easting", "northing"), crs = crs_utm44n, remove = FALSE)

message("\nSpatially joining all ", nrow(all_pts_sf), " points to district polygons...")
pts_with_district <- st_join(all_pts_sf, boundary_raw %>% select(district_label), join = st_within)

n_unmatched <- sum(is.na(pts_with_district$district_label))
message("Points with no district match (outside all polygons): ", n_unmatched)
if (n_unmatched > 0) {
  message("These will be excluded from district-stratified analysis. Check for")
  message("boundary/point CRS or edge-snapping issues if this number is large.")
}

full_matrix_district <- pts_with_district %>%
  st_drop_geometry() %>%
  filter(!is.na(district_label)) %>%
  mutate(district_label = str_to_title(str_trim(district_label)))

message("\nDistrict assignment complete: ", nrow(full_matrix_district), " points retained")
message("\nPoint counts by district x point_type:")
print(table(full_matrix_district$district_label, full_matrix_district$point_type))


# -------------------------------------------------------------
# SECTION 2: Sanity Check — Site District Assignment vs Script 01
# -------------------------------------------------------------
# Sites already carry an administrative district label from Script 01
# source records. This spatial-join label should match almost
# perfectly; large disagreement would indicate a boundary/CRS problem
# that must be resolved before trusting anything downstream.

sites_check <- full_matrix_district %>% filter(point_type == "site")
message("\n=== SANITY CHECK: spatial-join district vs Script 01 source district ===")
message("(Comparison requires sites_validated.csv — loading for cross-check)")

sites_validated <- read_csv(file.path(chapter_root, "01_Data_Processed/sites_validated.csv"),
                            show_col_types = FALSE) %>%
  mutate(point_id = paste0("SITE_", sprintf("%03d", site_id))) %>%
  select(point_id, district_source = district)

sites_compare <- sites_check %>%
  left_join(sites_validated, by = "point_id") %>%
  mutate(match = str_to_title(str_trim(district_label)) == str_to_title(str_trim(district_source)))

n_mismatch <- sum(!sites_compare$match, na.rm = TRUE)
message("Site district mismatches (spatial join vs Script 01 source): ", n_mismatch,
        " of ", nrow(sites_compare))
if (n_mismatch > 5) {
  message("*** WARNING: more than 5 mismatches. Inspect boundary CRS/geometry before")
  message("*** trusting the district-stratified results below.")
  print(sites_compare %>% filter(!match) %>% select(point_id, district_label, district_source))
}


# -------------------------------------------------------------
# SECTION 3: District-Stratified Mann-Whitney — Headline Structural Vars
# -------------------------------------------------------------

message("\n=== DISTRICT-STRATIFIED ANALYSIS: headline structural variables ===")

run_mw_stratified <- function(data, variable, period_filter, district_filter) {
  test_data <- data %>%
    filter(district_label == district_filter) %>%
    filter(point_type == "background" |
             (point_type == "site" & period == period_filter)) %>%
    mutate(group = ifelse(point_type == "site", "site", "background")) %>%
    select(group, value = all_of(variable)) %>%
    filter(!is.na(value))
  
  n_site <- sum(test_data$group == "site")
  n_bg   <- sum(test_data$group == "background")
  
  if (n_site < 5 || n_bg < 5) {
    return(tibble(variable = variable, period = period_filter, district = district_filter,
                  n_site = n_site, n_bg = n_bg, median_site = NA, median_bg = NA,
                  p_value = NA, direction = NA,
                  note = "insufficient n (<5) for reliable test"))
  }
  
  test_result <- tryCatch(
    wilcox.test(value ~ group, data = test_data, exact = FALSE, correct = TRUE),
    error = function(e) NULL
  )
  if (is.null(test_result)) {
    return(tibble(variable = variable, period = period_filter, district = district_filter,
                  n_site = n_site, n_bg = n_bg, median_site = NA, median_bg = NA,
                  p_value = NA, direction = NA, note = "test failed"))
  }
  
  med_site <- median(test_data$value[test_data$group == "site"], na.rm = TRUE)
  med_bg   <- median(test_data$value[test_data$group == "background"], na.rm = TRUE)
  
  tibble(variable = variable, period = period_filter, district = district_filter,
         n_site = n_site, n_bg = n_bg,
         median_site = round(med_site, 2), median_bg = round(med_bg, 2),
         p_value = round(test_result$p.value, 4),
         direction = ifelse(med_site > med_bg, "higher_at_sites", "lower_at_sites"),
         note = NA_character_)
}

headline_vars <- c("dist_fault", "dist_dyke", "dist_lineament", "dist_shear")

# Analysis A: LP within Chandrapur (primary confirmatory test)
message("\n--- Analysis A: LP within Chandrapur (n=45 of 48 LP sites) ---")
stratA <- bind_rows(lapply(headline_vars, run_mw_stratified,
                           data = full_matrix_district, period_filter = "LP",
                           district_filter = "Chandrapur"))
print(stratA)

# Analysis B: MP within Nagpur (primary confirmatory test)
message("\n--- Analysis B: MP within Nagpur (n=55 of 62 MP sites) ---")
stratB <- bind_rows(lapply(headline_vars, run_mw_stratified,
                           data = full_matrix_district, period_filter = "MP",
                           district_filter = "Nagpur"))
print(stratB)

# Analysis C: cross-district minority cells (small n, exploratory only)
message("\n--- Analysis C (exploratory, small n): LP within Nagpur (n=3) ---")
stratC <- bind_rows(lapply(headline_vars, run_mw_stratified,
                           data = full_matrix_district, period_filter = "LP",
                           district_filter = "Nagpur"))
print(stratC)

message("\n--- Analysis D (exploratory, small n): MP within Chandrapur (n=7) ---")
stratD <- bind_rows(lapply(headline_vars, run_mw_stratified,
                           data = full_matrix_district, period_filter = "MP",
                           district_filter = "Chandrapur"))
print(stratD)


# -------------------------------------------------------------
# SECTION 4: Pooled vs District-Stratified Comparison Table
# -------------------------------------------------------------
# This is the table the manuscript needs (reviewer's suggested
# format, Priority 1). Pooled values are pulled directly from the
# already-computed Script 06 mw_table for consistency — not
# recomputed here, so any discrepancy would indicate a real bug.

message("\n=== BUILDING POOLED VS DISTRICT-STRATIFIED COMPARISON TABLE ===")

pooled_lp <- mw_table %>% filter(period == "LP", variable %in% headline_vars) %>%
  select(variable, pooled_LP_p = p_value, pooled_LP_direction = direction,
         pooled_LP_median_site = median_site)
pooled_mp <- mw_table %>% filter(period == "MP", variable %in% headline_vars) %>%
  select(variable, pooled_MP_p = p_value, pooled_MP_direction = direction,
         pooled_MP_median_site = median_site)

comparison_table <- stratA %>%
  select(variable, chandrapur_LP_p = p_value, chandrapur_LP_direction = direction,
         chandrapur_LP_n_site = n_site) %>%
  left_join(pooled_lp, by = "variable") %>%
  left_join(
    stratB %>% select(variable, nagpur_MP_p = p_value, nagpur_MP_direction = direction,
                      nagpur_MP_n_site = n_site),
    by = "variable"
  ) %>%
  left_join(pooled_mp, by = "variable") %>%
  mutate(
    LP_direction_agrees = pooled_LP_direction == chandrapur_LP_direction,
    LP_still_significant_stratified = chandrapur_LP_p < alpha,
    MP_direction_agrees = pooled_MP_direction == nagpur_MP_direction,
    MP_still_significant_stratified = nagpur_MP_p < alpha
  )

message("\n=== POOLED vs DISTRICT-STRATIFIED: full comparison ===")
print(comparison_table, width = Inf)

n_lp_survives <- sum(comparison_table$LP_direction_agrees & comparison_table$LP_still_significant_stratified, na.rm = TRUE)
n_mp_survives <- sum(comparison_table$MP_direction_agrees & comparison_table$MP_still_significant_stratified, na.rm = TRUE)

message("\n=== VERDICT ===")
message("LP structural variables retaining same direction AND significance within")
message("Chandrapur alone: ", n_lp_survives, " of ", length(headline_vars))
message("MP structural variables retaining same direction AND significance within")
message("Nagpur alone: ", n_mp_survives, " of ", length(headline_vars))

if (n_lp_survives == length(headline_vars) && n_mp_survives == length(headline_vars)) {
  message("\nRESULT: the LP-MP structural reversal SURVIVES district stratification")
  message("in full. This substantially strengthens the manuscript's central")
  message("structural claim — the reversal is not solely a district-level artefact.")
} else if (n_lp_survives == 0 && n_mp_survives == 0) {
  message("\nRESULT: the LP-MP structural reversal DOES NOT SURVIVE district")
  message("stratification. This is a major finding in its own right — Section 7.7")
  message("anticipated this possibility. The manuscript's Section 6.2/7.1 structural")
  message("claims MUST be rewritten to state that the apparent period effect is")
  message("confounded with, and may be substantially attributable to, district-level")
  message("geographic variation rather than cultural-period behaviour specifically.")
} else {
  message("\nRESULT: PARTIAL survival — some structural variables hold within-district,")
  message("others do not. Report per-variable in the manuscript (use the comparison")
  message("table above directly), do not generalise to 'the structural finding' as")
  message("a whole. This is the most likely real-world outcome and is still a")
  message("legitimate, reportable result.")
}


# -------------------------------------------------------------
# SECTION 5: Save Outputs
# -------------------------------------------------------------

out_path <- file.path(chapter_root, "07_Statistics")

write_csv(full_matrix_district %>% select(point_id, point_type, period, district_label),
          file.path(out_path, "point_district_assignment.csv"))
write_csv(bind_rows(
  stratA %>% mutate(analysis = "A_LP_in_Chandrapur"),
  stratB %>% mutate(analysis = "B_MP_in_Nagpur"),
  stratC %>% mutate(analysis = "C_LP_in_Nagpur_exploratory"),
  stratD %>% mutate(analysis = "D_MP_in_Chandrapur_exploratory")
), file.path(out_path, "district_stratified_results_full.csv"))
write_csv(comparison_table, file.path(chapter_root, "11_Tables", "Table_district_stratified_comparison.csv"))

message("\nSaved: point_district_assignment.csv, district_stratified_results_full.csv,")
message("Table_district_stratified_comparison.csv [THIS is the manuscript Section 7.7")
message("Priority 1 table — insert into Section 6.2/7.1 rewrite depending on verdict]")

save(full_matrix_district, stratA, stratB, stratC, stratD, comparison_table,
     n_lp_survives, n_mp_survives,
     file = file.path(out_path, "district_stratified_analysis.RData"))
message("Saved: district_stratified_analysis.RData")


# -------------------------------------------------------------
# SECTION 6: Log Session
# -------------------------------------------------------------

log_file <- file.path(chapter_root, "13_Logs", paste0("script09_log_", Sys.Date(), ".txt"))
sink(log_file)
cat("SCRIPT 09 LOG — District-Stratified Sensitivity Analysis (Priority 1)\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("District field used:", district_field_name, "\n")
cat("Site district mismatches (spatial join vs Script 01 source):", n_mismatch, "\n\n")
cat("Pooled vs district-stratified comparison:\n")
print(comparison_table)
cat("\nLP survives (direction+sig) within Chandrapur:", n_lp_survives, "of", length(headline_vars), "\n")
cat("MP survives (direction+sig) within Nagpur:", n_mp_survives, "of", length(headline_vars), "\n")
sink()
message("Log saved: ", log_file)


# -------------------------------------------------------------
# SCRIPT 09 COMPLETE
# -------------------------------------------------------------

message("\n=============================================================")
message("Script 09 complete. Priority 1 outstanding limitation now tested.")
message("LP-in-Chandrapur survival: ", n_lp_survives, "/", length(headline_vars))
message("MP-in-Nagpur survival: ", n_mp_survives, "/", length(headline_vars))
message("Manuscript Sections 2, 6.2, 7.1, 7.7, and 8 (Conclusion) must be updated")
message("to state this result explicitly and definitively — replacing the current")
message("'has not yet been disentangled' language with the actual outcome.")
message("Remaining outstanding analyses (manuscript Section 7.7, Priorities 2-7)")
message("still pending: survey-envelope comparison, spatial block-size sensitivity,")
message("PCoA missingness diagnostic, PERMANOVA/PERMDISP, geochem buffer-scale")
message("sensitivity, zero-replacement sensitivity, manual classification check.")
message("=============================================================")