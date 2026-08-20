# =============================================================
# SCRIPT 12 (NEW): Survey-Envelope Sensitivity — Structural Geology
# =============================================================
# Project: Reading the Palaeolithic Landscape
# Chapter: Mapping the Past (Springer Nature, 2026)
# Author: Sushant Begade | RTMNU Nagpur
# ORCID: 0009-0003-0804-1763
# Created: August 2026
# =============================================================
# PURPOSE: Priority 3 of the outstanding-analyses list. Scripts 09-11
# established that of the chapter's structural findings, only LP
# fault and LP lineament proximity survive spatial block-bootstrap,
# district stratification, AND hold stable relative importance under
# district-transfer RF. This is the last unresolved test bearing on
# those two survivors: does the result hold against the 728-point
# survey-envelope background (3km per-site buffers, Script 02),
# representing a more realistic "landscape actually surveyed" null,
# or does it depend on comparing sites against the full unsurveyed
# 21,262 km2 boundary?
#
# SCOPE NOTE: this script extracts structural distances only for the
# envelope background points (never done before — Script 03 deferred
# this explicitly). Full 30-variable extraction (esp. geochemistry,
# which would require re-projecting onto Script 04's saved PCA
# rotation, not retained) is NOT done here — out of scope for the
# specific open question (structural geology), and a larger task
# left for future work if the full envelope comparison is ever needed.
# =============================================================

load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")
load(file.path(chapter_root, "02_Background_Points/background_pts.RData"))
load(file.path(chapter_root, "07_Statistics/mann_whitney_results.RData"))

library(sf)
library(tidyverse)

set.seed(rf_seed)

message("Loaded combined_envelope: ", nrow(combined_envelope), " points",
        " (should be 197 sites + 728 envelope background = 925)")
stopifnot("combined_envelope row count unexpected — check Script 02 output" =
            nrow(combined_envelope) == 925)


# -------------------------------------------------------------
# SECTION 1: Extract Structural Distances for Envelope Background
# -------------------------------------------------------------
# Sites already have dist_fault etc. from Script 04 (full_matrix_coords) —
# only the 728 envelope background points need extraction here, since
# Script 03/04 only ever processed the UNIFORM background (n=1000).

envelope_vect <- st_as_sf(combined_envelope, coords = c("easting","northing"),
                          crs = crs_utm44n, remove = FALSE)

safe_distance <- function(vector_path, points_sf, var_name) {
  tryCatch({
    vec <- st_read(vector_path, quiet = TRUE) %>% st_transform(crs = crs_utm44n)
    vec <- vec[!st_is_empty(vec), ]
    dists <- st_distance(points_sf, vec) %>% apply(1, min) %>% as.numeric()
    df <- data.frame(dists); names(df) <- var_name
    message("  Extracted: ", var_name, " | Range: ", round(min(dists),0), "-", round(max(dists),0), " m")
    return(df)
  }, error = function(e) {
    message("  ERROR: ", var_name, " - ", e$message)
    df <- data.frame(rep(NA, nrow(points_sf))); names(df) <- var_name
    return(df)
  })
}

message("\nExtracting structural distances for all 925 envelope-set points...")
dist_fault_e     <- safe_distance(paths$fault,     envelope_vect, "dist_fault")
dist_dyke_e      <- safe_distance(paths$dyke,      envelope_vect, "dist_dyke")
dist_lineament_e <- safe_distance(paths$lineament, envelope_vect, "dist_lineament")
dist_shear_e     <- safe_distance(paths$shear,     envelope_vect, "dist_shear")

envelope_structural <- bind_cols(
  combined_envelope %>% select(point_id, point_type, period, bg_source),
  dist_fault_e, dist_dyke_e, dist_lineament_e, dist_shear_e
)

message("\nEnvelope structural matrix: ", nrow(envelope_structural), " rows")


# -------------------------------------------------------------
# SECTION 2: Mann-Whitney — LP/MP vs Envelope Background
# -------------------------------------------------------------

run_mw_envelope <- function(data, variable, period_filter) {
  test_data <- data %>%
    filter(point_type == "background" | (point_type == "site" & period == period_filter)) %>%
    mutate(group = ifelse(point_type == "site", "site", "background")) %>%
    select(group, value = all_of(variable)) %>% filter(!is.na(value))
  
  n_site <- sum(test_data$group == "site"); n_bg <- sum(test_data$group == "background")
  if (n_site < 5 || n_bg < 5) {
    return(tibble(variable = variable, period = period_filter, n_site = n_site, n_bg = n_bg,
                  median_site = NA, median_bg = NA, p_value = NA, direction = NA,
                  note = "insufficient n"))
  }
  
  test_result <- wilcox.test(value ~ group, data = test_data, exact = FALSE, correct = TRUE)
  med_site <- median(test_data$value[test_data$group=="site"], na.rm=TRUE)
  med_bg   <- median(test_data$value[test_data$group=="background"], na.rm=TRUE)
  
  tibble(variable = variable, period = period_filter, n_site = n_site, n_bg = n_bg,
         median_site = round(med_site,1), median_bg = round(med_bg,1),
         p_value = round(test_result$p.value,4),
         direction = ifelse(med_site > med_bg, "higher_at_sites", "lower_at_sites"),
         note = NA_character_)
}

headline_vars <- c("dist_fault", "dist_dyke", "dist_lineament", "dist_shear")

message("\n=== LP vs survey-envelope background (n=728) ===")
lp_envelope <- bind_rows(lapply(headline_vars, run_mw_envelope, data = envelope_structural, period_filter = "LP"))
print(lp_envelope)

message("\n=== MP vs survey-envelope background (n=728) ===")
mp_envelope <- bind_rows(lapply(headline_vars, run_mw_envelope, data = envelope_structural, period_filter = "MP"))
print(mp_envelope)


# -------------------------------------------------------------
# SECTION 3: Three-Way Comparison — Uniform bg vs Envelope bg vs District-Stratified
# -------------------------------------------------------------

pooled_lp <- mw_table %>% filter(period == "LP", variable %in% headline_vars) %>%
  select(variable, uniform_bg_p = p_value, uniform_bg_direction = direction)
pooled_mp <- mw_table %>% filter(period == "MP", variable %in% headline_vars) %>%
  select(variable, uniform_bg_p = p_value, uniform_bg_direction = direction)

comparison_lp <- lp_envelope %>%
  select(variable, envelope_bg_p = p_value, envelope_bg_direction = direction) %>%
  left_join(pooled_lp, by = "variable") %>%
  mutate(direction_agrees = envelope_bg_direction == uniform_bg_direction,
         still_significant = envelope_bg_p < alpha,
         survives_envelope_check = direction_agrees & still_significant,
         period = "LP")

comparison_mp <- mp_envelope %>%
  select(variable, envelope_bg_p = p_value, envelope_bg_direction = direction) %>%
  left_join(pooled_mp, by = "variable") %>%
  mutate(direction_agrees = envelope_bg_direction == uniform_bg_direction,
         still_significant = envelope_bg_p < alpha,
         survives_envelope_check = direction_agrees & still_significant,
         period = "MP")

full_comparison <- bind_rows(comparison_lp, comparison_mp)

message("\n=== FULL COMPARISON: uniform-background vs envelope-background ===")
print(full_comparison)


# -------------------------------------------------------------
# SECTION 4: Focus Verdict — the Two Actual Survivors
# -------------------------------------------------------------

message("\n=== VERDICT: the two variables that survived every prior test ===")

lp_fault_check <- full_comparison %>% filter(period=="LP", variable=="dist_fault")
lp_lineament_check <- full_comparison %>% filter(period=="LP", variable=="dist_lineament")

message("\nLP dist_fault — survives envelope-background check: ",
        lp_fault_check$survives_envelope_check,
        " (p=", lp_fault_check$envelope_bg_p, ", direction=", lp_fault_check$envelope_bg_direction, ")")
message("LP dist_lineament — survives envelope-background check: ",
        lp_lineament_check$survives_envelope_check,
        " (p=", lp_lineament_check$envelope_bg_p, ", direction=", lp_lineament_check$envelope_bg_direction, ")")

both_survive <- isTRUE(lp_fault_check$survives_envelope_check) && isTRUE(lp_lineament_check$survives_envelope_check)

if (both_survive) {
  message("\nRESULT: BOTH surviving structural claims (LP fault, LP lineament) now hold")
  message("across FOUR independent tests: pooled MW, spatial block-bootstrap, district")
  message("stratification, AND survey-envelope background. This is the strongest possible")
  message("evidential basis this chapter can offer for a single finding. State this")
  message("explicitly in the manuscript as the chapter's headline, fully-triangulated result.")
} else {
  message("\nRESULT: at least one of the two remaining survivors does NOT hold against the")
  message("survey-envelope background. This means even the chapter's most robust apparent")
  message("finding may partly reflect where surveys were conducted rather than a genuine")
  message("environmental association. Report exactly which variable fails and why — this")
  message("would be a significant, honest downgrade of the chapter's central claim.")
}


# -------------------------------------------------------------
# SECTION 5: Save Outputs
# -------------------------------------------------------------

out_path <- file.path(chapter_root, "07_Statistics")

write_csv(envelope_structural, file.path(out_path, "envelope_structural_matrix.csv"))
write_csv(full_comparison, file.path(out_path, "envelope_vs_uniform_comparison.csv"))
write_csv(full_comparison, file.path(chapter_root, "11_Tables", "Table_survey_envelope_sensitivity.csv"))

message("\nSaved: envelope_structural_matrix.csv, envelope_vs_uniform_comparison.csv,")
message("Table_survey_envelope_sensitivity.csv [table-budget note: same 15-item cap")
message("issue as Scripts 10-11 — three new candidate tables now exist (district")
message("comparison, district-transfer RF, envelope sensitivity). Decide which of the")
message("current 4 manuscript tables to consolidate/replace before final formatting.]")

save(envelope_structural, lp_envelope, mp_envelope, full_comparison, both_survive,
     file = file.path(out_path, "envelope_sensitivity.RData"))
message("Saved: envelope_sensitivity.RData")


# -------------------------------------------------------------
# SECTION 6: Log Session
# -------------------------------------------------------------

log_file <- file.path(chapter_root, "13_Logs", paste0("script12_log_", Sys.Date(), ".txt"))
sink(log_file)
cat("SCRIPT 12 LOG — Survey-Envelope Sensitivity (Priority 3)\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Full comparison (uniform vs envelope background):\n")
print(full_comparison)
cat("\nLP fault survives:", lp_fault_check$survives_envelope_check, "\n")
cat("LP lineament survives:", lp_lineament_check$survives_envelope_check, "\n")
cat("Both survive:", both_survive, "\n")
sink()
message("Log saved: ", log_file)


# -------------------------------------------------------------
# SCRIPT 12 COMPLETE
# -------------------------------------------------------------

message("\n=============================================================")
message("Script 12 complete. Priority 3 (survey-envelope sensitivity) done.")
message("Both LP fault and LP lineament survive envelope-background check: ", both_survive)
message("\nWith Priorities 1-4 now complete (geochemistry, MULTI, district-transfer RF,")
message("survey-envelope), the manuscript's evidence hierarchy is substantially settled:")
message("LP fault + LP lineament = quadruple-triangulated (or fewer, per verdict above).")
message("Remaining Priorities 5-8 (block-size sensitivity, PCoA missingness+PERMANOVA,")
message("geochem buffer/zero-replacement sensitivity, manual classification check) are")
message("lower-stakes refinements, not central to the chapter's now-established headline")
message("finding. Recommend: proceed to final manuscript rewrite, listing Priorities 5-8")
message("honestly as future work rather than continuing to defer writing further.")
message("=============================================================")