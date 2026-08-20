# =============================================================
# SCRIPT 10 (NEW): District-Stratified Geochemistry + MULTI Analysis
# =============================================================
# Project: Reading the Palaeolithic Landscape
# Chapter: Mapping the Past (Springer Nature, 2026)
# Author: Sushant Begade | RTMNU Nagpur
# ORCID: 0009-0003-0804-1763
# Created: August 2026
# =============================================================
# PURPOSE: Priorities 1 and 2 of the outstanding-analyses list
# (manuscript Section 7.4). The structural-geology finding changed
# substantially under district stratification (Script 09) — there
# is no basis for assuming geochemistry or the MULTI-period result
# are independent of the same confound until actually tested.
#
# Priority 1: district-stratify the 6 major-oxide geochemistry PCs
#   (stream/horizon/regolith x PC1/PC2) for LP-in-Chandrapur and
#   MP-in-Nagpur, compared against the pooled result (Script 06).
# Priority 2: district-stratify the MULTI-period-locality result for
#   the variables flagged in the manuscript as inheriting the district
#   caveat (dist_fault, dist_lineament, dist_shear,
#   geochem_stream_major_PC1) plus the CHELSA variables.
#
# Reuses the district assignment already built and verified in
# Script 09 (0 mismatches for sites against source records) rather
# than rebuilding the spatial join — avoids risk of a second,
# possibly inconsistent, district assignment.
# =============================================================

load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")
load(file.path(chapter_root, "07_Statistics/mann_whitney_results.RData"))
load(file.path(chapter_root, "07_Statistics/district_stratified_analysis.RData"))

library(tidyverse)

set.seed(rf_seed)

message("Loaded full_matrix_district: ", nrow(full_matrix_district), " points")
message("District field already verified in Script 09 (0 site mismatches).")


# -------------------------------------------------------------
# SECTION 1: Reusable Stratified MW Test (site vs background, MULTI vs background)
# -------------------------------------------------------------

run_mw_stratified <- function(data, variable, group_type, group_filter, district_filter = NULL) {
  # group_type: "period" (LP/MP/UP) or "multi" (MULTI vs background)
  d <- data
  if (!is.null(district_filter)) d <- d %>% filter(district_label == district_filter)
  
  if (group_type == "period") {
    test_data <- d %>%
      filter(point_type == "background" | (point_type == "site" & period == group_filter)) %>%
      mutate(group = ifelse(point_type == "site", "site", "background"))
    label <- group_filter
  } else {
    test_data <- d %>%
      filter(point_type == "background" | (point_type == "site" & period == "MULTI")) %>%
      mutate(group = ifelse(point_type == "site", "site", "background"))
    label <- "MULTI"
  }
  
  test_data <- test_data %>% select(group, value = all_of(variable)) %>% filter(!is.na(value))
  n_site <- sum(test_data$group == "site")
  n_bg   <- sum(test_data$group == "background")
  
  district_note <- if (is.null(district_filter)) "pooled" else district_filter
  
  if (n_site < 5 || n_bg < 5) {
    return(tibble(variable = variable, group = label, district = district_note,
                  n_site = n_site, n_bg = n_bg, median_site = NA, median_bg = NA,
                  p_value = NA, direction = NA, note = "insufficient n (<5)"))
  }
  
  test_result <- tryCatch(
    wilcox.test(value ~ group, data = test_data, exact = FALSE, correct = TRUE),
    error = function(e) NULL
  )
  if (is.null(test_result)) {
    return(tibble(variable = variable, group = label, district = district_note,
                  n_site = n_site, n_bg = n_bg, median_site = NA, median_bg = NA,
                  p_value = NA, direction = NA, note = "test failed"))
  }
  
  med_site <- median(test_data$value[test_data$group == "site"], na.rm = TRUE)
  med_bg   <- median(test_data$value[test_data$group == "background"], na.rm = TRUE)
  
  tibble(variable = variable, group = label, district = district_note,
         n_site = n_site, n_bg = n_bg,
         median_site = round(med_site, 4), median_bg = round(med_bg, 4),
         p_value = round(test_result$p.value, 4),
         direction = ifelse(med_site > med_bg, "higher_at_sites", "lower_at_sites"),
         note = NA_character_)
}


# =============================================================
# PRIORITY 1: District-Stratified Geochemistry
# =============================================================

message("\n=== PRIORITY 1: District-stratified geochemistry ===")

geochem_vars <- c("geochem_stream_major_PC1", "geochem_stream_major_PC2",
                  "geochem_horizon_major_PC1", "geochem_horizon_major_PC2",
                  "geochem_regolith_major_PC1", "geochem_regolith_major_PC2")

message("\n--- LP within Chandrapur (n=45) ---")
geochem_lp_chandrapur <- bind_rows(lapply(geochem_vars, run_mw_stratified,
                                          data = full_matrix_district, group_type = "period",
                                          group_filter = "LP", district_filter = "Chandrapur"))
print(geochem_lp_chandrapur, width = Inf)

message("\n--- MP within Nagpur (n=55) ---")
geochem_mp_nagpur <- bind_rows(lapply(geochem_vars, run_mw_stratified,
                                      data = full_matrix_district, group_type = "period",
                                      group_filter = "MP", district_filter = "Nagpur"))
print(geochem_mp_nagpur, width = Inf)

# Pooled comparison, pulled from Script 06's already-computed mw_table for consistency
pooled_geochem <- mw_table %>% filter(variable %in% geochem_vars, period %in% c("LP","MP")) %>%
  select(variable, period, pooled_p = p_value, pooled_direction = direction)

geochem_comparison <- bind_rows(
  geochem_lp_chandrapur %>% mutate(period = "LP") %>%
    select(variable, period, stratified_p = p_value, stratified_direction = direction, n_site),
  geochem_mp_nagpur %>% mutate(period = "MP") %>%
    select(variable, period, stratified_p = p_value, stratified_direction = direction, n_site)
) %>%
  left_join(pooled_geochem, by = c("variable","period")) %>%
  mutate(
    direction_agrees = pooled_direction == stratified_direction,
    still_significant = stratified_p < alpha,
    survives = direction_agrees & still_significant
  )

message("\n=== GEOCHEMISTRY: pooled vs district-stratified ===")
print(geochem_comparison, width = Inf)

n_geochem_survive <- sum(geochem_comparison$survives, na.rm = TRUE)
n_geochem_total   <- sum(!is.na(geochem_comparison$survives))
message("\nGeochemistry variables surviving district stratification (direction+sig): ",
        n_geochem_survive, " of ", n_geochem_total, " tested")

if (n_geochem_survive == n_geochem_total) {
  message("RESULT: geochemistry finding SURVIVES district stratification in full —")
  message("this is independent evidence, not confounded the way structural geology was.")
} else if (n_geochem_survive == 0) {
  message("RESULT: geochemistry finding DOES NOT SURVIVE — same district-confound")
  message("pattern as structural geology. Section 7.2/Abstract must be rewritten to")
  message("state geochemistry differentiation is ALSO substantially district-driven.")
} else {
  message("RESULT: PARTIAL survival — report per-variable, do not generalise 'geochemistry'")
  message("as a whole the same way the structural finding could not be generalised.")
}


# =============================================================
# PRIORITY 2: District-Stratified MULTI-Period Locality Check
# =============================================================

message("\n\n=== PRIORITY 2: District-stratified MULTI-period localities ===")
message("MULTI sites (n=60) are NOT concentrated as heavily by district as LP/MP —")
message("checking distribution before proceeding:")

multi_district_dist <- full_matrix_district %>% filter(point_type == "site", period == "MULTI") %>%
  count(district_label)
print(multi_district_dist)

multi_vars_flagged <- c("dist_fault", "dist_lineament", "dist_shear",
                        "geochem_stream_major_PC1",
                        "chelsa_bio15_lgm", "chelsa_bio15_modern",
                        "chelsa_bio12_lgm", "chelsa_bio12_modern", "chelsa_bio01_lgm")

# Only proceed with per-district stratification for districts where MULTI n >= 15
# (below this, per Section 5.1's own n<5 gate philosophy, results are not
# reliably interpretable and are reported as insufficient rather than computed
# and potentially over-interpreted)
multi_district_counts <- multi_district_dist %>% filter(n >= 15) %>% pull(district_label)

if (length(multi_district_counts) == 0) {
  message("\nNo district has >=15 MULTI sites — district-stratified MULTI testing")
  message("is not reliably interpretable with this sample. Reporting pooled MULTI")
  message("result as-is with the existing caveat language; do not attempt to")
  message("resolve this further without more MULTI-period data.")
  multi_comparison <- tibble(variable = character(), note = "insufficient district n for MULTI stratification")
} else {
  message("\nProceeding with district-stratified MULTI testing for: ",
          paste(multi_district_counts, collapse = ", "))
  
  multi_strat_results <- list()
  for (d in multi_district_counts) {
    for (v in multi_vars_flagged) {
      multi_strat_results[[paste(d, v)]] <- run_mw_stratified(
        full_matrix_district, v, group_type = "multi", group_filter = "MULTI", district_filter = d
      )
    }
  }
  multi_strat_table <- bind_rows(multi_strat_results)
  
  pooled_multi <- multi_table %>% filter(variable %in% multi_vars_flagged) %>%
    select(variable, pooled_p = p_value, pooled_direction = direction)
  
  multi_comparison <- multi_strat_table %>%
    left_join(pooled_multi, by = "variable") %>%
    mutate(direction_agrees = pooled_direction == direction,
           still_significant = p_value < alpha,
           survives = direction_agrees & still_significant)
  
  message("\n=== MULTI-period: pooled vs district-stratified ===")
  print(multi_comparison, width = Inf)
}


# -------------------------------------------------------------
# SECTION: Save Outputs
# -------------------------------------------------------------

out_path <- file.path(chapter_root, "07_Statistics")

write_csv(geochem_comparison, file.path(out_path, "district_stratified_geochemistry.csv"))
write_csv(geochem_comparison, file.path(chapter_root, "11_Tables", "Table_district_stratified_geochemistry.csv"))
write_csv(multi_comparison, file.path(out_path, "district_stratified_multi.csv"))

message("\nSaved: district_stratified_geochemistry.csv, Table_district_stratified_geochemistry.csv,")
message("district_stratified_multi.csv")
message("NOTE ON FIGURE/TABLE BUDGET: adding this as a 5th table breaks the 15-item cap")
message("(11 figures + 4 tables already = 15). If this result is reported in the main")
message("text, fold it into Table 3's existing structure or replace a less essential")
message("table rather than adding a 5th — decide once you see whether geochem survives.")

save(geochem_comparison, multi_comparison, n_geochem_survive, n_geochem_total,
     file = file.path(out_path, "district_stratified_geochem_multi.RData"))
message("Saved: district_stratified_geochem_multi.RData")


# -------------------------------------------------------------
# SECTION: Log Session
# -------------------------------------------------------------

log_file <- file.path(chapter_root, "13_Logs", paste0("script10_log_", Sys.Date(), ".txt"))
sink(log_file)
cat("SCRIPT 10 LOG — District-Stratified Geochemistry + MULTI (Priorities 1-2)\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Geochemistry: pooled vs district-stratified\n"); print(geochem_comparison)
cat("\nSurvival:", n_geochem_survive, "of", n_geochem_total, "\n\n")
cat("MULTI: pooled vs district-stratified\n"); print(multi_comparison)
sink()
message("Log saved: ", log_file)


# -------------------------------------------------------------
# SCRIPT 10 COMPLETE
# -------------------------------------------------------------

message("\n=============================================================")
message("Script 10 complete.")
message("Priority 1 (geochemistry): ", n_geochem_survive, "/", n_geochem_total, " survive district stratification")
message("Priority 2 (MULTI): see multi_comparison above / district_stratified_multi.csv")
message("Remaining priorities per reviewer: 3=survey-envelope sensitivity,")
message("4=district-transfer RF, 5=block-size sensitivity, 6=PCoA missingness+PERMANOVA,")
message("7=geochem buffer-scale+zero-replacement sensitivity, 8=manual classification check.")
message("Per reviewer instruction: do NOT rewrite manuscript prose again until the full")
message("priority list is done. Next recommended: Priority 3 (survey-envelope) or")
message("Priority 4 (district-transfer RF) — both directly test whether the surviving")
message("LP fault/lineament result is genuinely robust or itself an artefact.")
message("=============================================================")