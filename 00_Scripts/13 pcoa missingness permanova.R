# =============================================================
# SCRIPT 13 (NEW): PCoA Missingness Diagnostic + PERMANOVA/PERMDISP
# =============================================================
# Project: Reading the Palaeolithic Landscape
# Chapter: Mapping the Past (Springer Nature, 2026)
# Author: Sushant Begade | RTMNU Nagpur
# ORCID: 0009-0003-0804-1763
# Created: August 2026
# =============================================================
# PURPOSE:
#   Task A (mandatory): is the 716-point PCoA complete-case subset
#     geographically/environmentally representative of the full
#     1,197-point sample, or geographically selective (given uneven
#     geochemical coverage by layer: 98.1%/59.6%/46.1%)?
#   Task B (strongly recommended): formal PERMANOVA test of whether
#     LP/MP/UP differ as multivariate distributions in Gower space
#     (not just at centroids), plus PERMDISP to separate genuine
#     centroid separation from differences in within-group dispersion.
#
# ALSO FLAGS A FOUND ERROR: the manuscript's PCoA section reports
# centroid n as the full period counts (LP=48, MP=62, UP=27), but
# Script 07's period_centroids table was computed on the 716-point
# complete-case subset, where actual per-period n was LP=30, MP=38,
# UP=17 (printed in the original Script 07 run). This script confirms
# those n's directly from pcoa_scores so the manuscript can be
# corrected with a verified number, not a recollection.
# =============================================================

load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")
load(file.path(chapter_root, "09_RandomForest/pcoa_rf_results.RData"))
load(file.path(chapter_root, "07_Statistics/district_stratified_analysis.RData"))

library(tidyverse)
library(vegan)

set.seed(rf_seed)

message("Loaded pcoa_scores: ", nrow(pcoa_scores), " rows (the complete-case PCoA subset)")
message("Loaded gower_dist: ", attr(gower_dist,"Size"), " x ", attr(gower_dist,"Size"))
message("Loaded full_matrix_district: ", nrow(full_matrix_district), " rows (full 1,197)")


# -------------------------------------------------------------
# SECTION 0: CONFIRM the centroid-n manuscript error
# -------------------------------------------------------------

message("\n=== CONFIRMING PCoA CENTROID n (manuscript correction needed) ===")
actual_centroid_n <- pcoa_scores %>% filter(point_type=="site", period %in% c("LP","MP","UP")) %>%
  count(period)
message("Actual per-period n used in Script 07 centroid computation:")
print(actual_centroid_n)
message("Manuscript currently states n=48 (LP), n=62 (MP), n=27 (UP) next to the")
message("centroid coordinates — these are the FULL period counts, not the complete-case")
message("counts actually used. CORRECT the manuscript's PCoA section to use the n's above.")


# =============================================================
# TASK A: Missingness / Representativeness Diagnostic
# =============================================================

message("\n=== TASK A: Missingness diagnostic (716 retained vs 481 excluded) ===")

full_matrix_district <- full_matrix_district %>%
  mutate(retained_in_pcoa = point_id %in% pcoa_scores$point_id)

n_retained <- sum(full_matrix_district$retained_in_pcoa)
n_excluded <- sum(!full_matrix_district$retained_in_pcoa)
message("Retained: ", n_retained, " | Excluded: ", n_excluded,
        " (should sum to 1197: ", n_retained + n_excluded, ")")

# --- Categorical comparisons: district, period, point_type ---

message("\n--- District: retained vs excluded ---")
district_tab <- table(full_matrix_district$district_label, full_matrix_district$retained_in_pcoa)
print(district_tab)
district_chisq <- chisq.test(district_tab)
message("Chi-square test (district x retained): p = ", round(district_chisq$p.value, 4))

message("\n--- Period: retained vs excluded ---")
period_tab <- table(full_matrix_district$period, full_matrix_district$retained_in_pcoa)
print(period_tab)
period_chisq <- chisq.test(period_tab)
message("Chi-square test (period x retained): p = ", round(period_chisq$p.value, 4))

message("\n--- Point type (site/background): retained vs excluded ---")
type_tab <- table(full_matrix_district$point_type, full_matrix_district$retained_in_pcoa)
print(type_tab)
type_chisq <- chisq.test(type_tab)
message("Chi-square test (point_type x retained): p = ", round(type_chisq$p.value, 4))

# --- Continuous comparisons: elevation, key structural distances ---

message("\n--- Continuous variables: retained vs excluded (Mann-Whitney) ---")
continuous_check_vars <- c("elevation", "dist_fault", "dist_lineament", "dist_shear",
                           "chelsa_bio12_modern")

continuous_diagnostic <- map_dfr(continuous_check_vars, function(v) {
  d <- full_matrix_district %>% select(retained_in_pcoa, value = all_of(v)) %>% filter(!is.na(value))
  test <- wilcox.test(value ~ retained_in_pcoa, data = d, exact = FALSE)
  tibble(variable = v,
         median_retained = round(median(d$value[d$retained_in_pcoa], na.rm=TRUE), 2),
         median_excluded = round(median(d$value[!d$retained_in_pcoa], na.rm=TRUE), 2),
         p_value = round(test$p.value, 4))
})
print(continuous_diagnostic)

n_biased_continuous <- sum(continuous_diagnostic$p_value < alpha)
n_biased_categorical <- sum(c(district_chisq$p.value, period_chisq$p.value, type_chisq$p.value) < alpha)

message("\n=== MISSINGNESS VERDICT ===")
message("Categorical variables with significant retained/excluded difference: ",
        n_biased_categorical, " of 3 (district, period, point_type)")
message("Continuous variables with significant retained/excluded difference: ",
        n_biased_continuous, " of ", length(continuous_check_vars))

if (n_biased_categorical == 0 && n_biased_continuous == 0) {
  message("\nRESULT: no significant retained/excluded differences detected. The 716-point")
  message("PCoA subset appears reasonably representative of the full 1,197-point sample")
  message("on the dimensions tested. This strengthens confidence in Section 6.3's PCoA")
  message("results — state this explicitly rather than leaving representativeness")
  message("'untested' as the manuscript currently does.")
} else {
  message("\nRESULT: significant retained/excluded differences found on ", n_biased_categorical +
            n_biased_continuous, " of ", 3 + length(continuous_check_vars), " dimensions tested.")
  message("The PCoA subset IS geographically/environmentally selective. Manuscript Section")
  message("6.3/7.3 must state this explicitly and downgrade confidence in the PCoA result")
  message("accordingly — specify exactly which dimensions are biased using the tables above.")
}


# =============================================================
# TASK B: PERMANOVA + PERMDISP
# =============================================================

message("\n\n=== TASK B: PERMANOVA + PERMDISP (site-only, LP/MP/UP) ===")

site_idx <- which(pcoa_scores$point_type == "site" & pcoa_scores$period %in% c("LP","MP","UP"))
message("Site subset for PERMANOVA/PERMDISP: n=", length(site_idx),
        " (LP=", sum(pcoa_scores$period[site_idx]=="LP"),
        ", MP=", sum(pcoa_scores$period[site_idx]=="MP"),
        ", UP=", sum(pcoa_scores$period[site_idx]=="UP"), ")")

gower_mat_full <- as.matrix(gower_dist)
gower_site_dist <- as.dist(gower_mat_full[site_idx, site_idx])
period_factor <- factor(pcoa_scores$period[site_idx], levels = c("LP","MP","UP"))

message("\nRunning PERMANOVA (adonis2, 999 permutations)...")
permanova_result <- adonis2(gower_site_dist ~ period_factor, permutations = 999)
print(permanova_result)

message("\nRunning PERMDISP (betadisper)...")
permdisp_model <- betadisper(gower_site_dist, period_factor)
permdisp_test <- permutest(permdisp_model, permutations = 999)
print(permdisp_test)

message("\nPERMDISP pairwise (Tukey HSD on distances to centroid):")
permdisp_tukey <- TukeyHSD(permdisp_model)
print(permdisp_tukey)

message("\n=== PERMANOVA/PERMDISP VERDICT ===")
permanova_p <- permanova_result$`Pr(>F)`[1]
permdisp_p  <- permdisp_test$tab$`Pr(>F)`[1]

message("PERMANOVA p = ", round(permanova_p, 4), " (tests whether period groups differ",
        " as multivariate distributions — centroid AND/OR dispersion)")
message("PERMDISP p = ", round(permdisp_p, 4), " (tests whether within-group dispersion",
        " differs by period — if significant, PERMANOVA result may partly reflect",
        " dispersion differences rather than pure centroid separation)")

if (permanova_p < alpha && permdisp_p >= alpha) {
  message("\nRESULT: periods differ significantly as multivariate distributions (PERMANOVA),")
  message("and this is NOT attributable to differing within-group dispersion (PERMDISP n.s.).")
  message("This is the cleanest possible outcome — genuine centroid separation. Report both")
  message("results together in the manuscript as formal confirmation of the descriptive")
  message("centroid displacements already reported (Section 6.3).")
} else if (permanova_p < alpha && permdisp_p < alpha) {
  message("\nRESULT: periods differ significantly (PERMANOVA), BUT dispersion also differs")
  message("significantly (PERMDISP). The PERMANOVA result may partly reflect unequal spread")
  message("rather than purely centroid separation — this caveat MUST be stated explicitly")
  message("wherever the PCoA centroid displacement is discussed in the manuscript.")
} else {
  message("\nRESULT: PERMANOVA does not reach significance. The descriptive centroid")
  message("displacements reported in Section 6.3 are NOT confirmed by a formal multivariate")
  message("test at this sample size. This should be stated plainly — the manuscript should")
  message("downgrade PCoA centroid language from 'differentiated' to 'descriptively distinct,")
  message("not confirmed by formal multivariate testing'.")
}


# -------------------------------------------------------------
# SECTION: Save Outputs
# -------------------------------------------------------------

out_path <- file.path(chapter_root, "07_Statistics")

write_csv(continuous_diagnostic, file.path(out_path, "pcoa_missingness_diagnostic.csv"))
write_csv(actual_centroid_n, file.path(out_path, "pcoa_centroid_n_corrected.csv"))

permanova_summary <- tibble(
  test = c("PERMANOVA (period)", "PERMDISP (period)"),
  statistic = c(permanova_result$F[1], permdisp_test$tab$F[1]),
  p_value = c(permanova_p, permdisp_p)
)
write_csv(permanova_summary, file.path(out_path, "permanova_permdisp_results.csv"))

message("\nSaved: pcoa_missingness_diagnostic.csv, pcoa_centroid_n_corrected.csv,")
message("permanova_permdisp_results.csv")
message("TABLE BUDGET NOTE: these are diagnostic/confirmatory results, best reported as")
message("in-text statistics (a sentence with F/p values) rather than new tables — the")
message("15-item cap is already fully allocated across the existing 4 tables + 11 figures.")

save(continuous_diagnostic, district_chisq, period_chisq, type_chisq,
     actual_centroid_n, permanova_result, permdisp_model, permdisp_test, permdisp_tukey,
     file = file.path(out_path, "pcoa_diagnostics_permanova.RData"))
message("Saved: pcoa_diagnostics_permanova.RData")


# -------------------------------------------------------------
# SECTION: Log Session
# -------------------------------------------------------------

log_file <- file.path(chapter_root, "13_Logs", paste0("script13_log_", Sys.Date(), ".txt"))
sink(log_file)
cat("SCRIPT 13 LOG — PCoA Missingness Diagnostic + PERMANOVA/PERMDISP\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("CORRECTED centroid n:\n"); print(actual_centroid_n)
cat("\nMissingness diagnostic:\n"); print(continuous_diagnostic)
cat("District chi-sq p:", round(district_chisq$p.value,4), "\n")
cat("Period chi-sq p:", round(period_chisq$p.value,4), "\n")
cat("Point-type chi-sq p:", round(type_chisq$p.value,4), "\n\n")
cat("PERMANOVA p:", round(permanova_p,4), "| PERMDISP p:", round(permdisp_p,4), "\n")
sink()
message("Log saved: ", log_file)


# -------------------------------------------------------------
# SCRIPT 13 COMPLETE
# -------------------------------------------------------------

message("\n=============================================================")
message("Script 13 complete.")
message("Task A (missingness): ", n_biased_categorical + n_biased_continuous,
        " of ", 3 + length(continuous_check_vars), " dimensions show significant bias")
message("Task B (PERMANOVA): p=", round(permanova_p,4), " | PERMDISP: p=", round(permdisp_p,4))
message("Manuscript correction required regardless of verdict: PCoA centroid n's must")
message("change from 48/62/27 to the actual complete-case n's confirmed above.")
message("Remaining: geochemical buffer-scale + zero-replacement sensitivity (task 4),")
message("manual classification review (task 1 — human review aid to follow separately).")
message("=============================================================")