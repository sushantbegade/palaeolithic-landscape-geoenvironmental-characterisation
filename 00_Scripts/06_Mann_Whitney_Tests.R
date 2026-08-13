# =============================================================
# SCRIPT 06 (REBUILT): Mann-Whitney U Tests + Significance Matrix
# =============================================================
# Project: Reading the Palaeolithic Landscape
# Chapter: Mapping the Past (Springer Nature, 2026)
# Author: Sushant Begade | RTMNU Nagpur
# ORCID: 0009-0003-0804-1763
# Rebuild date: August 2026
# =============================================================
# CHANGES FROM ORIGINAL SCRIPT 06:
#   1. analysis_vars REBUILT: burial_depth_m REMOVED (Script 05 found
#      the borehole data is GSI mineral-exploration DDH data with
#      <8% site coverage at any defensible interpolation distance —
#      not usable as an analysis variable). 6 geochem major-oxide
#      CLR-PCA variables ADDED (Script 04) — geochemistry now
#      genuinely enters the statistical results for the first time.
#      Net: 25 -> 30 analysis variables.
#   2. CHELSA-LGM variables RETAINED but flagged explicitly at every
#      relevant point as a spatially-structured covariate, NOT a
#      period-matched palaeoclimate reconstruction (LP/MP predate the
#      ~21ka TraCE21k window by >100ka). Script 07 will run RF with/
#      without these variables as an explicit sensitivity check.
#   3. NEW — Persistent Places analysis: MULTI-period sites (n=60,
#      Script 01) compared against background on the same variable
#      battery, reported separately from the LP/MP/UP period-
#      stratified tests. Turns the 60 excluded multi-period sites
#      into a positive finding (palimpsest/reoccupation locations)
#      instead of a silent data loss.
#   4. NEW — Spatial block bootstrap for the headline structural-
#      geology result (dist_fault/dyke/lineament/shear, LP & MP,
#      the chapter's strongest claim per the original draft). Standard
#      MW+Bonferroni does not correct for spatial pseudoreplication
#      (reviewer critical problem #2) — this provides a spatially-
#      aware robustness check on the specific claims the Discussion
#      leans on hardest, without rebuilding the entire inference
#      framework as full spatial logistic/GAM modelling (deferred to
#      Script 07 where RF spatial-block CV is implemented instead).
# =============================================================

load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")
load(file.path(chapter_root, "03_Extracted_Values/master_extraction_matrix.RData"))
load(file.path(chapter_root, "05_Structural_Geology/vector_extraction_matrix.RData"))
load(file.path(chapter_root, "06_Taphonomy/taphonomy_data.RData"))

library(tidyverse)
library(rstatix)
library(coin)
library(pheatmap)
library(ggplot2)
library(patchwork)

set.seed(rf_seed)

message("All data loaded.")
message("Master matrix: ", nrow(master_matrix), " rows")
message("Vector matrix (structural + geochem): ", nrow(vector_matrix), " rows")


# -------------------------------------------------------------
# SECTION 2: Assemble Full Analysis Matrix
# -------------------------------------------------------------

message("\nAssembling full analysis matrix...")

core_vars <- master_matrix %>%
  select(
    point_id, point_type, period, coord_precision, bg_source,
    elevation,
    NDVI, NDWI, MNDWI, NDBI, BSI, SAVI, MSAVI,
    soil_depth, soil_erosion, soil_productivity, soil_slope, soil_texture,
    chelsa_bio01_lgm, chelsa_bio12_lgm, chelsa_bio15_lgm,
    chelsa_bio01_modern, chelsa_bio12_modern, chelsa_bio15_modern
  )

struct_vars <- vector_matrix %>%
  select(point_id, dist_fault, dist_dyke, dist_lineament, dist_shear, dist_mineral)

geochem_vars <- vector_matrix %>%
  select(point_id, starts_with("geochem_"))

# burial_depth_m DELIBERATELY EXCLUDED — see Script 05 findings.
# Retained here ONLY as a descriptive join for reporting coverage/
# limitations in Methods/Discussion, NOT as an analysis_var.
burial_descriptive <- taphonomy_matrix %>%
  select(point_id, burial_depth_m, burial_risk_label)

full_matrix <- core_vars %>%
  left_join(struct_vars, by = "point_id") %>%
  left_join(geochem_vars, by = "point_id") %>%
  left_join(burial_descriptive, by = "point_id")

message("Full analysis matrix: ", nrow(full_matrix), " rows x ", ncol(full_matrix), " cols")

analysis_vars <- setdiff(
  names(full_matrix),
  c("point_id", "point_type", "period", "coord_precision", "bg_source",
    "burial_depth_m", "burial_risk_label")
)

message("\nAnalysis variables (n=", length(analysis_vars), "):")
message(paste(analysis_vars, collapse = ", "))
message("\nburial_depth_m / burial_risk_label retained for DESCRIPTIVE reporting")
message("only (Section 6.4 rewrite) — NOT included in analysis_vars.")
message("\nCHELSA-LGM vars (chelsa_*_lgm) are spatially-structured covariates,")
message("NOT period-matched palaeoclimate for LP/MP. Script 07 runs RF with/")
message("without these as an explicit sensitivity check — flag now, resolve there.")


# -------------------------------------------------------------
# SECTION 3: Mann-Whitney U Helper (unchanged logic)
# -------------------------------------------------------------

run_mw_test <- function(data, variable, group_filter_type, period_filter = NULL) {
  # group_filter_type: "period" (site period == period_filter vs background)
  #                     "multi" (site period == MULTI vs background)
  
  if (group_filter_type == "period") {
    test_data <- data %>%
      filter(point_type == "background" |
               (point_type == "site" & period == period_filter)) %>%
      mutate(group = ifelse(point_type == "site", "site", "background")) %>%
      select(group, value = all_of(variable)) %>%
      filter(!is.na(value))
    label <- period_filter
  } else {
    test_data <- data %>%
      filter(point_type == "background" |
               (point_type == "site" & period == "MULTI")) %>%
      mutate(group = ifelse(point_type == "site", "multi", "background")) %>%
      select(group, value = all_of(variable)) %>%
      filter(!is.na(value))
    label <- "MULTI"
  }
  
  n_site <- sum(test_data$group %in% c("site", "multi"))
  n_bg   <- sum(test_data$group == "background")
  
  if (n_site < 3 || n_bg < 3) {
    return(data.frame(variable = variable, period = label, p_value = NA, w_stat = NA,
                      median_site = NA, median_bg = NA, direction = NA,
                      significant = NA, n_site = n_site, n_bg = n_bg))
  }
  
  test_result <- tryCatch(
    wilcox.test(value ~ group, data = test_data, exact = FALSE,
                correct = TRUE, alternative = "two.sided"),
    error = function(e) NULL
  )
  if (is.null(test_result)) {
    return(data.frame(variable = variable, period = label, p_value = NA, w_stat = NA,
                      median_site = NA, median_bg = NA, direction = NA,
                      significant = NA, n_site = n_site, n_bg = n_bg))
  }
  
  med_site <- median(test_data$value[test_data$group %in% c("site","multi")], na.rm = TRUE)
  med_bg   <- median(test_data$value[test_data$group == "background"], na.rm = TRUE)
  direction <- ifelse(med_site > med_bg, "higher_at_sites", "lower_at_sites")
  
  data.frame(
    variable = variable, period = label,
    p_value = round(test_result$p.value, 4), w_stat = round(test_result$statistic, 0),
    median_site = round(med_site, 4), median_bg = round(med_bg, 4),
    direction = direction, significant = test_result$p.value < alpha,
    n_site = n_site, n_bg = n_bg
  )
}


# -------------------------------------------------------------
# SECTION 4: Run MW — All Variables x LP/MP/UP  (unchanged approach)
# -------------------------------------------------------------

message("\nRunning Mann-Whitney U tests: ", length(analysis_vars), " vars x ",
        length(periods), " periods = ", length(analysis_vars) * length(periods), " tests")

mw_results <- list()
for (per in periods) {
  for (var in analysis_vars) {
    mw_results[[paste(per, var, sep = "_")]] <- run_mw_test(full_matrix, var, "period", per)
  }
  message("  Done: ", per)
}
mw_table <- bind_rows(mw_results)

mw_table <- mw_table %>%
  mutate(p_bonferroni = p.adjust(p_value, method = "bonferroni"),
         sig_bonferroni = p_bonferroni < alpha)

message("\nTotal tests: ", nrow(mw_table))
message("Significant (p<", alpha, "): ", sum(mw_table$significant, na.rm = TRUE))
message("Significant (Bonferroni): ", sum(mw_table$sig_bonferroni, na.rm = TRUE))


# -------------------------------------------------------------
# SECTION 5: Persistent Places — MULTI vs Background  [NEW]
# -------------------------------------------------------------

message("\n=== PERSISTENT PLACES: MULTI-period sites (n=60) vs background ===")
message("These are sites attributed to more than one Palaeolithic period in")
message("source literature — repeat/persistent-use locations across the")
message("cultural sequence. Tested here on the same variable battery,")
message("reported separately from LP/MP/UP period-stratified results.")

multi_results <- list()
for (var in analysis_vars) {
  multi_results[[var]] <- run_mw_test(full_matrix, var, "multi")
}
multi_table <- bind_rows(multi_results) %>%
  mutate(p_bonferroni = p.adjust(p_value, method = "bonferroni"),
         sig_bonferroni = p_bonferroni < alpha)

message("\nPersistent Places — significant variables (p<", alpha, "):")
print(multi_table %>% filter(significant == TRUE) %>% arrange(p_value) %>%
        select(variable, p_value, direction, median_site, median_bg, n_site))

write_csv(multi_table, file.path(chapter_root, "07_Statistics", "persistent_places_MULTI_vs_background.csv"))
message("Saved: persistent_places_MULTI_vs_background.csv [NEW — potential")
message("Results subsection: 'Persistent Places: multi-period reoccupation']")


# -------------------------------------------------------------
# SECTION 6: Significance Matrices (unchanged logic)
# -------------------------------------------------------------

sig_matrix_p <- mw_table %>% select(variable, period, p_value) %>%
  pivot_wider(names_from = period, values_from = p_value) %>% column_to_rownames("variable")
sig_matrix_bin <- mw_table %>% select(variable, period, significant) %>%
  pivot_wider(names_from = period, values_from = significant) %>% column_to_rownames("variable")
dir_matrix <- mw_table %>% select(variable, period, direction) %>%
  pivot_wider(names_from = period, values_from = direction) %>% column_to_rownames("variable")

message("\n=== SIGNIFICANT VARIABLES PER PERIOD ===")
for (per in periods) {
  sig_vars <- mw_table %>% filter(period == per, significant == TRUE) %>%
    arrange(p_value) %>% select(variable, p_value, direction, median_site, median_bg)
  message("\n", per, " — ", nrow(sig_vars), " significant (p<", alpha, "):")
  print(sig_vars)
}


# -------------------------------------------------------------
# SECTION 7: Effect Sizes — Rank-Biserial Correlation (unchanged)
# -------------------------------------------------------------

compute_rbc <- function(data, variable, period_filter) {
  test_data <- data %>%
    filter(point_type == "background" |
             (point_type == "site" & period == period_filter)) %>%
    mutate(group = ifelse(point_type == "site", 1, 0)) %>%
    select(group, value = all_of(variable)) %>% filter(!is.na(value))
  n_site <- sum(test_data$group == 1); n_bg <- sum(test_data$group == 0)
  if (n_site < 3) return(NA)
  w <- wilcox.test(value ~ group, data = test_data, exact = FALSE)$statistic
  round(as.numeric(1 - (2 * w) / (n_site * n_bg)), 3)
}

effect_sizes <- expand.grid(variable = analysis_vars, period = periods, stringsAsFactors = FALSE) %>%
  rowwise() %>% mutate(rbc = compute_rbc(full_matrix, variable, period)) %>% ungroup()
mw_table <- mw_table %>% left_join(effect_sizes, by = c("variable", "period"))


# -------------------------------------------------------------
# SECTION 8: Period Summary (unchanged logic)
# -------------------------------------------------------------

period_summary <- mw_table %>% filter(!is.na(significant)) %>% group_by(period) %>%
  summarise(n_vars_tested = n(), n_significant = sum(significant, na.rm = TRUE),
            n_sig_bonf = sum(sig_bonferroni, na.rm = TRUE),
            pct_significant = round(100 * mean(significant, na.rm = TRUE), 1),
            top_variable = variable[which.min(p_value)],
            lowest_p = round(min(p_value, na.rm = TRUE), 4), .groups = "drop")
print(period_summary)


# =============================================================
# SECTION 9: Spatial Block Bootstrap — Headline Structural Claim [NEW]
# =============================================================
# Standard MW + Bonferroni does not correct for spatial pseudo-
# replication: neither sites nor background points are spatially
# independent (structural geology, soils, climate are all spatially
# autocorrelated). This provides a spatially-aware robustness check
# specifically for the structural-geology result the Discussion
# leans on hardest (dist_fault/dyke/lineament/shear, LP & MP).
# Block bootstrap: resample spatial BLOCKS with replacement (not
# individual points), preserving within-block spatial correlation,
# recompute the median site-background difference each iteration,
# build an empirical 95% CI. If the CI excludes zero, the effect
# survives accounting for spatial structure.

message("\n=== SPATIAL BLOCK BOOTSTRAP: structural geology, LP & MP ===")

spatial_block_size_m <- 10000  # 10km grid blocks — adjust if needed

full_matrix_coords <- full_matrix %>%
  left_join(master_matrix %>% select(point_id, easting, northing), by = "point_id") %>%
  mutate(
    block_x = floor(easting  / spatial_block_size_m),
    block_y = floor(northing / spatial_block_size_m),
    block_id = paste0(block_x, "_", block_y)
  )

n_blocks <- length(unique(full_matrix_coords$block_id))
message("Study area divided into ", n_blocks, " spatial blocks (", spatial_block_size_m/1000, "km grid).")

block_bootstrap_test <- function(data, variable, period_filter, n_boot = 1000) {
  test_data <- data %>%
    filter(point_type == "background" |
             (point_type == "site" & period == period_filter)) %>%
    mutate(group = ifelse(point_type == "site", "site", "background")) %>%
    select(group, value = all_of(variable), block_id) %>%
    filter(!is.na(value))
  
  obs_diff <- median(test_data$value[test_data$group == "site"], na.rm = TRUE) -
    median(test_data$value[test_data$group == "background"], na.rm = TRUE)
  
  blocks <- unique(test_data$block_id)
  n_b <- length(blocks)
  boot_diffs <- numeric(n_boot)
  
  for (i in seq_len(n_boot)) {
    sampled_blocks <- sample(blocks, size = n_b, replace = TRUE)
    boot_data <- bind_rows(lapply(sampled_blocks, function(b) test_data %>% filter(block_id == b)))
    site_vals <- boot_data$value[boot_data$group == "site"]
    bg_vals   <- boot_data$value[boot_data$group == "background"]
    if (length(site_vals) < 3 || length(bg_vals) < 3) { boot_diffs[i] <- NA; next }
    boot_diffs[i] <- median(site_vals, na.rm = TRUE) - median(bg_vals, na.rm = TRUE)
  }
  
  boot_diffs <- boot_diffs[!is.na(boot_diffs)]
  ci <- quantile(boot_diffs, c(0.025, 0.975), na.rm = TRUE)
  
  tibble(variable = variable, period = period_filter, n_blocks = n_b,
         observed_diff = round(obs_diff, 1),
         boot_ci_low = round(ci[1], 1), boot_ci_high = round(ci[2], 1),
         ci_excludes_zero = (ci[1] > 0 & ci[2] > 0) | (ci[1] < 0 & ci[2] < 0),
         n_valid_boot = length(boot_diffs))
}

headline_vars <- c("dist_fault", "dist_dyke", "dist_lineament", "dist_shear")
headline_periods <- c("LP", "MP")

message("\nRunning ", length(headline_vars) * length(headline_periods),
        " block-bootstrap tests (1000 iterations each — this takes a moment)...")

block_boot_results <- list()
for (per in headline_periods) {
  for (var in headline_vars) {
    message("  ", per, " - ", var, "...")
    block_boot_results[[paste(per, var)]] <- block_bootstrap_test(full_matrix_coords, var, per)
  }
}
block_boot_table <- bind_rows(block_boot_results)

message("\n=== SPATIAL BLOCK BOOTSTRAP RESULTS ===")
print(block_boot_table)

n_survive <- sum(block_boot_table$ci_excludes_zero)
message("\n", n_survive, " of ", nrow(block_boot_table),
        " headline structural results have a bootstrap CI excluding zero")
message("(i.e. survive spatial-block resampling). Compare against the",
        " original MW/Bonferroni significance — any result significant in")
message("MW but NOT surviving here should be downgraded in Discussion language")
message("from 'robust' to 'nominally significant, sensitive to spatial structure'.")

write_csv(block_boot_table, file.path(chapter_root, "07_Statistics", "spatial_block_bootstrap_structural.csv"))
message("Saved: spatial_block_bootstrap_structural.csv [REQUIRED before using",
        " 'analytically most robust' language in Discussion 7.2]")


# -------------------------------------------------------------
# SECTION 10: Save All Statistical Outputs
# -------------------------------------------------------------

out_path <- file.path(chapter_root, "07_Statistics")

write_csv(mw_table, file.path(out_path, "mann_whitney_results_full.csv"))
write_csv(sig_matrix_p %>% rownames_to_column("variable"), file.path(out_path, "significance_matrix_pvalues.csv"))
write_csv(sig_matrix_bin %>% rownames_to_column("variable"), file.path(out_path, "significance_matrix_binary.csv"))
write_csv(period_summary, file.path(chapter_root, "11_Tables", "Table_02_MW_Period_Summary.csv"))
message("\nSaved: mann_whitney_results_full.csv, significance matrices, Table_02")

save(mw_table, sig_matrix_p, sig_matrix_bin, dir_matrix, period_summary,
     multi_table, block_boot_table, full_matrix, full_matrix_coords, analysis_vars,
     file = file.path(out_path, "mann_whitney_results.RData"))
message("Saved: mann_whitney_results.RData")


# -------------------------------------------------------------
# SECTION 11: Log Session
# -------------------------------------------------------------

log_file <- file.path(chapter_root, "13_Logs", paste0("script06_log_", Sys.Date(), ".txt"))
sink(log_file)
cat("SCRIPT 06 LOG (REBUILT) — Mann-Whitney U Tests\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Analysis variables:", length(analysis_vars), "(burial_depth_m excluded,",
    "6 geochem major-oxide PCs added)\n")
cat("Total tests:", nrow(mw_table), "\n")
cat("Significant (p<0.05):", sum(mw_table$significant, na.rm=TRUE), "\n")
cat("Significant (Bonferroni):", sum(mw_table$sig_bonferroni, na.rm=TRUE), "\n\n")
cat("Period summary:\n"); print(period_summary)
cat("\nPersistent Places (MULTI vs background) — significant vars:\n")
print(multi_table %>% filter(significant==TRUE) %>% arrange(p_value))
cat("\nSpatial block bootstrap — headline structural claim:\n")
print(block_boot_table)
sink()
message("Log saved: ", log_file)


# -------------------------------------------------------------
# SCRIPT 06 COMPLETE
# -------------------------------------------------------------

message("\n=============================================================")
message("Script 06 (REBUILT) complete.")
message("Analysis vars: ", length(analysis_vars), " (was 25, now +6 geochem -1 burial)")
message("MW tests: ", nrow(mw_table), " | Persistent Places: ", nrow(multi_table))
message("Spatial block bootstrap: ", n_survive, "/", nrow(block_boot_table), " headline results survive")
message("Next: Script 07 — PCA + Random Forest")
message("  (FAMD/split ordinal-continuous vars, spatial block CV for RF,")
message("  permutation/grouped importance, CHELSA-LGM with/without sensitivity)")
message("=============================================================")