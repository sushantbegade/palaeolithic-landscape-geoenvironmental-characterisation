# =============================================================
# SCRIPT 07 (REBUILT): PCoA + Random Forest (Spatially Validated)
# =============================================================
# Project: Reading the Palaeolithic Landscape
# Chapter: Mapping the Past (Springer Nature, 2026)
# Author: Sushant Begade | RTMNU Nagpur
# ORCID: 0009-0003-0804-1763
# Rebuild date: August 2026
# =============================================================
# CHANGES FROM ORIGINAL SCRIPT 07:
#   1. PCA -> Gower distance + PCoA. Original ran Euclidean PCA on
#      25 z-scored variables including 5 ORDINAL soil-class variables
#      (soil_texture/slope/erosion/productivity coded 1-5) alongside
#      continuous variables — invalid (reviewer Problem #10). Gower
#      distance handles mixed continuous/ordinal data correctly;
#      classical PCoA (cmdscale) on the Gower matrix replaces prcomp.
#   2. Random Forest — TWO validation regimes now reported side by
#      side: naive OOB (as before, kept for comparability) AND
#      spatial block k-fold CV (5-fold, using the SAME block_id
#      scheme from Script 06's bootstrap) via manual fold loop +
#      pROC AUC. This directly answers reviewer Problem #4/#7 — no
#      longer claiming OOB accuracy as if it were spatially
#      independent performance.
#   3. Importance — primary metric switched to permutation-based
#      MeanDecreaseAccuracy (already computed by randomForest with
#      importance=TRUE, just wasn't the headline metric before).
#      NEW: grouped importance — correlated variable blocks (NDVI/
#      SAVI/MSAVI; NDWI/MNDWI; CHELSA-LGM; CHELSA-modern; structural;
#      pedology; geochem) permuted JOINTLY, answering reviewer
#      Problem #12 (duplicated environmental information inflating
#      individual-variable importance).
#   4. CHELSA-LGM sensitivity — RF run with all 30 vars AND with the
#      3 chelsa_*_lgm vars dropped (27 vars), rankings compared side
#      by side. Resolves the flag carried from Script 06.
#   5. burial_depth_m remains excluded (Script 05 finding).
# =============================================================

load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")
load(file.path(chapter_root, "07_Statistics/mann_whitney_results.RData"))
load(file.path(chapter_root, "03_Extracted_Values/master_extraction_matrix.RData"))

library(terra)
library(sf)
library(tidyverse)
library(cluster)      # Gower distance
library(randomForest)
library(pROC)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(viridis)

set.seed(rf_seed)

message("All data loaded. full_matrix_coords: ", nrow(full_matrix_coords), " rows")
message("Analysis variables: ", length(analysis_vars))


# -------------------------------------------------------------
# SECTION 2: Prepare Clean Matrix
# -------------------------------------------------------------

ORDINAL_VARS <- c("soil_depth", "soil_erosion", "soil_productivity",
                  "soil_slope", "soil_texture")
CONTINUOUS_VARS <- setdiff(analysis_vars, ORDINAL_VARS)

message("\nOrdinal variables (", length(ORDINAL_VARS), "): ", paste(ORDINAL_VARS, collapse=", "))
message("Continuous variables (", length(CONTINUOUS_VARS), ")")

model_data <- full_matrix_coords %>%
  select(point_id, point_type, period, block_id, all_of(analysis_vars)) %>%
  filter(!is.na(point_type))

na_counts <- rowSums(is.na(model_data %>% select(all_of(analysis_vars))))
message("\nNA distribution across rows: max ", max(na_counts), " missing of ", length(analysis_vars), " vars")
model_data <- model_data[na_counts <= 3, ]
message("Rows after NA filtering (<=3 missing): ", nrow(model_data))

model_imputed <- model_data
for (var in analysis_vars) {
  if (any(is.na(model_imputed[[var]]))) {
    med_val <- median(model_imputed[[var]], na.rm = TRUE)
    model_imputed[[var]][is.na(model_imputed[[var]])] <- med_val
  }
}
message("Complete cases after imputation: ", sum(complete.cases(model_imputed %>% select(all_of(analysis_vars)))))


# =============================================================
# PART A: PCoA (Gower distance) — replaces Euclidean PCA
# =============================================================

message("\n--- PART A: PCoA on Gower Distance (mixed continuous+ordinal) ---")

# NOTE: full Gower distance matrix on >1000 points is O(n^2) memory —
# for n~1100 this is ~1.2M cells, fine. If this ever runs on a much
# larger n, subsample or use cluster::daisy's memory-efficient path.

gower_input <- model_imputed %>% select(all_of(analysis_vars))
for (v in ORDINAL_VARS) gower_input[[v]] <- factor(gower_input[[v]], ordered = TRUE)

message("Computing Gower distance matrix (", nrow(gower_input), " x ", nrow(gower_input), ")...")
gower_dist <- daisy(gower_input, metric = "gower",
                    type = list(ordratio = ORDINAL_VARS))
message("Gower distance computed.")

message("Running classical PCoA (cmdscale)...")
pcoa_result <- cmdscale(gower_dist, k = 10, eig = TRUE)

eig_positive <- pcoa_result$eig[pcoa_result$eig > 0]
var_explained <- round(100 * eig_positive / sum(eig_positive), 2)
message("\nVariance explained by first 5 positive-eigenvalue axes:")
print(head(var_explained, 5))
message("Cumulative (first 2 axes): ", round(sum(head(var_explained,2)),1), "%")

pcoa_scores <- as.data.frame(pcoa_result$points) %>% rename_with(~paste0("PC", seq_along(.)))
pcoa_scores <- bind_cols(model_imputed %>% select(point_id, point_type, period), pcoa_scores)

sites_pcoa <- pcoa_scores %>% filter(point_type == "site")
bg_pcoa    <- pcoa_scores %>% filter(point_type == "background")

# Period centroids (LP/MP/UP only — MULTI excluded from centroid analysis
# for consistency with the period-stratified framework; reported
# separately in Persistent Places if desired)
period_centroids <- sites_pcoa %>% filter(period %in% c("LP","MP","UP")) %>%
  group_by(period) %>%
  summarise(PC1_mean = mean(PC1, na.rm=TRUE), PC2_mean = mean(PC2, na.rm=TRUE),
            PC1_sd = sd(PC1, na.rm=TRUE), PC2_sd = sd(PC2, na.rm=TRUE),
            n = n(), .groups = "drop")

message("\nPeriod centroids (PCoA axes 1-2):")
print(period_centroids)

lp_c <- period_centroids %>% filter(period=="LP")
mp_c <- period_centroids %>% filter(period=="MP")
up_c <- period_centroids %>% filter(period=="UP")

disp_lp_mp <- sqrt((mp_c$PC1_mean-lp_c$PC1_mean)^2 + (mp_c$PC2_mean-lp_c$PC2_mean)^2)
disp_mp_up <- sqrt((up_c$PC1_mean-mp_c$PC1_mean)^2 + (up_c$PC2_mean-mp_c$PC2_mean)^2)
disp_lp_up <- sqrt((up_c$PC1_mean-lp_c$PC1_mean)^2 + (up_c$PC2_mean-lp_c$PC2_mean)^2)

message("\nCentroid displacements (PCoA space, descriptive — not yet significance-tested):")
message("LP-MP: ", round(disp_lp_mp,3), " | MP-UP: ", round(disp_mp_up,3), " | LP-UP: ", round(disp_lp_up,3))

# -------------------------------------------------------------
# Bootstrap CI on centroid displacements [NEW — answers reviewer
# Problem #11: "do not call this diachronic niche displacement...
# perform bootstrap confidence regions"]
# -------------------------------------------------------------

message("\nBootstrapping centroid displacement CIs (1000 iterations)...")
boot_displacement <- function(data, period_a, period_b, n_boot = 1000) {
  a_pts <- data %>% filter(period == period_a)
  b_pts <- data %>% filter(period == period_b)
  boot_d <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    a_s <- a_pts[sample(nrow(a_pts), replace = TRUE), ]
    b_s <- b_pts[sample(nrow(b_pts), replace = TRUE), ]
    boot_d[i] <- sqrt((mean(b_s$PC1)-mean(a_s$PC1))^2 + (mean(b_s$PC2)-mean(a_s$PC2))^2)
  }
  quantile(boot_d, c(0.025, 0.975))
}

ci_lp_mp <- boot_displacement(sites_pcoa %>% filter(period %in% c("LP","MP")), "LP", "MP")
ci_mp_up <- boot_displacement(sites_pcoa %>% filter(period %in% c("MP","UP")), "MP", "UP")
ci_lp_up <- boot_displacement(sites_pcoa %>% filter(period %in% c("LP","UP")), "LP", "UP")

displacements <- tibble(
  transition = c("LP_to_MP", "MP_to_UP", "LP_to_UP"),
  displacement = round(c(disp_lp_mp, disp_mp_up, disp_lp_up), 4),
  boot_ci_low  = round(c(ci_lp_mp[1], ci_mp_up[1], ci_lp_up[1]), 4),
  boot_ci_high = round(c(ci_lp_mp[2], ci_mp_up[2], ci_lp_up[2]), 4),
  interpretation = c(
    "LP to MP multivariate environmental centroid separation (descriptive; bootstrap CI, not PERMANOVA)",
    "MP to UP multivariate environmental centroid separation (descriptive; UP n=27, caution)",
    "LP to UP multivariate environmental centroid separation (descriptive)"
  )
)
message("\nCentroid displacements with bootstrap 95% CI:")
print(displacements)
message("\nNOTE: UP n=27 (small sample per earlier coordinate-precision finding —")
message("78% of UP sites are GIS-digitised/low-precision, Script 01). UP centroid")
message("position should be reported with explicit caution regardless of CI width.")


# =============================================================
# PART B: Random Forest — Spatial Block CV + Grouped Importance
# =============================================================

message("\n--- PART B: Random Forest ---")

rf_data_full <- model_imputed %>%
  select(point_type, block_id, all_of(analysis_vars)) %>%
  mutate(response = factor(ifelse(point_type=="site","site","background"),
                           levels=c("site","background"))) %>%
  select(-point_type)   # FIX: point_type was leaking as a predictor — it IS
# the response 1:1, caused the 100% OOB / AUC=1 /
# zero grouped-importance result in the prior run.
# Must never appear in any model_data / fold data
# passed to randomForest() or predict() below.

n_sites <- sum(rf_data_full$response == "site")
n_bg    <- sum(rf_data_full$response == "background")
message("Sites: ", n_sites, " | Background: ", n_bg)

set.seed(rf_seed)
site_idx <- which(rf_data_full$response == "site")
bg_idx   <- which(rf_data_full$response == "background")
bg_sampled <- sample(bg_idx, size = n_sites, replace = FALSE)
rf_balanced <- rf_data_full[c(site_idx, bg_sampled), ]
message("Balanced RF dataset: ", nrow(rf_balanced), " rows")


# -------------------------------------------------------------
# Helper: fit RF, return model + naive OOB + permutation/grouped importance
# -------------------------------------------------------------

VAR_GROUPS <- list(
  vegetation      = c("NDVI","SAVI","MSAVI"),
  water           = c("NDWI","MNDWI"),
  surface         = c("NDBI","BSI"),
  climate_lgm     = c("chelsa_bio01_lgm","chelsa_bio12_lgm","chelsa_bio15_lgm"),
  climate_modern  = c("chelsa_bio01_modern","chelsa_bio12_modern","chelsa_bio15_modern"),
  structural      = c("dist_fault","dist_dyke","dist_lineament","dist_shear","dist_mineral"),
  pedology        = ORDINAL_VARS,
  geochem_stream  = c("geochem_stream_major_PC1","geochem_stream_major_PC2"),
  geochem_horizon = c("geochem_horizon_major_PC1","geochem_horizon_major_PC2"),
  geochem_regolith= c("geochem_regolith_major_PC1","geochem_regolith_major_PC2"),
  elevation       = c("elevation")
)

grouped_importance <- function(rf_model, data, vars, groups, n_reps = 25) {
  # FIX: single-shuffle permutation on n=244 rows is too noisy — accuracy
  # only moves in ~0.004 increments per misclassification, so several
  # groups showed identical or exact-zero acc_drop by chance in the first
  # run, not because they carry no signal. Repeat each group's permutation
  # n_reps times and average, standard practice for permutation importance
  # stability (Breiman 2001; Strobl et al 2008, both already in refs).
  baseline_pred <- predict(rf_model, data)
  baseline_acc  <- mean(baseline_pred == data$response)
  
  results <- tibble(group = character(), n_vars = integer(),
                    acc_drop_mean = numeric(), acc_drop_sd = numeric())
  for (g in names(groups)) {
    gv <- intersect(groups[[g]], vars)
    if (length(gv) == 0) next
    rep_drops <- numeric(n_reps)
    for (r in seq_len(n_reps)) {
      perm_data <- data
      perm_idx <- sample(nrow(perm_data))
      perm_data[gv] <- perm_data[gv][perm_idx, ]
      perm_pred <- predict(rf_model, perm_data)
      perm_acc  <- mean(perm_pred == data$response)
      rep_drops[r] <- baseline_acc - perm_acc
    }
    results <- bind_rows(results, tibble(group = g, n_vars = length(gv),
                                         acc_drop_mean = mean(rep_drops),
                                         acc_drop_sd = sd(rep_drops)))
  }
  message("  NOTE: acc_drop computed on TRAINING data (rf_formula_data), not")
  message("  strictly held-out OOB rows — an optimistic-leaning estimate, same")
  message("  caveat as applies to naive OOB elsewhere in this script. Relative")
  message("  ranking across groups is still informative even if absolute")
  message("  magnitudes run slightly high.")
  results %>% arrange(desc(acc_drop_mean))
}


# -------------------------------------------------------------
# Fit primary RF model (all 30 vars) — naive OOB
# -------------------------------------------------------------

message("\nFitting primary RF (all ", length(analysis_vars), " vars, ntree=", rf_ntree, ")...")

rf_formula_data <- rf_balanced %>% select(-block_id)
rf_model <- randomForest(response ~ ., data = rf_formula_data, ntree = rf_ntree,
                         mtry = floor(sqrt(length(analysis_vars))),
                         importance = TRUE, keep.forest = TRUE)

oob_acc <- round((1 - rf_model$err.rate[rf_ntree, "OOB"]) * 100, 1)
message("Naive OOB accuracy: ", oob_acc, "%")
message("Confusion matrix:"); print(rf_model$confusion)


# -------------------------------------------------------------
# Spatial Block k-fold CV [NEW]
# -------------------------------------------------------------

message("\n=== SPATIAL BLOCK 5-FOLD CROSS-VALIDATION ===")
message("Using the SAME block_id scheme as Script 06's bootstrap (10km grid).")

k_folds <- 5
unique_blocks <- unique(rf_balanced$block_id)
set.seed(rf_seed)
block_fold_assign <- sample(rep(1:k_folds, length.out = length(unique_blocks)))
names(block_fold_assign) <- unique_blocks

rf_balanced$fold <- block_fold_assign[rf_balanced$block_id]

spatial_cv_preds <- tibble()
for (k in 1:k_folds) {
  train_data <- rf_balanced %>% filter(fold != k) %>% select(-block_id, -fold)
  test_data  <- rf_balanced %>% filter(fold == k) %>% select(-block_id, -fold)
  
  if (nrow(test_data) < 5 || length(unique(train_data$response)) < 2) {
    message("  Fold ", k, ": skipped (insufficient data)")
    next
  }
  
  fold_model <- randomForest(response ~ ., data = train_data, ntree = rf_ntree,
                             mtry = floor(sqrt(length(analysis_vars))))
  fold_prob  <- predict(fold_model, test_data, type = "prob")[, "site"]
  fold_pred  <- predict(fold_model, test_data)
  
  spatial_cv_preds <- bind_rows(spatial_cv_preds, tibble(
    fold = k, observed = test_data$response, predicted = fold_pred, prob_site = fold_prob
  ))
  message("  Fold ", k, ": n_test=", nrow(test_data), " done")
}

spatial_cv_acc <- mean(spatial_cv_preds$observed == spatial_cv_preds$predicted)
spatial_cv_roc <- roc(spatial_cv_preds$observed, spatial_cv_preds$prob_site,
                      levels = c("background", "site"), quiet = TRUE)
spatial_cv_auc <- as.numeric(auc(spatial_cv_roc))

message("\n=== NAIVE OOB vs SPATIAL BLOCK CV ===")
message("Naive OOB accuracy:        ", oob_acc, "%")
message("Spatial block CV accuracy: ", round(100*spatial_cv_acc, 1), "%")
message("Spatial block CV AUC:      ", round(spatial_cv_auc, 3))
if ((oob_acc/100 - spatial_cv_acc) > 0.05) {
  message("\n*** OOB overstates performance by >5 percentage points under spatial")
  message("*** dependence. REPORT SPATIAL BLOCK CV FIGURES IN RESULTS, not just OOB.")
} else {
  message("\nOOB and spatial-CV accuracy are reasonably close — spatial dependence")
  message("inflation appears modest for this model, but report both regardless.")
}

validation_comparison <- tibble(
  method = c("Naive OOB", "Spatial Block 5-fold CV"),
  accuracy_pct = c(oob_acc, round(100*spatial_cv_acc,1)),
  auc = c(NA, round(spatial_cv_auc, 3))
)


# -------------------------------------------------------------
# Grouped Importance [NEW] + Standard Permutation Importance
# -------------------------------------------------------------

message("\n=== IMPORTANCE: permutation (MeanDecreaseAccuracy) as primary ===")

importance_df <- as.data.frame(importance(rf_model)) %>%
  rownames_to_column("variable") %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  mutate(
    rank = row_number(),
    variable_clean = gsub("_", " ", variable),
    category = case_when(
      str_detect(variable, "NDVI|NDWI|MNDWI|NDBI|BSI|SAVI|MSAVI") ~ "Spectral",
      str_detect(variable, "soil_") ~ "Pedogenic",
      str_detect(variable, "chelsa") ~ "Climate",
      str_detect(variable, "dist_") ~ "Structural Geology",
      str_detect(variable, "elevation") ~ "Topographic",
      str_detect(variable, "geochem_") ~ "Geochemistry",
      TRUE ~ "Other"
    )
  )

message("\nTop 15 by permutation importance (MeanDecreaseAccuracy):")
print(importance_df %>% select(rank, variable, MeanDecreaseAccuracy, MeanDecreaseGini, category) %>% head(15))

message("\nRunning grouped-variable permutation importance...")
grp_importance <- grouped_importance(rf_model, rf_formula_data, analysis_vars, VAR_GROUPS)
message("\nGrouped importance (accuracy drop when whole group permuted jointly):")
print(grp_importance)
message("\nCompare to individual vegetation vars (NDVI/SAVI/MSAVI) in importance_df —")
message("if grouped 'vegetation' importance is much lower than 3x any single")
message("vegetation var's individual importance, this confirms redundancy")
message("inflation (reviewer Problem #12) — report grouped importance as the")
message("primary evidence for spectral/climate variable contribution in Discussion.")


# -------------------------------------------------------------
# CHELSA-LGM Sensitivity [NEW]
# -------------------------------------------------------------

message("\n=== CHELSA-LGM SENSITIVITY: RF with vs without LGM variables ===")

vars_no_lgm <- setdiff(analysis_vars, c("chelsa_bio01_lgm","chelsa_bio12_lgm","chelsa_bio15_lgm"))
rf_data_no_lgm <- rf_balanced %>% select(response, all_of(vars_no_lgm))

rf_model_no_lgm <- randomForest(response ~ ., data = rf_data_no_lgm, ntree = rf_ntree,
                                mtry = floor(sqrt(length(vars_no_lgm))), importance = TRUE)

oob_acc_no_lgm <- round((1 - rf_model_no_lgm$err.rate[rf_ntree,"OOB"]) * 100, 1)

importance_no_lgm <- as.data.frame(importance(rf_model_no_lgm)) %>%
  rownames_to_column("variable") %>% arrange(desc(MeanDecreaseAccuracy)) %>%
  mutate(rank_no_lgm = row_number())

lgm_comparison <- importance_df %>% select(variable, rank_with_lgm = rank) %>%
  left_join(importance_no_lgm %>% select(variable, rank_no_lgm), by = "variable") %>%
  filter(!is.na(rank_no_lgm)) %>%
  mutate(rank_shift = rank_with_lgm - rank_no_lgm) %>%
  arrange(rank_with_lgm)

message("With-LGM OOB: ", oob_acc, "% | Without-LGM OOB: ", oob_acc_no_lgm, "%")
message("\nTop-10 variable rank comparison (with vs without CHELSA-LGM):")
print(lgm_comparison %>% head(10))
message("\nIf top structural/pedogenic variables hold similar rank either way,")
message("the chapter's core findings do NOT depend on the temporally-mismatched")
message("CHELSA-LGM variables — state this explicitly in Methods 5.5/Discussion 7.4.")


# -------------------------------------------------------------
# Geoenvironmental Suitability Surface (unchanged approach, uses
# primary rf_model; caveat language strengthened)
# -------------------------------------------------------------

message("\nGenerating geoenvironmental suitability surface (primary RF model)...")
message("NOTE: this surface reflects the NAIVE-OOB-fitted model. Given the OOB")
message("vs spatial-CV gap reported above, treat absolute suitability VALUES as")
message("relative/ordinal (higher vs lower) rather than calibrated probabilities —")
message("already the chapter's stated convention, now with a quantified basis for it.")

boundary <- st_read(paths$boundary, quiet = TRUE) %>% st_transform(crs_utm44n)
study_area_vect <- vect(st_union(boundary))
dem_ref <- rast(paths$dem) %>% terra::project("EPSG:32644") %>%
  crop(study_area_vect) %>% mask(study_area_vect) %>% aggregate(fact = 8)

raster_vars <- list(
  elevation = paths$dem, NDVI = paths$NDVI, NDWI = paths$NDWI, MNDWI = paths$MNDWI,
  NDBI = paths$NDBI, BSI = paths$BSI, SAVI = paths$SAVI, MSAVI = paths$MSAVI,
  soil_depth = paths$soil_depth, soil_erosion = paths$soil_erosion,
  soil_productivity = paths$soil_productivity, soil_slope = paths$soil_slope,
  soil_texture = paths$soil_texture,
  chelsa_bio01_lgm = paths$chelsa_bio01_lgm, chelsa_bio12_lgm = paths$chelsa_bio12_lgm,
  chelsa_bio15_lgm = paths$chelsa_bio15_lgm,
  chelsa_bio01_modern = paths$chelsa_bio01_modern, chelsa_bio12_modern = paths$chelsa_bio12_modern,
  chelsa_bio15_modern = paths$chelsa_bio15_modern
)

raster_stack <- list()
for (name in names(raster_vars)) {
  tryCatch({
    r <- rast(raster_vars[[name]]) %>% terra::project("EPSG:32644") %>%
      crop(study_area_vect) %>% mask(study_area_vect) %>% resample(dem_ref, method="bilinear")
    raster_stack[[name]] <- r
  }, error = function(e) message("  SKIP: ", name, " - ", e$message))
}
stack_rast <- rast(raster_stack); names(stack_rast) <- names(raster_stack)
pred_df <- as.data.frame(stack_rast, na.rm = FALSE)

# structural + geochem vars have no continuous raster surface — fill median
# (unchanged limitation from original; flagged explicitly)
missing_vars <- setdiff(analysis_vars, names(pred_df))
message("\nVariables filled with sample median for prediction surface (no raster",
        " available): ", paste(missing_vars, collapse=", "))
for (mv in missing_vars) pred_df[[mv]] <- median(model_imputed[[mv]], na.rm = TRUE)

for (v in names(pred_df)) {
  if (any(is.na(pred_df[[v]]))) pred_df[[v]][is.na(pred_df[[v]])] <- median(pred_df[[v]], na.rm=TRUE)
}

pred_input <- pred_df %>% select(all_of(analysis_vars))
rf_probs <- predict(rf_model, pred_input, type = "prob")
site_probs <- rf_probs[, "site"]

attract_rast <- dem_ref
values(attract_rast) <- site_probs
attract_rast <- mask(attract_rast, study_area_vect)
message("Suitability surface range: ", round(minmax(attract_rast)[1],3), " to ",
        round(minmax(attract_rast)[2],3))

writeRaster(attract_rast, file.path(chapter_root, "08_PCA", "geoenvironmental_suitability.tif"),
            overwrite = TRUE)


# -------------------------------------------------------------
# Save All Outputs
# -------------------------------------------------------------

message("\nSaving outputs...")

write_csv(pcoa_scores, file.path(chapter_root, "08_PCA", "pcoa_scores.csv"))
write_csv(period_centroids, file.path(chapter_root, "08_PCA", "period_centroids.csv"))
write_csv(displacements, file.path(chapter_root, "08_PCA", "centroid_displacement_vectors_bootstrapCI.csv"))
write_csv(tibble(axis = paste0("PCoA",seq_along(var_explained)), variance_pct = var_explained),
          file.path(chapter_root, "08_PCA", "pcoa_variance_explained.csv"))

write_csv(importance_df, file.path(chapter_root, "09_RandomForest", "rf_variable_importance.csv"))
write_csv(importance_df %>% head(15), file.path(chapter_root, "11_Tables", "Table_03_RF_Variable_Importance.csv"))
write_csv(grp_importance, file.path(chapter_root, "09_RandomForest", "rf_grouped_importance.csv"))
write_csv(validation_comparison, file.path(chapter_root, "11_Tables", "Table_RF_validation_comparison.csv"))
write_csv(lgm_comparison, file.path(chapter_root, "09_RandomForest", "chelsa_lgm_sensitivity.csv"))
write_csv(displacements, file.path(chapter_root, "11_Tables", "Table_04_Centroid_Displacements.csv"))

save(pcoa_result, pcoa_scores, period_centroids, displacements, gower_dist,
     rf_model, rf_model_no_lgm, importance_df, grp_importance,
     validation_comparison, lgm_comparison, spatial_cv_preds, spatial_cv_auc,
     attract_rast,
     file = file.path(chapter_root, "09_RandomForest", "pcoa_rf_results.RData"))
message("Saved: all PCoA/RF outputs")


# -------------------------------------------------------------
# Log + Final Summary
# -------------------------------------------------------------

log_file <- file.path(chapter_root, "13_Logs", paste0("script07_log_", Sys.Date(), ".txt"))
sink(log_file)
cat("SCRIPT 07 LOG (REBUILT) — PCoA + Random Forest\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("PCoA variance (axis 1-2):", round(sum(head(var_explained,2)),1), "%\n")
cat("Centroid displacements + bootstrap CI:\n"); print(displacements)
cat("\nNaive OOB:", oob_acc, "% | Spatial Block CV:", round(100*spatial_cv_acc,1),
    "% | AUC:", round(spatial_cv_auc,3), "\n")
cat("\nTop 15 permutation importance:\n")
print(importance_df %>% head(15) %>% select(rank, variable, MeanDecreaseAccuracy, category))
cat("\nGrouped importance:\n"); print(grp_importance)
cat("\nCHELSA-LGM sensitivity (with vs without OOB):", oob_acc, "vs", oob_acc_no_lgm, "\n")
sink()
message("Log saved: ", log_file)

message("\n=============================================================")
message("Script 07 (REBUILT) complete.")
message("PCoA (Gower, mixed types) replaces invalid Euclidean PCA.")
message("RF: naive OOB=", oob_acc, "% | spatial block CV=", round(100*spatial_cv_acc,1),
        "% | AUC=", round(spatial_cv_auc,3))
message("Primary importance metric: permutation (MeanDecreaseAccuracy), +grouped importance")
message("CHELSA-LGM sensitivity checked: with=", oob_acc, "% without=", oob_acc_no_lgm, "%")
message("Next: Script 08 — Publication Figures (rebuild all from corrected outputs)")
message("=============================================================")