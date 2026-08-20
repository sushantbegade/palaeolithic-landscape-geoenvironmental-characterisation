# =============================================================
# SCRIPT 11 (NEW): District-Transfer Random Forest Validation
# =============================================================
# Project: Reading the Palaeolithic Landscape
# Chapter: Mapping the Past (Springer Nature, 2026)
# Author: Sushant Begade | RTMNU Nagpur
# ORCID: 0009-0003-0804-1763
# Created: August 2026
# =============================================================
# PURPOSE: Priority 4 of the outstanding-analyses list. Scripts 09-10
# established that cultural period is confounded with district for
# structural geology (severely) and geochemistry (partially). The
# Script 07 Random Forest (spatial block CV: 62.3% acc, AUC=0.685)
# was validated only against spatially-held-out blocks WITHIN the
# same pooled, cross-district sample — it has never been tested for
# whether it learned genuine, spatially transferable site/background
# environmental structure, or substantially learned to distinguish
# Chandrapur from Nagpur (which would masquerade as predictive skill
# without representing real environmental preference).
#
# Model 1: train on Chandrapur, predict on Nagpur (held out entirely)
# Model 2: train on Nagpur, predict on Chandrapur (held out entirely)
#
# If AUC collapses toward 0.5 on the held-out district, the model is
# substantially district-specific. If it holds close to the pooled
# spatial-CV AUC (0.685), this is evidence of genuinely transferable
# environmental signal, independent of the district confound
# affecting individual variables (Scripts 09-10).
# =============================================================

load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")
load(file.path(chapter_root, "07_Statistics/mann_whitney_results.RData"))
load(file.path(chapter_root, "07_Statistics/district_stratified_analysis.RData"))
load(file.path(chapter_root, "09_RandomForest/pcoa_rf_results.RData"))

library(tidyverse)
library(randomForest)
library(pROC)

set.seed(rf_seed)

message("Loaded full_matrix_district: ", nrow(full_matrix_district), " points")
message("Reference: Script 07 pooled spatial block CV — accuracy=", round(100*mean(spatial_cv_preds$observed==spatial_cv_preds$predicted),1),
        "%, AUC=", round(spatial_cv_auc,3))


# -------------------------------------------------------------
# SECTION 1: Rebuild Imputed Model Data (with district_label retained)
# -------------------------------------------------------------
# Same deterministic NA-filter/imputation logic as Script 07 Section 2,
# reproduced here (not reloaded) because Script 07's model_imputed did
# not retain district_label. Same seed, same rule -> identical rows.

model_data <- full_matrix_district %>%
  select(point_id, point_type, period, district_label, all_of(analysis_vars)) %>%
  filter(!is.na(point_type))

na_counts <- rowSums(is.na(model_data %>% select(all_of(analysis_vars))))
model_data <- model_data[na_counts <= 3, ]

model_imputed_d <- model_data
for (var in analysis_vars) {
  if (any(is.na(model_imputed_d[[var]]))) {
    med_val <- median(model_imputed_d[[var]], na.rm = TRUE)
    model_imputed_d[[var]][is.na(model_imputed_d[[var]])] <- med_val
  }
}
message("model_imputed_d rows: ", nrow(model_imputed_d), " (should match Script 07's 716)")
stopifnot("Row count mismatch vs Script 07 — investigate before trusting district-transfer RF" =
            nrow(model_imputed_d) == 716)

message("\nSite/background counts by district in the modelling sample:")
print(table(model_imputed_d$district_label, model_imputed_d$point_type))


# -------------------------------------------------------------
# SECTION 2: Helper — Fit on One District, Test on the Other
# -------------------------------------------------------------

run_district_transfer <- function(train_district, test_district, data, vars, seed) {
  message("\n--- Train: ", train_district, " | Test: ", test_district, " ---")
  
  train_raw <- data %>% filter(district_label == train_district) %>%
    mutate(response = factor(ifelse(point_type=="site","site","background"),
                             levels=c("site","background")))
  test_raw <- data %>% filter(district_label == test_district) %>%
    mutate(response = factor(ifelse(point_type=="site","site","background"),
                             levels=c("site","background")))
  
  # Balance TRAINING set only (standard practice, as Script 07). Test set is
  # left at its natural imbalanced ratio — AUC is threshold/prevalence-robust,
  # unlike raw accuracy, so this is the fairer transfer-performance metric.
  set.seed(seed)
  site_idx <- which(train_raw$response == "site")
  bg_idx   <- which(train_raw$response == "background")
  n_site_train <- length(site_idx)
  if (n_site_train < 10) {
    message("  Fewer than 10 site rows available for training in ", train_district,
            " — result will be unstable, interpret with extreme caution.")
  }
  bg_sampled <- sample(bg_idx, size = min(n_site_train, length(bg_idx)), replace = FALSE)
  train_balanced <- train_raw[c(site_idx, bg_sampled), ] %>%
    select(response, all_of(vars))
  
  message("  Training set (balanced): ", nrow(train_balanced), " rows (",
          sum(train_balanced$response=="site"), " site + ",
          sum(train_balanced$response=="background"), " background)")
  message("  Test set (natural imbalance, held out entirely): ", nrow(test_raw), " rows (",
          sum(test_raw$response=="site"), " site + ",
          sum(test_raw$response=="background"), " background)")
  
  fit <- randomForest(response ~ ., data = train_balanced, ntree = rf_ntree,
                      mtry = floor(sqrt(length(vars))), importance = TRUE)
  
  test_input <- test_raw %>% select(all_of(vars))
  pred_prob  <- predict(fit, test_input, type = "prob")[, "site"]
  pred_class <- predict(fit, test_input)
  
  acc <- mean(pred_class == test_raw$response)
  roc_obj <- tryCatch(
    roc(test_raw$response, pred_prob, levels = c("background","site"), quiet = TRUE),
    error = function(e) NULL
  )
  auc_val <- if (!is.null(roc_obj)) as.numeric(auc(roc_obj)) else NA
  
  message("  Transfer accuracy: ", round(100*acc,1), "% | Transfer AUC: ", round(auc_val,3))
  
  importance_transfer <- as.data.frame(importance(fit)) %>%
    rownames_to_column("variable") %>% arrange(desc(MeanDecreaseAccuracy)) %>%
    mutate(rank = row_number())
  
  list(model = fit, accuracy = acc, auc = auc_val,
       n_train = nrow(train_balanced), n_test = nrow(test_raw),
       importance = importance_transfer,
       confusion = table(observed = test_raw$response, predicted = pred_class))
}


# -------------------------------------------------------------
# SECTION 3: Run Both Transfer Directions
# -------------------------------------------------------------

message("\n=== DISTRICT-TRANSFER VALIDATION ===")

model1 <- run_district_transfer("Chandrapur", "Nagpur", model_imputed_d, analysis_vars, rf_seed)
message("  Confusion matrix (Chandrapur -> Nagpur):")
print(model1$confusion)

model2 <- run_district_transfer("Nagpur", "Chandrapur", model_imputed_d, analysis_vars, rf_seed)
message("  Confusion matrix (Nagpur -> Chandrapur):")
print(model2$confusion)


# -------------------------------------------------------------
# SECTION 4: Compare Against Pooled Spatial Block CV
# -------------------------------------------------------------

message("\n=== COMPARISON: pooled spatial-CV vs district-transfer ===")

comparison <- tibble(
  method = c("Pooled naive OOB (Script 07)",
             "Pooled spatial block 5-fold CV (Script 07)",
             "District-transfer: Chandrapur -> Nagpur",
             "District-transfer: Nagpur -> Chandrapur"),
  accuracy_pct = c(round((1-rf_model$err.rate[rf_ntree,"OOB"])*100,1),
                   round(100*mean(spatial_cv_preds$observed==spatial_cv_preds$predicted),1),
                   round(100*model1$accuracy,1),
                   round(100*model2$accuracy,1)),
  auc = c(NA, round(spatial_cv_auc,3), round(model1$auc,3), round(model2$auc,3))
)
print(comparison)

auc_drop_1 <- spatial_cv_auc - model1$auc
auc_drop_2 <- spatial_cv_auc - model2$auc
mean_transfer_auc <- mean(c(model1$auc, model2$auc), na.rm = TRUE)

message("\n=== VERDICT ===")
message("Mean district-transfer AUC: ", round(mean_transfer_auc,3),
        " vs pooled spatial-CV AUC: ", round(spatial_cv_auc,3))

if (mean_transfer_auc >= (spatial_cv_auc - 0.05)) {
  message("\nRESULT: district-transfer performance HOLDS UP (within 0.05 AUC of pooled).")
  message("This is evidence the RF captured genuinely transferable site/background")
  message("environmental structure, not merely district discrimination — a meaningfully")
  message("positive finding that should be reported prominently, especially given the")
  message("fragility found for individual structural variables (Scripts 09-10).")
} else if (mean_transfer_auc <= 0.55) {
  message("\nRESULT: district-transfer performance COLLAPSES toward chance (AUC<=0.55).")
  message("The pooled RF's apparent skill is substantially attributable to learning")
  message("district identity rather than transferable environmental preference. This")
  message("must be stated plainly in the manuscript — the suitability surface (Fig. 10)")
  message("should be reframed as descriptive of the sampled districts only, NOT as")
  message("evidence of a general Palaeolithic landscape suitability model.")
} else {
  message("\nRESULT: PARTIAL transfer — performance drops meaningfully but does not")
  message("collapse to chance. Report both directions' AUC explicitly rather than an")
  message("average; do not claim general transferability, but do not discard the model")
  message("either. State this as a moderate-confidence, district-limited result.")
}

message("\nDirection asymmetry check: Chandrapur->Nagpur AUC=", round(model1$auc,3),
        " vs Nagpur->Chandrapur AUC=", round(model2$auc,3))
if (abs(model1$auc - model2$auc) > 0.10) {
  message("Large asymmetry between transfer directions (>0.10 AUC) — one district's")
  message("environmental structure predicts the other much better than the reverse.")
  message("This itself is worth a sentence: transferability is not symmetric here.")
}


# -------------------------------------------------------------
# SECTION 5: Variable Importance Shift Under Transfer
# -------------------------------------------------------------

message("\n=== Variable importance: pooled vs district-transfer models ===")

pooled_importance <- importance_df %>% select(variable, pooled_rank = rank)

importance_shift <- pooled_importance %>%
  left_join(model1$importance %>% select(variable, chandrapur_to_nagpur_rank = rank), by = "variable") %>%
  left_join(model2$importance %>% select(variable, nagpur_to_chandrapur_rank = rank), by = "variable") %>%
  arrange(pooled_rank)

message("\nTop 10 pooled variables and their rank under each transfer direction:")
print(importance_shift %>% head(10))


# -------------------------------------------------------------
# SECTION 6: Save Outputs
# -------------------------------------------------------------

out_path <- file.path(chapter_root, "09_RandomForest")

write_csv(comparison, file.path(chapter_root, "11_Tables", "Table_district_transfer_RF.csv"))
write_csv(importance_shift, file.path(out_path, "district_transfer_importance_shift.csv"))

message("\nSaved: Table_district_transfer_RF.csv [table-budget note: same 15-item cap")
message("issue as Script 10 — decide which existing table this replaces or is folded into]")
message("Saved: district_transfer_importance_shift.csv")

save(model1, model2, comparison, importance_shift, mean_transfer_auc,
     file = file.path(out_path, "district_transfer_rf.RData"))
message("Saved: district_transfer_rf.RData")


# -------------------------------------------------------------
# SECTION 7: Log Session
# -------------------------------------------------------------

log_file <- file.path(chapter_root, "13_Logs", paste0("script11_log_", Sys.Date(), ".txt"))
sink(log_file)
cat("SCRIPT 11 LOG — District-Transfer Random Forest (Priority 4)\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Pooled spatial-CV: accuracy=", round(100*mean(spatial_cv_preds$observed==spatial_cv_preds$predicted),1),
    "%, AUC=", round(spatial_cv_auc,3), "\n")
cat("Chandrapur->Nagpur: accuracy=", round(100*model1$accuracy,1), "%, AUC=", round(model1$auc,3), "\n")
cat("Nagpur->Chandrapur: accuracy=", round(100*model2$accuracy,1), "%, AUC=", round(model2$auc,3), "\n\n")
cat("Comparison table:\n"); print(comparison)
cat("\nImportance shift (top 10 pooled):\n"); print(importance_shift %>% head(10))
sink()
message("Log saved: ", log_file)


# -------------------------------------------------------------
# SCRIPT 11 COMPLETE
# -------------------------------------------------------------

message("\n=============================================================")
message("Script 11 complete.")
message("Chandrapur->Nagpur AUC: ", round(model1$auc,3), " | Nagpur->Chandrapur AUC: ", round(model2$auc,3))
message("Mean transfer AUC: ", round(mean_transfer_auc,3), " vs pooled spatial-CV AUC: ", round(spatial_cv_auc,3))
message("This determines whether Fig. 10 (suitability surface) and Section 6.3/7.1's")
message("RF discussion can claim general transferable predictive skill, or must be")
message("reframed as district-limited/descriptive only.")
message("Remaining priorities: 3=survey-envelope sensitivity, 5=block-size sensitivity,")
message("6=PCoA missingness+PERMANOVA, 7=geochem buffer-scale+zero-replacement,")
message("8=manual classification check. Per reviewer instruction: continue analyses")
message("before the next manuscript rewrite.")
message("=============================================================")