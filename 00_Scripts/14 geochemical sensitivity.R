# =============================================================
# SCRIPT 14 (NEW): Geochemical Buffer-Scale + Zero-Replacement Sensitivity
# =============================================================
# Project: Reading the Palaeolithic Landscape
# Chapter: Mapping the Past (Springer Nature, 2026)
# Author: Sushant Begade | RTMNU Nagpur
# ORCID: 0009-0003-0804-1763
# Created: August 2026
# =============================================================
# PURPOSE: Recommended task 4. The 5km buffer scale was chosen on
# coverage grounds (Script 04/Section 5.3), not independently
# established; the zero-replacement multiplier (half column minimum)
# was flagged as a pragmatic, untested choice. This tests both.
#
# SCOPE: focused on stream sediment major oxides — the best-coverage
# layer (98.1% at 5km) and the one carrying the chapter's cited
# geochem_stream_major_PC1/PC2 results (Sections 6.1, 6.2, 6.3 of the
# manuscript — the RF importance ranking, the district-stratified
# geochemistry check, and the district-transfer importance shift all
# reference this specific layer). Horizon/regolith layers have far
# lower coverage (59.6%/46.1%) and were not central findings — testing
# them is a lower-value use of time and is left as further work if
# needed.
# =============================================================

load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")
load(file.path(chapter_root, "05_Structural_Geology/vector_extraction_matrix.RData"))
load(file.path(chapter_root, "07_Statistics/mann_whitney_results.RData"))

library(tidyverse)

set.seed(rf_seed)

message("Loaded geoch_matrix: ", nrow(geoch_matrix), " rows x ", ncol(geoch_matrix), " cols")
message("(contains raw buffer means at 1km, 2km, AND 5km — no new extraction needed)")

MAJOR_OXIDES <- c("sio2","al2o3","fe2o3","tio2","cao","mgo","mno","na2o","k2o","p2o5")


# -------------------------------------------------------------
# SECTION 1: Reusable CLR-PCA Function (parameterised by buffer + zero-mult)
# -------------------------------------------------------------

clr_transform <- function(mat) {
  log_mat <- log(mat)
  gm <- rowMeans(log_mat, na.rm = TRUE)
  sweep(log_mat, 1, gm, "-")
}

replace_zeros <- function(mat, mult = 0.5, min_positive_frac = 0.5) {
  dropped <- character(0); keep <- rep(TRUE, ncol(mat))
  for (j in seq_len(ncol(mat))) {
    col <- mat[, j]
    frac_pos <- mean(col > 0, na.rm = TRUE)
    if (frac_pos < min_positive_frac) { keep[j] <- FALSE; dropped <- c(dropped, colnames(mat)[j]); next }
    pos_min <- min(col[col > 0], na.rm = TRUE)
    mat[!is.na(col) & col <= 0, j] <- pos_min * mult
  }
  mat[, keep, drop = FALSE]
}

run_stream_clr_pca <- function(buffer_km, zero_mult, winsorize = TRUE) {
  buf_label <- paste0("_stream_", buffer_km, "km")
  layer_cols <- names(geoch_matrix)[str_ends(names(geoch_matrix), buf_label)]
  layer_cols <- setdiff(layer_cols, paste0("n_pts", buf_label))
  element_base <- str_remove(layer_cols, buf_label)
  layer_cols <- layer_cols[element_base %in% MAJOR_OXIDES]
  
  sub <- geoch_matrix %>% select(point_id, all_of(layer_cols))
  complete_idx <- complete.cases(sub %>% select(-point_id))
  n_complete <- sum(complete_idx)
  if (n_complete < 30) return(NULL)
  
  mat <- as.matrix(sub[complete_idx, layer_cols])
  if (winsorize) {
    for (j in seq_len(ncol(mat))) {
      qs <- quantile(mat[,j], c(0.01,0.99), na.rm=TRUE)
      mat[,j] <- pmin(pmax(mat[,j], qs[1]), qs[2])
    }
  }
  mat <- replace_zeros(mat, mult = zero_mult)
  if (ncol(mat) < 3) return(NULL)
  clr_mat <- clr_transform(mat)
  
  pca_fit <- prcomp(clr_mat, center = TRUE, scale. = FALSE)
  var_pct <- round(100 * summary(pca_fit)$importance[2, 1:2], 1)
  scores <- as.data.frame(pca_fit$x[, 1:2])
  names(scores) <- c("PC1", "PC2")
  scores$point_id <- sub$point_id[complete_idx]
  
  list(scores = scores, loadings = pca_fit$rotation[,1:2], var_pct = var_pct, n = n_complete)
}


# =============================================================
# PART A: Buffer-Scale Sensitivity (1km / 2km / 5km, mult=0.5 fixed)
# =============================================================

message("\n=== PART A: Buffer-scale sensitivity (stream sediment major oxides) ===")

buffer_results <- list(
  "1km" = run_stream_clr_pca(1, zero_mult = 0.5),
  "2km" = run_stream_clr_pca(2, zero_mult = 0.5),
  "5km" = run_stream_clr_pca(5, zero_mult = 0.5)
)

for (nm in names(buffer_results)) {
  r <- buffer_results[[nm]]
  if (is.null(r)) { message(nm, ": insufficient data, skipped"); next }
  message(nm, ": n=", r$n, " | PC1 var=", r$var_pct[1], "% | PC2 var=", r$var_pct[2], "%")
  message("  Top-3 PC1 loadings: ", paste(names(sort(abs(r$loadings[,1]), decreasing=TRUE))[1:3], collapse=", "))
}

# Correlate PC1 scores across buffer scales (common point_ids only)
message("\n--- PC1 score correlation across buffer scales (common points) ---")
common_ids_12 <- intersect(buffer_results[["1km"]]$scores$point_id, buffer_results[["2km"]]$scores$point_id)
common_ids_15 <- intersect(buffer_results[["1km"]]$scores$point_id, buffer_results[["5km"]]$scores$point_id)
common_ids_25 <- intersect(buffer_results[["2km"]]$scores$point_id, buffer_results[["5km"]]$scores$point_id)

cor_1v2 <- cor(buffer_results[["1km"]]$scores$PC1[match(common_ids_12, buffer_results[["1km"]]$scores$point_id)],
               buffer_results[["2km"]]$scores$PC1[match(common_ids_12, buffer_results[["2km"]]$scores$point_id)])
cor_1v5 <- cor(buffer_results[["1km"]]$scores$PC1[match(common_ids_15, buffer_results[["1km"]]$scores$point_id)],
               buffer_results[["5km"]]$scores$PC1[match(common_ids_15, buffer_results[["5km"]]$scores$point_id)])
cor_2v5 <- cor(buffer_results[["2km"]]$scores$PC1[match(common_ids_25, buffer_results[["2km"]]$scores$point_id)],
               buffer_results[["5km"]]$scores$PC1[match(common_ids_25, buffer_results[["5km"]]$scores$point_id)])

message("1km vs 2km PC1 correlation: r=", round(cor_1v2,3), " (n=", length(common_ids_12), ")")
message("1km vs 5km PC1 correlation: r=", round(cor_1v5,3), " (n=", length(common_ids_15), ")")
message("2km vs 5km PC1 correlation: r=", round(cor_2v5,3), " (n=", length(common_ids_25), ")")

# Rerun MW test for LP/MP using 1km and 2km PC1 in place of the primary 5km PC1
message("\n--- MW test stability: does LP/MP significance for stream_PC1 hold at 1km/2km? ---")

master_pts <- mw_table %>% distinct() # placeholder, need point-level data — pull from full_matrix_coords
load(file.path(chapter_root, "07_Statistics/mann_whitney_results.RData")) # reload for full_matrix_coords

run_mw_scale_check <- function(scores_df, buffer_label) {
  d <- full_matrix_coords %>% select(point_id, point_type, period) %>%
    inner_join(scores_df %>% select(point_id, PC1), by = "point_id")
  out <- list()
  for (per in c("LP","MP")) {
    test_data <- d %>% filter(point_type=="background" | (point_type=="site" & period==per)) %>%
      mutate(group = ifelse(point_type=="site","site","background")) %>%
      select(group, value=PC1) %>% filter(!is.na(value))
    if (sum(test_data$group=="site") < 5) next
    test_result <- wilcox.test(value ~ group, data = test_data, exact=FALSE)
    med_site <- median(test_data$value[test_data$group=="site"])
    med_bg <- median(test_data$value[test_data$group=="background"])
    out[[per]] <- tibble(buffer = buffer_label, period = per,
                         p_value = round(test_result$p.value,4),
                         direction = ifelse(med_site>med_bg,"higher_at_sites","lower_at_sites"))
  }
  bind_rows(out)
}

scale_stability <- bind_rows(
  run_mw_scale_check(buffer_results[["1km"]]$scores, "1km"),
  run_mw_scale_check(buffer_results[["2km"]]$scores, "2km"),
  run_mw_scale_check(buffer_results[["5km"]]$scores, "5km (primary)")
)
print(scale_stability)


# =============================================================
# PART B: Zero-Replacement Multiplier Sensitivity (5km, mult=0.1/0.5/0.9)
# =============================================================

message("\n\n=== PART B: Zero-replacement sensitivity (5km stream sediment) ===")

mult_results <- list(
  "mult_0.1" = run_stream_clr_pca(5, zero_mult = 0.1),
  "mult_0.5" = run_stream_clr_pca(5, zero_mult = 0.5),   # primary, as used in Script 04
  "mult_0.9" = run_stream_clr_pca(5, zero_mult = 0.9)
)

for (nm in names(mult_results)) {
  r <- mult_results[[nm]]
  message(nm, ": PC1 var=", r$var_pct[1], "% | Top-3 loadings: ",
          paste(names(sort(abs(r$loadings[,1]), decreasing=TRUE))[1:3], collapse=", "))
}

loadings_cor <- tibble(
  comparison = c("mult_0.1 vs mult_0.5", "mult_0.5 vs mult_0.9", "mult_0.1 vs mult_0.9"),
  pc1_loadings_correlation = c(
    cor(mult_results[["mult_0.1"]]$loadings[,1], mult_results[["mult_0.5"]]$loadings[,1]),
    cor(mult_results[["mult_0.5"]]$loadings[,1], mult_results[["mult_0.9"]]$loadings[,1]),
    cor(mult_results[["mult_0.1"]]$loadings[,1], mult_results[["mult_0.9"]]$loadings[,1])
  )
)
message("\nPC1 loadings vector correlation across zero-replacement multipliers:")
print(loadings_cor)


# -------------------------------------------------------------
# SECTION: Combined Verdict
# -------------------------------------------------------------

message("\n=== COMBINED SENSITIVITY VERDICT ===")

buffer_stable <- all(c(cor_1v2, cor_1v5, cor_2v5) > 0.7, na.rm=TRUE)
mult_stable <- all(loadings_cor$pc1_loadings_correlation > 0.9, na.rm=TRUE)

message("Buffer-scale PC1 correlation >0.7 across all pairs: ", buffer_stable)
message("Zero-replacement PC1 loadings correlation >0.9 across all pairs: ", mult_stable)

if (buffer_stable && mult_stable) {
  message("\nRESULT: geochemistry PC1 is STABLE to both buffer-scale choice and zero-")
  message("replacement method. The 5km/half-minimum choices used throughout the chapter")
  message("were pragmatic but not consequential — state this explicitly in Methods 5.3,")
  message("this materially strengthens confidence in the geochemistry results.")
} else {
  message("\nRESULT: geochemistry PC1 shows meaningful sensitivity to ",
          ifelse(!buffer_stable, "buffer scale ", ""),
          ifelse(!mult_stable, "and/or zero-replacement method ", ""),
          "choice. Report this explicitly as a limitation on the geochemistry findings —")
  message("do not present geochem_stream_major_PC1 results as scale/method-independent.")
}


# -------------------------------------------------------------
# SECTION: Save Outputs
# -------------------------------------------------------------

out_path <- file.path(chapter_root, "07_Statistics")
write_csv(scale_stability, file.path(out_path, "geochem_buffer_scale_sensitivity.csv"))
write_csv(loadings_cor, file.path(out_path, "geochem_zero_replacement_sensitivity.csv"))
message("\nSaved: geochem_buffer_scale_sensitivity.csv, geochem_zero_replacement_sensitivity.csv")

save(buffer_results, mult_results, scale_stability, loadings_cor, buffer_stable, mult_stable,
     file = file.path(out_path, "geochem_sensitivity_full.RData"))
message("Saved: geochem_sensitivity_full.RData")

log_file <- file.path(chapter_root, "13_Logs", paste0("script14_log_", Sys.Date(), ".txt"))
sink(log_file)
cat("SCRIPT 14 LOG — Geochemical Sensitivity\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Buffer-scale PC1 correlations: 1v2=", round(cor_1v2,3), " 1v5=", round(cor_1v5,3),
    " 2v5=", round(cor_2v5,3), "\n")
cat("MW stability across scales:\n"); print(scale_stability)
cat("\nZero-replacement loadings correlations:\n"); print(loadings_cor)
sink()
message("Log saved: ", log_file)

message("\n=============================================================")
message("Script 14 complete. Buffer-stable: ", buffer_stable, " | Zero-mult-stable: ", mult_stable)
message("All 4 recommended/mandatory-supporting analyses now done except Task 1")
message("(manual classification verification — human review aid to follow).")
message("=============================================================")