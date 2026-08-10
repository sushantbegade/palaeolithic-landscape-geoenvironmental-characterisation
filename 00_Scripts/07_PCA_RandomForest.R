# =============================================================
# SCRIPT 07: PCA + Random Forest Variable Importance
# =============================================================
# Project: Reading the Palaeolithic Landscape
# Chapter: Mapping the Past (Springer Nature, 2026)
# Author: Sushant Begade | RTMNU Nagpur
# ORCID: 0009-0003-0804-1763
# Date: August 2026
# =============================================================
# PURPOSE: Multivariate integration of all geoenvironmental
# variables through:
# PART A — Principal Component Analysis (PCA)
#   - Site vs background discrimination
#   - LP / MP / UP period centroids in PCA space
#   - Diachronic centroid displacement vectors
#   - Geoenvironmental attractiveness surface
# PART B — Random Forest Classification
#   - Site vs background classification
#   - Variable importance ranking (Mean Decrease Gini)
#   - OOB error + confusion matrix
# =============================================================


# -------------------------------------------------------------
# SECTION 1: Load Global Parameters + All Data
# -------------------------------------------------------------

load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")

load(file.path(chapter_root,
               "07_Statistics/mann_whitney_results.RData"))
load(file.path(chapter_root,
               "06_Taphonomy/taphonomy_data.RData"))
load(file.path(chapter_root,
               "03_Extracted_Values/master_extraction_matrix.RData"))

# Redirect terra temp files to E drive
terraOptions(tempdir = "E:/R_terra_temp")

# Create folder if doesn't exist
if (!dir.exists("E:/R_terra_temp")) {
  dir.create("E:/R_terra_temp", recursive = TRUE)
}

message("Terra temp dir set to E drive: ", terraOptions()$tempdir)

library(terra)
library(sf)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(randomForest)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(viridis)

set.seed(rf_seed)

message("All data loaded.")
message("Full matrix: ", nrow(full_matrix), " rows × ",
        ncol(full_matrix), " columns")
message("Analysis variables: ", length(analysis_vars))


# -------------------------------------------------------------
# SECTION 2: Prepare Clean Matrix for PCA + RF
# -------------------------------------------------------------

message("\nPreparing clean matrix for multivariate analysis...")

# Select analysis variables only
model_data <- full_matrix %>%
  select(point_id, point_type, period,
         all_of(analysis_vars)) %>%
  filter(!is.na(point_type))

# Remove rows with excessive NAs
na_counts <- rowSums(is.na(model_data %>%
                             select(all_of(analysis_vars))))
model_data <- model_data[na_counts <= 3, ]

message("Rows after NA filtering: ", nrow(model_data))

# Impute remaining NAs with column medians
model_imputed <- model_data
for (var in analysis_vars) {
  if (any(is.na(model_imputed[[var]]))) {
    med_val <- median(model_imputed[[var]], na.rm = TRUE)
    model_imputed[[var]][is.na(model_imputed[[var]])] <- med_val
    message("  Imputed NAs in: ", var, " with median = ",
            round(med_val, 3))
  }
}

# Separate sites and background
sites_data  <- model_imputed %>% filter(point_type == "site")
bg_data     <- model_imputed %>% filter(point_type == "background")

message("Sites: ",      nrow(sites_data))
message("Background: ", nrow(bg_data))
message("Complete cases: ", sum(complete.cases(
  model_imputed %>% select(all_of(analysis_vars))
)))


# =============================================================
# PART A: PRINCIPAL COMPONENT ANALYSIS
# =============================================================

message("\n--- PART A: Principal Component Analysis ---")


# -------------------------------------------------------------
# SECTION 3: Run PCA
# -------------------------------------------------------------

message("\nRunning PCA...")

# Scale and center
pca_input <- model_imputed %>%
  select(all_of(analysis_vars)) %>%
  scale(center = TRUE, scale = TRUE) %>%
  as.data.frame()

# Run PCA using FactoMineR
pca_result <- PCA(
  pca_input,
  ncp   = 10,
  graph = FALSE
)

# Variance explained
eig_vals <- get_eigenvalue(pca_result)
message("\nVariance explained by first 5 PCs:")
print(round(eig_vals[1:5, ], 2))

# Variables contributing most to PC1 + PC2
var_contrib <- get_pca_var(pca_result)
message("\nTop 10 variables on PC1:")
print(sort(var_contrib$contrib[, 1], decreasing = TRUE)[1:10])

message("\nTop 10 variables on PC2:")
print(sort(var_contrib$contrib[, 2], decreasing = TRUE)[1:10])


# -------------------------------------------------------------
# SECTION 4: Extract PCA Scores
# -------------------------------------------------------------

message("\nExtracting PCA scores...")

pca_scores <- as.data.frame(pca_result$ind$coord) %>%
  rename_with(~paste0("PC", seq_along(.)))

# Add metadata
pca_scores <- bind_cols(
  model_imputed %>% select(point_id, point_type, period),
  pca_scores
)

# Separate
sites_pca  <- pca_scores %>% filter(point_type == "site")
bg_pca     <- pca_scores %>% filter(point_type == "background")

message("PCA scores extracted for ", nrow(pca_scores), " points.")


# -------------------------------------------------------------
# SECTION 5: Period Centroids in PCA Space
# -------------------------------------------------------------

message("\nCalculating period centroids in PCA space...")

# Centroids per period
period_centroids <- sites_pca %>%
  group_by(period) %>%
  summarise(
    PC1_mean = mean(PC1, na.rm = TRUE),
    PC2_mean = mean(PC2, na.rm = TRUE),
    PC1_sd   = sd(PC1, na.rm = TRUE),
    PC2_sd   = sd(PC2, na.rm = TRUE),
    n        = n(),
    .groups  = "drop"
  )

# Background centroid
bg_centroid <- bg_pca %>%
  summarise(
    period   = "Background",
    PC1_mean = mean(PC1, na.rm = TRUE),
    PC2_mean = mean(PC2, na.rm = TRUE),
    PC1_sd   = sd(PC1, na.rm = TRUE),
    PC2_sd   = sd(PC2, na.rm = TRUE),
    n        = n()
  )

all_centroids <- bind_rows(period_centroids, bg_centroid)

message("\nPeriod centroids in PC1-PC2 space:")
print(all_centroids)

# Centroid displacement vectors — diachronic niche shift
message("\nCentroid displacement vectors:")

# LP to MP shift
lp_centroid <- period_centroids %>% filter(period == "LP")
mp_centroid <- period_centroids %>% filter(period == "MP")
up_centroid <- period_centroids %>% filter(period == "UP")

disp_lp_mp <- sqrt(
  (mp_centroid$PC1_mean - lp_centroid$PC1_mean)^2 +
    (mp_centroid$PC2_mean - lp_centroid$PC2_mean)^2
)
disp_mp_up <- sqrt(
  (up_centroid$PC1_mean - mp_centroid$PC1_mean)^2 +
    (up_centroid$PC2_mean - mp_centroid$PC2_mean)^2
)
disp_lp_up <- sqrt(
  (up_centroid$PC1_mean - lp_centroid$PC1_mean)^2 +
    (up_centroid$PC2_mean - lp_centroid$PC2_mean)^2
)

message("LP → MP displacement: ", round(disp_lp_mp, 3), " PC units")
message("MP → UP displacement: ", round(disp_mp_up, 3), " PC units")
message("LP → UP displacement: ", round(disp_lp_up, 3), " PC units")

displacements <- data.frame(
  transition     = c("LP_to_MP", "MP_to_UP", "LP_to_UP"),
  displacement   = round(c(disp_lp_mp, disp_mp_up, disp_lp_up), 4),
  interpretation = c(
    "LP to MP geoenvironmental niche shift",
    "MP to UP geoenvironmental niche shift",
    "Total LP to UP niche displacement"
  )
)

print(displacements)


# -------------------------------------------------------------
# SECTION 6: Figure — PCA Biplot (Fig. 8)
# -------------------------------------------------------------

message("\nGenerating PCA biplot (Fig. 8)...")

# Variance labels
pc1_var <- round(eig_vals[1, "variance.percent"], 1)
pc2_var <- round(eig_vals[2, "variance.percent"], 1)

# Period colours
period_colors <- c(
  "LP"         = "#2166AC",
  "MP"         = "#F4A582",
  "UP"         = "#D6604D",
  "background" = "grey80"
)

# Plot data
plot_pca <- pca_scores %>%
  mutate(
    group = case_when(
      point_type == "background" ~ "background",
      TRUE                       ~ period
    ),
    alpha_val = ifelse(point_type == "background", 0.15, 0.6),
    size_val  = ifelse(point_type == "background", 0.3, 1.0)
  )

fig8 <- ggplot() +
  # Background points
  geom_point(
    data = plot_pca %>% filter(group == "background"),
    aes(x = PC1, y = PC2),
    color = "grey75", alpha = 0.2, size = 0.3
  ) +
  # Site points by period
  geom_point(
    data = plot_pca %>% filter(group != "background"),
    aes(x = PC1, y = PC2, color = group),
    alpha = 0.7, size = 1.2
  ) +
  # Period centroids
  geom_point(
    data = period_centroids,
    aes(x = PC1_mean, y = PC2_mean, color = period),
    size = 5, shape = 18
  ) +
  # Centroid labels
  geom_label(
    data = period_centroids,
    aes(x = PC1_mean, y = PC2_mean,
        label = period, color = period),
    size = 3, fontface = "bold",
    nudge_y = 0.3, show.legend = FALSE
  ) +
  # Centroid displacement arrows
  annotate(
    "segment",
    x    = lp_centroid$PC1_mean,
    y    = lp_centroid$PC2_mean,
    xend = mp_centroid$PC1_mean,
    yend = mp_centroid$PC2_mean,
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    color = "grey30", lwd = 0.8, linetype = "dashed"
  ) +
  annotate(
    "segment",
    x    = mp_centroid$PC1_mean,
    y    = mp_centroid$PC2_mean,
    xend = up_centroid$PC1_mean,
    yend = up_centroid$PC2_mean,
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    color = "grey30", lwd = 0.8, linetype = "dashed"
  ) +
  # Zero lines
  geom_hline(yintercept = 0, color = "grey50",
             linetype = "dotted", lwd = 0.3) +
  geom_vline(xintercept = 0, color = "grey50",
             linetype = "dotted", lwd = 0.3) +
  scale_color_manual(values = period_colors) +
  labs(
    title    = "PCA Biplot — Geoenvironmental Space",
    subtitle = "Site and background locations with LP/MP/UP period centroids",
    x        = paste0("PC1 (", pc1_var, "% variance)"),
    y        = paste0("PC2 (", pc2_var, "% variance)"),
    color    = "Period",
    caption  = "Dashed arrows = diachronic centroid displacement vectors"
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title    = element_text(size = 9, face = "bold"),
    plot.subtitle = element_text(size = 8),
    legend.position = "right"
  )

ggsave(
  file.path(chapter_root, "10_Figures",
            "Fig08_PCA_Biplot.png"),
  fig8,
  width  = fig_width / 25.4,
  height = 130 / 25.4,
  dpi    = fig_dpi,
  units  = "in"
)
message("Fig08 saved: PCA_Biplot.png")


# =============================================================
# PART B: RANDOM FOREST CLASSIFICATION
# =============================================================

message("\n--- PART B: Random Forest Classification ---")


# -------------------------------------------------------------
# SECTION 7: Prepare RF Data
# -------------------------------------------------------------

message("\nPreparing Random Forest data...")

rf_data <- model_imputed %>%
  select(point_type, all_of(analysis_vars)) %>%
  mutate(
    response = factor(
      ifelse(point_type == "site", "site", "background"),
      levels = c("site", "background")
    )
  ) %>%
  select(-point_type)

# Balance classes — downsample background to match site count
n_sites <- sum(rf_data$response == "site")
n_bg    <- sum(rf_data$response == "background")

message("Sites: ", n_sites, " | Background: ", n_bg)

# Downsample background
set.seed(rf_seed)
bg_idx   <- which(rf_data$response == "background")
site_idx <- which(rf_data$response == "site")
bg_sampled <- sample(bg_idx, size = n_sites, replace = FALSE)

rf_balanced <- rf_data[c(site_idx, bg_sampled), ]

message("Balanced RF dataset: ", nrow(rf_balanced),
        " rows (", n_sites, " sites + ", n_sites, " background)")


# -------------------------------------------------------------
# SECTION 8: Run Random Forest
# -------------------------------------------------------------

message("\nRunning Random Forest (", rf_ntree, " trees)...")

rf_model <- randomForest(
  response ~ .,
  data       = rf_balanced,
  ntree      = rf_ntree,
  mtry       = floor(sqrt(length(analysis_vars))),
  importance = TRUE,
  keep.forest = TRUE
)

message("\nRandom Forest complete.")
message("OOB Error Rate: ",
        round(rf_model$err.rate[rf_ntree, "OOB"] * 100, 2), "%")
message("\nConfusion Matrix:")
print(rf_model$confusion)


# -------------------------------------------------------------
# SECTION 9: Variable Importance
# -------------------------------------------------------------

message("\nExtracting variable importance...")

# Mean Decrease Gini
importance_df <- as.data.frame(importance(rf_model)) %>%
  rownames_to_column("variable") %>%
  arrange(desc(MeanDecreaseGini)) %>%
  mutate(
    rank          = row_number(),
    variable_clean = gsub("_", " ", variable),
    variable_clean = gsub("chelsa ", "CHELSA ", variable_clean),
    variable_clean = gsub("soil ", "Soil ", variable_clean),
    variable_clean = gsub("dist ", "Dist ", variable_clean),
    variable_clean = gsub("burial depth m", "Burial Depth", variable_clean),
    variable_clean = gsub("elevation", "Elevation", variable_clean),
    category = case_when(
      str_detect(variable, "NDVI|NDWI|MNDWI|NDBI|BSI|SAVI|MSAVI|S2_") ~
        "Spectral",
      str_detect(variable, "soil_") ~ "Pedogenic",
      str_detect(variable, "chelsa|wc_") ~ "Climate",
      str_detect(variable, "dist_") ~ "Structural Geology",
      str_detect(variable, "elevation") ~ "Topographic",
      str_detect(variable, "burial") ~ "Taphonomy",
      TRUE ~ "Other"
    )
  )

message("\nTop 15 variables by Mean Decrease Gini:")
print(importance_df %>%
        select(rank, variable, MeanDecreaseGini, category) %>%
        head(15))


# -------------------------------------------------------------
# SECTION 10: Figure — RF Variable Importance (Fig. 9)
# -------------------------------------------------------------

message("\nGenerating RF variable importance plot (Fig. 9)...")

category_colors <- c(
  "Spectral"         = "#1F78B4",
  "Pedogenic"        = "#33A02C",
  "Climate"          = "#FF7F00",
  "Structural Geology" = "#E31A1C",
  "Topographic"      = "#6A3D9A",
  "Taphonomy"        = "#B15928"
)

fig9 <- importance_df %>%
  head(20) %>%
  mutate(variable_clean = fct_reorder(variable_clean,
                                      MeanDecreaseGini)) %>%
  ggplot(aes(x = MeanDecreaseGini,
             y = variable_clean,
             fill = category)) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = paste0("#", rank)),
    hjust = -0.2, size = 2.5
  ) +
  scale_fill_manual(values = category_colors) +
  labs(
    title    = "Random Forest Variable Importance",
    subtitle = paste0("Top 20 variables | OOB accuracy: ",
                      round((1 - rf_model$err.rate[rf_ntree, "OOB"])
                            * 100, 1), "%"),
    x        = "Mean Decrease Gini",
    y        = NULL,
    fill     = "Variable Category"
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title      = element_text(size = 9, face = "bold"),
    plot.subtitle   = element_text(size = 8),
    legend.position = "bottom",
    legend.text     = element_text(size = 7),
    axis.text.y     = element_text(size = 7)
  ) +
  expand_limits(x = max(importance_df$MeanDecreaseGini[1:20]) * 1.1)

ggsave(
  file.path(chapter_root, "10_Figures",
            "Fig09_RF_Variable_Importance.png"),
  fig9,
  width  = fig_width / 25.4,
  height = 130 / 25.4,
  dpi    = fig_dpi,
  units  = "in"
)
message("Fig09 saved: RF_Variable_Importance.png")


# -------------------------------------------------------------
# SECTION 11: Geoenvironmental Attractiveness Surface (Fig. 11)
# -------------------------------------------------------------

message("\nGenerating geoenvironmental attractiveness surface (Fig. 11)...")

# Load all rasters used in analysis
message("Loading rasters for prediction surface...")

raster_vars <- list(
  elevation         = paths$dem,
  NDVI              = paths$NDVI,
  NDWI              = paths$NDWI,
  MNDWI             = paths$MNDWI,
  NDBI              = paths$NDBI,
  BSI               = paths$BSI,
  SAVI              = paths$SAVI,
  MSAVI             = paths$MSAVI,
  soil_depth        = paths$soil_depth,
  soil_erosion      = paths$soil_erosion,
  soil_productivity = paths$soil_productivity,
  soil_slope        = paths$soil_slope,
  soil_texture      = paths$soil_texture,
  chelsa_bio01_lgm  = paths$chelsa_bio01_lgm,
  chelsa_bio12_lgm  = paths$chelsa_bio12_lgm,
  chelsa_bio15_lgm  = paths$chelsa_bio15_lgm,
  chelsa_bio01_modern = paths$chelsa_bio01_modern,
  chelsa_bio12_modern = paths$chelsa_bio12_modern,
  chelsa_bio15_modern = paths$chelsa_bio15_modern
)

# Load + reproject + crop each raster
boundary <- st_read(paths$boundary, quiet = TRUE) %>%
  st_transform(crs_utm44n)
study_area_vect <- vect(st_union(boundary))

# Reference raster — DEM at 250m for tractable prediction
dem_ref <- rast(paths$dem) %>%
  terra::project("EPSG:32644") %>%
  crop(study_area_vect) %>%
  mask(study_area_vect) %>%
  aggregate(fact = 8)  # ~250m resolution for prediction

message("Reference raster resolution: ", res(dem_ref)[1], "m")

# Build raster stack
raster_stack <- list()
for (name in names(raster_vars)) {
  tryCatch({
    r <- rast(raster_vars[[name]]) %>%
      terra::project("EPSG:32644") %>%
      crop(study_area_vect) %>%
      mask(study_area_vect) %>%
      resample(dem_ref, method = "bilinear")
    raster_stack[[name]] <- r
    message("  Loaded: ", name)
  }, error = function(e) {
    message("  SKIP: ", name, " — ", e$message)
  })
}

# Convert stack to dataframe for prediction
stack_rast <- rast(raster_stack)
names(stack_rast) <- names(raster_stack)

pred_df <- as.data.frame(stack_rast, na.rm = FALSE)
message("Prediction grid: ", nrow(pred_df), " pixels")

# Add missing structural variables with median values
# (structural geology not available as continuous raster)
for (sv in c("dist_fault", "dist_dyke", "dist_lineament",
             "dist_shear", "dist_mineral", "burial_depth_m")) {
  if (sv %in% analysis_vars) {
    med_val <- median(model_imputed[[sv]], na.rm = TRUE)
    pred_df[[sv]] <- med_val
  }
}

# Impute NAs in prediction grid
for (v in names(pred_df)) {
  if (any(is.na(pred_df[[v]]))) {
    pred_df[[v]][is.na(pred_df[[v]])] <-
      median(pred_df[[v]], na.rm = TRUE)
  }
}

# Keep only model variables
pred_vars_available <- intersect(analysis_vars, names(pred_df))
missing_vars <- setdiff(analysis_vars, names(pred_df))
if (length(missing_vars) > 0) {
  message("Missing vars for prediction (using median): ",
          paste(missing_vars, collapse = ", "))
  for (mv in missing_vars) {
    pred_df[[mv]] <- median(model_imputed[[mv]], na.rm = TRUE)
  }
}

pred_input <- pred_df %>% select(all_of(analysis_vars))

# Predict probability of "site"
message("Predicting geoenvironmental attractiveness...")
rf_probs <- predict(rf_model, pred_input, type = "prob")
site_probs <- rf_probs[, "site"]

# Map back to raster
attract_vals <- rep(NA, nrow(as.data.frame(dem_ref, na.rm = FALSE)))
non_na_idx <- which(!is.na(as.data.frame(dem_ref,
                                         na.rm = FALSE)[, 1]))

# Align predictions with non-NA pixels
if (length(site_probs) == length(non_na_idx)) {
  attract_vals[non_na_idx] <- site_probs
} else {
  attract_vals <- site_probs
}

attract_rast <- dem_ref
values(attract_rast) <- site_probs

attract_rast <- mask(attract_rast, study_area_vect)

message("Attractiveness surface range: ",
        round(minmax(attract_rast)[1], 3), " to ",
        round(minmax(attract_rast)[2], 3))

# Save raster
writeRaster(
  attract_rast,
  file.path(chapter_root, "08_PCA",
            "geoenvironmental_attractiveness.tif"),
  overwrite = TRUE
)
message("Attractiveness surface saved.")

# -------------------------------------------------------------
# Fig 11 — Simple base plot, no sf dependency
# -------------------------------------------------------------

message("Generating Fig11...")

# Site coordinates directly from model_imputed
sites_coords <- full_matrix %>%
  filter(point_type == "site") %>%
  select(point_id, period) %>%
  left_join(
    master_matrix %>% select(point_id, easting, northing),
    by = "point_id"
  ) %>%
  filter(!is.na(easting), !is.na(northing))

message("Site coords for plot: ", nrow(sites_coords))

png(
  file.path(chapter_root, "10_Figures",
            "Fig11_Geoenvironmental_Attractiveness.png"),
  width  = fig_width,
  height = fig_height,
  units  = "mm",
  res    = fig_dpi
)

plot(attract_rast,
     main = "Geoenvironmental Attractiveness Surface",
     col  = viridis(100),
     axes = TRUE)

points(
  x   = sites_coords$easting,
  y   = sites_coords$northing,
  pch = 16,
  cex = 0.5,
  col = "red"
)

dev.off()
message("Fig11 saved: Geoenvironmental_Attractiveness.png")


# -------------------------------------------------------------
# SECTION 12: Save All Outputs
# -------------------------------------------------------------

message("\nSaving outputs...")

# PCA outputs
write_csv(
  pca_scores,
  file.path(chapter_root, "08_PCA", "pca_scores.csv")
)
write_csv(
  all_centroids,
  file.path(chapter_root, "08_PCA", "period_centroids.csv")
)
write_csv(
  displacements,
  file.path(chapter_root, "08_PCA",
            "centroid_displacement_vectors.csv")
)
write_csv(
  as.data.frame(eig_vals) %>% rownames_to_column("PC"),
  file.path(chapter_root, "08_PCA", "pca_eigenvalues.csv")
)
message("Saved: PCA outputs")

# RF outputs
write_csv(
  importance_df,
  file.path(chapter_root, "09_RandomForest",
            "rf_variable_importance.csv")
)
write_csv(
  importance_df %>% head(15),
  file.path(chapter_root, "11_Tables",
            "Table_03_RF_Variable_Importance.csv")
)
message("Saved: RF variable importance")

# Displacement table
write_csv(
  displacements,
  file.path(chapter_root, "11_Tables",
            "Table_04_Centroid_Displacements.csv")
)
message("Saved: Table_04_Centroid_Displacements.csv")

# RData
save(
  pca_result,
  pca_scores,
  period_centroids,
  all_centroids,
  displacements,
  rf_model,
  importance_df,
  attract_rast,
  file = file.path(chapter_root, "09_RandomForest",
                   "pca_rf_results.RData")
)
message("Saved: pca_rf_results.RData")


# -------------------------------------------------------------
# SECTION 13: Final Summary
# -------------------------------------------------------------

message("\n=== FINAL ANALYTICAL SUMMARY ===")

message("\nPCA:")
message("  PC1 variance: ", round(eig_vals[1, "variance.percent"], 1), "%")
message("  PC2 variance: ", round(eig_vals[2, "variance.percent"], 1), "%")
message("  PC1+PC2 total: ",
        round(sum(eig_vals[1:2, "variance.percent"]), 1), "%")

message("\nPeriod centroids (PC1, PC2):")
for (i in 1:nrow(period_centroids)) {
  message("  ", period_centroids$period[i], ": (",
          round(period_centroids$PC1_mean[i], 3), ", ",
          round(period_centroids$PC2_mean[i], 3), ")")
}

message("\nDiachronic displacements:")
message("  LP → MP: ", round(disp_lp_mp, 3), " PC units")
message("  MP → UP: ", round(disp_mp_up, 3), " PC units")

message("\nRandom Forest:")
message("  OOB accuracy: ",
        round((1 - rf_model$err.rate[rf_ntree, "OOB"]) * 100, 1), "%")
message("  Top variable: ", importance_df$variable[1])
message("  Top category: ", importance_df$category[1])


# -------------------------------------------------------------
# SECTION 14: Log Session
# -------------------------------------------------------------

log_file <- file.path(chapter_root, "13_Logs",
                      paste0("script07_log_", Sys.Date(), ".txt"))
sink(log_file)
cat("SCRIPT 07 LOG — PCA + Random Forest\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("PCA — PC1 variance:", round(eig_vals[1, "variance.percent"], 1), "%\n")
cat("PCA — PC2 variance:", round(eig_vals[2, "variance.percent"], 1), "%\n")
cat("\nPeriod centroids:\n")
print(all_centroids)
cat("\nDisplacement vectors:\n")
print(displacements)
cat("\nRF OOB accuracy:",
    round((1 - rf_model$err.rate[rf_ntree, "OOB"]) * 100, 1), "%\n")
cat("\nTop 15 variable importance:\n")
print(importance_df %>% head(15) %>%
        select(rank, variable, MeanDecreaseGini, category))
sink()

message("Log saved: ", log_file)


# -------------------------------------------------------------
# SCRIPT 07 COMPLETE
# -------------------------------------------------------------

message("\n=============================================================")
message("Script 07 complete. All analyses done.")
message("PCA biplot:           Fig08")
message("RF importance:        Fig09")
message("Attractiveness surface: Fig11")
message("Next: Run Script 08 — Final Figures Production")
message("=============================================================")