# =============================================================
# SCRIPT 06: Mann-Whitney U Tests + Significance Matrix
# =============================================================
# Project: Reading the Palaeolithic Landscape
# Chapter: Mapping the Past (Springer Nature, 2026)
# Author: Sushant Begade | RTMNU Nagpur
# ORCID: 0009-0003-0804-1763
# Date: August 2026
# =============================================================
# PURPOSE: Enhancement 3 — Period-specific significance matrix.
# Run separate Mann-Whitney U tests per variable per cultural
# period (LP, MP, UP) comparing site vs background locations.
# Generate period × variable significance matrix.
# Identify which environmental dimensions mattered in
# LP vs MP vs UP.
# All variables tested:
# - Sentinel-2A bands + indices
# - Soil properties
# - DEM (elevation)
# - Climate variables
# - Structural geological proximity
# - Burial depth (taphonomy)
# =============================================================


# -------------------------------------------------------------
# SECTION 1: Load Global Parameters + All Data
# -------------------------------------------------------------

load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")

load(file.path(chapter_root,
               "03_Extracted_Values/master_extraction_matrix.RData"))
load(file.path(chapter_root,
               "05_Structural_Geology/vector_extraction_matrix.RData"))
load(file.path(chapter_root,
               "06_Taphonomy/taphonomy_data.RData"))

library(tidyverse)
library(rstatix)
library(coin)
library(pheatmap)
library(ggplot2)
library(patchwork)

set.seed(rf_seed)

message("All data loaded.")
message("Master matrix: ", nrow(master_matrix), " rows × ",
        ncol(master_matrix), " columns")
message("Vector matrix: ", nrow(vector_matrix), " rows × ",
        ncol(vector_matrix), " columns")


# -------------------------------------------------------------
# SECTION 2: Assemble Full Analysis Matrix
# -------------------------------------------------------------

message("\nAssembling full analysis matrix...")

# Core variables from master matrix
core_vars <- master_matrix %>%
  select(
    point_id, point_type, period,
    # DEM
    elevation,
    # Sentinel-2A indices (primary — not raw bands)
    NDVI, NDWI, MNDWI, NDBI, BSI, SAVI, MSAVI,
    # Soil properties
    soil_depth, soil_erosion, soil_productivity,
    soil_slope, soil_texture,
    # CHELSA LGM palaeoclimate
    chelsa_bio01_lgm, chelsa_bio12_lgm, chelsa_bio15_lgm,
    # CHELSA modern
    chelsa_bio01_modern, chelsa_bio12_modern, chelsa_bio15_modern
  )

# Structural geology from vector matrix
struct_vars <- structural_matrix %>%
  select(point_id, dist_fault, dist_dyke,
         dist_lineament, dist_shear, dist_mineral)

# Taphonomy
taph_vars <- taphonomy_matrix %>%
  select(point_id, burial_depth_m)

# Join all
full_matrix <- core_vars %>%
  left_join(struct_vars, by = "point_id") %>%
  left_join(taph_vars,   by = "point_id")

message("Full analysis matrix: ", nrow(full_matrix),
        " rows × ", ncol(full_matrix), " columns")

# Define analysis variables (exclude ID + metadata columns)
analysis_vars <- names(full_matrix)[
  !names(full_matrix) %in%
    c("point_id", "point_type", "period")
]

message("Variables to test: ", length(analysis_vars))
message(paste(analysis_vars, collapse = ", "))


# -------------------------------------------------------------
# SECTION 3: Helper Function — Mann-Whitney U Test
# -------------------------------------------------------------

run_mw_test <- function(data, variable, period_filter) {
  
  # Filter to period sites + all background
  test_data <- data %>%
    filter(point_type == "background" |
             (point_type == "site" & period == period_filter)) %>%
    mutate(group = ifelse(point_type == "site", "site", "background")) %>%
    select(group, value = all_of(variable)) %>%
    filter(!is.na(value))
  
  n_sites <- sum(test_data$group == "site")
  n_bg    <- sum(test_data$group == "background")
  
  # Need minimum observations
  if (n_sites < 3 || n_bg < 3) {
    return(data.frame(
      variable = variable,
      period   = period_filter,
      p_value  = NA,
      w_stat   = NA,
      median_site = NA,
      median_bg   = NA,
      direction   = NA,
      significant = NA
    ))
  }
  
  # Run Wilcoxon rank-sum (Mann-Whitney U)
  test_result <- tryCatch({
    wilcox.test(
      value ~ group,
      data        = test_data,
      exact       = FALSE,
      correct     = TRUE,
      alternative = "two.sided"
    )
  }, error = function(e) NULL)
  
  if (is.null(test_result)) {
    return(data.frame(
      variable    = variable,
      period      = period_filter,
      p_value     = NA,
      w_stat      = NA,
      median_site = NA,
      median_bg   = NA,
      direction   = NA,
      significant = NA
    ))
  }
  
  # Medians
  med_site <- median(
    test_data$value[test_data$group == "site"], na.rm = TRUE
  )
  med_bg <- median(
    test_data$value[test_data$group == "background"], na.rm = TRUE
  )
  
  # Direction
  direction <- ifelse(med_site > med_bg, "higher_at_sites",
                      "lower_at_sites")
  
  data.frame(
    variable    = variable,
    period      = period_filter,
    p_value     = round(test_result$p.value, 4),
    w_stat      = round(test_result$statistic, 0),
    median_site = round(med_site, 4),
    median_bg   = round(med_bg, 4),
    direction   = direction,
    significant = test_result$p.value < alpha
  )
}


# -------------------------------------------------------------
# SECTION 4: Run Mann-Whitney U — All Variables × All Periods
# -------------------------------------------------------------

message("\nRunning Mann-Whitney U tests...")
message("Variables: ", length(analysis_vars))
message("Periods: ", paste(periods, collapse = ", "))
message("Total tests: ", length(analysis_vars) * length(periods))

mw_results <- list()

for (per in periods) {
  message("\nPeriod: ", per)
  for (var in analysis_vars) {
    result <- run_mw_test(full_matrix, var, per)
    mw_results[[paste(per, var, sep = "_")]] <- result
  }
  message("  Done: ", per)
}

# Combine all results
mw_table <- bind_rows(mw_results)

message("\nMann-Whitney U tests complete.")
message("Total tests run: ", nrow(mw_table))
message("Significant results (p < ", alpha, "): ",
        sum(mw_table$significant, na.rm = TRUE))


# -------------------------------------------------------------
# SECTION 5: Apply Bonferroni Correction
# -------------------------------------------------------------

message("\nApplying Bonferroni correction...")

mw_table <- mw_table %>%
  mutate(
    p_bonferroni = p.adjust(p_value, method = "bonferroni"),
    sig_bonferroni = p_bonferroni < alpha
  )

message("Significant after Bonferroni (p < ", alpha, "): ",
        sum(mw_table$sig_bonferroni, na.rm = TRUE))


# -------------------------------------------------------------
# SECTION 6: Build Period × Variable Significance Matrix
# -------------------------------------------------------------

message("\nBuilding significance matrix...")

# Significance matrix — p_value
sig_matrix_p <- mw_table %>%
  select(variable, period, p_value) %>%
  pivot_wider(names_from  = period,
              values_from = p_value) %>%
  column_to_rownames("variable")

# Significance matrix — binary (significant / not)
sig_matrix_bin <- mw_table %>%
  select(variable, period, significant) %>%
  pivot_wider(names_from  = period,
              values_from = significant) %>%
  column_to_rownames("variable")

# Direction matrix
dir_matrix <- mw_table %>%
  select(variable, period, direction) %>%
  pivot_wider(names_from  = period,
              values_from = direction) %>%
  column_to_rownames("variable")

message("Significance matrix dimensions: ",
        nrow(sig_matrix_p), " variables × ",
        ncol(sig_matrix_p), " periods")

# Print full results
message("\n=== SIGNIFICANCE MATRIX (p-values) ===")
print(round(sig_matrix_p, 4))

message("\n=== SIGNIFICANT VARIABLES PER PERIOD ===")
for (per in periods) {
  sig_vars <- mw_table %>%
    filter(period == per, significant == TRUE) %>%
    arrange(p_value) %>%
    select(variable, p_value, direction, median_site, median_bg)
  message("\n", per, " — ", nrow(sig_vars),
          " significant variables (p < ", alpha, "):")
  print(sig_vars)
}


# -------------------------------------------------------------
# SECTION 7: Compute Effect Sizes (Rank-Biserial Correlation)
# -------------------------------------------------------------

message("\nComputing effect sizes...")

compute_rbc <- function(data, variable, period_filter) {
  test_data <- data %>%
    filter(point_type == "background" |
             (point_type == "site" & period == period_filter)) %>%
    mutate(group = ifelse(point_type == "site", 1, 0)) %>%
    select(group, value = all_of(variable)) %>%
    filter(!is.na(value))
  
  n_site <- sum(test_data$group == 1)
  n_bg   <- sum(test_data$group == 0)
  
  if (n_site < 3) return(NA)
  
  w <- wilcox.test(
    value ~ group, data = test_data,
    exact = FALSE
  )$statistic
  
  # Rank-biserial correlation
  rbc <- 1 - (2 * w) / (n_site * n_bg)
  round(as.numeric(rbc), 3)
}

effect_sizes <- expand.grid(
  variable = analysis_vars,
  period   = periods,
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    rbc = compute_rbc(full_matrix, variable, period)
  ) %>%
  ungroup()

# Add to main table
mw_table <- mw_table %>%
  left_join(effect_sizes, by = c("variable", "period"))

message("Effect sizes computed.")


# -------------------------------------------------------------
# SECTION 8: Summary Statistics per Period
# -------------------------------------------------------------

message("\n=== PERIOD SUMMARY ===")

period_summary <- mw_table %>%
  filter(!is.na(significant)) %>%
  group_by(period) %>%
  summarise(
    n_vars_tested   = n(),
    n_significant   = sum(significant, na.rm = TRUE),
    n_sig_bonf      = sum(sig_bonferroni, na.rm = TRUE),
    pct_significant = round(100 * mean(significant, na.rm = TRUE), 1),
    top_variable    = variable[which.min(p_value)],
    lowest_p        = round(min(p_value, na.rm = TRUE), 4),
    .groups = "drop"
  )

print(period_summary)


# -------------------------------------------------------------
# SECTION 9: Save Statistical Outputs
# -------------------------------------------------------------

message("\nSaving outputs...")

out_path <- file.path(chapter_root, "07_Statistics")

# Full MW results table
write_csv(
  mw_table,
  file.path(out_path, "mann_whitney_results_full.csv")
)
message("Saved: mann_whitney_results_full.csv")

# Significance matrix p-values
write_csv(
  sig_matrix_p %>% rownames_to_column("variable"),
  file.path(out_path, "significance_matrix_pvalues.csv")
)
message("Saved: significance_matrix_pvalues.csv")

# Significance matrix binary
write_csv(
  sig_matrix_bin %>% rownames_to_column("variable"),
  file.path(out_path, "significance_matrix_binary.csv")
)
message("Saved: significance_matrix_binary.csv")

# Period summary table
write_csv(
  period_summary,
  file.path(chapter_root, "11_Tables",
            "Table_02_MW_Period_Summary.csv")
)
message("Saved: Table_02_MW_Period_Summary.csv")

# RData
save(
  mw_table,
  sig_matrix_p,
  sig_matrix_bin,
  dir_matrix,
  period_summary,
  full_matrix,
  analysis_vars,
  file = file.path(out_path, "mann_whitney_results.RData")
)
message("Saved: mann_whitney_results.RData")


# -------------------------------------------------------------
# SECTION 10: Figure — Significance Heatmap (Fig. 10)
# -------------------------------------------------------------

message("\nGenerating significance heatmap (Fig. 10)...")

# Prepare -log10(p) matrix for heatmap
heat_data <- sig_matrix_p %>%
  mutate(across(everything(), ~-log10(pmax(., 1e-10))))

# Significance threshold line at -log10(0.05) = 1.301
sig_threshold <- -log10(alpha)

# Clean variable labels
rownames(heat_data) <- gsub("_", " ", rownames(heat_data))
rownames(heat_data) <- gsub("chelsa ", "CHELSA ", rownames(heat_data))
rownames(heat_data) <- gsub("soil ", "Soil ", rownames(heat_data))
rownames(heat_data) <- gsub("dist ", "Dist ", rownames(heat_data))
rownames(heat_data) <- gsub("burial depth m", "Burial Depth", rownames(heat_data))
rownames(heat_data) <- gsub("elevation", "Elevation", rownames(heat_data))

png(
  file.path(chapter_root, "10_Figures",
            "Fig10_Significance_Matrix_Heatmap.png"),
  width  = 120,
  height = fig_width,
  units  = "mm",
  res    = fig_dpi
)

pheatmap(
  as.matrix(heat_data),
  color            = colorRampPalette(
    c("white", "#FFF7BC", "#FD8D3C", "#BD0026"))(100),
  cluster_rows     = TRUE,
  cluster_cols     = FALSE,
  display_numbers  = FALSE,
  border_color     = "grey80",
  fontsize         = 7,
  fontsize_row     = 6,
  fontsize_col     = 9,
  main             = expression(paste(
    "Geoenvironmental Significance Matrix (", -log[10](p), ")"
  )),
  angle_col        = 0,
  legend_breaks    = c(0, 1.301, 2, 3, 4),
  legend_labels    = c("0", "p=0.05", "p=0.01", "p=0.001", "p=0.0001"),
  annotation_col   = data.frame(
    Period = colnames(heat_data),
    row.names = colnames(heat_data)
  )
)

dev.off()
message("Fig10 saved: Significance_Matrix_Heatmap.png")


# -------------------------------------------------------------
# SECTION 11: Figure — Box Plots (Fig. 3)
# -------------------------------------------------------------

message("\nGenerating spectral index box plots (Fig. 3)...")

spectral_vars <- c("NDVI", "NDWI", "MNDWI", "NDBI", "BSI",
                   "SAVI", "MSAVI")

# Prepare long format for plotting
plot_data <- full_matrix %>%
  filter(point_type %in% c("site", "background")) %>%
  mutate(
    group = case_when(
      point_type == "background" ~ "Background",
      period == "LP" ~ "LP Sites",
      period == "MP" ~ "MP Sites",
      period == "UP" ~ "UP Sites"
    ),
    group = factor(group,
                   levels = c("Background", "LP Sites",
                              "MP Sites", "UP Sites"))
  ) %>%
  select(group, all_of(spectral_vars)) %>%
  pivot_longer(
    cols      = all_of(spectral_vars),
    names_to  = "index",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

# Plot
fig3 <- ggplot(plot_data,
               aes(x = group, y = value, fill = group)) +
  geom_boxplot(
    outlier.size  = 0.3,
    outlier.alpha = 0.3,
    lwd           = 0.3
  ) +
  facet_wrap(~index, scales = "free_y", ncol = 4) +
  scale_fill_manual(
    values = c(
      "Background" = "grey70",
      "LP Sites"   = "#2166AC",
      "MP Sites"   = "#F4A582",
      "UP Sites"   = "#D6604D"
    )
  ) +
  labs(
    title    = "Sentinel-2A Spectral Indices at Palaeolithic Site vs Background Locations",
    subtitle = "Lower Palaeolithic (LP), Middle Palaeolithic (MP), Upper Palaeolithic (UP)",
    x        = NULL,
    y        = "Index Value",
    fill     = "Group",
    caption  = "Mann-Whitney U tests; p < 0.05 threshold"
  ) +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1,
                                    size = 6),
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", size = 7),
    legend.position  = "bottom",
    plot.title       = element_text(size = 8, face = "bold"),
    plot.subtitle    = element_text(size = 7)
  )

ggsave(
  file.path(chapter_root, "10_Figures",
            "Fig03_Spectral_Boxplots.png"),
  fig3,
  width  = fig_width / 25.4,
  height = 140 / 25.4,
  dpi    = fig_dpi,
  units  = "in"
)
message("Fig03 saved: Spectral_Boxplots.png")


# -------------------------------------------------------------
# SECTION 12: Figure — Structural Geology Box Plots
# -------------------------------------------------------------

message("\nGenerating structural geology box plots...")

struct_vars_plot <- c("dist_fault", "dist_dyke",
                      "dist_lineament", "dist_shear",
                      "dist_mineral")

struct_plot_data <- full_matrix %>%
  filter(point_type %in% c("site", "background")) %>%
  mutate(
    group = case_when(
      point_type == "background" ~ "Background",
      period == "LP" ~ "LP Sites",
      period == "MP" ~ "MP Sites",
      period == "UP" ~ "UP Sites"
    ),
    group = factor(group,
                   levels = c("Background", "LP Sites",
                              "MP Sites", "UP Sites"))
  ) %>%
  select(group, all_of(struct_vars_plot)) %>%
  mutate(across(all_of(struct_vars_plot), ~. / 1000)) %>%
  pivot_longer(
    cols      = all_of(struct_vars_plot),
    names_to  = "feature",
    values_to = "distance_km"
  ) %>%
  mutate(
    feature = recode(feature,
                     "dist_fault"     = "Fault",
                     "dist_dyke"      = "Dyke",
                     "dist_lineament" = "Lineament",
                     "dist_shear"     = "Shear Zone",
                     "dist_mineral"   = "Mineral Deposit"
    )
  ) %>%
  filter(!is.na(distance_km))

fig5 <- ggplot(struct_plot_data,
               aes(x = group, y = distance_km, fill = group)) +
  geom_boxplot(
    outlier.size  = 0.3,
    outlier.alpha = 0.3,
    lwd           = 0.3
  ) +
  facet_wrap(~feature, scales = "free_y", ncol = 3) +
  scale_fill_manual(
    values = c(
      "Background" = "grey70",
      "LP Sites"   = "#2166AC",
      "MP Sites"   = "#F4A582",
      "UP Sites"   = "#D6604D"
    )
  ) +
  labs(
    title    = "Proximity to Structural Geological Features",
    subtitle = "Distance (km) from site and background locations",
    x        = NULL,
    y        = "Distance (km)",
    fill     = "Group"
  ) +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1,
                                    size = 6),
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", size = 7),
    legend.position  = "bottom",
    plot.title       = element_text(size = 8, face = "bold")
  )

ggsave(
  file.path(chapter_root, "10_Figures",
            "Fig05_Structural_Geology_Boxplots.png"),
  fig5,
  width  = fig_width / 25.4,
  height = 120 / 25.4,
  dpi    = fig_dpi,
  units  = "in"
)
message("Fig05 saved: Structural_Geology_Boxplots.png")


# -------------------------------------------------------------
# SECTION 13: Log Session
# -------------------------------------------------------------

log_file <- file.path(chapter_root, "13_Logs",
                      paste0("script06_log_", Sys.Date(), ".txt"))
sink(log_file)
cat("SCRIPT 06 LOG — Mann-Whitney U Tests\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Variables tested:", length(analysis_vars), "\n")
cat("Periods tested:", paste(periods, collapse = ", "), "\n")
cat("Total tests:", nrow(mw_table), "\n")
cat("Significant (p < 0.05):",
    sum(mw_table$significant, na.rm = TRUE), "\n")
cat("Significant (Bonferroni):",
    sum(mw_table$sig_bonferroni, na.rm = TRUE), "\n\n")
cat("Period summary:\n")
print(period_summary)
cat("\nFull significance matrix (p-values):\n")
print(round(sig_matrix_p, 4))
sink()

message("Log saved: ", log_file)


# -------------------------------------------------------------
# SCRIPT 06 COMPLETE
# -------------------------------------------------------------

message("\n=============================================================")
message("Script 06 complete.")
message("Mann-Whitney U tests: ", nrow(mw_table), " total")
message("Significant (p<0.05): ",
        sum(mw_table$significant, na.rm = TRUE))
message("Significant (Bonferroni): ",
        sum(mw_table$sig_bonferroni, na.rm = TRUE))
message("Figures saved: Fig03, Fig05, Fig10")
message("Next: Run Script 07 — PCA + Random Forest")
message("=============================================================")