# =============================================================
# SCRIPT 00: Environment Setup and Project Architecture
# =============================================================
# Project: Reading the Palaeolithic Landscape
# Chapter: Mapping the Past (Springer Nature, 2026)
# Author: Sushant Begade
# ORCID: 0009-0003-0804-1763
# Institution: RTMNU Nagpur
# Date: August 2026
# GitHub: [your GitHub repo URL]
# Zenodo: [your Zenodo DOI — add after upload]
# =============================================================
# PURPOSE: Install all required packages, establish project
# folder architecture, set global parameters, verify all
# input data files exist, and initialise reproducibility
# infrastructure (session info + seed).
# =============================================================


# -------------------------------------------------------------
# SECTION 1: Set Global Seed for Reproducibility
# -------------------------------------------------------------

set.seed(42)
# All random processes (background point generation,
# Random Forest) use this seed throughout all scripts.


# -------------------------------------------------------------
# SECTION 2: Define Master Project Path
# -------------------------------------------------------------

master_path <- "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset"

# Verify master path exists
if (!dir.exists(master_path)) {
  stop("ERROR: Master data path not found. Check E drive connection.")
} else {
  message("Master data path confirmed: ", master_path)
  
  library(terra)
  # Redirect terra temp to E drive
  dir.create("E:/R_terra_temp", showWarnings = FALSE)
  terraOptions(tempdir = "E:/R_terra_temp")
}


# -------------------------------------------------------------
# SECTION 3: Create Chapter Project Folder Architecture
# -------------------------------------------------------------

chapter_root <- file.path(master_path, "Reading_the_Palaeolithic_Landscape_Chapter")

# Define all subfolders
folders <- c(
  "00_Scripts",           # All R scripts (push to GitHub)
  "01_Data_Processed",    # Cleaned + processed intermediate data
  "02_Background_Points", # Generated pseudo-absence points
  "03_Extracted_Values",  # Variable extraction outputs (site + background)
  "04_Geochemistry",      # Buffer join outputs
  "05_Structural_Geology",# Proximity distance outputs
  "06_Taphonomy",         # Borehole interpolation + burial depth surface
  "07_Statistics",        # Mann-Whitney U results + significance matrices
  "08_PCA",               # PCA outputs + biplots
  "09_RandomForest",      # RF model + variable importance
  "10_Figures",           # All final figures (300 dpi minimum)
  "11_Tables",            # All final tables
  "12_Chapter_Draft",     # Chapter writing outputs
  "13_Logs"               # Session info + run logs
)

# Create all folders
for (f in folders) {
  path <- file.path(chapter_root, f)
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("Created: ", path)
  } else {
    message("Already exists: ", path)
  }
}

message("\nProject folder architecture complete.")


# -------------------------------------------------------------
# SECTION 4: Install and Load Required Packages
# -------------------------------------------------------------

# List of all packages needed across all scripts
required_packages <- c(
  # Spatial data handling
  "terra",          # Raster processing (replaces raster package)
  "sf",             # Vector data (shapefiles)
  "sp",             # Legacy spatial support
  
  # Geospatial analysis
  "gstat",          # IDW interpolation (borehole taphonomy)
  "spatstat",       # Spatial point pattern analysis
  "spdep",          # Spatial autocorrelation
  
  # Statistical analysis
  "coin",           # Exact Mann-Whitney U tests
  "rstatix",        # Tidy statistical tests
  "corrplot",       # Correlation matrices
  
  # Machine learning
  "randomForest",   # Random Forest classification
  "caret",          # Model training + validation framework
  "pROC",           # ROC + AUC evaluation
  
  # Multivariate analysis
  "FactoMineR",     # PCA
  "factoextra",     # PCA visualisation
  "vegan",          # Multivariate ecology methods
  
  # Data manipulation
  "tidyverse",      # dplyr, ggplot2, tidyr, readr, purrr
  "readxl",         # Read Excel files
  "janitor",        # Clean column names
  
  # Visualisation
  "ggplot2",        # Core plotting
  "ggpubr",         # Publication-ready plots
  "patchwork",      # Multi-panel figure assembly
  "viridis",        # Colour-blind-safe palettes
  "RColorBrewer",   # Colour palettes
  "scales",         # Axis formatting
  "pheatmap",       # Heatmaps (significance matrix)
  "cowplot",        # Plot arrangement
  
  # Reproducibility
  "renv",           # Package version locking
  "here"            # Relative path management
)

# Install missing packages only
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    message("Installed: ", pkg)
  } else {
    message("Already installed: ", pkg)
  }
}

message("\nChecking and installing packages...")
invisible(lapply(required_packages, install_if_missing))

# Load all packages
message("\nLoading packages...")
invisible(lapply(required_packages, function(pkg) {
  library(pkg, character.only = TRUE)
  message("Loaded: ", pkg)
}))

message("\nAll packages installed and loaded successfully.")


# -------------------------------------------------------------
# SECTION 5: Define All Input Data Paths
# -------------------------------------------------------------

paths <- list(
  
  # Site data
  sites_csv  = file.path(master_path, "Nagpur Chandrapur Palaeolithic Site Raw Data.csv"),
  sites_xlsx = file.path(master_path, "Nagpur Chandrapur Palaeolithic Site Raw Data.xlsx"),
  
  # District Boundary
  boundary   = file.path(master_path, "District Boundary/Nagpur_and_Chandrapur_Boundry.shp"),
  
  # DEM
  dem        = file.path(master_path, "DEM/DEM.tif"),
  
  # Sentinel-2A Bands
  B01 = file.path(master_path, "Sentinel 2A Multispectral Data/B01_Nagpur_Chandrapur.tif"),
  B02 = file.path(master_path, "Sentinel 2A Multispectral Data/B02_Nagpur_Chandrapur.tif"),
  B03 = file.path(master_path, "Sentinel 2A Multispectral Data/B03_Nagpur_Chandrapur.tif"),
  B04 = file.path(master_path, "Sentinel 2A Multispectral Data/B04_Nagpur_Chandrapur.tif"),
  B05 = file.path(master_path, "Sentinel 2A Multispectral Data/B05_Nagpur_Chandrapur.tif"),
  B06 = file.path(master_path, "Sentinel 2A Multispectral Data/B06_Nagpur_Chandrapur.tif"),
  B07 = file.path(master_path, "Sentinel 2A Multispectral Data/B07_Nagpur_Chandrapur.tif"),
  B08 = file.path(master_path, "Sentinel 2A Multispectral Data/B08_Nagpur_Chandrapur.tif"),
  B09 = file.path(master_path, "Sentinel 2A Multispectral Data/B09_Nagpur_Chandrapur.tif"),
  B11 = file.path(master_path, "Sentinel 2A Multispectral Data/B11_Nagpur_Chandrapur.tif"),
  B12 = file.path(master_path, "Sentinel 2A Multispectral Data/B12_Nagpur_Chandrapur.tif"),
  B8A = file.path(master_path, "Sentinel 2A Multispectral Data/B8A_Nagpur_Chandrapur.tif"),
  
  # Sentinel-2A Derived Indices
  NDVI  = file.path(master_path, "Sentinel 2A Multispectral Data/NDVI.tif"),
  NDWI  = file.path(master_path, "Sentinel 2A Multispectral Data/NDWI.tif"),
  MNDWI = file.path(master_path, "Sentinel 2A Multispectral Data/MNDWI.tif"),
  NDBI  = file.path(master_path, "Sentinel 2A Multispectral Data/NDBI.tif"),
  BSI   = file.path(master_path, "Sentinel 2A Multispectral Data/BSI.tif"),
  SAVI  = file.path(master_path, "Sentinel 2A Multispectral Data/SAVI.tif"),
  MSAVI = file.path(master_path, "Sentinel 2A Multispectral Data/MSAVI.tif"),
  
  # Soil Data
  soil_depth      = file.path(master_path, "Soil Data/Soil Depth.tif"),
  soil_erosion    = file.path(master_path, "Soil Data/Soil Erosion.tif"),
  soil_productivity = file.path(master_path, "Soil Data/Soil Productivity.tif"),
  soil_slope      = file.path(master_path, "Soil Data/Soil Slope.tif"),
  soil_texture    = file.path(master_path, "Soil Data/Soil Texture.tif"),
  
  # Geochemistry
  geoch_horizon  = file.path(master_path, "Geochemistry Data/Soil Horizon/Soil Horizon Data.shp"),
  geoch_regolith = file.path(master_path, "Geochemistry Data/Soil Regolith/Soil Regolith Data.shp"),
  geoch_stream   = file.path(master_path, "Geochemistry Data/Stream Sediment/Stream Sediment Data.shp"),
  
  # Structural Geology
  fault     = file.path(master_path, "Fault/Fault.shp"),
  dyke      = file.path(master_path, "Dyke/Dyke.shp"),
  lineament = file.path(master_path, "Linement/Linement.shp"),   # Note: folder typo in source
  shear     = file.path(master_path, "Shear Zone/Shear Zone.shp"),
  mineral   = file.path(master_path, "Mineral Deposit/Mineral Deposit.shp"),
  
  # Supporting Layers
  rivers      = file.path(master_path, "Rivers/Rivers.shp"),
  waterbody   = file.path(master_path, "Waterbody/Waterbody.shp"),
  geology     = file.path(master_path, "Geology/Geology.shp"),
  lithology   = file.path(master_path, "Lithology/Lithology.shp"),
  geomorphology = file.path(master_path, "Geomorphology/Geomorphology.shp"),
  
  # Borehole
  borehole_shp  = file.path(master_path, "Drilling Borehole Data/Drilling Borehole Data.shp"),
  borehole_xlsx = file.path(master_path, "Drilling Borehole Data/Drilling Borehole Data.xlsx"),
  
  # Climate — CHELSA Modern
  chelsa_bio01_modern = file.path(master_path, "Climate Data/CHELSA_bio01_1981-2010_V.2.1.tif"),
  chelsa_bio12_modern = file.path(master_path, "Climate Data/CHELSA_bio12_1981-2010_V.2.1.tif"),
  chelsa_bio15_modern = file.path(master_path, "Climate Data/CHELSA_bio15_1981-2010_V.2.1.tif"),
  
  # Climate — CHELSA TraCE21k Palaeoclimate (~20ka LGM)
  chelsa_bio01_lgm = file.path(master_path, "Climate Data/CHELSA_TraCE21k_bio01_-200_V.1.0.tif"),
  chelsa_bio12_lgm = file.path(master_path, "Climate Data/CHELSA_TraCE21k_bio12_-200_V.1.0.tif"),
  chelsa_bio15_lgm = file.path(master_path, "Climate Data/CHELSA_TraCE21k_bio15_-200_V.1.0.tif"),
  
  # Climate — WorldClim 2.1
  wc_bio01 = file.path(master_path, "Climate Data/wc2.1_30s_bio_1.tif"),
  wc_bio12 = file.path(master_path, "Climate Data/wc2.1_30s_bio_12.tif"),
  wc_bio15 = file.path(master_path, "Climate Data/wc2.1_30s_bio_15.tif")
)


# -------------------------------------------------------------
# SECTION 6: Verify All Input Files Exist
# -------------------------------------------------------------

message("\nVerifying all input data files...")

missing_files <- c()
for (name in names(paths)) {
  if (!file.exists(paths[[name]])) {
    missing_files <- c(missing_files, paste0(name, ": ", paths[[name]]))
    message("MISSING: ", name, " -> ", paths[[name]])
  } else {
    message("OK: ", name)
  }
}

if (length(missing_files) == 0) {
  message("\nAll input files verified. Proceed to Script 01.")
} else {
  message("\nWARNING: ", length(missing_files), " file(s) missing:")
  for (f in missing_files) message("  - ", f)
  message("Resolve before running subsequent scripts.")
}


# -------------------------------------------------------------
# SECTION 7: Define Global CRS
# -------------------------------------------------------------

# WGS84 geographic CRS for input data
crs_wgs84 <- "EPSG:4326"

# UTM Zone 44N for projected analysis (distance calculations)
# Appropriate for Nagpur-Chandrapur region (73-78E longitude)
crs_utm44n <- "EPSG:32644"

message("\nCRS defined:")
message("  Geographic: ", crs_wgs84)
message("  Projected:  ", crs_utm44n)


# -------------------------------------------------------------
# SECTION 8: Define Global Parameters
# -------------------------------------------------------------

# Background points
n_background <- 1000       # Number of pseudo-absence points
min_dist_bg  <- 2000       # Minimum separation (metres) — avoid spatial autocorrelation

# Buffer distances for geochemical catchment analysis
buffer_distances <- c(1000, 2000, 5000)  # metres

# Cultural period labels
periods <- c("LP", "MP", "UP")           # Lower, Middle, Upper Palaeolithic

# Statistical threshold
alpha <- 0.05              # Significance level for Mann-Whitney U

# Random Forest parameters
rf_ntree <- 1000           # Number of trees
rf_seed  <- 42             # Reproducibility

# Figure export settings
fig_dpi    <- 300          # Minimum per Springer Nature guidelines
fig_width  <- 180          # mm (single column = 90mm, double = 180mm)
fig_height <- 120          # mm default; adjust per figure

message("\nGlobal parameters set.")


# -------------------------------------------------------------
# SECTION 9: Save Paths and Parameters as RData
# -------------------------------------------------------------

# Save to chapter root for loading in all subsequent scripts
save(
  paths,
  master_path,
  chapter_root,
  crs_wgs84,
  crs_utm44n,
  n_background,
  min_dist_bg,
  buffer_distances,
  periods,
  alpha,
  rf_ntree,
  rf_seed,
  fig_dpi,
  fig_width,
  fig_height,
  file = file.path(chapter_root, "00_Scripts/global_params.RData")
)

message("\nGlobal parameters saved to: ",
        file.path(chapter_root, "00_Scripts/global_params.RData"))
message("Load in all subsequent scripts with:")
message('  load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/',
        'Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")')


# -------------------------------------------------------------
# SECTION 10: Initialise renv for Package Version Locking
# -------------------------------------------------------------

# Run once to lock package versions for reproducibility
# renv::init()       # Uncomment on first run
# renv::snapshot()   # Uncomment after all packages confirmed working

message("\nrenv ready. After confirming all packages work:")
message("  Run: renv::init()")
message("  Run: renv::snapshot()")
message("  Push renv.lock to GitHub with all scripts.")


# -------------------------------------------------------------
# SECTION 11: Log Session Information
# -------------------------------------------------------------

log_file <- file.path(chapter_root, "13_Logs",
                      paste0("session_info_", Sys.Date(), ".txt"))

sink(log_file)
cat("=============================================================\n")
cat("SESSION LOG: Reading the Palaeolithic Landscape Chapter\n")
cat("Author: Sushant Begade | RTMNU Nagpur\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Script: 00_Environment_Setup.R\n")
cat("=============================================================\n\n")
print(sessionInfo())
sink()

message("\nSession info logged to: ", log_file)

# -------------------------------------------------------------
# SCRIPT 00 COMPLETE
# -------------------------------------------------------------

message("\n=============================================================")
message("Script 00 complete. Environment ready.")
message("Next: Run Script 01 — Site Data Loading + Validation")
message("=============================================================")