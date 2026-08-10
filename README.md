# Reading the Palaeolithic Landscape: Multi-Source Geoenvironmental Characterisation

[![DOI](https://zenodo.org/badge/DOI.svg)](https://doi.org/10.5281/zenodo.21866496)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-%3E%3D4.5.0-blue)](https://www.r-project.org/)

## Overview

This repository contains all R scripts used in the geoenvironmental analysis presented in:

> Begade, S. (2026). Reading the Palaeolithic Landscape: A Multi-Source
> Geoenvironmental Approach to Site Characterisation in the Deccan Trap
> Region of Central India Using GIS and Remote Sensing. In Dhanaraj K.,
> Das S., & Anand S. (Eds.), *Mapping the Past: GIS and Remote Sensing
> in Landscape Archaeology*. Springer Nature.

The study presents the **first multi-source geoenvironmental characterisation**
of 197 georeferenced Palaeolithic sites spanning the Lower, Middle, and Upper
Palaeolithic in the Nagpur and Chandrapur Districts of the Wainganga-Wardha
Basin, Vidarbha, Central India.

---

## Study Area

- **Districts:** Nagpur and Chandrapur, Maharashtra, India
- **Basin:** Wainganga-Wardha, Vidarbha
- **Geological context:** Deccan Trap volcanic province
- **Sites:** 197 georeferenced Palaeolithic sites (LP=75, MP=95, UP=27)
- **Coordinate system:** WGS84 / UTM Zone 44N (EPSG:32644)

---

## Repository Structure

```
├── 00_Environment_Setup.R # Package installation + project folder architecture
├── 01_Site_Data_Loading_Validation.R # Site data cleaning, validation, period classification
├── 02_Background_Point_Generation.R # Pseudo-absence background point generation (n=1000)
├── 03_Raster_Variable_Extraction.R # Raster value extraction at site + background points
├── 04_Vector_Variable_Extraction.R # Structural geology proximity + geochemistry buffer joins
├── 05_Taphonomic_Bias_Correction.R # Borehole IDW interpolation + burial depth surface
├── 06_Mann_Whitney_Tests.R # Period-specific significance matrix (75 tests)
├── 07_PCA_RandomForest.R # PCA biplot + Random Forest variable importance
├── 08_Publication_Figures.R # All 12 publication-ready figures (300 dpi)
├── global_params.RData # Global parameters + file paths
├── renv.lock # Package version lockfile (reproducibility)
└── README.md # This file
```
---

## Data Sources

All input datasets are archived on Zenodo (DOI: https://doi.org/10.5281/zenodo.21866544).

| Layer | Source | Type | Resolution |
|---|---|---|---|
| 197 Palaeolithic sites | Field survey + literature | Vector points | GPS ±5–10m |
| Sentinel-2A bands B01–B12 | ESA Copernicus | Raster | 10–60m |
| Derived indices (NDVI, NDWI, MNDWI, NDBI, BSI, SAVI, MSAVI) | Computed | Raster | 10m |
| DEM | SRTM / ALOS | Raster | 30m |
| Stream sediment geochemistry | GSI | Vector points | — |
| Soil horizon geochemistry | GSI | Vector points | — |
| Soil regolith geochemistry | GSI | Vector points | — |
| Soil depth, erosion, productivity, slope, texture | NBSS&LUP | Raster | 250m |
| Fault, dyke, lineament, shear zone | GSI | Vector | 1:250,000 |
| Mineral deposit | GSI | Vector | 1:250,000 |
| Drilling borehole data | CGWB | Vector points | — |
| CHELSA TraCE21k bio01, bio12, bio15 (~20ka LGM) | CHELSA V1.0 | Raster | 1km |
| CHELSA modern bio01, bio12, bio15 (1981–2010) | CHELSA V2.1 | Raster | 1km |
| WorldClim 2.1 bio01, bio12, bio15 | WorldClim | Raster | 30s |
| District boundary, geology, lithology, geomorphology, rivers, waterbody | GSI / Survey of India | Vector | 1:250,000 |
---

## Analytical Workflow
Script 00 → Environment + project folder setup
Script 01 → Site data loading + validation + period classification
Script 02 → Background point generation (n=1,000, min 2km separation)
Script 03 → Raster variable extraction (39 variables × 1,197 points)
Script 04 → Structural geology proximity + geochemistry buffer joins
Script 05 → Taphonomic bias correction (IDW borehole interpolation)
Script 06 → Mann-Whitney U tests (25 variables × 3 periods = 75 tests)
+ Bonferroni correction + effect sizes
Script 07 → PCA (25 variables) + Random Forest (1,000 trees)
+ Geoenvironmental attractiveness surface
Script 08 → 12 publication-ready figures (300 dpi, Springer Nature format)
---

## Analytical Methods

### Background Sampling
- n=1,000 pseudo-absence background points
- Spatial filtering: minimum 2km separation (seed=42)
- Constrained within district boundary

### Mann-Whitney U Tests
- Non-parametric, two-sided
- Run separately per variable per cultural period (LP, MP, UP)
- Bonferroni correction applied
- Rank-biserial correlation effect sizes computed
- Significance threshold: α=0.05

### Principal Component Analysis
- 25 standardised geoenvironmental variables
- FactoMineR implementation
- PC1=29.3%, PC2=21.5% variance explained
- Period centroids + diachronic displacement vectors

### Random Forest
- 1,000 trees, mtry=5
- Balanced classes (n=197 sites + 197 background)
- OOB accuracy: 66.8%
- Variable importance: Mean Decrease Gini

### Taphonomic Bias Correction
- IDW interpolation (power=2) of borehole depth data (n=93)
- Burial risk classification: Low (<36.9m) / Medium (36.9–135.5m) / High (>135.5m)

---

## Key Results

| Finding | Detail |
|---|---|
| LP → MP niche displacement | 3.62 PC units |
| MP → UP niche displacement | 3.58 PC units |
| LP → UP proximity | 1.32 PC units |
| Top RF variable | Distance to Shear Zone |
| Top RF category | Structural Geology (4 of top 5) |
| LP significant variables | 19/25 (76%) |
| MP significant variables | 18/25 (72%) |
| UP significant variables | 6/25 (24%) |
| LP high burial risk | 65% of sites |
| MP high burial risk | 6% of sites |

---

## Requirements

**R version:** ≥ 4.5.0

**Key packages:**
```r
terra, sf, tidyverse, FactoMineR, factoextra,
randomForest, caret, gstat, pheatmap, ggplot2,
patchwork, ggrepel, ggnewscale, viridis, renv
```

**Restore exact package versions:**
```r
renv::restore()
```

---

## Usage

### 1. Clone repository
```bash
git clone https://github.com/YOUR_USERNAME/palaeolithic-landscape-geoenvironmental-characterisation.git
```

### 2. Download data from Zenodo
Download and extract the data archive from:
[https://doi.org/PASTE_ZENODO_DOI_HERE](https://doi.org/10.5281/zenodo.21866544)

### 3. Update file paths
In `00_Environment_Setup.R`, update:
```r
master_path <- "YOUR_LOCAL_PATH/Nagpur-Chandrapur Enhanced Dataset"
chapter_root <- "YOUR_LOCAL_PATH/Reading_the_Palaeolithic_Landscape_Chapter"
```

### 4. Restore packages
```r
renv::restore()
```

### 5. Run scripts in order
```r
source("00_Environment_Setup.R")
source("01_Site_Data_Loading_Validation.R")
source("02_Background_Point_Generation.R")
source("03_Raster_Variable_Extraction.R")
source("04_Vector_Variable_Extraction.R")
source("05_Taphonomic_Bias_Correction.R")
source("06_Mann_Whitney_Tests.R")
source("07_PCA_RandomForest.R")
source("08_Publication_Figures.R")
```

---

## Reproducibility

- All random processes use `set.seed(42)`
- Package versions locked via `renv.lock`
- Intermediate outputs saved as `.RData` files between scripts
- Full session logs saved to `13_Logs/`

---

## Citation

If you use these scripts, please cite:

**Chapter:**
> Begade, S. (2026). Reading the Palaeolithic Landscape: A Multi-Source
> Geoenvironmental Approach to Site Characterisation in the Deccan Trap
> Region of Central India Using GIS and Remote Sensing. In Dhanaraj K.,
> Das S., & Anand S. (Eds.), *Mapping the Past: GIS and Remote Sensing
> in Landscape Archaeology*. Springer Nature.

**Code:**
> Begade, S. (2026). R scripts for multi-source geoenvironmental
> characterisation of Palaeolithic sites, Wainganga-Wardha Basin,
> Central India [Software]. GitHub.
> https://github.com/YOUR_USERNAME/palaeolithic-landscape-geoenvironmental-characterisation
> DOI: https://doi.org/10.5281/zenodo.21866496

**Data:**
> Begade, S. (2026). Multi-source geoenvironmental dataset for
> Palaeolithic landscape characterisation, Nagpur-Chandrapur Districts,
> Central India [Dataset]. Zenodo.
> DOI: https://doi.org/10.5281/zenodo.21866544

---

## Author

**Sushant Begade**
Ph.D. Research Scholar
Department of Ancient Indian History, Culture and Archaeology
Rashtrasant Tukadoji Maharaj Nagpur University
Nagpur, Maharashtra, India
ORCID: [0009-0003-0804-1763](https://orcid.org/0009-0003-0804-1763)

---

## License

MIT License — see [LICENSE](LICENSE) file for details.

---

## Acknowledgements

Data sources: GSI (Geological Survey of India), CGWB (Central Ground Water Board), NBSS&LUP (National Bureau of Soil Survey and Land Use Planning), ESA Copernicus (Sentinel-2A), CHELSA (Karger et al.), WorldClim 2.1, SRTM/ALOS.