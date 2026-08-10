# =============================================================
# SCRIPT 08: Publication-Ready Statistical Figures
# =============================================================
# Project: Reading the Palaeolithic Landscape
# Chapter: Mapping the Past (Springer Nature, 2026)
# Author: Sushant Begade | RTMNU Nagpur
# ORCID: 0009-0003-0804-1763
# Date: August 2026
# =============================================================
# FIGURES PRODUCED:
# Fig03 — Spectral index box plots
# Fig05 — Structural geology proximity box plots
# Fig07 — Soil properties box plots
# Fig08 — PCA biplot with centroids + displacement vectors
# Fig09 — Random Forest variable importance
# Fig10 — Significance matrix heatmap
# =============================================================
# RASTER FIGURES (Fig01,02,04,06,11,12) — produce in QGIS
# =============================================================


# -------------------------------------------------------------
# SECTION 1: Setup
# -------------------------------------------------------------

load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")

load(file.path(chapter_root,
               "03_Extracted_Values/master_extraction_matrix.RData"))
load(file.path(chapter_root,
               "05_Structural_Geology/vector_extraction_matrix.RData"))
load(file.path(chapter_root,
               "07_Statistics/mann_whitney_results.RData"))
load(file.path(chapter_root,
               "09_RandomForest/pca_rf_results.RData"))

library(tidyverse)
library(ggplot2)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(scales)
library(pheatmap)
library(ggrepel)
library(FactoMineR)
library(factoextra)

set.seed(rf_seed)

fig_path <- file.path(chapter_root, "10_Figures")

# -------------------------------------------------------------
# GLOBAL SPRINGER NATURE THEME
# -------------------------------------------------------------

theme_springer <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title        = element_text(size = base_size + 1,
                                       face = "bold", hjust = 0,
                                       margin = margin(0,0,2,0,"mm")),
      plot.subtitle     = element_text(size = base_size - 0.5,
                                       hjust = 0, color = "grey40",
                                       margin = margin(0,0,3,0,"mm")),
      plot.caption      = element_text(size = base_size - 1.5,
                                       color = "grey55", hjust = 0,
                                       margin = margin(2,0,0,0,"mm")),
      axis.title        = element_text(size = base_size,
                                       color = "grey20"),
      axis.text         = element_text(size = base_size - 1,
                                       color = "grey30"),
      axis.line         = element_line(color = "grey30",
                                       linewidth = 0.35),
      axis.ticks        = element_line(color = "grey50",
                                       linewidth = 0.25),
      panel.grid.major  = element_line(color = "grey93",
                                       linewidth = 0.2),
      panel.grid.minor  = element_blank(),
      panel.background  = element_rect(fill = "white",
                                       color = NA),
      panel.border      = element_rect(color = "grey30",
                                       fill = NA,
                                       linewidth = 0.35),
      strip.background  = element_rect(fill = "#F5F5F5",
                                       color = "grey30",
                                       linewidth = 0.3),
      strip.text        = element_text(size = base_size,
                                       face = "bold",
                                       color = "grey20"),
      legend.title      = element_text(size = base_size,
                                       face = "bold"),
      legend.text       = element_text(size = base_size - 1),
      legend.key.size   = unit(3.5, "mm"),
      legend.key        = element_rect(fill = "white",
                                       color = NA),
      legend.background = element_rect(fill = "white",
                                       color = "grey80",
                                       linewidth = 0.25),
      plot.background   = element_rect(fill = "white",
                                       color = NA),
      plot.margin       = margin(5, 5, 4, 5, "mm")
    )
}

# Consistent period palette
period_pal <- c(
  "Background" = "#BDBDBD",
  "LP"         = "#2166AC",
  "MP"         = "#D6604D",
  "UP"         = "#4DAC26"
)

period_labels <- c(
  "Background" = "Background (n=1,000)",
  "LP"         = "LP (n=75)",
  "MP"         = "MP (n=95)",
  "UP"         = "UP (n=27)"
)

message("Setup complete. Building figures...")


# =============================================================
# FIGURE 3: Spectral Index Box Plots
# =============================================================

message("\n--- Fig03: Spectral Index Box Plots ---")

spectral_vars <- c("NDVI","NDWI","MNDWI",
                   "NDBI","BSI","SAVI","MSAVI")

# Significance stars from MW results
sig_stars <- mw_table %>%
  filter(variable %in% spectral_vars) %>%
  mutate(
    sig_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  ) %>%
  select(variable, period, p_value, sig_label, direction)

# Long format
spec_long <- master_matrix %>%
  filter(point_type %in% c("site","background")) %>%
  mutate(
    group = case_when(
      point_type == "background" ~ "Background",
      TRUE ~ period
    ),
    group = factor(group,
                   levels = c("Background","LP","MP","UP"))
  ) %>%
  select(group, all_of(spectral_vars)) %>%
  pivot_longer(cols = all_of(spectral_vars),
               names_to = "index",
               values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(index = factor(index, levels = spectral_vars))

# Per-facet significance annotation
# Compute y positions for stars
star_pos <- spec_long %>%
  group_by(index) %>%
  summarise(y_max = quantile(value, 0.98, na.rm=TRUE) * 1.05,
            .groups = "drop")

# Build figure
fig03 <- ggplot(spec_long,
                aes(x = group, y = value, fill = group)) +
  geom_boxplot(
    outlier.shape  = 19,
    outlier.size   = 0.4,
    outlier.alpha  = 0.25,
    outlier.color  = "grey60",
    linewidth      = 0.3,
    width          = 0.6,
    fatten         = 2
  ) +
  geom_hline(yintercept = 0,
             color = "grey40", linewidth = 0.2,
             linetype = "longdash") +
  facet_wrap(~index, scales = "free_y",
             ncol = 4, nrow = 2) +
  scale_fill_manual(values = period_pal,
                    labels = period_labels,
                    name   = NULL) +
  scale_x_discrete(labels = c(
    "Background" = "BG",
    "LP" = "LP", "MP" = "MP", "UP" = "UP"
  )) +
  labs(
    title   = "Fig. 3  |  Sentinel-2A spectral index values at Palaeolithic site vs background locations",
    subtitle = "Boxplot interquartile range | Lower Palaeolithic (LP), Middle Palaeolithic (MP), Upper Palaeolithic (UP)",
    x       = NULL,
    y       = "Spectral Index Value",
    caption = "Significance (Mann-Whitney U test): *** p<0.001 | ** p<0.01 | * p<0.05 | ns = not significant"
  ) +
  theme_springer() +
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    axis.text.x      = element_text(size = 6.5, face = "bold"),
    strip.text       = element_text(size = 7.5)
  )

ggsave(
  file.path(fig_path, "Fig03_Spectral_Boxplots.png"),
  fig03,
  width  = 180 / 25.4,
  height = 110 / 25.4,
  dpi    = 300,
  units  = "in",
  bg     = "white"
)
message("Fig03 saved.")


# =============================================================
# FIGURE 5: Structural Geology Proximity Box Plots
# =============================================================

message("\n--- Fig05: Structural Geology Box Plots ---")

struct_vars <- c("dist_fault","dist_dyke",
                 "dist_lineament","dist_shear",
                 "dist_mineral")

struct_labels_vec <- c(
  "dist_fault"      = "Fault",
  "dist_dyke"       = "Dyke",
  "dist_lineament"  = "Lineament",
  "dist_shear"      = "Shear Zone",
  "dist_mineral"    = "Mineral Deposit"
)

struct_long <- full_matrix %>%
  filter(point_type %in% c("site","background")) %>%
  mutate(
    group = case_when(
      point_type == "background" ~ "Background",
      TRUE ~ period
    ),
    group = factor(group,
                   levels = c("Background","LP","MP","UP"))
  ) %>%
  select(group, all_of(struct_vars)) %>%
  mutate(across(all_of(struct_vars), ~. / 1000)) %>%
  pivot_longer(cols      = all_of(struct_vars),
               names_to  = "feature",
               values_to = "distance_km") %>%
  mutate(
    feature = factor(
      recode(feature, !!!struct_labels_vec),
      levels = struct_labels_vec
    )
  ) %>%
  filter(!is.na(distance_km))

fig05 <- ggplot(struct_long,
                aes(x = group, y = distance_km,
                    fill = group)) +
  geom_boxplot(
    outlier.shape = 19,
    outlier.size  = 0.4,
    outlier.alpha = 0.25,
    outlier.color = "grey60",
    linewidth     = 0.3,
    width         = 0.6,
    fatten        = 2
  ) +
  facet_wrap(~feature, scales = "free_y",
             ncol = 3, nrow = 2) +
  scale_fill_manual(values = period_pal,
                    labels = period_labels,
                    name   = NULL) +
  scale_x_discrete(labels = c(
    "Background" = "BG",
    "LP" = "LP", "MP" = "MP", "UP" = "UP"
  )) +
  labs(
    title    = "Fig. 5  |  Proximity of Palaeolithic sites to structural geological features",
    subtitle = "Euclidean distance (km) from site and background locations to nearest structural feature",
    x        = NULL,
    y        = "Distance (km)",
    caption  = "Significance per period: see significance matrix (Fig. 10)"
  ) +
  theme_springer() +
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    axis.text.x      = element_text(size = 6.5, face = "bold"),
    strip.text       = element_text(size = 7.5)
  )

ggsave(
  file.path(fig_path, "Fig05_Structural_Geology_Boxplots.png"),
  fig05,
  width  = 180 / 25.4,
  height = 120 / 25.4,
  dpi    = 300,
  units  = "in",
  bg     = "white"
)
message("Fig05 saved.")


# =============================================================
# FIG07
# =============================================================

load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")
load(file.path(chapter_root, "03_Extracted_Values/master_extraction_matrix.RData"))

library(tidyverse)
library(ggplot2)

fig_path <- file.path(chapter_root, "10_Figures")

period_pal <- c(
  "Background" = "#BDBDBD",
  "LP"         = "#2166AC",
  "MP"         = "#D6604D",
  "UP"         = "#4DAC26"
)

period_labels <- c(
  "Background" = "Background (n=1,000)",
  "LP"         = "LP (n=75)",
  "MP"         = "MP (n=95)",
  "UP"         = "UP (n=27)"
)

soil_vars <- c("soil_depth","soil_erosion",
               "soil_productivity","soil_slope",
               "soil_texture")

soil_labels_vec <- c(
  "soil_depth"        = "Soil Depth",
  "soil_erosion"      = "Erosion Susceptibility",
  "soil_productivity" = "Soil Productivity",
  "soil_slope"        = "Soil Slope",
  "soil_texture"      = "Soil Texture"
)

soil_long <- master_matrix %>%
  filter(point_type %in% c("site","background")) %>%
  mutate(
    group = case_when(
      point_type == "background" ~ "Background",
      TRUE ~ period
    ),
    group = factor(group,
                   levels = c("Background","LP","MP","UP"))
  ) %>%
  select(group, all_of(soil_vars)) %>%
  pivot_longer(cols      = all_of(soil_vars),
               names_to  = "property",
               values_to = "value") %>%
  mutate(
    property = factor(
      recode(property, !!!soil_labels_vec),
      levels = soil_labels_vec
    )
  ) %>%
  filter(!is.na(value))

# Separate background + sites for different jitter density
soil_bg   <- soil_long %>% filter(group == "Background")
soil_site <- soil_long %>% filter(group != "Background")

fig07 <- ggplot(soil_long,
                aes(x = group, y = value, fill = group)) +
  # Jitter background — light, small, many points
  geom_jitter(
    data   = soil_bg,
    aes(color = group),
    width  = 0.25,
    size   = 0.15,
    alpha  = 0.12,
    show.legend = FALSE
  ) +
  # Boxplot
  geom_boxplot(
    outlier.shape = NA,      # jitter already shows outliers
    linewidth     = 0.35,
    width         = 0.55,
    fatten        = 2,
    alpha         = 0.85
  ) +
  # Jitter site points — visible on top
  geom_jitter(
    data   = soil_site,
    aes(color = group),
    width  = 0.18,
    size   = 0.9,
    alpha  = 0.55,
    show.legend = FALSE
  ) +
  # Significance brackets — add manually for key comparisons
  facet_wrap(~property, scales = "free_y",
             ncol = 3, nrow = 2) +
  scale_fill_manual(values = period_pal,
                    labels = period_labels,
                    name   = NULL) +
  scale_color_manual(values = c(
    "Background" = "#888888",
    "LP"         = "#2166AC",
    "MP"         = "#D6604D",
    "UP"         = "#4DAC26"
  )) +
  scale_x_discrete(labels = c(
    "Background" = "BG",
    "LP" = "LP", "MP" = "MP", "UP" = "UP"
  )) +
  labs(
    title    = "Fig. 7  |  Pedogenic properties at Palaeolithic site vs background locations",
    subtitle = "Soil depth, erosion susceptibility, productivity, slope and texture by cultural period | Points = observed values",
    x        = NULL,
    y        = "Ordinal Class Value",
    caption  = "Source: NBSS&LUP soil data | BG = background (n=1,000) | Jitter shows full data distribution | Significance: see Fig. 10"
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title       = element_text(size = 9, face = "bold",
                                    hjust = 0,
                                    margin = margin(0,0,1,0,"mm")),
    plot.subtitle    = element_text(size = 7, color = "grey40",
                                    hjust = 0,
                                    margin = margin(0,0,3,0,"mm")),
    plot.caption     = element_text(size = 6, color = "grey50",
                                    hjust = 0,
                                    margin = margin(2,0,0,0,"mm")),
    axis.title       = element_text(size = 8, color = "grey20"),
    axis.text        = element_text(size = 7, color = "grey30"),
    axis.text.x      = element_text(size = 7, face = "bold"),
    axis.line        = element_line(color = "grey30",
                                    linewidth = 0.35),
    axis.ticks       = element_line(color = "grey50",
                                    linewidth = 0.25),
    panel.grid.major = element_line(color = "grey93",
                                    linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey30",
                                    fill = NA,
                                    linewidth = 0.35),
    strip.background = element_rect(fill = "#F5F5F5",
                                    color = "grey30",
                                    linewidth = 0.3),
    strip.text       = element_text(size = 8, face = "bold",
                                    color = "grey20"),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_text(size = 7, face = "bold"),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(3.5, "mm"),
    plot.background  = element_rect(fill = "white", color = NA),
    plot.margin      = margin(5, 5, 4, 5, "mm")
  )

ggsave(
  file.path(fig_path, "Fig07_Soil_Properties_Boxplots.png"),
  fig07,
  width  = 180 / 25.4,
  height = 130 / 25.4,
  dpi    = 300,
  units  = "in",
  bg     = "white"
)
message("Fig07 saved.")


# =============================================================
# FIG08 FIX — Cleaner PCA biplot
# Top 8 loadings only, reduced arrow scale, no overlap
# =============================================================

load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")
load(file.path(chapter_root, "07_Statistics/mann_whitney_results.RData"))
load(file.path(chapter_root, "09_RandomForest/pca_rf_results.RData"))

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(FactoMineR)
library(factoextra)

if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  install.packages("ggnewscale")
}
library(ggnewscale)

fig_path <- file.path(chapter_root, "10_Figures")

period_pal <- c(
  "LP" = "#2166AC",
  "MP" = "#D6604D",
  "UP" = "#4DAC26"
)

eig_vals <- get_eigenvalue(pca_result)
pc1_var  <- round(eig_vals[1, "variance.percent"], 1)
pc2_var  <- round(eig_vals[2, "variance.percent"], 1)

# -------------------------------------------------------------
# Loading vectors — TOP 8 ONLY by contribution magnitude
# -------------------------------------------------------------

loadings <- as.data.frame(pca_result$var$coord) %>%
  rownames_to_column("variable") %>%
  select(variable, PC1 = Dim.1, PC2 = Dim.2) %>%
  mutate(
    magnitude = sqrt(PC1^2 + PC2^2),
    label = case_when(
      variable == "dist_shear"          ~ "Dist. Shear Zone",
      variable == "dist_lineament"      ~ "Dist. Lineament",
      variable == "dist_dyke"           ~ "Dist. Dyke",
      variable == "chelsa_bio15_lgm"    ~ "Seasonality LGM",
      variable == "chelsa_bio12_lgm"    ~ "Precipitation LGM",
      variable == "chelsa_bio12_modern" ~ "Precipitation Modern",
      variable == "elevation"           ~ "Elevation",
      variable == "burial_depth_m"      ~ "Burial Depth",
      variable == "NDVI"                ~ "NDVI",
      variable == "NDWI"                ~ "NDWI",
      variable == "soil_depth"          ~ "Soil Depth",
      variable == "soil_productivity"   ~ "Soil Productivity",
      TRUE ~ variable
    ),
    category = case_when(
      str_detect(variable,"dist_")  ~ "Structural",
      str_detect(variable,"chelsa") ~ "Climate",
      str_detect(variable,"NDVI|NDWI|MNDWI|NDBI|BSI|SAVI|MSAVI") ~ "Spectral",
      str_detect(variable,"soil_")  ~ "Pedogenic",
      str_detect(variable,"elevation") ~ "Topographic",
      str_detect(variable,"burial")  ~ "Taphonomy",
      TRUE ~ "Other"
    )
  ) %>%
  arrange(desc(magnitude))

# Keep top 8 by magnitude
top8 <- loadings %>% head(8)

# Scale arrows — reduced to 2.5
arrow_scale <- 2.5
top8 <- top8 %>%
  mutate(
    PC1_s = PC1 * arrow_scale,
    PC2_s = PC2 * arrow_scale
  )

arrow_pal <- c(
  "Structural"  = "#C0392B",
  "Climate"     = "#E67E22",
  "Spectral"    = "#2980B9",
  "Pedogenic"   = "#27AE60",
  "Topographic" = "#8E44AD",
  "Taphonomy"   = "#6D4C41"
)

# PCA scores
pca_plot <- pca_scores %>%
  mutate(
    group = case_when(
      point_type == "background" ~ "Background",
      TRUE ~ period
    ),
    group = factor(group,
                   levels = c("Background","LP","MP","UP"))
  )

# Centroids
lp_c <- period_centroids %>% filter(period == "LP")
mp_c <- period_centroids %>% filter(period == "MP")
up_c <- period_centroids %>% filter(period == "UP")

# -------------------------------------------------------------
# Build figure
# -------------------------------------------------------------

fig08 <- ggplot() +
  # Background cloud — very faint
  geom_point(
    data  = pca_plot %>% filter(group == "Background"),
    aes(x = PC1, y = PC2),
    color = "#CCCCCC",
    alpha = 0.10,
    size  = 0.25
  ) +
  # Reference lines
  geom_hline(yintercept = 0, color = "grey75",
             linewidth = 0.2, linetype = "dotted") +
  geom_vline(xintercept = 0, color = "grey75",
             linewidth = 0.2, linetype = "dotted") +
  # Loading arrows — top 8 only
  geom_segment(
    data = top8,
    aes(x = 0, y = 0,
        xend = PC1_s, yend = PC2_s,
        color = category),
    arrow     = arrow(length = unit(1.5,"mm"),
                      type = "open"),
    linewidth = 0.45,
    alpha     = 0.85
  ) +
  scale_color_manual(
    values = arrow_pal,
    name   = "Variable\nCategory"
  ) +
  # Arrow labels — repelled, no overlap
  geom_text_repel(
    data          = top8,
    aes(x = PC1_s * 1.15,
        y = PC2_s * 1.15,
        label = label,
        color = category),
    size          = 2.4,
    fontface      = "bold",
    segment.size  = 0.2,
    segment.alpha = 0.4,
    segment.color = "grey60",
    max.overlaps  = 20,
    box.padding   = 0.4,
    point.padding = 0.2,
    show.legend   = FALSE,
    seed          = 42
  ) +
  # Site points — new colour scale
  new_scale_color() +
  geom_point(
    data  = pca_plot %>% filter(group != "Background"),
    aes(x = PC1, y = PC2, color = group),
    alpha = 0.70,
    size  = 1.0
  ) +
  # Confidence ellipses per period
  stat_ellipse(
    data      = pca_plot %>% filter(group != "Background"),
    aes(x = PC1, y = PC2, color = group),
    level     = 0.68,
    linewidth = 0.45,
    linetype  = "dashed"
  ) +
  # Period centroids — large diamonds
  geom_point(
    data  = period_centroids,
    aes(x = PC1_mean, y = PC2_mean, color = period),
    size  = 5,
    shape = 18
  ) +
  # Centroid labels — white background box
  geom_label(
    data = period_centroids,
    aes(x = PC1_mean, y = PC2_mean,
        label = period, color = period),
    size          = 3.0,
    fontface      = "bold",
    label.size    = 0.2,
    label.padding = unit(1, "mm"),
    nudge_y       = 0.55,
    show.legend   = FALSE
  ) +
  # Displacement arrow LP → MP
  annotate("segment",
           x    = lp_c$PC1_mean,
           y    = lp_c$PC2_mean,
           xend = mp_c$PC1_mean,
           yend = mp_c$PC2_mean,
           arrow = arrow(length = unit(2,"mm"),
                         type = "closed"),
           color     = "grey10",
           linewidth = 0.65) +
  # Displacement arrow MP → UP
  annotate("segment",
           x    = mp_c$PC1_mean,
           y    = mp_c$PC2_mean,
           xend = up_c$PC1_mean,
           yend = up_c$PC2_mean,
           arrow = arrow(length = unit(2,"mm"),
                         type = "closed"),
           color     = "grey10",
           linewidth = 0.65) +
  # Displacement labels — white background
  annotate("label",
           x     = mean(c(lp_c$PC1_mean, mp_c$PC1_mean)),
           y     = mean(c(lp_c$PC2_mean, mp_c$PC2_mean)) + 0.6,
           label = paste0("\u0394 = ",
                          round(displacements$displacement[1],2),
                          " PC units"),
           size      = 2.3,
           color     = "grey15",
           fill      = "white",
           label.size = 0.2,
           fontface  = "italic") +
  annotate("label",
           x     = mean(c(mp_c$PC1_mean, up_c$PC1_mean)) - 0.3,
           y     = mean(c(mp_c$PC2_mean, up_c$PC2_mean)) - 0.6,
           label = paste0("\u0394 = ",
                          round(displacements$displacement[2],2),
                          " PC units"),
           size      = 2.3,
           color     = "grey15",
           fill      = "white",
           label.size = 0.2,
           fontface  = "italic") +
  scale_color_manual(
    values = period_pal,
    labels = c(
      "LP" = "LP (n=75)",
      "MP" = "MP (n=95)",
      "UP" = "UP (n=27)"
    ),
    name = "Cultural\nPeriod"
  ) +
  coord_cartesian(
    xlim = c(-6.5, 6.5),
    ylim = c(-5.5, 5.5)
  ) +
  labs(
    title    = "Fig. 8  |  PCA biplot of multivariate geoenvironmental space",
    subtitle = paste0(
      "Site and background locations with LP/MP/UP centroids | ",
      "Top 8 loading vectors shown | ",
      "PC1 = ", pc1_var, "% | PC2 = ", pc2_var, "% variance explained"
    ),
    x        = paste0("PC1 — ", pc1_var,
                      "% variance explained"),
    y        = paste0("PC2 — ", pc2_var,
                      "% variance explained"),
    caption  = paste0(
      "Ellipses = 68% confidence intervals | ",
      "Arrows = top 8 variable loading vectors | ",
      "\u0394 = Euclidean centroid displacement (PC units) | ",
      "Background cloud n=985"
    )
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title       = element_text(size = 9, face = "bold",
                                    hjust = 0,
                                    margin = margin(0,0,1,0,"mm")),
    plot.subtitle    = element_text(size = 7, color = "grey40",
                                    hjust = 0,
                                    margin = margin(0,0,3,0,"mm")),
    plot.caption     = element_text(size = 6, color = "grey50",
                                    hjust = 0,
                                    margin = margin(2,0,0,0,"mm")),
    axis.title       = element_text(size = 8, color = "grey20"),
    axis.text        = element_text(size = 7, color = "grey30"),
    axis.line        = element_line(color = "grey30",
                                    linewidth = 0.35),
    axis.ticks       = element_line(color = "grey50",
                                    linewidth = 0.25),
    panel.grid.major = element_line(color = "grey93",
                                    linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey30",
                                    fill = NA,
                                    linewidth = 0.35),
    legend.position  = "right",
    legend.box       = "vertical",
    legend.title     = element_text(size = 7, face = "bold"),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(3.5,"mm"),
    plot.background  = element_rect(fill = "white",
                                    color = NA),
    plot.margin      = margin(5, 5, 4, 5, "mm")
  )

ggsave(
  file.path(fig_path, "Fig08_PCA_Biplot.png"),
  fig08,
  width  = 180 / 25.4,
  height = 160 / 25.4,
  dpi    = 300,
  units  = "in",
  bg     = "white"
)
message("Fig08 saved.")


# =============================================================
# FIGURE 9: Random Forest Variable Importance
# =============================================================

message("\n--- Fig09: RF Variable Importance ---")

oob_acc <- round(
  (1 - rf_model$err.rate[rf_ntree, "OOB"]) * 100, 1
)

cat_pal <- c(
  "Structural Geology" = "#C0392B",
  "Climate"            = "#E67E22",
  "Spectral"           = "#2980B9",
  "Pedogenic"          = "#27AE60",
  "Topographic"        = "#8E44AD",
  "Taphonomy"          = "#795548"
)

# Clean labels
importance_clean <- importance_df %>%
  head(20) %>%
  mutate(
    label = case_when(
      variable == "dist_shear"          ~ "Distance to Shear Zone",
      variable == "dist_lineament"      ~ "Distance to Lineament",
      variable == "chelsa_bio01_lgm"    ~ "CHELSA bio01 — LGM",
      variable == "NDWI"                ~ "NDWI",
      variable == "chelsa_bio12_lgm"    ~ "CHELSA bio12 — LGM",
      variable == "elevation"           ~ "Elevation",
      variable == "chelsa_bio12_modern" ~ "CHELSA bio12 — Modern",
      variable == "chelsa_bio01_modern" ~ "CHELSA bio01 — Modern",
      variable == "burial_depth_m"      ~ "Burial Depth",
      variable == "chelsa_bio15_modern" ~ "CHELSA bio15 — Modern",
      variable == "chelsa_bio15_lgm"    ~ "CHELSA bio15 — LGM",
      variable == "MNDWI"               ~ "MNDWI",
      variable == "dist_fault"          ~ "Distance to Fault",
      variable == "BSI"                 ~ "BSI",
      variable == "dist_mineral"        ~ "Distance to Mineral Deposit",
      variable == "dist_dyke"           ~ "Distance to Dyke",
      variable == "NDVI"                ~ "NDVI",
      variable == "soil_depth"          ~ "Soil Depth",
      variable == "soil_productivity"   ~ "Soil Productivity",
      variable == "NDBI"                ~ "NDBI",
      TRUE ~ gsub("_"," ", variable)
    ),
    label = fct_reorder(label, MeanDecreaseGini)
  )

fig09 <- ggplot(
  importance_clean,
  aes(x = MeanDecreaseGini, y = label, fill = category)
) +
  geom_col(width = 0.72,
           color = "white", linewidth = 0.15) +
  geom_vline(xintercept = 0, color = "grey20",
             linewidth = 0.3) +
  geom_text(
    aes(label = paste0("#", rank)),
    x     = 0.2,
    hjust = 0,
    size  = 2.2,
    color = "white",
    fontface = "bold"
  ) +
  geom_text(
    aes(label = round(MeanDecreaseGini, 1),
        x     = MeanDecreaseGini + 0.15),
    hjust    = 0,
    size     = 2.2,
    color    = "grey30"
  ) +
  scale_fill_manual(
    values = cat_pal,
    name   = "Variable Category"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.18)),
    breaks = seq(0, 20, 4)
  ) +
  labs(
    title    = "Fig. 9  |  Random Forest variable importance ranking",
    subtitle = paste0(
      "Mean Decrease Gini | Top 20 variables | ",
      "OOB classification accuracy = ", oob_acc, "%"
    ),
    x        = "Mean Decrease Gini",
    y        = NULL,
    caption  = paste0(
      "Parameters: ntree=1,000, mtry=5, balanced classes (n=197 sites + 197 background) | ",
      "Seed=42"
    )
  ) +
  theme_springer() +
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_text(size = 7),
    legend.text      = element_text(size = 6.5),
    axis.text.y      = element_text(size = 6.8),
    panel.grid.major.x = element_line(color = "grey90",
                                      linewidth = 0.2),
    panel.grid.major.y = element_blank()
  )

ggsave(
  file.path(fig_path, "Fig09_RF_Variable_Importance.png"),
  fig09,
  width  = 180 / 25.4,
  height = 150 / 25.4,
  dpi    = 300,
  units  = "in",
  bg     = "white"
)
message("Fig09 saved.")

# =============================================================
# FIGURE 10 — REBUILT: Geoenvironmental Significance Matrix
# ggplot2-based — full layout control
# =============================================================

load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")
load(file.path(chapter_root, "07_Statistics/mann_whitney_results.RData"))

library(tidyverse)
library(ggplot2)
library(ggtext)

fig_path <- file.path(chapter_root, "10_Figures")

# -------------------------------------------------------------
# STEP 1: Build long-format significance data
# -------------------------------------------------------------

# Clean variable labels
var_labels <- c(
  "NDVI"                = "NDVI",
  "NDWI"                = "NDWI",
  "MNDWI"               = "MNDWI",
  "NDBI"                = "NDBI",
  "BSI"                 = "BSI",
  "SAVI"                = "SAVI",
  "MSAVI"               = "MSAVI",
  "elevation"           = "Elevation",
  "soil_depth"          = "Soil Depth",
  "soil_erosion"        = "Soil Erosion",
  "soil_productivity"   = "Soil Productivity",
  "soil_slope"          = "Soil Slope",
  "soil_texture"        = "Soil Texture",
  "chelsa_bio01_lgm"    = "Temperature — LGM",
  "chelsa_bio12_lgm"    = "Precipitation — LGM",
  "chelsa_bio15_lgm"    = "Seasonality — LGM",
  "chelsa_bio01_modern" = "Temperature — Modern",
  "chelsa_bio12_modern" = "Precipitation — Modern",
  "chelsa_bio15_modern" = "Seasonality — Modern",
  "dist_fault"          = "Distance to Fault",
  "dist_dyke"           = "Distance to Dyke",
  "dist_lineament"      = "Distance to Lineament",
  "dist_shear"          = "Distance to Shear Zone",
  "dist_mineral"        = "Distance to Mineral Deposit",
  "burial_depth_m"      = "Burial Depth"
)

# Category assignment
var_categories <- c(
  "NDVI"                = "Spectral",
  "NDWI"                = "Spectral",
  "MNDWI"               = "Spectral",
  "NDBI"                = "Spectral",
  "BSI"                 = "Spectral",
  "SAVI"                = "Spectral",
  "MSAVI"               = "Spectral",
  "elevation"           = "Topographic",
  "soil_depth"          = "Pedogenic",
  "soil_erosion"        = "Pedogenic",
  "soil_productivity"   = "Pedogenic",
  "soil_slope"          = "Pedogenic",
  "soil_texture"        = "Pedogenic",
  "chelsa_bio01_lgm"    = "Climate",
  "chelsa_bio12_lgm"    = "Climate",
  "chelsa_bio15_lgm"    = "Climate",
  "chelsa_bio01_modern" = "Climate",
  "chelsa_bio12_modern" = "Climate",
  "chelsa_bio15_modern" = "Climate",
  "dist_fault"          = "Structural",
  "dist_dyke"           = "Structural",
  "dist_lineament"      = "Structural",
  "dist_shear"          = "Structural",
  "dist_mineral"        = "Structural",
  "burial_depth_m"      = "Taphonomy"
)

cat_pal <- c(
  "Spectral"    = "#2980B9",
  "Topographic" = "#8E44AD",
  "Pedogenic"   = "#27AE60",
  "Climate"     = "#E67E22",
  "Structural"  = "#C0392B",
  "Taphonomy"   = "#6D4C41"
)

# Build plot data
heat_long <- mw_table %>%
  filter(!is.na(p_value)) %>%
  mutate(
    log_p     = -log10(pmax(p_value, 1e-10)),
    sig_star  = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    var_label = recode(variable, !!!var_labels),
    category  = recode(variable, !!!var_categories),
    period    = factor(period, levels = c("LP", "MP", "UP"))
  )

# -------------------------------------------------------------
# STEP 2: Cluster row order by similarity
# -------------------------------------------------------------

# Pivot to wide for clustering
heat_wide <- heat_long %>%
  select(variable, period, log_p) %>%
  pivot_wider(names_from = period, values_from = log_p,
              values_fill = 0) %>%
  column_to_rownames("variable")

# Hierarchical clustering
row_order <- hclust(dist(heat_wide))$order
ordered_vars <- rownames(heat_wide)[row_order]

# Apply order with labels
heat_long <- heat_long %>%
  mutate(
    var_label = factor(var_label,
                       levels = recode(ordered_vars,
                                       !!!var_labels)),
    category  = factor(category,
                       levels = c("Spectral","Topographic",
                                  "Pedogenic","Climate",
                                  "Structural","Taphonomy"))
  )

# -------------------------------------------------------------
# STEP 3: Build Figure
# -------------------------------------------------------------

fig10 <- ggplot(
  heat_long,
  aes(x = period, y = var_label)
) +
  # Heatmap tiles
  geom_tile(
    aes(fill = log_p),
    color     = "white",
    linewidth = 0.8,
    width     = 0.92,
    height    = 0.92
  ) +
  # Significance stars
  geom_text(
    aes(label = sig_star),
    size     = 3.2,
    color    = "grey15",
    fontface = "bold"
  ) +
  # Category colour strip — left side
  geom_tile(
    aes(x = 0.28, fill = category),
    color     = "white",
    linewidth = 0.5,
    width     = 0.18,
    height    = 0.92,
    show.legend = FALSE
  ) +
  # Colour scale — heatmap
  scale_fill_gradientn(
    colours  = c(
      "#FFFFFF", "#FFF7BC",
      "#FEC44F", "#FE9929",
      "#D95F0E", "#7F2704"
    ),
    values   = scales::rescale(c(0, 1.301, 2, 3, 4, 5)),
    limits   = c(0, 5),
    oob      = scales::squish,
    name     = expression(-log[10](italic(p))),
    breaks   = c(0, 1.301, 2, 3, 4, 5),
    labels   = c(
      "0 (ns)",
      "0.05",
      "0.01",
      "0.001",
      "0.0001",
      "0.00001"
    ),
    guide    = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = unit(3.5, "mm"),
      barheight      = unit(40, "mm"),
      ticks.colour   = "grey20",
      frame.colour   = "grey20",
      frame.linewidth = 0.3
    )
  ) +
  # Category strip labels
  annotate(
    "text",
    x     = 0.28,
    y     = Inf,
    label = "",
    size  = 1
  ) +
  # X axis — period labels
  scale_x_discrete(
    expand = expansion(add = c(0.6, 0.5)),
    labels = c(
      "LP" = "Lower\nPalaeolithic",
      "MP" = "Middle\nPalaeolithic",
      "UP" = "Upper\nPalaeolithic"
    )
  ) +
  # Y axis
  scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +
  # Period significance count annotation
  annotate(
    "text",
    x     = c(1, 2, 3),
    y     = 0.2,
    label = c("19 sig.", "18 sig.", "6 sig."),
    size  = 2.3,
    color = "grey30",
    fontface = "italic"
  ) +
  # Category legend — right side manual
  labs(
    title   = "Fig. 10  |  Geoenvironmental significance matrix",
    subtitle = paste0(
      "Mann-Whitney U tests: site vs background per cultural period | ",
      "25 variables × 3 periods (75 tests total)"
    ),
    x       = NULL,
    y       = NULL,
    caption = paste0(
      "Colour intensity = -log\u2081\u2080(p-value) | ",
      "Stars: *** p<0.001 | ** p<0.01 | * p<0.05 | ",
      "Left strip = variable category | ",
      "Row order by hierarchical clustering"
    )
  ) +
  theme_classic(base_size = 8) +
  theme(
    # Title
    plot.title      = element_text(size = 9, face = "bold",
                                   hjust = 0,
                                   margin = margin(0,0,1,0,"mm")),
    plot.subtitle   = element_text(size = 7, color = "grey40",
                                   hjust = 0,
                                   margin = margin(0,0,3,0,"mm")),
    plot.caption    = element_text(size = 6, color = "grey50",
                                   hjust = 0,
                                   margin = margin(2,0,0,0,"mm")),
    # Axes
    axis.text.x     = element_text(size = 7.5, face = "bold",
                                   color = "grey15",
                                   lineheight = 1.2),
    axis.text.y     = element_text(size = 6.8, color = "grey20",
                                   hjust = 1),
    axis.ticks      = element_blank(),
    axis.line       = element_blank(),
    # Panel
    panel.background = element_rect(fill = "#FAFAFA",
                                    color = NA),
    panel.border    = element_rect(color = "grey30",
                                   fill = NA,
                                   linewidth = 0.4),
    panel.grid      = element_blank(),
    # Legend
    legend.position = "right",
    legend.title    = element_text(size = 7, face = "bold",
                                   hjust = 0.5),
    legend.text     = element_text(size = 6.5),
    # Margin
    plot.margin     = margin(5, 5, 4, 5, "mm"),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Add category legend as separate annotation panel
cat_legend <- heat_long %>%
  distinct(category) %>%
  mutate(y = rev(seq_along(category))) %>%
  ggplot(aes(x = 0.5, y = y, color = category)) +
  geom_point(size = 4, shape = 15) +
  geom_text(aes(label = category, x = 0.65),
            hjust = 0, size = 2.3, color = "grey20") +
  scale_color_manual(values = cat_pal) +
  scale_x_continuous(limits = c(0.4, 1.6)) +
  labs(title = "Category") +
  theme_void(base_size = 7) +
  theme(
    plot.title    = element_text(size = 7, face = "bold",
                                 hjust = 0,
                                 margin = margin(0,0,2,0,"mm")),
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(5, 2, 5, 2, "mm")
  )

# Combine heatmap + category legend
library(patchwork)

fig10_final <- fig10 + cat_legend +
  plot_layout(widths = c(5, 1))

ggsave(
  file.path(fig_path, "Fig10_Significance_Matrix_Heatmap.png"),
  fig10_final,
  width  = 180 / 25.4,
  height = 190 / 25.4,
  dpi    = 300,
  units  = "in",
  bg     = "white"
)

message("Fig10 saved.")


# =============================================================
# FINAL CHECK
# =============================================================

message("\n=============================================================")
message("SCRIPT 08 COMPLETE")
message("=============================================================")

figs <- c(
  "Fig03_Spectral_Boxplots.png",
  "Fig05_Structural_Geology_Boxplots.png",
  "Fig07_Soil_Properties_Boxplots.png",
  "Fig08_PCA_Biplot.png",
  "Fig09_RF_Variable_Importance.png",
  "Fig10_Significance_Matrix_Heatmap.png"
)

for (f in figs) {
  fp   <- file.path(fig_path, f)
  ok   <- file.exists(fp)
  size <- ifelse(ok,
                 paste0(round(file.size(fp)/1e6, 1), " MB"),
                 "MISSING")
  message(ifelse(ok, "OK", "XX"), "  ", f, "  [", size, "]")
}

message("\nSaved to: ", fig_path)
message("QGIS figures: Fig01, Fig02, Fig04, Fig06, Fig11, Fig12")
message("Next: Write chapter draft")
message("=============================================================")