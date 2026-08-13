# =============================================================
# SCRIPT 08 (REBUILT): Publication-Ready Statistical Figures
# =============================================================
# Project: Reading the Palaeolithic Landscape
# Chapter: Mapping the Past (Springer Nature, 2026)
# Author: Sushant Begade | RTMNU Nagpur
# ORCID: 0009-0003-0804-1763
# Rebuild date: August 2026
# =============================================================
# CHANGES FROM ORIGINAL SCRIPT 08:
#   1. All source data now pulled from REBUILT Scripts 01-07 outputs.
#   2. Fig08 (biplot): PCoA (Gower) replaces Euclidean PCA. cmdscale
#      has no native variable loadings, so vegan::envfit() fits
#      vector correlations of continuous variables onto the PCoA
#      ordination — the standard approach for PCoA biplots. Centroid
#      displacement labels now show bootstrap 95% CI, not a bare
#      descriptive number.
#   3. Fig09 (RF importance): primary metric switched to permutation
#      MeanDecreaseAccuracy. NEW second panel: grouped importance
#      (answers reviewer Problem #12 directly in the figure, not just
#      in text).
#   4. Fig10 (significance matrix): rebuilt for 30 variables x 3
#      periods (90 tests, was 25x3=75). NEW: structural-geology cells
#      in LP/MP columns get a marker showing whether that result
#      survived the Script 06 spatial block bootstrap — distinguishes
#      "significant" from "significant AND spatially robust" directly
#      on the figure.
# FIGURES PRODUCED (11-figure final chapter set, Fig 6 removed entirely):
# Fig01 — Study area map [QGIS]
# Fig02 — Sentinel-2A false colour composite [QGIS]
# Fig03 — Spectral index box plots
# Fig04 — Geochemistry spatial distribution [QGIS]
# Fig05 — Structural geology proximity box plots
# Fig06 — Soil properties box plots
# Fig07 — PCoA biplot with centroids + bootstrap-CI displacement vectors
# Fig08 — Random Forest variable importance (individual + grouped)
# Fig09 — Significance matrix heatmap (+ spatial bootstrap markers)
# Fig10 — Geoenvironmental suitability surface [QGIS]
# Fig11 — CHELSA palaeoclimate spatial context [QGIS]

load("E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset/Reading_the_Palaeolithic_Landscape_Chapter/00_Scripts/global_params.RData")

load(file.path(chapter_root, "03_Extracted_Values/master_extraction_matrix.RData"))
load(file.path(chapter_root, "05_Structural_Geology/vector_extraction_matrix.RData"))
load(file.path(chapter_root, "07_Statistics/mann_whitney_results.RData"))
load(file.path(chapter_root, "09_RandomForest/pcoa_rf_results.RData"))

library(tidyverse)
library(ggplot2)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(scales)
library(ggrepel)
library(vegan)

set.seed(rf_seed)

fig_path <- file.path(chapter_root, "10_Figures")

theme_springer <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title    = element_text(size = base_size+1, face="bold", hjust=0, margin=margin(0,0,2,0,"mm")),
      plot.subtitle = element_text(size = base_size-0.5, hjust=0, color="grey40", margin=margin(0,0,3,0,"mm")),
      plot.caption  = element_text(size = base_size-1.5, color="grey55", hjust=0, margin=margin(2,0,0,0,"mm")),
      axis.title    = element_text(size = base_size, color="grey20"),
      axis.text     = element_text(size = base_size-1, color="grey30"),
      axis.line     = element_line(color="grey30", linewidth=0.35),
      axis.ticks    = element_line(color="grey50", linewidth=0.25),
      panel.grid.major = element_line(color="grey93", linewidth=0.2),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill="white", color=NA),
      panel.border  = element_rect(color="grey30", fill=NA, linewidth=0.35),
      strip.background = element_rect(fill="#F5F5F5", color="grey30", linewidth=0.3),
      strip.text    = element_text(size = base_size, face="bold", color="grey20"),
      legend.title  = element_text(size = base_size, face="bold"),
      legend.text   = element_text(size = base_size-1),
      legend.key.size = unit(3.5,"mm"),
      legend.key    = element_rect(fill="white", color=NA),
      legend.background = element_rect(fill="white", color="grey80", linewidth=0.25),
      plot.background = element_rect(fill="white", color=NA),
      plot.margin   = margin(5,5,4,5,"mm")
    )
}

period_pal <- c("Background"="#BDBDBD", "LP"="#2166AC", "MP"="#D6604D", "UP"="#4DAC26")
period_labels <- c("Background"="Background (n=1,000)", "LP"="LP (n=48)",
                   "MP"="MP (n=62)", "UP"="UP (n=27)")
# n's updated: post period-classification-bug-fix (Script 01), LP=48 MP=62 UP=27,
# MULTI=60 handled separately (Persistent Places, not in these boxplots)

message("Setup complete. Building figures from rebuilt Script 01-07 outputs...")

# -------------------------------------------------------------
# Recompute model_imputed / CONTINUOUS_VARS / ORDINAL_VARS
# -------------------------------------------------------------
# These were local to Script 07's session, not saved into
# pcoa_rf_results.RData. Deterministic given the same seed + same
# NA-filter/imputation logic as Script 07 Section 2 — recomputing
# here reproduces identical values without re-running the expensive
# Gower distance / RF fitting steps.

ORDINAL_VARS <- c("soil_depth", "soil_erosion", "soil_productivity",
                  "soil_slope", "soil_texture")
CONTINUOUS_VARS <- setdiff(analysis_vars, ORDINAL_VARS)

model_data_fig <- full_matrix_coords %>%
  select(point_id, point_type, period, block_id, all_of(analysis_vars)) %>%
  filter(!is.na(point_type))
na_counts_fig <- rowSums(is.na(model_data_fig %>% select(all_of(analysis_vars))))
model_data_fig <- model_data_fig[na_counts_fig <= 3, ]

model_imputed <- model_data_fig
for (var in analysis_vars) {
  if (any(is.na(model_imputed[[var]]))) {
    med_val <- median(model_imputed[[var]], na.rm = TRUE)
    model_imputed[[var]][is.na(model_imputed[[var]])] <- med_val
  }
}
message("Recomputed model_imputed: ", nrow(model_imputed), " rows (should match Script 07's 716)")
stopifnot("model_imputed row count mismatch vs Script 07 log (716) — investigate before trusting envfit" =
            nrow(model_imputed) == 716)


# =============================================================
# FIGURE 3: Spectral Index Box Plots  (30-var significance now)
# =============================================================

message("\n--- Fig03: Spectral Index Box Plots ---")

spectral_vars <- c("NDVI","NDWI","MNDWI","NDBI","BSI","SAVI","MSAVI")

spec_long <- full_matrix_coords %>%
  filter(point_type %in% c("site","background"), period %in% c("LP","MP","UP","background")) %>%
  mutate(group = ifelse(point_type=="background","Background",period),
         group = factor(group, levels=c("Background","LP","MP","UP"))) %>%
  select(group, all_of(spectral_vars)) %>%
  pivot_longer(all_of(spectral_vars), names_to="index", values_to="value") %>%
  filter(!is.na(value)) %>% mutate(index=factor(index, levels=spectral_vars))

fig03 <- ggplot(spec_long, aes(x=group, y=value, fill=group)) +
  geom_boxplot(outlier.shape=19, outlier.size=0.4, outlier.alpha=0.25,
               outlier.color="grey60", linewidth=0.3, width=0.6, median.linewidth=2) +
  geom_hline(yintercept=0, color="grey40", linewidth=0.2, linetype="longdash") +
  facet_wrap(~index, scales="free_y", ncol=4, nrow=2) +
  scale_fill_manual(values=period_pal, labels=period_labels, name=NULL) +
  scale_x_discrete(labels=c("Background"="BG","LP"="LP","MP"="MP","UP"="UP")) +
  labs(title="Fig. 3  |  Sentinel-2A spectral index values at Palaeolithic site vs background locations",
       subtitle="Boxplot interquartile range | LP (n=48), MP (n=62), UP (n=27) — post coordinate-precision-corrected classification",
       x=NULL, y="Spectral Index Value",
       caption="Significance (Mann-Whitney U, 30-variable battery): see Fig. 9") +
  theme_springer() + theme(legend.position="bottom", legend.direction="horizontal",
                           axis.text.x=element_text(size=6.5, face="bold"), strip.text=element_text(size=7.5))

ggsave(file.path(fig_path,"Fig03_Spectral_Boxplots.png"), fig03,
       width=180/25.4, height=110/25.4, dpi=300, units="in", bg="white")
message("Fig03 saved.")


# =============================================================
# FIGURE 5: Structural Geology Proximity Box Plots
# =============================================================

message("\n--- Fig05: Structural Geology Box Plots ---")

struct_vars <- c("dist_fault","dist_dyke","dist_lineament","dist_shear","dist_mineral")
struct_labels_vec <- c(dist_fault="Fault", dist_dyke="Dyke", dist_lineament="Lineament",
                       dist_shear="Shear Zone", dist_mineral="Mineral Deposit")

struct_long <- full_matrix_coords %>%
  filter(point_type %in% c("site","background"), period %in% c("LP","MP","UP","background")) %>%
  mutate(group = ifelse(point_type=="background","Background",period),
         group = factor(group, levels=c("Background","LP","MP","UP"))) %>%
  select(group, all_of(struct_vars)) %>%
  mutate(across(all_of(struct_vars), ~./1000)) %>%
  pivot_longer(all_of(struct_vars), names_to="feature", values_to="distance_km") %>%
  mutate(feature = factor(recode(feature, !!!struct_labels_vec), levels=struct_labels_vec)) %>%
  filter(!is.na(distance_km))

fig05 <- ggplot(struct_long, aes(x=group, y=distance_km, fill=group)) +
  geom_boxplot(outlier.shape=19, outlier.size=0.4, outlier.alpha=0.25,
               outlier.color="grey60", linewidth=0.3, width=0.6, median.linewidth=2) +
  facet_wrap(~feature, scales="free_y", ncol=3, nrow=2) +
  scale_fill_manual(values=period_pal, labels=period_labels, name=NULL) +
  scale_x_discrete(labels=c("Background"="BG","LP"="LP","MP"="MP","UP"="UP")) +
  labs(title="Fig. 5  |  Proximity of Palaeolithic sites to structural geological features",
       subtitle="Euclidean distance (km) to nearest structural feature by cultural period",
       x=NULL, y="Distance (km)",
       caption="LP structural associations survive spatial block bootstrap (4/4 variables); MP partially (2/4) — see Fig. 9 markers") +
  theme_springer() + theme(legend.position="bottom", legend.direction="horizontal",
                           axis.text.x=element_text(size=6.5, face="bold"), strip.text=element_text(size=7.5))

ggsave(file.path(fig_path,"Fig05_Structural_Geology_Boxplots.png"), fig05,
       width=180/25.4, height=120/25.4, dpi=300, units="in", bg="white")
message("Fig05 saved.")


# =============================================================
# FIGURE 7: Soil Properties Box Plots
# =============================================================

message("\n--- Fig07: Soil Properties Box Plots ---")

soil_vars <- c("soil_depth","soil_erosion","soil_productivity","soil_slope","soil_texture")
soil_labels_vec <- c(soil_depth="Soil Depth", soil_erosion="Erosion Susceptibility",
                     soil_productivity="Soil Productivity", soil_slope="Soil Slope",
                     soil_texture="Soil Texture")

soil_long <- full_matrix_coords %>%
  filter(point_type %in% c("site","background"), period %in% c("LP","MP","UP","background")) %>%
  mutate(group = ifelse(point_type=="background","Background",period),
         group = factor(group, levels=c("Background","LP","MP","UP"))) %>%
  select(group, all_of(soil_vars)) %>%
  pivot_longer(all_of(soil_vars), names_to="property", values_to="value") %>%
  mutate(property = factor(recode(property, !!!soil_labels_vec), levels=soil_labels_vec)) %>%
  filter(!is.na(value))

soil_bg <- soil_long %>% filter(group=="Background")
soil_site <- soil_long %>% filter(group!="Background")

fig07 <- ggplot(soil_long, aes(x=group, y=value, fill=group)) +
  geom_jitter(data=soil_bg, aes(color=group), width=0.25, size=0.15, alpha=0.12, show.legend=FALSE) +
  geom_boxplot(outlier.shape=NA, linewidth=0.35, width=0.55, median.linewidth=2, alpha=0.85) +
  geom_jitter(data=soil_site, aes(color=group), width=0.18, size=0.9, alpha=0.55, show.legend=FALSE) +
  facet_wrap(~property, scales="free_y", ncol=3, nrow=2) +
  scale_fill_manual(values=period_pal, labels=period_labels, name=NULL) +
  scale_color_manual(values=c(Background="#888888", LP="#2166AC", MP="#D6604D", UP="#4DAC26")) +
  scale_x_discrete(labels=c(Background="BG", LP="LP", MP="MP", UP="UP")) +
  labs(title="Fig. 6  |  Pedogenic properties at Palaeolithic site vs background locations",
       subtitle="Ordinal soil-class variables — Wilcoxon rank-sum valid on ordinal data (see Methods; PCoA also treats these correctly via Gower distance)",
       x=NULL, y="Ordinal Class Value",
       caption="Source: NBSS&LUP soil data | Jitter shows full data distribution | Significance: see Fig. 9") +
  theme_springer() + theme(legend.position="bottom", legend.direction="horizontal",
                           axis.text.x=element_text(size=6.5, face="bold"), strip.text=element_text(size=7.5))

ggsave(file.path(fig_path,"Fig06_Soil_Properties_Boxplots.png"), fig07,
       width=180/25.4, height=130/25.4, dpi=300, units="in", bg="white")
message("Fig07 saved.")


# =============================================================
# FIGURE 8: PCoA Biplot (Gower) — replaces PCA biplot
# =============================================================

message("\n--- Fig08: PCoA Biplot (Gower distance) ---")

var_explained_pct <- round(100 * (pcoa_result$eig[pcoa_result$eig>0]) /
                             sum(pcoa_result$eig[pcoa_result$eig>0]), 1)
pc1_var <- var_explained_pct[1]; pc2_var <- var_explained_pct[2]

# envfit: fit continuous-variable vectors onto the PCoA ordination —
# cmdscale has no native loadings, this is the standard PCoA-biplot
# substitute. Fit on the SAME rows used to build pcoa_scores.
message("Fitting environmental vectors (vegan::envfit) onto PCoA axes...")
env_data_for_fit <- model_imputed %>% select(all_of(CONTINUOUS_VARS))
ord_matrix <- pcoa_scores %>% select(PC1, PC2) %>% as.matrix()

env_fit <- envfit(ord_matrix, env_data_for_fit, permutations = 999)
env_vectors <- as.data.frame(env_fit$vectors$arrows * sqrt(env_fit$vectors$r)) %>%
  rownames_to_column("variable") %>%
  mutate(r2 = env_fit$vectors$r, pval = env_fit$vectors$pvals)

top8_env <- env_vectors %>% arrange(desc(r2)) %>% head(8) %>%
  mutate(category = case_when(
    str_detect(variable,"dist_") ~ "Structural",
    str_detect(variable,"chelsa") ~ "Climate",
    str_detect(variable,"NDVI|NDWI|MNDWI|NDBI|BSI|SAVI|MSAVI") ~ "Spectral",
    str_detect(variable,"geochem") ~ "Geochemistry",
    str_detect(variable,"elevation") ~ "Topographic",
    TRUE ~ "Other"
  ))

message("Top 8 environmental vectors by fit R2:")
print(top8_env %>% select(variable, r2, pval, category))

arrow_scale <- 1.5
top8_env <- top8_env %>% mutate(PC1_s = PC1*arrow_scale, PC2_s = PC2*arrow_scale)

# FIX: PCoA (Gower) axis units are tiny (roughly +/-0.1 to +/-1.5) compared
# to the old Euclidean PCA's +/-6 range that arrow_scale=1.5 was tuned for.
# Arrows at raw envfit length were dwarfing the actual point/centroid cloud,
# forcing autoscale to squash all real data into an unreadable clump at the
# center (confirmed in first render — LP/UP centroid labels were invisible,
# overlapping each other near the origin). Rescale arrows to a FIXED
# fraction of the actual site/background data spread instead of a hardcoded
# constant tuned for a different ordination method.
data_spread <- max(abs(c(pcoa_plot_data$PC1, pcoa_plot_data$PC2)), na.rm = TRUE)
raw_arrow_len <- sqrt(top8_env$PC1^2 + top8_env$PC2^2)
target_arrow_len <- data_spread * 0.9   # arrows reach ~90% of data spread
arrow_scale <- target_arrow_len / max(raw_arrow_len)
top8_env <- top8_env %>% mutate(PC1_s = PC1*arrow_scale, PC2_s = PC2*arrow_scale)
message("Data spread: ", round(data_spread,3), " | arrow_scale recalculated: ", round(arrow_scale,3))

arrow_pal <- c(Structural="#C0392B", Climate="#E67E22", Spectral="#2980B9",
               Geochemistry="#16A085", Topographic="#8E44AD", Other="#7F8C8D")

pcoa_plot_data <- pcoa_scores %>%
  mutate(group = ifelse(point_type=="background","Background",period),
         group = factor(group, levels=c("Background","LP","MP","UP")))

lp_c <- period_centroids %>% filter(period=="LP")
mp_c <- period_centroids %>% filter(period=="MP")
up_c <- period_centroids %>% filter(period=="UP")

fig08 <- ggplot() +
  geom_point(data=pcoa_plot_data %>% filter(group=="Background"),
             aes(x=PC1,y=PC2), color="#CCCCCC", alpha=0.15, size=0.3) +
  geom_hline(yintercept=0, color="grey75", linewidth=0.2, linetype="dotted") +
  geom_vline(xintercept=0, color="grey75", linewidth=0.2, linetype="dotted") +
  geom_segment(data=top8_env, aes(x=0,y=0,xend=PC1_s,yend=PC2_s,color=category),
               arrow=arrow(length=unit(1.5,"mm"), type="open"), linewidth=0.45, alpha=0.85) +
  scale_color_manual(values=arrow_pal, name="Variable\nCategory") +
  geom_text_repel(data=top8_env, aes(x=PC1_s*1.15,y=PC2_s*1.15,label=variable,color=category),
                  size=2.2, fontface="bold", segment.size=0.2, segment.alpha=0.4,
                  segment.color="grey60", max.overlaps=20, box.padding=0.4,
                  point.padding=0.2, show.legend=FALSE, seed=42) +
  ggnewscale::new_scale_color() +
  geom_point(data=pcoa_plot_data %>% filter(group %in% c("LP","MP","UP")),
             aes(x=PC1,y=PC2,color=group), alpha=0.7, size=1.1) +
  stat_ellipse(data=pcoa_plot_data %>% filter(group %in% c("LP","MP","UP")),
               aes(x=PC1,y=PC2,color=group), level=0.68, linewidth=0.45, linetype="dashed") +
  geom_point(data=period_centroids, aes(x=PC1_mean,y=PC2_mean,color=period), size=5, shape=18) +
  geom_label_repel(data=period_centroids, aes(x=PC1_mean,y=PC2_mean,label=period,color=period),
                   size=3, fontface="bold", show.legend=FALSE, seed=42,
                   box.padding=0.6, min.segment.length=0, segment.size=0.3) +
  annotate("segment", x=lp_c$PC1_mean, y=lp_c$PC2_mean, xend=mp_c$PC1_mean, yend=mp_c$PC2_mean,
           arrow=arrow(length=unit(2,"mm"), type="closed"), color="grey10", linewidth=0.65) +
  annotate("segment", x=mp_c$PC1_mean, y=mp_c$PC2_mean, xend=up_c$PC1_mean, yend=up_c$PC2_mean,
           arrow=arrow(length=unit(2,"mm"), type="closed"), color="grey10", linewidth=0.65) +
  geom_label_repel(
    data = tibble(
      x = c(mean(c(lp_c$PC1_mean,mp_c$PC1_mean)), mean(c(mp_c$PC1_mean,up_c$PC1_mean))),
      y = c(mean(c(lp_c$PC2_mean,mp_c$PC2_mean)), mean(c(mp_c$PC2_mean,up_c$PC2_mean))),
      lab = c(
        paste0("d=",displacements$displacement[1]," [",displacements$boot_ci_low[1],"-",displacements$boot_ci_high[1],"]"),
        paste0("d=",displacements$displacement[2]," [",displacements$boot_ci_low[2],"-",displacements$boot_ci_high[2],"]")
      )
    ),
    aes(x=x,y=y,label=lab), size=2.0, color="grey15", fill="white",
    fontface="italic", box.padding=0.8, min.segment.length=0, segment.size=0.25,
    segment.color="grey50", seed=43
  ) +
  scale_color_manual(values=period_pal, labels=period_labels, name="Cultural\nPeriod") +
  coord_cartesian(xlim = c(-data_spread*1.1, data_spread*1.1),
                  ylim = c(-data_spread*1.6, data_spread*1.3)) +  # extra room bottom for climate arrows
  labs(title="Fig. 7  |  PCoA (Gower distance) of multivariate geoenvironmental space",
       subtitle=paste0("Mixed continuous+ordinal variables (n=30) | PCoA axis 1=",pc1_var,
                       "% | axis 2=",pc2_var,"% variance | Top 8 environmental vectors (vegan::envfit)"),
       x=paste0("PCoA Axis 1 (",pc1_var,"%)"), y=paste0("PCoA Axis 2 (",pc2_var,"%)"),
       caption="Centroid displacement labels show bootstrap 95% CI (1000 iter), not a bare descriptive distance | Ellipses = 68% CI | UP (n=27) — interpret cautiously, low coordinate precision (Fig. S1)") +
  theme_classic(base_size=8) +
  theme(plot.title=element_text(size=9,face="bold"), plot.subtitle=element_text(size=7,color="grey40"),
        plot.caption=element_text(size=5.5,color="grey50"), legend.position="right")

ggsave(file.path(fig_path,"Fig07_PCoA_Biplot.png"), fig08,
       width=180/25.4, height=160/25.4, dpi=300, units="in", bg="white")
message("Fig08 saved.")


# =============================================================
# FIGURE 9: RF Variable Importance — permutation + grouped (2 panels)
# =============================================================

message("\n--- Fig09: RF Variable Importance (permutation + grouped) ---")

cat_pal <- c(Structural="#C0392B", "Structural Geology"="#C0392B", Climate="#E67E22",
             Spectral="#2980B9", Pedogenic="#27AE60", Topographic="#8E44AD",
             Geochemistry="#16A085")

importance_clean <- importance_df %>% head(20) %>%
  mutate(label = fct_reorder(gsub("_"," ",variable), MeanDecreaseAccuracy))

panel_a <- ggplot(importance_clean, aes(x=MeanDecreaseAccuracy, y=label, fill=category)) +
  geom_col(width=0.72, color="white", linewidth=0.15) +
  geom_text(aes(label=paste0("#",rank)), x=0.3, hjust=0, size=2.0, color="white", fontface="bold") +
  scale_fill_manual(values=cat_pal, name="Category") +
  scale_x_continuous(expand=expansion(mult=c(0,0.15))) +
  labs(title="A. Individual variable importance",
       subtitle=paste0("Permutation (Mean Decrease Accuracy) | Top 20 of 30 | Naive OOB=",
                       round((1-rf_model$err.rate[rf_ntree,"OOB"])*100,1),
                       "% | Spatial block CV=", round(100*mean(spatial_cv_preds$observed==spatial_cv_preds$predicted),1),
                       "%, AUC=",round(spatial_cv_auc,3)),
       x="Mean Decrease Accuracy", y=NULL) +
  theme_springer() + theme(legend.position="bottom", axis.text.y=element_text(size=6.5))

grp_plot_data <- grp_importance %>%
  mutate(group_clean = fct_reorder(str_to_title(gsub("_"," ",group)), acc_drop_mean))

panel_b <- ggplot(grp_plot_data, aes(x=acc_drop_mean, y=group_clean)) +
  geom_col(fill="#34495E", width=0.65) +
  geom_errorbar(aes(xmin=pmax(0,acc_drop_mean-acc_drop_sd), xmax=acc_drop_mean+acc_drop_sd),
                orientation = "y", width=0.2, color="grey30") +
  labs(title="B. Grouped importance (correlated variables permuted jointly)",
       subtitle="Answers redundancy inflation: vegetation (NDVI/SAVI/MSAVI) ranks higher grouped than any single index alone",
       x="Mean accuracy drop (25 reps, +/-SD)", y=NULL) +
  theme_springer()

fig09 <- panel_a / panel_b + plot_layout(heights = c(1.3, 1))

ggsave(file.path(fig_path,"Fig08_RF_Variable_Importance.png"), fig09,
       width=180/25.4, height=200/25.4, dpi=300, units="in", bg="white")
message("Fig09 saved.")


# =============================================================
# FIGURE 10: Significance Matrix + Spatial Bootstrap Markers
# =============================================================

message("\n--- Fig10: Significance Matrix Heatmap (30 vars, bootstrap markers) ---")

var_categories <- setNames(
  c(rep("Spectral",7), "Topographic", rep("Pedogenic",5), rep("Climate",6),
    rep("Structural",5), rep("Geochemistry",6)),
  c("NDVI","NDWI","MNDWI","NDBI","BSI","SAVI","MSAVI","elevation",
    "soil_depth","soil_erosion","soil_productivity","soil_slope","soil_texture",
    "chelsa_bio01_lgm","chelsa_bio12_lgm","chelsa_bio15_lgm",
    "chelsa_bio01_modern","chelsa_bio12_modern","chelsa_bio15_modern",
    "dist_fault","dist_dyke","dist_lineament","dist_shear","dist_mineral",
    "geochem_stream_major_PC1","geochem_stream_major_PC2",
    "geochem_horizon_major_PC1","geochem_horizon_major_PC2",
    "geochem_regolith_major_PC1","geochem_regolith_major_PC2")
)
cat_pal2 <- c(Spectral="#2980B9", Topographic="#8E44AD", Pedogenic="#27AE60",
              Climate="#E67E22", Structural="#C0392B", Geochemistry="#16A085")

heat_long <- mw_table %>% filter(!is.na(p_value)) %>%
  mutate(log_p = -log10(pmax(p_value,1e-10)),
         sig_star = case_when(p_value<0.001~"***", p_value<0.01~"**", p_value<0.05~"*", TRUE~""),
         category = recode(variable, !!!var_categories),
         period = factor(period, levels=c("LP","MP","UP"))) %>%
  # NEW: mark structural cells with spatial-bootstrap-robustness status
  left_join(block_boot_table %>% select(variable, period, ci_excludes_zero),
            by = c("variable","period")) %>%
  mutate(bootstrap_marker = case_when(
    !is.na(ci_excludes_zero) & ci_excludes_zero  ~ "robust",
    !is.na(ci_excludes_zero) & !ci_excludes_zero ~ "not_robust",
    TRUE ~ NA_character_
  ))

heat_wide <- heat_long %>% select(variable, period, log_p) %>%
  pivot_wider(names_from=period, values_from=log_p, values_fill=0) %>% column_to_rownames("variable")
row_order <- hclust(dist(heat_wide))$order
ordered_vars <- rownames(heat_wide)[row_order]

heat_long <- heat_long %>%
  mutate(var_label = factor(variable, levels=ordered_vars),
         category = factor(category, levels=c("Spectral","Topographic","Pedogenic",
                                              "Climate","Structural","Geochemistry")))

fig10 <- ggplot(heat_long, aes(x=period, y=var_label)) +
  geom_tile(aes(fill=log_p), color="white", linewidth=0.8, width=0.92, height=0.92) +
  geom_text(aes(label=sig_star), size=3.0, color="grey15", fontface="bold", nudge_y=0.08) +
  # bootstrap-robustness marker: diamond for robust, cross for not-robust.
  # geom_point() has no nudge_y argument (caused a warning in the first
  # render, silently ignored) — offset the y-position in the data instead.
  geom_point(data = heat_long %>% filter(bootstrap_marker=="robust") %>%
               mutate(var_label_num = as.numeric(var_label) - 0.28),
             aes(x=period,y=var_label_num), shape=18, size=2.2, color="black") +
  geom_point(data = heat_long %>% filter(bootstrap_marker=="not_robust") %>%
               mutate(var_label_num = as.numeric(var_label) - 0.28),
             aes(x=period,y=var_label_num), shape=4, size=2.0, color="black") +
  scale_fill_gradientn(colours=c("#FFFFFF","#FFF7BC","#FEC44F","#FE9929","#D95F0E","#7F2704"),
                       values=scales::rescale(c(0,1.301,2,3,4,5)), limits=c(0,5), oob=scales::squish,
                       name=expression(-log[10](italic(p))),
                       breaks=c(0,1.301,2,3,4,5), labels=c("0 (ns)","0.05","0.01","0.001","0.0001","0.00001")) +
  # FIX: aes(fill=category) below is a DISCRETE scale sharing the "fill"
  # aesthetic with the CONTINUOUS scale_fill_gradientn above — ggplot2
  # 4.x rejects this outright (error confirmed in first render: "Discrete
  # value supplied to a continuous scale"). new_scale_fill() creates an
  # independent second fill scale for the category strip, same pattern
  # already used for the two color scales in Fig08.
  ggnewscale::new_scale_fill() +
  geom_tile(aes(x=0.28, fill=category), color="white", linewidth=0.5,
            width=0.18, height=0.92, show.legend=FALSE) +
  scale_fill_manual(values = cat_pal2) +
  scale_x_discrete(expand=expansion(add=c(0.6,0.5)),
                   labels=c(LP="Lower\nPalaeolithic", MP="Middle\nPalaeolithic", UP="Upper\nPalaeolithic")) +
  scale_y_discrete(expand=expansion(add=c(0.5,0.5))) +
  labs(title="Fig. 9  |  Geoenvironmental significance matrix (30 variables x 3 periods, 90 tests)",
       subtitle="Mann-Whitney U tests, site vs background per cultural period",
       x=NULL, y=NULL,
       caption="Diamond = survives spatial block bootstrap (Fig. 5/Discussion 7.2) | Cross = significant but NOT spatially robust | Stars: MW significance | Left strip = category | Rows: hierarchical clustering") +
  theme_classic(base_size=8) +
  theme(plot.title=element_text(size=9,face="bold"), plot.subtitle=element_text(size=7,color="grey40"),
        plot.caption=element_text(size=5.5,color="grey50"), axis.text.x=element_text(size=7.5,face="bold"),
        axis.text.y=element_text(size=6.2), axis.ticks=element_blank(), axis.line=element_blank(),
        panel.background=element_rect(fill="#FAFAFA",color=NA),
        panel.border=element_rect(color="grey30",fill=NA,linewidth=0.4), panel.grid=element_blank(),
        legend.position="right")

cat_legend <- tibble(category=names(cat_pal2)) %>% mutate(y=rev(seq_along(category))) %>%
  ggplot(aes(x=0.5,y=y,color=category)) + geom_point(size=4,shape=15) +
  geom_text(aes(label=category,x=0.65), hjust=0, size=2.3, color="grey20") +
  scale_color_manual(values=cat_pal2) + scale_x_continuous(limits=c(0.4,1.8)) +
  labs(title="Category") + theme_void(base_size=7) +
  theme(plot.title=element_text(size=7,face="bold"), legend.position="none")

fig10_final <- fig10 + cat_legend + plot_layout(widths=c(5,1))

ggsave(file.path(fig_path,"Fig09_Significance_Matrix_Heatmap.png"), fig10_final,
       width=180/25.4, height=210/25.4, dpi=300, units="in", bg="white")
message("Fig10 saved.")


# =============================================================
# FIGURE 11: Geoenvironmental Suitability Surface
# =============================================================

message("\n--- Fig11: Geoenvironmental Suitability Surface ---")

sites_coords_plot <- full_matrix_coords %>%
  filter(point_type=="site", period %in% c("LP","MP","UP")) %>%
  select(period, easting, northing)

png(file.path(fig_path,"Fig11_Geoenvironmental_Suitability.png"),
    width=fig_width, height=fig_height, units="mm", res=fig_dpi)
plot(attract_rast, main="Geoenvironmental Suitability (spatially-validated RF)", col=viridis(100), axes=TRUE)
points(sites_coords_plot$easting, sites_coords_plot$northing, pch=16, cex=0.5, col="red")
mtext(paste0("Naive OOB=", round((1-rf_model$err.rate[rf_ntree,"OOB"])*100,1),
             "% | Spatial block CV=", round(100*mean(spatial_cv_preds$observed==spatial_cv_preds$predicted),1),
             "%, AUC=", round(spatial_cv_auc,3), " | Relative/ordinal scores, not calibrated probabilities"),
      side=1, line=3, cex=0.55, col="grey40")
dev.off()
message("Fig11 saved.")


# =============================================================
# FINAL CHECK
# =============================================================

message("\n=============================================================")
message("SCRIPT 08 (REBUILT) COMPLETE")
message("=============================================================")

figs <- c("Fig03_Spectral_Boxplots.png","Fig05_Structural_Geology_Boxplots.png",
          "Fig06_Soil_Properties_Boxplots.png","Fig07_PCoA_Biplot.png",
          "Fig08_RF_Variable_Importance.png","Fig09_Significance_Matrix_Heatmap.png",
          "Fig11_Geoenvironmental_Suitability.png")

for (f in figs) {
  fp <- file.path(fig_path, f); ok <- file.exists(fp)
  size <- ifelse(ok, paste0(round(file.size(fp)/1e6,1)," MB"), "MISSING")
  message(ifelse(ok,"OK","XX"), "  ", f, "  [", size, "]")
}

message("\nRemoved from figure set entirely: (a) burial-depth surface map —")
message("Script 05 found the borehole data unusable, see Methods 5.6; (b) the")
message("original 'structural features map' slot (old Fig 6) — dropped per")
message("explicit instruction, chapter now holds at 11 figures total. All")
message("R-generated files renumbered accordingly: Fig06=soil, Fig07=PCoA,")
message("Fig08=RF importance, Fig09=significance matrix. QGIS figures needed:")
message("Fig01 (study area), Fig02 (S2A composite), Fig04 (geochemistry map),")
message("Fig10 (suitability surface), Fig11 (CHELSA context).")
message("Persistent Places (MULTI vs background) reported as Table only, not")
message("a dedicated figure — conserves figure budget.")
message("Saved to: ", fig_path)
message("=============================================================")