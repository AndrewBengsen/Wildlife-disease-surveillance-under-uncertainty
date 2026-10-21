## =============================================================================
## Reproducible workflow for: Bengsen et al., "Wildlife disease surveillance
## under uncertainty: an adaptive search-theoretic framework for early
## detection of transboundary animal diseases"

## PURPOSE
## This script demonstrates the process for creating a spatial representation
## of expected risk of a foot and mouth disease virus incursion in a feral pig
## population. Several base layers representing different aspects of risk are 
## scaled and combined to provide an overall risk rating from 0:1 for each of
## 1873 grid cells across the state of NSW, Australia. The numeric ratings are 
## then categorised into one of four risk classes. The process is deliberately
## simple and aims to reduce the risk of over-fitting a beautiful model to 
## data and processes that are all subject to high uncertainty. 

## DATA
## We do not have the authority to share the base feral pig, goat and deer data.
## For privacy reasons, we are not sharing point locations of piggeries,
## saleyards and feedlots. Instead, we have pre-loaded the processed output
## of these layers in `input/LS_grid_preload2.RDS` (object: LS_grid). We
## provide the processing code for the feral pig and point-hazard layers
## below for transparency; the code for deer and goats follows the same
## logic as for the feral pigs below and is not reproduced here.
##
## Graphical representations of the pig/deer/goat distribution data are
## available at:
## https://www.dpird.nsw.gov.au/dpi/biosecurity/invasive-plants-and-animals/pest-animals/distribution-maps

## REQUIRED INPUTS (all under ./input/)
## - LLS_Layer.shp                    Local Land Services region boundaries,
##                                    used to clip layers
## - LS_grid_preload2.RDS              Pre-processed hex grid (see DATA note)
##                                    Columns: cell, LS_pig, LS_goats, LS_deer,
##                                    LS_points, geometry
## - 5_Ct_2020_Da_clip_4326.tif       Cattle density, Gridded Livestock of the
##                                    World (FAO 2024), clipped to the LLS
##                                    region and reprojected to EPSG:4326
## - 5_Sh_2020_Da_clip_4326.tif       Sheep density, same source/processing
##
## GLW source: https://data.amerigeoss.org/dataset/9d1e149b-d63f-4213-978b-317a8eb42d02

## SOFTWARE
## Produced using R v4.3.1 (R Core Team 2023).
## Full session info is provided at the tail of this document

library(tidyverse)          
library(sf)                 
library(terra)              
library(raster)
library(here)                
library(scico)

## Dependencies -----------------------------------------------------------

terra_extract_mean <- function(r, polys, na.rm = TRUE) {
  out <- terra::extract(r, terra::vect(polys), fun = mean, na.rm = na.rm, ID = FALSE)
  as.numeric(out[[1]])
}

range01 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

## Load background data ========================================================

LS_grid <- readRDS(here("input", "LS_grid_preload2.RDS")) |>
  st_transform(4326)

# Cropping layer, for standardisation
lls <- st_read(here("input", "LLS_Layer.shp")) 

## Introduction risk weights (sum to 1) ========================================
LS_pig_w <- 0.51   # Feral pig density
LS_cattle_w <- 0.10   # Cattle density: primary host
LS_sheep_w <- 0.08   # Sheep density: maintenance host
LS_goat_w <- 0.08   # Goat density: maintenance host
LS_deer_w <- 0.08   # Deer density: maintenance host
LS_point_w <- 0.15   # Piggeries, saleyards and feedlots

## Feral pig density layer =====================================================

# The feral pig density layer was adapted from Crittle and Millyn (2023) but is
# not available to be reproduced here. An alternative version that can be
# georeferenced can be found at:
# https://www.dpird.nsw.gov.au/__data/assets/image/0006/1649733/Feral-pig-relative-abundance-2023-and-distribution-change-2020-2023.png

# Here, we've pre-loaded the output of the following process into LS_grid.
# This block is retained for transparency and is not run.
if(1==2){
  feralpigs <- raster(here("input", "Pig_2023_geotiff_3857.tif")) 
  feralpigs <- projectRaster(feralpigs, crs=4326)
  feralpigs <- crop(feralpigs, lls)
  
  # Rescale
  feralpigs[which(feralpigs[] == 0)] <- 1      # High density
  feralpigs[which(feralpigs[] == 115)] <- 0.66 # Med density
  feralpigs[which(feralpigs[] == 190)] <- 0.33 # Low density
  feralpigs[which(feralpigs[] == 245)] <- 0.33 # Assign present as low density
  feralpigs[which(feralpigs[] > 1)] <- 0       # Pigs not known to occur     
  
  feralpigs <- rast(feralpigs)  
  
  # Fit to grid using mean raster value for each grid cell
  LS_grid <- LS_grid |>
    mutate(LS_pig = terra_extract_mean(feralpigs, LS_grid))
}


## Expected cattle density =====================================================
# Adapted from Gridded Livestock of the World database (FAO 2024)
# https://data.apps.fao.org/catalog//iso/9d1e149b-d63f-4213-978b-317a8eb42d02

cattle <- rast(here("input", "5_Ct_2020_Da_clip_4326.tif"))

cattle[cattle > 100] <- 0  # convert urban/sentinel errors to 0

# Extract mean value for each risk map grid cell
cattle_grid <- terra_extract_mean(cattle, LS_grid)

# Scale 0:1
LS_grid <- LS_grid |>
  mutate(LS_cattle = range01(cattle_grid))


## Expected sheep density ======================================================
# Adapted from Gridded Livestock of the World database (FAO 2024)
# https://data.apps.fao.org/catalog//iso/9d1e149b-d63f-4213-978b-317a8eb42d02

sheep <- rast(here("input", "5_Sh_2020_Da_clip_4326.tif"))

sheep[sheep == 373.09647] <- 0    # convert urban error (Broken Hill) to 0

sheep_grid <- terra_extract_mean(sheep, LS_grid)

LS_grid <- LS_grid |>
  mutate(LS_sheep = range01(sheep_grid))

## Point hazards: piggeries, saleyards and feedlots ============================
# As noted in the header, we do not provide point locations for these features
# for privacy reasons. We provide processing code to combine the three point
# layers into a count of point hazards per cell, for transparency. This block
# is retained for documentation and is not run.

if(1 == 2){
 feedlots  <- read.csv(here("input", "NSW_feedlots.csv"))
 piggeries <- read.csv(here("input", "NSW_piggeries.csv"))
 saleyards <- read.csv(here("input", "saleyards_nsw_4326_2023-05.csv"))

 risk_points <- bind_rows(feedlots[, c("category", "lat", "lon")],
                          piggeries[, c("category", "lat", "lon")],
                          saleyards[, c("category", "lat", "lon")]) |>
    st_as_sf(coords = c("lon", "lat"), crs = st_crs(4326))

 points_3308  <- st_transform(risk_points, 3308)
 LS_grid_3308 <- st_transform(LS_grid, 3308)

 LS_grid$risk_points <- lengths(st_intersects(LS_grid_3308, points_3308))
 LS_grid$LS_points   <- range01(LS_grid$risk_points)
}


## Combine layers ==============================================================

LSw <- st_drop_geometry(LS_grid) |>
  transmute(LS_pig = LS_pig * LS_pig_w,
            LS_cattle = LS_cattle * LS_cattle_w,
            LS_sheep = LS_sheep * LS_sheep_w,
            LS_goats = LS_goats * LS_goat_w,
            LS_deer = LS_deer * LS_deer_w,
            LS_points = LS_points * LS_point_w)

# Sum total weighted risk for each cell and rescale to 0-1
LS_grid <- LS_grid |>
  mutate(LS_tot = range01(rowSums(LSw, na.rm = TRUE)))

LS_grid_cat <- LS_grid

## Categorise each cell by risk weight =========================================
# Target: mean(high) : mean(moderate) : mean(low) : mean(very low)
#         approx. 9 : 7 : 3 : 1
# Solve for the three quantile cut-points of LS_tot that come closest to
# producing group means in the ratio above, using a bounded optimiser over 
# the cut quantiles.

risk_weights <- c(high = 9, moderate = 7, low = 3, very_low = 1)

category_breaks <- function(values, grp) {
  cut(values, breaks = grp, labels = c("very low", "low", "moderate", "high"),
      right = FALSE, include.lowest = TRUE)
}

objective <- function(p, values) {
  # ensure p[1] > p[2] > p[3], so that cuts[1] > cuts[2] > cuts[3]
  p <- sort(p, decreasing = TRUE) 
  
  # convert probabilities to actual cut values
  cuts <- quantile(values, probs = p, na.rm = T)
  
  # Reject this parameter set if quantile() failed
  # or if > 1 cut point lands on the same value
  # The large penalty pushes optim() away without killing it
  if (is.null(cuts) || anyDuplicated(cuts)) return(10^6) 
  
  # Four breakpoints for cut()
  grp_breaks <- c(-Inf, cuts[3], cuts[2], cuts[1], Inf)
  
  # Assign each cell's LS_tot value to a category
  grp <- category_breaks(values, grp_breaks)
  
  # Mean LS_tot in each category, to check performance
  means <- tapply(values, grp, mean)
  
  # Penalise any empty or broken categories to avoid failure
  if (any(!is.finite(means))) return(10^6)
  
  # Extract reference mean from "very low" category
  vh <- means["very low"]

  # Penalise any zero or non-finite mean for vh that will cause a failure
  if (vh == 0) return(10^6)
  if (!is.finite(vh)) return(10^6)  
  
  # Express target weights for each category as a ratio relative to "very low"
  target_ratio <- risk_weights[c("high", "moderate", "low")] / risk_weights["very_low"]
  
  # Objective: sum of squares between achieved and target mean ratios
  # optim() will search for the cut quantiles (p) that minimise the error
  # and come closest to the target 9:7:3:1 relationship
  # err == 0 is a perfect match
  err <- sum((c(means["high"], means["moderate"], means["low"]) / vh - target_ratio)^2)
  if (!is.finite(err)) return(10^6)   
  err
}

# Nelder-Mead has no bounds, so optimise over unconstrained values (theta) and
# then map them back to bounded quantile probabilities (p) within the objective.

lower <- c(0.30, 0.05, 0.001)
upper <- c(0.99, 0.80, 0.500)

# Optimiser starting values
props <- risk_weights / sum(risk_weights)
p1 <- 1 - props["high"]
p2 <- 1 - (props["high"] + props["moderate"])
p3 <- 1 - (props["high"] + props["moderate"] + props["low"])
par_init <- c(p1, p2, p3)

## Optimise cut-point quantiles

# p (bounded) -> theta (unconstrained): rescale to (0,1), then logit
p_to_theta <- function(p) qlogis((p - lower) / (upper - lower))

# theta (unconstrained) -> p (bounded): inverse logit, then rescale to bounds
theta_to_p <- function(theta) lower + (upper - lower) * plogis(theta)

# Wrapper so optim() searches over theta; objective() still receives
# probabilities (and sorts them internally, so ordering is always valid)
objective_theta <- function(theta, values) objective(theta_to_p(theta), values)

# Starting values on the theta scale
theta_init <- p_to_theta(as.numeric(par_init))

fit_nm <- function(theta0, values) {
  optim(par     = theta0,
        fn      = objective_theta,
        values  = values,
        method  = "Nelder-Mead",
        control = list(maxit = 2000, reltol = 1e-10))
}

# The objective is piecewise constant (cut points jump between data values),
# so a single Nelder-Mead run can stall on a plateau. Run from the original
# starting point plus several jittered starts, and keep the best fit.
set.seed(5678)
n_starts <- 5
starts <- rbind(theta_init,
                t(replicate(n_starts - 1, theta_init + rnorm(3, sd = 1))))

fits <- lapply(seq_len(nrow(starts)),
               function(i) fit_nm(starts[i, ], LS_grid$LS_tot))

errs <- vapply(fits, function(f) f$value, numeric(1))
summary(errs)                       # spread across starts; 1e6 = penalty plateau
best <- fits[[which.min(errs)]]
best$value                          # 0 would be a perfect match to 9:7:3:1

# Back-transform to quantile probabilities, ordered so p_opt[1] > p_opt[2] > p_opt[3]
p_opt <- sort(theta_to_p(best$par), decreasing = TRUE)
p_opt

cuts <- quantile(LS_grid_cat$LS_tot, probs = p_opt, na.rm = TRUE)

# It's a tricky optimisation: check performance
grp_check <- category_breaks(LS_grid$LS_tot, c(-Inf, cuts[3], cuts[2], cuts[1], Inf))
achieved_means <- tapply(LS_grid$LS_tot, grp_check, mean)

# Achieved risk ratios relative to 'very low' (target 1 : 3 : 7 : 9) 
round(achieved_means / achieved_means["very low"], 2)

# Final categories
LS_grid <- LS_grid_cat |>
  mutate(
    risk_category = category_breaks(LS_tot, c(-Inf, cuts[3], cuts[2], cuts[1], Inf)),
    risk_category = factor(risk_category, 
                           levels = c("high", "moderate", "low", "very low")))

## Plots ========================================================================

LS_grid |> 
  st_drop_geometry() |>
  group_by(risk_category) |>
  summarise(n = n()) |>
  mutate(p = n/1873) |>
  knitr::kable(digits=2)

ggplot() +
  geom_sf(aes(fill = LS_tot), col = NA, data = LS_grid, show.legend = TRUE) +
  scico::scale_fill_scico(palette = "lajolla", direction = -1, name = "Risk") +
  theme_void() +
  labs(title = "Total expected introduction risk (continuous, 0:1)")


ggplot() +
  geom_sf(aes(fill = risk_category), col = NA, data = LS_grid, show.legend = TRUE) +
  scale_fill_manual(values = c("red4", "orange3", "khaki", "grey"), name = "Risk") +
  theme_void() +
  labs(title = "Total expected introduction risk (categorical)")

# Save for use in `BST_DBSCAN_synthetic_example_v1-1.R`
st_write(LS_grid, here("input", "FMDV_risk_4326_v1-3.shp"), append = FALSE)


## Sensitivity analysis ========================================================

# How much does the final risk map change if a layer weight is wrong? Each of
# the five secondary layers is varied in turn across nine weights (0.1 to 0.5),
# with all other weights held at their base values. The feral pig weight is held
# constant (see the pig-vs-total comparison at the end). For each scenario the
# total risk is rescaled to 0:1 and compared with the original map using RMSE,
# Spearman correlation and the coefficient of variation (CV).

LS_grid <- LS_grid |>
  arrange(cell)

# Base weights
base_w <- c(LS_pig = LS_pig_w, LS_cattle = LS_cattle_w, LS_sheep = LS_sheep_w,
            LS_goats = LS_goat_w, LS_deer = LS_deer_w, LS_points = LS_point_w)

# Matrix of scaled layers (one column per layer, same order as base_w)
layer_mat <- LS_grid |>
  st_drop_geometry() |>
  dplyr::select(all_of(names(base_w))) |>
  as.matrix()

# Layers to vary (column name = display label). Order sets the panel order below
layer_labels <- c(LS_goats = "Goats", LS_cattle = "Cattle", LS_deer = "Deer",
                  LS_sheep = "Sheep", LS_points = "Point hazards")

weight_seq <- round(seq(0.1, 0.5, length.out = 9), 2)

# Total risk (0:1) for a vector of layer weights.
# Same calculation as in 'Combine layers' above
risk_total <- function(w) {
  range01(rowSums(sweep(layer_mat, 2, w, "*"), na.rm = TRUE))
}

# One scenario: `weight` on `layer`, base weights on everything else
sens_scenario <- function(layer, weight) {
  tibble(cell = LS_grid$cell,
         LS_ref = LS_grid$LS_tot,    # original risk map
         LS_tot = risk_total(replace(base_w, layer, weight)))
}

# Every layer x weight scenario, one row per cell per scenario
sens <- expand_grid(layer = names(layer_labels), weight = weight_seq) |>
  mutate(res = map2(layer, weight, sens_scenario)) |>
  unnest(res)

# Compare each scenario map with the original map
sens_tab <- sens |>
  group_by(layer, weight) |>
  summarise(RMSE = sqrt(mean((LS_tot - LS_ref)^2, na.rm = TRUE)),
            rho = cor(LS_tot, LS_ref, method = "spearman", use = "complete.obs"),
            cv = sd(LS_tot, na.rm = TRUE) / mean(LS_tot, na.rm = TRUE),
            .groups = "drop") |>
  mutate(layer = factor(layer, levels = names(layer_labels), labels = layer_labels))

print(sens_tab, n = Inf)

# Layer and weight with the greatest impact on the risk map (largest RMSE)
slice_max(sens_tab, RMSE)

# Agreement with the original map at the highest weight, for each layer
filter(sens_tab, weight == max(weight))

## Sensitivity plots -----------------------------------------------------------

# Maps of total risk at each weight for one layer (`lyr` = layer column name)
plot_sens_map <- function(lyr) {
  LS_grid |>
    dplyr::select(cell) |>
    inner_join(filter(sens, layer == lyr), by = "cell") |>
    ggplot() +
    geom_sf(aes(fill = LS_tot), col = NA, show.legend = FALSE) +
    scico::scale_fill_scico(palette = "lajolla", direction = -1) +
    facet_wrap(~weight) +
    theme_void(base_size = 10)
}

sens_bar <- ggplot(sens_tab) +
  geom_col(aes(RMSE, factor(weight))) +
  scale_x_continuous(limits = c(0, 0.15)) +
  lemon::facet_rep_wrap(~layer) +
  labs(y = "Weight") +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(size = 7),
        panel.spacing.x = unit(10, "pt"),
        axis.text = element_text(colour = 1))
sens_bar

# RMSE panel (a) alongside the risk maps for one layer (b)
combine_sens <- function(lyr) {
  cowplot::plot_grid(sens_bar, plot_sens_map(lyr), ncol = 2,
                     labels = c("a)", "b)"), rel_widths = c(0.5, 0.5),
                     align = "v", axis = "b", rel_heights = c(1, 0.95))
}

dir.create(here("Figures"), showWarnings = FALSE)
ggsave(here("Figures", "combined_sensitivity.png"), combine_sens("LS_cattle"),
       height = 8, width = 18, units = "cm")

## Is the total risk map just the feral pig layer? -----------------------------
# Feral pigs carry the largest weight, so check how closely the total risk
# tracks the pig layer alone

range(LS_grid$LS_pig, na.rm = TRUE)
range(LS_grid$LS_tot, na.rm = TRUE)

df <- st_drop_geometry(LS_grid) |>
  filter(!is.na(LS_pig))

# Very highly correlated under current weightings (0.98, 0.94)
cor(df$LS_pig, df$LS_tot, method = "pearson")
cor(df$LS_pig, df$LS_tot, method = "spearman")

# Lin's concordance correlation coefficient (CCC) is also high (0.96)
DescTools::CCC(df$LS_pig, df$LS_tot)$rho.c


## Session info ================================================================

# sessioninfo::session_info()
# ─ Session info ───────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 (2023-06-16 ucrt)
# os       Windows 11 x64 (build 26200)
# system   x86_64, mingw32
# ui       RStudio
# language (EN)
# collate  English_Australia.utf8
# ctype    English_Australia.utf8
# tz       Australia/Sydney
# date     2026-09-29
# rstudio  2025.05.0+496 Mariposa Orchid (desktop)
# pandoc   NA
# 
# ─ Packages ───────────────────────────────────────────────────────────────────
# package     * version date (UTC) lib source
# class         7.3-22  2023-05-03 [1] CRAN (R 4.3.1)
# classInt      0.4-10  2023-09-05 [1] CRAN (R 4.3.2)
# cli           3.6.2   2023-12-11 [1] CRAN (R 4.3.2)
# codetools     0.2-19  2023-02-01 [1] CRAN (R 4.3.1)
# colorspace    2.1-0   2023-01-23 [1] CRAN (R 4.3.2)
# DBI           1.2.2   2024-02-16 [1] CRAN (R 4.3.2)
# dplyr       * 1.1.4   2023-11-17 [1] CRAN (R 4.3.2)
# e1071         1.7-14  2023-12-06 [1] CRAN (R 4.3.2)
# farver        2.1.1   2022-07-06 [1] CRAN (R 4.3.2)
# forcats     * 1.0.0   2023-01-29 [1] CRAN (R 4.3.2)
# generics      0.1.4   2025-05-09 [1] CRAN (R 4.3.1)
# ggplot2     * 3.5.2   2025-04-09 [1] CRAN (R 4.3.1)
# glue          1.7.0   2024-01-09 [1] CRAN (R 4.3.2)
# gtable        0.3.4   2023-08-21 [1] CRAN (R 4.3.2)
# here        * 1.0.1   2020-12-13 [1] CRAN (R 4.3.3)
# hms           1.1.3   2023-03-21 [1] CRAN (R 4.3.2)
# KernSmooth    2.23-21 2023-05-03 [1] CRAN (R 4.3.1)
# labeling      0.4.3   2023-08-29 [1] CRAN (R 4.3.1)
# lifecycle     1.0.5   2026-01-08 [1] CRAN (R 4.3.1)
# lubridate   * 1.9.3   2023-09-27 [1] CRAN (R 4.3.2)
# magrittr      2.0.3   2022-03-30 [1] CRAN (R 4.3.2)
# munsell       0.5.0   2018-06-12 [1] CRAN (R 4.3.2)
# pillar        1.11.1  2025-09-17 [1] CRAN (R 4.3.1)
# pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.3.2)
# pkgload       1.3.4   2024-01-16 [1] CRAN (R 4.3.2)
# proxy         0.4-27  2022-06-09 [1] CRAN (R 4.3.2)
# purrr       * 1.0.2   2023-08-10 [1] CRAN (R 4.3.2)
# R6            2.6.1   2025-02-15 [1] CRAN (R 4.3.3)
# Rcpp          1.0.12  2024-01-09 [1] CRAN (R 4.3.2)
# readr       * 2.1.5   2024-01-10 [1] CRAN (R 4.3.2)
# rlang         1.1.3   2024-01-10 [1] CRAN (R 4.3.2)
# rprojroot     2.0.4   2023-11-05 [1] CRAN (R 4.3.2)
# rstudioapi    0.15.0  2023-07-07 [1] CRAN (R 4.3.2)
# s2            1.1.6   2023-12-19 [1] CRAN (R 4.3.2)
# scales        1.3.0   2023-11-28 [1] CRAN (R 4.3.2)
# scico       * 1.5.0   2023-08-14 [1] CRAN (R 4.3.3)
# sessioninfo   1.2.2   2021-12-06 [1] CRAN (R 4.3.2)
# sf          * 1.0-15  2023-12-18 [1] CRAN (R 4.3.2)
# stringi       1.8.3   2023-12-11 [1] CRAN (R 4.3.2)
# stringr     * 1.5.1   2023-11-14 [1] CRAN (R 4.3.2)
# terra       * 1.8-29  2025-02-26 [1] CRAN (R 4.3.3)
# tibble      * 3.2.1   2023-03-20 [1] CRAN (R 4.3.2)
# tidyr       * 1.3.1   2024-01-24 [1] CRAN (R 4.3.2)
# tidyselect    1.2.1   2024-03-11 [1] CRAN (R 4.3.3)
# tidyverse   * 2.0.0   2023-02-22 [1] CRAN (R 4.3.2)
# timechange    0.3.0   2024-01-18 [1] CRAN (R 4.3.2)
# tzdb          0.4.0   2023-05-12 [1] CRAN (R 4.3.2)
# units         0.8-5   2023-11-28 [1] CRAN (R 4.3.2)
# utf8          1.2.4   2023-10-22 [1] CRAN (R 4.3.2)
# vctrs         0.6.5   2023-12-01 [1] CRAN (R 4.3.2)
# withr         3.0.2   2024-10-28 [1] CRAN (R 4.3.3)
# wk            0.9.1   2023-11-29 [1] CRAN (R 4.3.2)
# 
# [1] C:/Program Files/R/R-4.3.1/library