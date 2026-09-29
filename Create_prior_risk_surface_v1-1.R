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
## of these layers in `input/LS_grid_preload.RDS` (object: LS_grid). We
## provide the processing code for the feral pig and point-hazard layers
## below for transparency; the code for deer and goats follows the same
## logic as cattle/sheep below and is not reproduced here.
##
## Graphical representations of the pig/deer/goat distribution data are
## available at:
## https://www.dpird.nsw.gov.au/dpi/biosecurity/invasive-plants-and-animals/pest-animals/distribution-maps

## REQUIRED INPUTS (all under ./input/)
## - LLS_Layer.shp                    Local Land Services region boundaries
## - Buffer_NSW_50km_4326.shp         NSW + 50 km buffer, used to clip layers
## - LS_grid_preload.RDS              Pre-processed hex grid (see DATA note)
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

LS_grid <- readRDS(here("input", "LS_grid_preload.RDS")) |>
  st_transform(4326)

## Introduction risk weights (sum to 1) ========================================
LS_pig_w    <- 0.51   # Feral pig density
LS_cattle_w <- 0.10   # Cattle density: primary host
LS_sheep_w  <- 0.08   # Sheep density: maintenance host
LS_goat_w   <- 0.08   # Goat density: maintenance host
LS_deer_w   <- 0.08   # Deer density: maintenance host
LS_point_w  <- 0.15   # Piggeries, saleyards and feedlots

## Feral pig density layer =====================================================

# The feral pig density layer was adapted from Crittle and Millyn (2023) but is
# not available to be reproduced here. An alternative version that can be
# georeferenced can be found at:
# https://www.dpird.nsw.gov.au/__data/assets/image/0006/1649733/Feral-pig-relative-abundance-2023-and-distribution-change-2020-2023.png

# Here, we've pre-loaded the output of the following process into LS_grid.
# This block is retained for transparency and is not run.
if(1==2){
 feralpigs <- rast(here("input", "Pig_2023_geotiff_3857.tif"))
 feralpigs <- project(feralpigs, "EPSG:4326")
 feralpigs <- crop(feralpigs, lls)
 feralpigs <- project(feralpigs, "EPSG:4326", method = "near")

 # Recode class codes to a 0-1 scale (high -> 1, not-known-to-occur -> 0)
 feralpigs[feralpigs == 0]   <- 1     # High density
 feralpigs[feralpigs == 115] <- 0.66  # Medium density
 feralpigs[feralpigs == 190] <- 0.33  # Low density
 feralpigs[feralpigs == 245] <- 0.33  # Present, assigned low density
 feralpigs[feralpigs == 225] <- 0     # Not known to occur
 feralpigs[feralpigs > 1]    <- 0     # Remove any other/unexpected values

 # Fit to grid using mean raster value for each grid cell
 LS_grid <- LS_grid %>%
   mutate(LS_pig = terra_extract_mean(feralpigs, LS_grid))
}


## Expected cattle density =====================================================
# Adapted from Gridded Livestock of the World database (FAO 2024)
# https://data.amerigeoss.org/dataset/9d1e149b-d63f-4213-978b-317a8eb42d02

cattle <- rast(here("input", "5_Ct_2020_Da_clip_4326.tif"))

cattle[cattle > 6500] <- 0  # convert urban/sentinel errors to 0

cattle_range <- minmax(cattle, compute = TRUE)
cattle <- (cattle - cattle_range[1]) / (cattle_range[2] - cattle_range[1])

# Fit to grid using mean raster value for each grid cell
LS_grid <- LS_grid %>%
  mutate(LS_cattle = terra_extract_mean(cattle, LS_grid))


## Expected sheep density ======================================================
# Adapted from Gridded Livestock of the World database (FAO 2024)
# https://data.amerigeoss.org/dataset/9d1e149b-d63f-4213-978b-317a8eb42d02

sheep <- rast(here("input", "5_Sh_2020_Da_clip_4326.tif"))

sheep[sheep > 6500] <- 0
sheep_range <- minmax(sheep, compute = TRUE)
sheep <- (sheep - sheep_range[1]) / (sheep_range[2] - sheep_range[1])

LS_grid <- LS_grid %>%
  mutate(LS_sheep = terra_extract_mean(sheep, LS_grid))

## Point hazards: piggeries, saleyards and feedlots ============================
# As noted in the header, we do not provide point locations for these features
# for privacy reasons. We provide processing code to combine the three point
# layers into a count of point hazards per cell, for transparency. This block
# is retained for documentation and is not run.

if(1 == 2){
 feedlots  <- read.csv(here("input", "NSW_feedlots.csv"))
 piggeries <- read.csv(here("input", "NSW_piggeries.csv"))
 saleyards <- read.csv(here("input", "NSW_saleyards.csv"))

 risk_points <- bind_rows(feedlots[, c("category", "lat", "lon")],
                          piggeries[, c("category", "lat", "lon")],
                          saleyards[, c("category", "lat", "lon")]) %>%
    st_as_sf(coords = c("lon", "lat"), crs = st_crs(4326))

 points_3308  <- st_transform(risk_points, 3308)
 LS_grid_3308 <- st_transform(LS_grid, 3308)

 LS_grid$risk_points <- lengths(st_intersects(LS_grid_3308, points_3308))
 LS_grid$LS_points   <- range01(LS_grid$risk_points)
}


## Combine layers ==============================================================

# Check completeness of each cell - 
# some boundary cells may be missing raster risk layers

component_cols <- c("LS_pig", "LS_cattle", "LS_sheep", "LS_goats", "LS_deer", "LS_points")
n_incomplete <- sum(!complete.cases(st_drop_geometry(LS_grid)[component_cols]))
if (n_incomplete > 0) {
  message(n_incomplete, " of ", nrow(LS_grid),
          " grid cells have at least one missing risk component and will be",
          " scored using only their available layers (rowSums na.rm = TRUE).")
}

LSw <- st_drop_geometry(LS_grid) %>%
  transmute(LS_pig = LS_pig * LS_pig_w,
            LS_cattle = LS_cattle * LS_cattle_w,
            LS_sheep = LS_sheep * LS_sheep_w,
            LS_goats = LS_goats * LS_goat_w,
            LS_deer = LS_deer * LS_deer_w,
            LS_points = LS_points * LS_point_w)

# Sum total weighted risk for each cell and rescale to 0-1
LS_grid <- LS_grid %>%
  mutate(LS_tot = range01(rowSums(LSw, na.rm = TRUE)))

LS_grid_cat <- LS_grid

## Categorise each cell by risk weight =========================================
# Target: mean(high) : mean(moderate) : mean(low) : mean(very low)
#         approx. 9 : 7 : 3 : 1
# Solve for the three quantile cut-points of LS_tot that come closest to
# producing group means in the ratio above, using an unconstrained-but-bounded
# optimiser over the cut quantiles.

risk_weights <- c(high = 9, moderate = 7, low = 3, very_low = 1)

category_breaks <- function(values, grp) {
  cut(values, breaks = grp, labels = c("very low", "low", "moderate", "high"),
      right = FALSE, include.lowest = TRUE)
}

objective <- function(p, values) {
  p <- sort(p, decreasing = TRUE)
  cuts <- tryCatch(quantile(values, probs = p, na.rm = T), error = function(e) NULL)
  if (is.null(cuts) || anyDuplicated(cuts)) return(1e6) 
  grp_breaks <- c(-Inf, cuts[3], cuts[2], cuts[1], Inf)
  grp <- category_breaks(values, grp_breaks)
  means <- tapply(values, grp, mean)
  
  if (any(is.na(means)) || any(!is.finite(means))) return(1e6)
  vh <- means["very low"]
  if (!is.finite(vh) || vh == 0) return(1e6)  
  
  target_ratio <- risk_weights[c("high", "moderate", "low")] / risk_weights["very_low"]
  err <- sum((c(means["high"], means["moderate"], means["low"]) / vh - target_ratio) ^ 2)
  err
}

# Optimiser starting values
props <- risk_weights / sum(risk_weights)

p1 <- 1 - props["high"]
p2 <- 1 - (props["high"] + props["moderate"])
p3 <- 1 - (props["high"] + props["moderate"] + props["low"])
par_init <- c(p1, p2, p3)

res <- optim(par = par_init,
             fn = objective,
             values = LS_grid$LS_tot,
             method = "L-BFGS-B",
             lower = c(0.30, 0.05, 0.001),
             upper = c(0.99, 0.80, 0.500))

p_opt <- sort(res$par, decreasing = TRUE)
cuts  <- quantile(LS_grid_cat$LS_tot, probs = p_opt, na.rm = TRUE)

grp_check <- category_breaks(LS_grid$LS_tot, c(-Inf, cuts[3], cuts[2], cuts[1], Inf))
achieved_means <- tapply(LS_grid$LS_tot, grp_check, mean)
message("Achieved risk ratios relative to 'very low' (target 9 / 7 / 3 / 1):")
print(round(achieved_means / achieved_means["very low"], 2))

# Final categories
LS_grid <- LS_grid_cat %>%
  mutate(
    risk_category = category_breaks(LS_tot, c(-Inf, cuts[3], cuts[2], cuts[1], Inf)),
    risk_category = factor(risk_category, 
                           levels = c("high", "moderate", "low", "very low")))

LS_grid <- LS_grid_cat |>
  mutate(risk_category = case_when(
    LS_tot >= cuts[1] ~ "high",
    LS_tot >= cuts[2] ~ "moderate",
    LS_tot >= cuts[3] ~ "low",
    TRUE ~ "very low")) |>
  mutate(risk_category = factor(
    risk_category, levels = c("high", "moderate", "low", "very low")))

LS_grid |>
  st_drop_geometry() |>
  group_by(risk_category) |>
  summarise(n = n())

## Plots ========================================================================

ggplot() +
  geom_sf(aes(fill = LS_tot), col = NA, data = LS_grid, show.legend = TRUE) +
  scico::scale_fill_scico(palette = "lajolla", direction = -1, name = "Risk") +
  theme_void() +
  labs(title = "Total expected introduction risk (continuous, 0-1)")


ggplot() +
  geom_sf(aes(fill = risk_category), col = NA, data = LS_grid, show.legend = TRUE) +
  scale_fill_manual(values = c("red4", "orange3", "khaki", "grey"), name = "Risk") +
  theme_void() +
  labs(title = "Total expected introduction risk (categorical)")

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