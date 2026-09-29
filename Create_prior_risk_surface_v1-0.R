## =============================================================================
## Reproducible workflow for: Bengsen et al., "Wildlife disease surveillance
## under uncertainty: an adaptive search-theoretic framework for early
## detection of transboundary animal diseases"

## DATA
## We do not have the authority to share the base feral pig, goat and deer data.
## For privacy reasons, we are not sharing point locations of piggeries and 
## feedlots. Instead, we have pre-loaded the processed output in the `LS_grid` 
## We provide the processing code for the feral pig layer below. 
## The code for deer and goats follows the same logic and is not reproduced here
## Graphical representations of these data are available at:
## https://www.dpird.nsw.gov.au/dpi/biosecurity/invasive-plants-and-animals/pest-animals/distribution-maps



## REQUIRED INPUTS

## SOFTWARE
## Produced using R v4.3.1 (R Core Team 2024). 
## Full session info at the tail of this workflow

library(tidyverse)          
library(sf)                 
library(terra)              
library(here)                
library(scico)


## Dependencies

terra_extract_mean <- function(r, polys, na.rm = TRUE) {
  out <- terra::extract(r, terra::vect(polys), fun = mean, na.rm = na.rm, ID = FALSE)
  as.numeric(out[[1]])
}

range01 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

## Load background data =====================================================

LS_grid <- readRDS("input/LS_grid_preload.RDS") |>
  st_transform(4326)


## Introduction risk weights (sum to 1) ========================================
LS_pig_w <- 0.51         # Feral pig density
LS_cattle_w <- 0.1       # Cattle density  - primary host
LS_sheep_w <- 0.08        # Sheep density - maintenance host
LS_goat_w <- 0.08         # Goat density  - maintenance host
LS_deer_w <- 0.08         # Deer density  - maintenance host
LS_point_w <- 0.15        # Piggeries, saleyards and feedlots


## Feral pig density layer =====================================================

# The feral pig density layer was adapted from Crittle and Millyn (2023) but is 
# not available to be reproduced here. An alternative version that can be 
# georeferenced can be found at:
# https://www.dpird.nsw.gov.au/__data/assets/image/0006/1649733/Feral-pig-relative-abundance-2023-and-distribution-change-2020-2023.png

# Here, we have pre-loaded the output of the following process in LS_grid

#feralpigs <- raster("input/Pig_2023_geotiff_3857.tif") 
#feralpigs <- projectRaster(feralpigs, crs=4326)
#feralpigs <- crop(feralpigs, lls)

# Scale band 1 
# feralpigs[which(feralpigs[] == 0)] <- 1      # High density
# feralpigs[which(feralpigs[] == 115)] <- 0.66 # Med density
# feralpigs[which(feralpigs[] == 190)] <- 0.33 # Low density
# feralpigs[which(feralpigs[] == 245)] <- 0.33 # Assign present as low density
# feralpigs[which(feralpigs[] == 225)] <- 0    # Pigs not known to occur
# feralpigs[which(feralpigs[] > 1)] <- 0       # remove other values
# 
# # Fit to grid using mean raster value for each grid cell
# LS_grid <- LS_grid |>
#   mutate(LS_pig = raster_extract(x = feralpigs, y = LS_grid, fun = mean, na.rm=T))


## Expected cattle density =====================================================
# Adapted from Gridded Livestock of the World database (FAO 2024)
# https://data.amerigeoss.org/dataset/9d1e149b-d63f-4213-978b-317a8eb42d02

cattle <- terra::rast(here("input", "5_Ct_2020_Da_clip_4326.tif"))

# Rescale to 0–1
cattle[cattle > 6500] <- 0   # convert urban errors to 0
cattle <- (cattle - minmax(cattle)[1]) / (minmax(cattle)[2] - minmax(cattle)[1])

# Fit to grid using mean raster value for each grid cell
LS_grid <- LS_grid |>
  mutate(LS_cattle = terra_extract_mean(cattle, LS_grid))


## Expected sheep density ======================================================
# Adapted from Gridded Livestock of the World database (FAO 2024)
# https://data.amerigeoss.org/dataset/9d1e149b-d63f-4213-978b-317a8eb42d02

sheep <- terra::rast(here("input", "5_Sh_2020_Da_clip_4326.tif"))

sheep[sheep > 6500] <- 0
sheep_range <- terra::minmax(sheep, compute = TRUE)
sheep <- (sheep - sheep_range[1]) / (sheep_range[2] - sheep_range[1])

LS_grid <- LS_grid %>%
  mutate(LS_sheep = terra_extract_mean(sheep, LS_grid))

# Rescale to 0–1
sheep[sheep > 6500] <- 0   # convert urban errors to 0
sheep <- (sheep - minmax(sheep)[1]) / (minmax(sheep)[2] - minmax(sheep)[1])

# Fit to grid using mean raster value for each grid cell
LS_grid <- LS_grid |>
  mutate(LS_sheep = terra_extract_mean(sheep, LS_grid))


## Point hazards: piggeries, saleyards and feedlots ============================
# As noted in the header, we do not provide point locations for these features
# for privacy reasons. We provide processing code to combine the three separate
# point features into a count of point hazards per cell, for transparency.

# feedlots  <- read.csv("input/NSW_feedlots.csv")
# piggeries <- read.csv(...)
# saleyards <- read.csv(...) 

# risk_points <- rbind(feedlots[,c("category", "lat", "lon")], 
#                      piggeries[,c("category", "lat", "lon")], 
#                      saleyards[,c("category", "lat", "lon")]) |>
#   st_as_sf(coords = c("lon", "lat"), crs = st_crs(4326))

# Calculate point hazard (piggery, saleyard, feedlot) density per cell and
# count points per polygon 

# points_3308 <- st_transform(risk_points, 3308)
# LS_grid$risk_points <- lengths(st_intersects(LS_grid, points_3308))
# LS_grid$LS_points = range01(LS_grid$risk_points)


## Combine layers ==============================================================

# adjust by weight
LSw <- data.frame(LS_grid$LS_pig * LS_pig_w,
                  LS_grid$LS_cattle * LS_cattle_w,
                  LS_grid$LS_sheep * LS_sheep_w,
                  LS_grid$LS_goats * LS_goat_w,
                  LS_grid$LS_deer * LS_deer_w,
                  LS_grid$LS_points * LS_point_w)

# Sum total weighted risk for each cell and rescale from 0:1
LS_grid <- LS_grid |>
  mutate(LS_tot = range01(rowSums(LSw, na.rm=TRUE)))

# Categorise each cell by risk weight, such that 
# - high risk =       c. 9 times very low risk
# - moderate risk =   c. 7 times very low risk
# - low risk =        c. 3 times very low risk

LS_grid_cat <- LS_grid |>
  arrange(desc(LS_tot))

risk_weights <- c(high = 9, moderate = 7, low = 3, very_low = 1)

objective <- function(p, values) {
  # enforce ordering
  p <- sort(p, decreasing = TRUE)
  cuts <- quantile(values, probs = p)
  grp <- cut(values,
             breaks = c(-Inf, cuts[3], cuts[2], cuts[1], Inf),
             labels = c("very low","low","moderate","high"))
  means <- tapply(values, grp, mean)
  if (any(is.na(means))) return(1e6)
  vh <- means["very low"]
  err <- sum((
    c(means["high"]/vh,
      means["moderate"]/vh,
      means["low"]/vh)- risk_weights[c("high","moderate","low")])^2)
  return(err)
}


props <- risk_weights / sum(risk_weights)

# Convert to cumulative thresholds (from top)
p1 <- 1 - props["high"]
p2 <- 1 - (props["high"] + props["moderate"])
p3 <- 1 - (props["high"] + props["moderate"] + props["low"])
par_init <- c(p1, p2, p3)

res <- optim(par = par_init, 
             fn = objective,
             values = LS_grid$LS_tot,
             method = "L-BFGS-B",
             lower = c(0.3, 0.05, 0.001),
             upper = c(0.99, 0.8, 0.5))

p_opt <- sort(res$par, decreasing = TRUE)

cuts <- quantile(LS_grid_cat$LS_tot, probs = p_opt)

# Final categories
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

# Plot total expected risk (scaled 0:1)
ggplot() +
  geom_sf(aes(fill=LS_tot), col=NA, data=LS_grid, show.legend = TRUE) +
  scico::scale_fill_scico(palette = "lajolla", direction = -1, name = "Risk") +
  theme_void() +
  labs(title = "Total expected introduction risk, numeric")

ggplot() +
  geom_sf(aes(fill=risk_category), col=NA, data=LS_grid, show.legend = TRUE) +
  scale_fill_manual(values = c("red4", "orange3", "khaki","grey"), name="Risk") +
  theme_void() +
  labs(title = "Total introduction risk, categorical")

## Session info ================================================================
#sessioninfo::session_info() # When finalised
