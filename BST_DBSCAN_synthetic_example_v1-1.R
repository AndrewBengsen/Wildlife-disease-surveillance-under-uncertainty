## =============================================================================
## Reproducible workflow for: Bengsen et al., "Wildlife disease surveillance
## under uncertainty: an adaptive search-theoretic framework for early
## detection of transboundary animal diseases"

## PURPOSE
## This script demonstrates the full analytical workflow described in the
## Methods: constructing a spatial FMDv risk prior, converting it to expected
## search value, selecting a spatially dispersed set of high-value cells,
## clustering realised sampling effort with DBSCAN, updating search values
## via Bayesian search theory, and bootstrapping surveillance system
## sensitivity (SSe).

## DATA
## The real surveillance data cannot be shared because sampling was conducted
## on private properties. This script therefore runs on SYNTHETIC data that
## mimic the structure (but not the content) of the real dataset. Absolute
## results (risk values, cluster locations, SSe estimates) are illustrative
## only and should not be interpreted as the paper's actual findings. The
## code structure and parameter choices otherwise match the Methods.

## REQUIRED INPUTS 
##   - FMDV_risk_3308_v3   : polygon layer (.shp) of the hexagonal
##                           risk grid, with fields LS_pig, LS_tt_c, LS_tot,
##                           cell (see Table 1 of the manuscript for the
##                           component risk layers this is built from)
##   - synthetic_pigs_3308.RDS : sf point object of synthetic pig sample
##                           locations, EPSG:3308, with fields x, y, operation
##   - error_df.csv        : lookup table of expected detection probability
##                           (detect_prob) by sample size (n), from the
##                           two-stage sampling model (see Supplementary Fig S2)

## KEY PARAMETERS 
##   pstar_among/pstar_within : design prevalence assumptions (0.05 / 0.1)
##   se.elisa/sp.elisa        : NSP ELISA sensitivity/specificity (0.97/0.97)
##   risk_pref                : relative risk weighting by category (9:7:3:1)
##   eps/minpts               : DBSCAN neighbourhood radius (10 km) and
##                              minimum cluster size (5 pigs)
##   a, b                     : weights balancing search value vs spatial
##                              dispersion in cell selection (1, 5)

## SOFTWARE
## Produced using R v4.3.1 (R Core Team 2024). 
## Full session info at the tail of this workflow

library(tidyverse)          # v2.0.0
library(sf)                 # v1.0-15
library(dbscan)             # v1.2.2 
library(epiR)               # v2.0.91
library(scico)              # v1.5.0

## Dependencies ================================================================

# Scale 0:1
range01 <- function(x){(x-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))}

# Categorise risk or search value
classify_risk <- function(x, cuts){
  factor(case_when(x >= cuts[1] ~ "high",
                   x >= cuts[2] ~ "moderate",
                   x >= cuts[3] ~ "low",
                   TRUE         ~ "very low"),
         levels = c("high", "moderate", "low", "very low"))}

# update priors
update.priors <- function(prior, detect_prob){
  detect_prob[is.na(detect_prob)] <- 0
  posterior <- prior * (1 - detect_prob)
  posterior / sum(posterior)
}

# Calculate realised distribution of sampling effort among risk classes
effort_wt_fun <- function(cluster_current, risk_cat_prior){
  c(length(which(cluster_current  >=0 & risk_cat_prior == "high")),
    length(which(cluster_current  >=0 & risk_cat_prior == "moderate")),
    length(which(cluster_current  >=0 & risk_cat_prior == "low")),
    length(which(cluster_current  >=0 & risk_cat_prior == "very low"))) /
    length(which(cluster_current  >=0))
}

## Load the prior risk map =====================================================

# LS_tot provides a risk rating for each cell, in the range 0:1
dat <- st_read(dsn = getwd(), layer = "FMDV_risk_3308_v3") |>
  mutate(risk_cat = factor(LS_tt_c, 
                           levels = c("very low", "low", "moderate", "high"))) |>
  # mutate(risk_prior2 = recode(risk_cat, 
  #                             "very low" = 0.0001, 
  #                             "low" = 0.001, 
  #                             "moderate" = 0.005, 
  #                             "high" = 0.01)) |>
  dplyr::select(cell, pig_dens = LS_pig, LS_tot, risk_cat) 

# Create a risk prior that sums to 1 across the state
dat$risk_prior <- dat$LS_tot / sum(dat$LS_tot)

# What proportion of the state falls into each risk category?
# (for surveillance system sensitivity (SSe) estimation)
proportion_by_risk <- dat |>
  st_drop_geometry() |>
  group_by(risk_cat) |>
  summarise(n = n()) |>
  mutate(p=n/nrow(dat)) 
proportion_by_risk

## Assign likelihood ===========================================================

# All cells start with expected SSe = 0.95, based on:
 #  sample size = 30
 #  minimum expected prevalence / design prevalence
    pstar_among <- 0.05
    pstar_within <- 0.1
 # NSP ELISA sensitivity and specificity
    se.elisa <- 0.97   
    sp.elisa <- 0.97   
 # risk weight by site risk category
    risk_pref <- c(high = 9, moderate = 7, low = 3, very_low = 1)
 # constant risk weight among individuals within sites
    ppr.u <- 1

# prior proportional representation of each risk category across NSW 
proportional_rep <- proportion_by_risk$p

SSe <- 0.95 
dat$likelihood <- SSe 

## Assign belief ===============================================================

# Categorise expected search value based on risk weightings

df <- arrange(dat, desc(risk_prior))

objective <- function(p, values) {
  p <- sort(p, decreasing = TRUE)
  cuts <- quantile(values, probs = p)
  grp <- cut(values,
             breaks = c(-Inf, cuts[3], cuts[2], cuts[1], Inf),
             labels = c("very low","low","moderate","high"))
  means <- tapply(values, grp, mean)
  if (any(is.na(means))) return(1e6)
  vh <- means["very low"]
  err <- sum((c(means["high"]/vh,
                means["moderate"]/vh,
                means["low"]/vh) - c(9, 7, 3))^2)
  return(err)
}

props <- risk_pref / sum(risk_pref)

# Convert to cumulative thresholds (from top)
p1 <- 1 - props["high"]
p2 <- 1 - (props["high"] + props["moderate"])
p3 <- 1 - (props["high"] + props["moderate"] + props["low"])

par_init <- c(p1, p2, p3)

res <- optim(par = par_init, 
             fn = objective,
             values = df$risk_prior,
             method = "L-BFGS-B",
             lower = c(0.3, 0.05, 0.001),
             upper = c(0.99, 0.8, 0.5))

# NB: the optimiser fails to converge here, because it's a difficult problem,
# but it returns the best values that it had found up to that point. 
# However, vv << oo below, suggesting that a sensible solution had been found. 
oo <- objective(par_init, df$risk_prior)
vv <- res$value

p_opt <- sort(res$par, decreasing = TRUE)

cuts <- quantile(df$risk_prior, probs = p_opt)

dat <- dat |>
  mutate(search_value = round(risk_prior * likelihood, 5)) |>
  dplyr::select(cell, risk_prior, risk_cat, likelihood, search_value)

dat <- dat |>
  mutate(search_cat = case_when(
    risk_prior >= cuts[1] * likelihood ~ "high",
    risk_prior >= cuts[2] * likelihood ~ "moderate",
    risk_prior >= cuts[3] * likelihood ~ "low",
    TRUE ~ "very low")) |>
  mutate(search_cat = factor(search_cat, 
                             levels = c("high", "moderate", "low", "very low")))

## Select cells to balance risk rating and spatial dispersion ==================

# Centroids (x,y) for distance calculations
LS_grid <- dat |>
  mutate(centroid = st_centroid(geometry)) |>
  mutate(x = st_coordinates(centroid)[,1],
         y = st_coordinates(centroid)[,2])

# Set a risk threshold for selection 
prior_search_quantiles <- quantile(LS_grid$search_value, 
                                   probs = c(0.25, 0.5, 0.6, 0.75, 0.8, 1))
hist(LS_grid$search_value)

# Top quartile of cells, by search value, are available for selection
high_risk_grid_prior <- LS_grid |>
  dplyr::filter(search_value > prior_search_quantiles[4])

# How many cells in the top quartile?
nrow(high_risk_grid_prior) 

# Run the loop to select cells, weighted by (a) risk and (b) spatial dispersion 
a <- 1 
b <- 5 

  # Rescale search values across all remaining cells
  remaining_prior <- high_risk_grid_prior |>
    mutate(search_value_scaled = range01(search_value))
  
  # Calculate all pairwise distances between centroids
  all_dists <- st_distance(remaining_prior$centroid, remaining_prior$centroid) |> 
    as.numeric()
  
  # Remove infinite or NA values
  all_dists <- all_dists[is.finite(all_dists)]
  
  # Get min and max for normalisation
  min_dist_val <- min(all_dists)
  max_dist_val <- max(all_dists)
  
  # Selection loop
  selected_prior<- remaining_prior[0, ]  # empty sf object with same structure
  N <- 10  
  
  for (i in 1:N) {
    scores <- sapply(1:nrow(remaining_prior), function(j) {
      cell <- remaining_prior[j, ]
      
      if (nrow(selected_prior) == 0) {
        return(a * cell$search_value_scaled)
      }
      
      # Calculate distance to nearest selected cell
      dists <- st_distance(cell$centroid, selected_prior$centroid) |> 
        as.numeric()
      min_dist <- min(dists[is.finite(dists)])
      
      # Normalise the distance
      min_dist_norm <- (min_dist - min_dist_val) / (max_dist_val - min_dist_val)
      
      # Calculate the score
      score <- a * cell$search_value_scaled + b * min_dist_norm
      return(score)
    })
    cat("Selected cell", i, "of", N, "\n")
    # Select the best cell
    best_index <- which.max(scores)
    selected_prior <- rbind(selected_prior, remaining_prior[best_index, ])
    remaining_prior <- remaining_prior[-best_index, ]
  }
  
#Assign a rank. Highest number = greatest search value
selected_prior$rank <- c(nrow(selected_prior):1) 
selected_prior <- selected_prior |>
  dplyr::select(-centroid)


# Plot selected cells
max_search_val <- max(LS_grid$search_value, na.rm=T)

ggplot() +
  geom_sf(data = LS_grid, aes(fill = search_value), colour=NA, show.legend=T) +
  geom_sf(data = selected_prior, color = "white", size = 1, fill = NA) +
  scale_fill_scico(palette = "lajolla", direction = -1, 
                          name="Expected search value", 
                          guide = guide_colourbar(title.position = "top"),
                          breaks = c(0, max_search_val),
                          labels = c("Low", "High")) +
  labs(title = "Cells selected for surveuillance, prior to initial activities") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10))

## Cluster cells by realised search effort, using DBSCAN =======================

# Load and examine synthetic pig samples created based on spatial distributions
# of samples within operations in the real data between July 2024 and July 2025.
# The real data can not be made public for privacy reasons.

# sf object, projected crs = EPSG:3308
pig_sp <- readRDS("synthetic_pigs_3308.RDS")

ggplot() +
  geom_sf(data = LS_grid, aes(), fill="grey", colour=NA, show.legend=F) +
  geom_sf(data = pig_sp, aes(), size = 1, show.legend=F) +
  scale_fill_scico(palette = "lajolla", direction = -1, 
                          name="Expected search value", 
                          guide = guide_colourbar(title.position = "top"),
                          breaks = c(0, max_search_val),
                          labels = c("Low", "High")) +
  labs(title = "Synthetic pig locations") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10))

# Set clustering parameters
eps <- 10     # cluster neighbourhood radius (kms)
minpts <- 5   # minimum number of points to qualify as a cluster

nrow(pig_sp)                      # N animals sampled   
length(unique(pig_sp$operation))  # N sampling operations
pxy <- data.frame(x = pig_sp$x, y = pig_sp$y)

# DBSCAN. eps km scaled converted to m to match the projected coordinates
Dbscan_p <- dbscan(pxy, eps = 1000*eps, minPts = minpts)

# Assign cluster id to each sample 
pig_sp$cluster <- Dbscan_p$cluster

length(unique(pig_sp$cluster))-1    # N clusters (-1 to exclude non-cluster 0's)
length(which(pig_sp$cluster == 0))  # N pigs outside of a cluster

# Create a minuimum convex polygon hull around each sample 
# with an outer buffer of 0.5\*eps 
hulls <- pig_sp |>
  filter(!cluster %in% c(0)) |> 
  group_by(cluster) |>
  summarise(geometry = st_combine(geometry)) |>
  st_convex_hull() |>
  st_buffer(dist=1000*eps/2) # sets a buffer around the MCP of length half eps

# Add sample size to each hull
hulls$n <- lengths(st_intersects(hulls, pig_sp))

# Assign cluster values and cluster sample sizes to the risk map cells
cl_grid <- st_join(dat, hulls, left = T) 

# Note that this can induce duplication of some cells that intersect > 1 cluster
# Although it doesn't with these synthetic data
# Retain the duplicate row with the highest samples (lowest new_search_value)
nrow(cl_grid)
cl_grid <- cl_grid |>
  group_by(cell) |>
  slice_max(coalesce(n, -Inf), n = 1, with_ties = FALSE) |>
  ungroup()
nrow(cl_grid)

# Plot synthetic sample clusters
ggplot() +
  geom_sf(aes(fill=search_value), col=NA, show.legend = F, data=dat) +
  scale_fill_scico(palette = "lajolla", alpha=0.9, direction = -1) +
  geom_sf(aes(), fill=NA, colour="white", data=hulls, show.legend = F) +
  labs(title = "Synthetic sample clusters") +
  theme_void()

# Realised distribution of sampling effort among search value classes
# high value : very low value
effort_weights_1  <- effort_wt_fun(cl_grid$cluster, cl_grid$search_cat)

## Update search values (beliefs) based on realised search effort ==============

# Load a table showing how detection probability is expected to decline with
# with diminishing sample sizes
# (based on previously stated sample size assumptions)
error_df <- read.csv("error_df.csv") |>
  arrange(desc(n)) 

# New fields to hold updated search value and identify whether 
# a cell has been searched
cl_grid$searched <- new_search_value <- 0
cl_grid$searched[which(is.na(cl_grid$cluster)==F)] <- 1

# Match error downgrading of likelihood (pD) to sample size (n)
cl_grid <- left_join(cl_grid, error_df, by="n")

# Add detection probability = 0.97 for samples > 35
cl_grid$detect_prob[which(cl_grid$n > 35)] <- 0.97

# Calculate detection probability for each cell
cl_grid$detect_prob[is.na(cl_grid$detect_prob)] <- 0

cl_grid$new_search_value <- update.priors(prior = cl_grid$risk_prior,
                                          detect_prob = cl_grid$detect_prob)

# Check that the new search value sums to 1 across the state 
sum(cl_grid$new_search_value)

# Categorise new search value categories for cells
risk_cuts <- quantile(dat$risk_prior, 
                      probs = sort(p_opt, decreasing = TRUE),
                      na.rm = TRUE)

cl_grid$new_search_cat <- classify_risk(cl_grid$new_search_value, risk_cuts)
table(cl_grid$new_search_cat)

# Plot updated search values
ggplot(cl_grid) + 
  geom_sf(aes(fill = new_search_value), col=NA, show.legend = F) +
  scale_fill_scico(palette = "lajolla", direction = -1) +
  theme_void(base_size = 12) +
  labs(title = "Updated search values following 12 months of sampling",
       caption = "Based on synthetic data")

## Run the second cell selection loop, based on updated search values ==========
# Following the previous selection procedure

# Centroids for distance calculations
LS_grid_1 <- cl_grid |>
  mutate(centroid = st_centroid(geometry)) |>
  mutate(x = st_coordinates(centroid)[,1],
         y = st_coordinates(centroid)[,2])

# Set a risk threshold for selection 
quantile(LS_grid_1$new_search_value, probs = c(0.25, 0.5, 0.6, 0.75, 0.8, 1))
hist(LS_grid_1$new_search_value)

threshold <- quantile(LS_grid_1$new_search_value, probs = 0.75, na.rm = TRUE)

high_risk_grid_1 <- LS_grid_1 |>
  filter(new_search_value > threshold) |> 
  mutate(rank = NA) 

# Identify cells that have already been effectively searched
searched_grid <- LS_grid_1 |>
   dplyr::filter(detect_prob > 0.94) 
 searched_buffer <- st_buffer(searched_grid, 12500) |>
   st_union()

# Run the loop
remaining_1 <- high_risk_grid_1 |>
  mutate(search_value_scaled = range01(new_search_value))
  
  # Remove cells that intersect the buffer around already searched cells
  remaining_intersects <- st_intersects(remaining_1, searched_buffer, sparse=F)
  remaining_1 <- remaining_1[!remaining_intersects,] 
  
  # Calculate all pairwise distances between centroids
  all_dists <- st_distance(remaining_1$centroid, remaining_1$centroid) |> 
    as.numeric()
  
  # Remove infinite or NA values
  all_dists <- all_dists[is.finite(all_dists)]
  
  # Get min and max for normalisation
  min_dist_val <- min(all_dists)
  max_dist_val <- max(all_dists)
  
  # Selection loop
  selected_1 <- remaining_1[0, ]  # empty sf object with same structure
  N <- 10  
  
  for (i in 1:N) {
    scores <- sapply(1:nrow(remaining_1), function(j) {
      cell <- remaining_1[j, ]
      
      if (nrow(selected_1) == 0) {
        return(a * cell$search_value_scaled)
      }
      
      # Compute distance to nearest selected cell
      dists <- st_distance(cell$centroid, selected_1$centroid) |> 
        as.numeric()
      min_dist <- min(dists[is.finite(dists)])
      
      # Normalize the distance
      min_dist_norm <- (min_dist - min_dist_val) / (max_dist_val - min_dist_val)
      
      # Compute score
      score <- a * cell$search_value_scaled + b * min_dist_norm
      return(score)
    })
    cat("Selected cell", i, "of", N, "\n")
    # Select the best cell
    best_index <- which.max(scores)
    selected_1 <- rbind(selected_1, remaining_1[best_index, ])
    remaining_1 <- remaining_1[-best_index, ]
  }
  
# Assign a rank. Highest number = greatest search value
selected_1$rank <- c(nrow(selected_1):1) 
selected_1 <- selected_1 |>
  dplyr::select(-centroid)

# Plot selected cells
ggplot() +
   geom_sf(data = LS_grid_1, aes(fill = new_search_value), colour=NA, 
           show.legend=T) +
   geom_sf(data = selected_1, colour = "white", size = 1, fill = NA) +
   scale_fill_scico(palette = "lajolla", direction = -1, 
                           name="Expected search value") +
   labs(title = "Selected cells after intial update") +
   theme_void() +
   theme(legend.position = "bottom",
         legend.text = element_blank())

## Given the sampling effort, what was the realised SSe? =======================

# Assign a dominant risk_cat (from the fixed original map) to each cluster
assign_cluster_risk <- function(hulls, dat) {
  ints <- st_intersects(hulls, dat)
  rg_lab <- sapply(ints, function(idx) {
    if (length(idx) == 0) return(NA_character_)
    tab <- table(dat$risk_cat[idx])
    names(tab)[which.max(tab)]          # modal risk_cat among intersecting cells
  })
  factor(rg_lab, levels = c("high","moderate","low","very low"))
}

# Bootstrap confidence intervals for SSe using the deterministic rsu.sep.rb2st()  
# Resample the clusters (rows of hulls) with replacement, 
# recalculate se.p each time, take percentiles

boot_realised_se_sys <- function(hulls, dat, B = 2000, ...) {
  rg_fac <- assign_cluster_risk(hulls, dat)
  keep   <- !is.na(rg_fac) & !is.na(hulls$n)
  hulls_k <- hulls[keep, ]
  rg_k    <- rg_fac[keep]
  
  se_boot <- numeric(B)
  n_clusters <- nrow(hulls_k)
  
  for (b in seq_len(B)) {
    idx <- sample.int(n_clusters, replace = TRUE)
    n_mat <- matrix(hulls_k$n[idx], ncol = 1)
    ppr.u.mat <- matrix(1, nrow = length(idx), ncol = 1)
    
    out <- rsu.sep.rb2st(
      H = NA, N = NA, n = n_mat,
      pstar.c = pstar_among, pstar.u = pstar_within,
      rg = as.integer(rg_k[idx]),
      rr.c = risk_pref, rr.u = ppr.u,
      ppr.c = proportional_rep, ppr.u = ppr.u.mat,
      se.u = se.elisa
    )
    se_boot[b] <- out$se.p
  }
  quantile(se_boot, c(0.025, 0.5, 0.975), na.rm = TRUE)
}

set.seed(666)
realised_SSe_CI <- round(boot_realised_se_sys(hulls,  dat), 3)

paste0("Realised SSe = ", realised_SSe_CI[2], " (95% CI = ", realised_SSe_CI[1],
       ", ", realised_SSe_CI[3], ")")

## Session info ================================================================
# setting  value
# version  R version 4.3.1 (2023-06-16 ucrt)
# os       Windows 11 x64 (build 26200)
# system   x86_64, mingw32
# ui       RStudio
# language (EN)
# collate  English_Australia.utf8
# ctype    English_Australia.utf8
# tz       Australia/Sydney
# date     2026-09-25
# rstudio  2025.05.0+496 Mariposa Orchid (desktop)
# pandoc   3.4 @ C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)

# Packages ─────────────────────────────────────────────────────────────────────
# package           * version date (UTC) lib source
# askpass             1.2.0   2023-09-03 [1] CRAN (R 4.3.2)
# BiasedUrn           2.0.12  2024-06-16 [1] CRAN (R 4.3.3)
# class               7.3-22  2023-05-03 [1] CRAN (R 4.3.1)
# classInt            0.4-10  2023-09-05 [1] CRAN (R 4.3.2)
# cli                 3.6.2   2023-12-11 [1] CRAN (R 4.3.2)
# colorspace          2.1-0   2023-01-23 [1] CRAN (R 4.3.2)
# crayon              1.5.2   2022-09-29 [1] CRAN (R 4.3.2)
# crul                1.4.0   2023-05-17 [1] CRAN (R 4.3.3)
# curl                5.2.0   2023-12-08 [1] CRAN (R 4.3.2)
# data.table          1.15.0  2024-01-30 [1] CRAN (R 4.3.2)
# DBI                 1.2.2   2024-02-16 [1] CRAN (R 4.3.2)
# dbscan            * 1.2.2   2025-01-26 [1] CRAN (R 4.3.3)
# digest              0.6.34  2024-01-11 [1] CRAN (R 4.3.2)
# dplyr             * 1.1.4   2023-11-17 [1] CRAN (R 4.3.2)
# e1071               1.7-14  2023-12-06 [1] CRAN (R 4.3.2)
# ellipsis            0.3.2   2021-04-29 [1] CRAN (R 4.3.2)
# epiR              * 2.0.91  2026-02-23 [1] CRAN (R 4.3.1)
# evaluate            0.23    2023-11-01 [1] CRAN (R 4.3.2)
# farver              2.1.1   2022-07-06 [1] CRAN (R 4.3.2)
# fastmap             1.1.1   2023-02-24 [1] CRAN (R 4.3.2)
# flextable           0.9.4   2023-10-22 [1] CRAN (R 4.3.3)
# fontBitstreamVera   0.1.1   2017-02-01 [1] CRAN (R 4.3.1)
# fontLiberation      0.1.0   2016-10-15 [1] CRAN (R 4.3.1)
# fontquiver          0.2.1   2017-02-01 [1] CRAN (R 4.3.3)
# forcats           * 1.0.0   2023-01-29 [1] CRAN (R 4.3.2)
# gdtools             0.3.6   2024-02-22 [1] CRAN (R 4.3.3)
# generics            0.1.4   2025-05-09 [1] CRAN (R 4.3.1)
# gfonts              0.2.0   2023-01-08 [1] CRAN (R 4.3.3)
# ggplot2           * 3.5.2   2025-04-09 [1] CRAN (R 4.3.1)
# glue                1.7.0   2024-01-09 [1] CRAN (R 4.3.2)
# gtable              0.3.4   2023-08-21 [1] CRAN (R 4.3.2)
# hms                 1.1.3   2023-03-21 [1] CRAN (R 4.3.2)
# htmltools           0.5.7   2023-11-03 [1] CRAN (R 4.3.2)
# httpcode            0.3.0   2020-04-10 [1] CRAN (R 4.3.3)
# httpuv              1.6.14  2024-01-26 [1] CRAN (R 4.3.2)
# jsonlite            1.8.8   2023-12-04 [1] CRAN (R 4.3.2)
# KernSmooth          2.23-21 2023-05-03 [1] CRAN (R 4.3.1)
# knitr               1.45    2023-10-30 [1] CRAN (R 4.3.2)
# labeling            0.4.3   2023-08-29 [1] CRAN (R 4.3.1)
# later               1.3.2   2023-12-06 [1] CRAN (R 4.3.2)
# lattice             0.21-8  2023-04-05 [1] CRAN (R 4.3.1)
# lifecycle           1.0.5   2026-01-08 [1] CRAN (R 4.3.1)
# lubridate         * 1.9.3   2023-09-27 [1] CRAN (R 4.3.2)
# magrittr            2.0.3   2022-03-30 [1] CRAN (R 4.3.2)
# Matrix              1.6-5   2024-01-11 [1] CRAN (R 4.3.2)
# mime                0.12    2021-09-28 [1] CRAN (R 4.3.1)
# munsell             0.5.0   2018-06-12 [1] CRAN (R 4.3.2)
# officer             0.6.5   2024-02-24 [1] CRAN (R 4.3.3)
# openssl             2.1.1   2023-09-25 [1] CRAN (R 4.3.2)
# pander              0.6.5   2022-03-18 [1] CRAN (R 4.3.3)
# pillar              1.11.1  2025-09-17 [1] CRAN (R 4.3.1)
# pkgconfig           2.0.3   2019-09-22 [1] CRAN (R 4.3.2)
# promises            1.2.1   2023-08-10 [1] CRAN (R 4.3.2)
# proxy               0.4-27  2022-06-09 [1] CRAN (R 4.3.2)
# purrr             * 1.0.2   2023-08-10 [1] CRAN (R 4.3.2)
# R6                  2.6.1   2025-02-15 [1] CRAN (R 4.3.3)
# ragg                1.2.7   2023-12-11 [1] CRAN (R 4.3.2)
# Rcpp                1.0.12  2024-01-09 [1] CRAN (R 4.3.2)
# readr             * 2.1.5   2024-01-10 [1] CRAN (R 4.3.2)
# rlang               1.1.3   2024-01-10 [1] CRAN (R 4.3.2)
# rmarkdown           2.25    2023-09-18 [1] CRAN (R 4.3.2)
# rsconnect           1.3.1   2024-06-04 [1] CRAN (R 4.3.1)
# rstudioapi          0.15.0  2023-07-07 [1] CRAN (R 4.3.2)
# scales              1.3.0   2023-11-28 [1] CRAN (R 4.3.2)
# scico             * 1.5.0   2023-08-14 [1] CRAN (R 4.3.3)
# sessioninfo         1.2.2   2021-12-06 [1] CRAN (R 4.3.2)
# sf                * 1.0-15  2023-12-18 [1] CRAN (R 4.3.2)
# shiny               1.8.0   2023-11-17 [1] CRAN (R 4.3.2)
# stringi             1.8.3   2023-12-11 [1] CRAN (R 4.3.2)
# stringr           * 1.5.1   2023-11-14 [1] CRAN (R 4.3.2)
# survival          * 3.7-0   2024-06-05 [1] CRAN (R 4.3.3)
# systemfonts         1.0.5   2023-10-09 [1] CRAN (R 4.3.2)
# textshaping         0.3.7   2023-10-09 [1] CRAN (R 4.3.2)
# tibble            * 3.2.1   2023-03-20 [1] CRAN (R 4.3.2)
# tidyr             * 1.3.1   2024-01-24 [1] CRAN (R 4.3.2)
# tidyselect          1.2.1   2024-03-11 [1] CRAN (R 4.3.3)
# tidyverse         * 2.0.0   2023-02-22 [1] CRAN (R 4.3.2)
# timechange          0.3.0   2024-01-18 [1] CRAN (R 4.3.2)
# tzdb                0.4.0   2023-05-12 [1] CRAN (R 4.3.2)
# units               0.8-5   2023-11-28 [1] CRAN (R 4.3.2)
# utf8                1.2.4   2023-10-22 [1] CRAN (R 4.3.2)
# uuid                1.2-0   2024-01-14 [1] CRAN (R 4.3.2)
# vctrs               0.6.5   2023-12-01 [1] CRAN (R 4.3.2)
# withr               3.0.2   2024-10-28 [1] CRAN (R 4.3.3)
# xfun                0.42    2024-02-08 [1] CRAN (R 4.3.2)
# xml2                1.3.6   2023-12-04 [1] CRAN (R 4.3.2)
# xtable              1.8-4   2019-04-21 [1] CRAN (R 4.3.2)
# zip                 2.3.1   2024-01-27 [1] CRAN (R 4.3.2)
# zoo                 1.8-12  2023-04-13 [1] CRAN (R 4.3.2)

# [1] C:/Program Files/R/R-4.3.1/library