# Wildlife-disease-surveillance-under-uncertainty

Reproducible workflow for: Bengsen et al., "Wildlife disease surveillance under uncertainty: an adaptive search-theoretic framework for early detection of transboundary animal diseases"

## PURPOSE
This script demonstrates the full analytical workflow described in the Methods: constructing a spatial FMDv risk prior, converting it to expected search value, selecting a spatially dispersed set of high-value cells, clustering realised sampling effort with DBSCAN, updating search values via Bayesian search theory, and bootstrapping surveillance system sensitivity (SSe).

## DATA
The real surveillance data cannot be shared because sampling was conducted on private properties. This script therefore runs on SYNTHETIC data that mimic the structure (but not the content) of the real dataset. Absolute results (risk values, cluster locations, SSe estimates) are illustrative only and should not be interpreted as the paper's actual findings. The code structure and parameter choices otherwise match the Methods.

## REQUIRED INPUTS 
   - FMDV_risk_3308_v3   : polygon layer (.shp) of the hexagonal
                           risk grid, with fields LS_pig, LS_tt_c, LS_tot,
                           cell (see Table 1 of the manuscript for the
                           component risk layers this is built from)
   - synthetic_pigs_3308.RDS : sf point object of synthetic pig sample
                           locations, EPSG:3308, with fields x, y, operation
   - error_df.csv        : lookup table of expected detection probability
                           (detect_prob) by sample size (n), from the
                           two-stage sampling model (see Supplementary
                           Figure S2)

## KEY PARAMETERS 
   pstar_among/pstar_within : design prevalence assumptions (0.05 / 0.1)
   se.elisa/sp.elisa        : NSP ELISA sensitivity/specificity (0.97/0.97)
   risk_pref                : relative risk weighting by category (9:7:3:1)
   eps/minpts               : DBSCAN neighbourhood radius (10 km) and
                              minimum cluster size (5 pigs)
   a, b                     : weights balancing search value vs. spatial
                              dispersion in cell selection (1, 5)

## SOFTWARE
Produced using R v4.3.1 (R Core Team 2024). See BST_DBSCA_synthetic_example_v1-1.R for full session info. 
