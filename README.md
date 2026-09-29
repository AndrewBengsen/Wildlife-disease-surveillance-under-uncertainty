# Wildlife-disease-surveillance-under-uncertainty

Reproducible workflow for: Bengsen et al., [Wildlife disease surveillance under uncertainty: an adaptive search-theoretic framework for early detection of transboundary animal diseases](https://doi.org/10.64898/2026.08.31.748179)

## PURPOSE
This script demonstrates the full workflow described in the Methods:  

`Create_prior_risk_surface_v1-1.R` combines spatial layers representing different aspects of the risk of FMDv incursion in feral pig populations.  

`BST_DBSCA_synthetic_example_v1-1.R` converts the risk map to a map expected search values, before selecting a spatially dispersed set of high-value cells, clustering realised sampling effort with DBSCAN, and updating search values via Bayesian search theory.

## DATA
The real surveillance data cannot be shared because sampling was conducted on private properties. This script therefore runs on SYNTHETIC data that mimic the structure (but not the content) of the real dataset. Absolute results (risk values, cluster locations, SSe estimates) are illustrative only and should not be interpreted as the paper's actual findings. The code structure and parameter choices otherwise match the Methods.

## SOFTWARE
Produced using R v4.3.1 (R Core Team 2023).  
See `BST_DBSCA_synthetic_example_v1-1.R` and `Create_prior_risk_surface_v1-1.R` for full session info. 

## LICENCE
Source code and datasets are licensed under the MIT License (see `LICENCE.md`).  
