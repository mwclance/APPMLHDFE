# The Tails of Gravity

[The Tails of Gravity: Using Expectiles to Quantify the Trade-Margins Effects of Economic Integration Agreements](https://www.surrey.ac.uk/sites/default/files/2025-02/DP01-25.pdf)

[Github repository](https://github.com/mwclance/APPMLHDFE)

## Authors
- Jeffrey H. Bergstrand
- Matthew W. Clance
- J.M.C. Santos Silva

## Overview
This repository introduces and implements the **Poisson-based expectile regressions** described in the article *"The Tails of Gravity: Using Expectiles to Quantify the Trade-Margins Effects of Economic Integration Agreements"*. The estimation method is called asymmetric Poisson pseudo maximum likelihood (APPML) estimation and is used **to analyze** heterogeneous effects in international trade, but it is also useful in many other applications.

The R script (`appml_r2`) has been adapted from the existing *Stata* command `appmlhdfe.ado`. 

**Stata:**  
Matthew Clance & J.M.C. Santos Silva, 2025. "APPMLHDFE: Stata module to estimate asymmetric Poisson regression with high-dimensional fixed effects," *Statistical Software Components* S459414, Boston College Department of Economics.

**Data:**  
Bergstrand, Jeffrey; Clance, Matthew; Santos Silva, Joao (2025), “The Tails of Gravity”, Mendeley Data, V1, doi: 10.17632/n67gft8fvm.1

### What Are Expectile Regressions?
Expectile regressions extend the flexibility of traditional regression models by estimating effects across the entire conditional distribution of a dependent variable. Unlike quantile regressions, expectiles are **global measures of location**, providing robust insights into how covariates influence not just the mean but also the tails of the distribution.

### Why Poisson-Based Expectiles?
The **Poisson-based expectile regression** approach combines the advantages of Poisson pseudo-maximum likelihood (PPML) with the ability to model heterogeneity across the distribution:
- Handles **non-negative dependent variables** (e.g., trade flows) with many zeros.
- Allows for high-dimensional panel data.
- Yields easily interpretable **elasticities or semi-elasticities** across different conditional expectiles.

### Key Features
- **Accommodates zeros**: Unlike quantile regression, which struggles with mass points at zero, Poisson-based expectiles maintain smooth functional relationships.
- **Heterogeneous effects**: Captures how relationships vary across different parts of the distribution, from lower to upper tails.
- **Easy implementation**: Simple extensions of PPML make expectile estimation computationally efficient.

---

## Repository Contents
- **R Scripts**: 
  - **appml_r2.R**: Function to estimate APPML is **essentially** a wrapper program using the `fixest` package's `feglm` function. Note that you can estimate either a single expectile or multiple expectiles. When estimating multiple expectiles, using starting values from sequential (previous) estimates will **speed up** the process.
  - **estimation_r2.R**: Replicates Figure 1 and Table 4 in *"The Tails of Gravity: Using Expectiles to Quantify the Trade-Margins Effects of Economic Integration Agreements"*.
  - **Sim_Chi2_updated.R**: Replicates the small simulation in *"The Tails of Gravity: Using Expectiles to Quantify the Trade-Margins Effects of Economic Integration Agreements"*.
  - **example/**: Example of APPML using a simulated dataset.


<head>
<meta name="google-site-verification" content="EzJYapKM-pCnhFPWOoB59ZSu2TeuSjfNA0rkmsIA88o" />
</head>



