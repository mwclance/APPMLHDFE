# The Tails of Gravity: Using Expectiles to Quantify the Trade-Margins Effects of Economic Integration Agreements

## Authors
- Jeffrey H. Bergstrand
- Matthew W. Clance
- J.M.C. Santos Silva

## Overview
This repository introduces and implements **Poisson-based expectile regressions** described in *The Tails of Gravity: Using Expectiles to Quantify the Trade-Margins Effects of Economic Integration Agreements* called asymmetric Poisson pseudo maximum likelihood (APPML). The estimation is used to analyze heterogeneous effects in international trade but it is also useful in many other applications. 

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
- **R Scripts**: Implementation of Poisson-based expectile regressions with sample datasets.


