# Estimation of the FGM Copula Dependence Parameter Using Classical and Informative Bayesian Methods

This repository contains the results of a simulation study focused on estimating the dependence parameter of the FGM copula. The corresponding density function is:

$$c_{\varphi}(u,v)= (1+\varphi(1-2u)(1-2v)),$$

where $\varphi$ is the dependence parameter, defined in the interval $(0,1)$.

## Table of Contents

* [Overview](#overview)
* [Installation](#installation)
* [Usage](#usage)
* [Included Functions](#included-functions)
* [Results](#results)
* [Contributions](#contributions)

## Overview
This repository presents a simulation study on the estimation of the dependence parameter of the FGM copula. Classical estimation methods, such as maximum likelihood and moment-based estimators using Kendall's Tau, Spearman's Rho, and Blomqvist’s Beta, are compared with informative Bayesian approaches employing Triangular, four-parameter Beta, and Uniform prior distributions.

Since the resulting posterior distributions—such as the one induced by the four-parameter Beta prior—do not have a closed form, the Metropolis-Hastings algorithm was applied in two variants: the standard form (using a uniform instrumental distribution) and the random walk version (using the four-parameter Beta distribution as instrumental). For example, the posterior is given by:

$$\dfrac{ (1+\varphi)^{\tilde\alpha-1} ( 1-\varphi)^{\tilde\beta-1}\prod\limits_{i = 1}^n c_{\varphi}(\hat u_i,\hat v_i)}{\int_{-1}^{1} (1+\varphi)^{\tilde\alpha-1} ( 1-\varphi)^{\tilde\beta-1}\cdot\prod\limits_{i = 1}^n c_{\varphi}(\hat u_i,\hat v_i)  d\varphi}$$

Convergence monitoring for this algorithm was conducted using graphical tools such as the histogram overlaid with the fitted density, trace plots, cumulative averages, and autocorrelation functions. An example of these results, based on a Beta prior and a chain of 5,000 iterations, is available here: [Beta Prior Results](https://github.com/LLerzy/FGM_Estimates/blob/main/Graphics/ResultsBetaBeta.png).

## Installation

To run the code, you need to have `R` installed with the following packages:

- `copula`
- `xlsx`
- `readxl`
- `ggplot2` 
- `gridExtra` 
- `correlation` 
- `betafunctions` 
- `parallel` 

You can install the required packages using the following command:

``` r
install.packages(c("copula","xlsx","readxl","ggplot2", "gridExtra", "correlation","betafunctions","parallel"))
```

## Usage

To run the designed algorithms, follow these steps:

1. Clone this repository:

```bash
git clone https://github.com/LLerzy/FGM_Estimates.git
```

2. Open the scripts in an R-compatible development environment. The main files are:

* [`FunctionsRequired.R`](https://github.com/LLerzy/FGM_Estimates/blob/main/FunctionsRequired.R): main functions for parameter estimation.
* [`SimulationStudy.R`](https://github.com/LLerzy/FGM_Estimates/blob/main/SimulationStudy.R): simulation study scheme.
* [`DengueDataSet.R`](https://github.com/LLerzy/FGM_Estimates/blob/main/DengueData/DengueDataSet.R): application using dengue data.

3. To load the functions and run the scripts:

```bash
Rscript FunctionsRequired.R
Rscript SimulationStudy.R
Rscript DengueData/DengueDataSet.R
```

## Included Functions

The functions defined in the `FunctionsRequired.R` script include:

-   `Density_FGM`: FGM copula density function.
-   `Der_log_lik_FGM`: Derivatives of the log-likelihood function.
-   `rfgm / rnfgm`: Generation of FGM random samples (1 or n observations).
-   `me`: Moment-based estimators (Pearson, Kendall, Spearman).
-   `metropolis_hastings`: MCMC sampling with flexible priors/instrumentals.
-   `Graphs`: MCMC convergence diagnostics.
-   `Est_A_Samp`: Estimation from a single sample (point, interval, credibility).
-   `Sim_Est_A_Samp`: Parameter estimation with filtering of invalid samples.
-   `Sim_Est_N_Samp`: Repeated simulations for performance evaluation.
-   `Mtovar`: Prior elicitation for the beta distribution (Tovar’s method).
-   `GeneralGraph`: Aggregated plot generation for simulation results.
-   `Marginal`: Empirical marginal distribution estimation.

**Note:** It is recommended to review and adapt each function to the specific needs of your analysis. Understanding each routine is essential to avoid misinterpretation of the results.

## Results

The graphical results from the simulation study can be found in the [`Graphics`](https://github.com/LLerzy/FGM_Estimates/tree/main/Graphics) folder.

Numerical results from the simulations are available in the [`Results`](https://github.com/LLerzy/FGM_Estimates/tree/main/Results) folder.

The results of the method applied to the dengue dataset can be found in the [`DengueData`](https://github.com/LLerzy/FGM_Estimates/tree/main/DengueData) folder, particularly in the `DengueData6.RData` file, which stores all the results obtained.

## Contributions

Contributions are welcome! Please submit a pull request or open an issue if you have any suggestions or improvements.
