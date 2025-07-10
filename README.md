# Estimación del parámetro de dependencia FGM usando métodos clásicos y métodos bayesianos informativos.

Este repositorio contiene los resultados de un estudio de simulación sobre la estimación del parámetro de dependencia FGM. La función de densidad de la copula FGM es:

$$c_{\varphi}(u,v)= (1+\varphi(1-2u)(1-2v)),$$

donde $\varphi$ es el parámetro de dependencia, cuyo espacio de valores es el intervalo $(0,1)$.

## Table of Contents

-   [Overview](#overview)
-   [Installation](#installation)
-   [Usage](#usage)
-   [Included Functions](#included-functions)
-   [Contributions](#contributions)
## Overview

Los métodos de estimación clásicos utilizados fue el método de máxima verosimilitud y métodos de momentos generados por el coeficiente de concordancia Tau de Kendall, Rho de Spearman y Blomqvist’s Beta,
mientras que desde el enfoque Bayesiano fueron utilizadas las distribuciones a priori Triangular, Beta de cuatro parámetros y Uniforme. Debido a que las integrales obtenidas para la distribución posterior no tienen forma cerrada,
como por ejemplo la producida por la prior beta de cuatro parámetros:

$$\dfrac{ (1+\varphi)^{\tilde\alpha-1} ( 1-\varphi)^{\tilde\beta-1}\prod\limits_{i = 1}^n c_{\varphi}(\hat u_i,\hat v_i)}{\int_{-1}^{1} (1+\varphi)^{\tilde\alpha-1} ( 1-\varphi)^{\tilde\beta-1}\cdot\prod\limits_{i = 1}^n c_{\varphi}(\hat u_i,\hat v_i)  d\varphi}$$

fue necesario utilizar el método de Metropolis-Hasting en su forma estándar (usando la distribución uniforme como instrumental) y en sú forma de camianta aleatoria (usando la distribución beta de cuatro parámetros como instrumental).

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

    ``` bash
    git clone [https://github.com/LLerzy/Estimation-Parameter-Beta/tree/main.git](https://github.com/LLerzy/FGM_Estimates.git)
    ```

2. Open the R script in any compatible Integrated Development Environment (IDE):

    -   For the script that contains the main functions in the shape parameter estimation process, refer to [`FunctionsRequired.R`](https://github.com/LLerzy/FGM_Estimates/blob/main/FunctionsRequired.R).
    -   El script que contiene el esquema de simulación utilizado para obtener los resultados de simulación puede ser consultado en: [`SimulationStudy.R`](https://github.com/LLerzy/FGM_Estimates/blob/main/SimulationStudy.R).
    -   EL script que contiene el código para obtener la estimación del parámetro de dependencia asociada al conjunto de dengue puede ser consultado en: [`DengueDataSet.R`](https://github.com/LLerzy/FGM_Estimates/blob/main/DengueData/DengueDataSet.R).

3. You can run the script and load all the functions with the following command:

    ``` bash
    Rscript FunctionsRequired.R
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

### Notes:

-   It is recommended to review and adapt each function according to the specific needs of each analysis.
-   Make sure you understand each function before using it to ensure accurate results and avoid potential errors.

## Contributions

Contributions are welcome! Please submit a pull request or open an issue if you have any suggestions or improvements.
