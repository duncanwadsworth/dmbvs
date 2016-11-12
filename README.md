dmbvs
===

---

### Description

Bayesian variable selection for Dirichlet-Multinomial regression. For $n$ 
samples, let $Y$ be an $n \times q$ matrix with $q$ taxa. Let $X$ be an 
$n \times p$ matrix with $p$ covariates. By applying spike-and-slab priors to 
the Dirichlet-Multinomial regression coefficients, we obtain concise 
taxa/covariate associations. The methodology is further described in the 
manuscript:

> Wadsworth, W. D., Argiento, R., Guindani, M., Galloway-Pena, J., Samuel, S. A., & 
> Vannucci, M. (2016). An Integrative Bayesian Dirichlet-Multinomial Regression 
> Model for the Analysis of Taxonomic Abundances in Microbiome data.

### Contents

The repository contains 

1. C code implementing an MCMC sampler for Dirichlet-Multinomial Bayesian 
Variable Selection (dmbvs) using spike-and-slab priors, 

2. R code to wrap and run the sampler from within R, 

3. a function for simulating data, 

4. and a "start-to-finish" script (`example_analysis_script.R`) demonstrating 
usage of the code. This script gives reasonable default settings for the 
hyperparameters and MCMC parameters for the example simulated data. Settings may 
change for other data.

### Usage

* The R code requires the `dirmult` and `MASS` packages. Please ensure those are 
installed first. 

* The C code relies on the [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/). 
GSL must be installed and modifications may need to be made to the 'library' 
and 'include' paths in the Makefile.

* The main C file (`dmbvs.c`), as well as the Makefile to be used for compilation 
of the code, can be found under the directory `code/`.

* On Linux and Mac the C code may be compiled with the Makefile from R using:

```{r}
# must be in package's root directory
setwd("dmbvs")
system("cd code; make")
```

* If compilation has been successful there will be an executable called `dmbvs.x` 
the `code/` directory. Data may be simulated and the MCMC code run using:

```{r}
source(file.path("code", "wrapper.R"))
source(file.path("code", "helper_functions.R"))
simdata = simulate_dirichlet_multinomial_regression(n_obs = 100, n_vars = 50,
                                                    n_taxa = 50, n_relevant_vars = 5,
                                                    n_relevant_taxa = 5)
results = dmbvs(XX = simdata$XX[,-1], YY = simdata$YY, 
                intercept_variance = 10, slab_variance = 10, 
                bb_alpha = 0.02, bb_beta = 1.98, GG = 1100L, thin = 10L, burn = 100L,
                exec = file.path(".", "code", "dmbvs.x"), output_location = ".")
```

* Tested on:
    * Mac with R version 3.2.2 and GSL version 2.1 with the Apple LLVM version 
    7.0.2 compiler
    * Red Hat Linux 6.6 with R version 3.1.2 and GSL version 2.1 and gcc 
    version 5.2.0

