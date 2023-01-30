# FMMInferences

This repository contains source codes to estimate 3DFMM model parameters and their confidence intervals. It also includes a code example in the `3DFMMUse.R` file.
The code is developed in the programming language R.

## `fitMultiFMM` function description.

`fitMultiFMM` performs the estimation of a 3DFMM model and it's confidence intervals based on the asymptotic parameter covariance matrix proposed in [1].

### Arguments

* `vDataMatrix`: double matrix. Each column contains the signal in a concrete direction.
* `nBack`: number of FMM components to be fitted.
* `maxIter`: maximum number of iterations for the backfitting algorithm.
* `lengthAlphaGrid`: number of values in the equally spaced grid for $\alpha$. By default it is set to 48.
* `lengthOmegaGrid`: number of values in the equally spaced grid for $\omega$. By default it is set to 24.
* `confidenceLevel`: confidence level for the parameter intervals. By default it is set to 0.95.
* `parallelize`: TRUE to use parallelized procedure to fit FMM model (recommended).

If `vDataMatrix` columns are named, the names will be used for label the directions in the resulting plot and the confidence intervals for the parameters. 

### Return values

An R list with two elements:
* `paramsPerWave`: an `R` list. Each element of the list contains the estimated values of the FMM parameters in a concrete direction. 
* `Confints`: confidence intervals for each parameter. The FMM parameters are named **parameter_component_direction**.
* It also returns a plot of the signal prediction and the fitted waves in the given directions.



## References

[1] Seber, G. A. and Wild C. J. (2003). Nonlinear regression. *New Jersey: Jhon Wiley & Sons, 62(63):1238.*
