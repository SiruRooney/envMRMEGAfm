# The weighted environment-adjusted linear regression.

Fit the linear regression with weights using GWAS data across cohorts.
The contribution of each GWAS is weighted by the estimated inverse
variance of the reference allele effect at the corresponding variant.

## Usage

``` r
lr_env_w(Y, X, W, n_pc)
```

## Arguments

- Y:

  BETA values of the genetic variants across all cohorts.

- X:

  PCs plus environment covariates.

- W:

  The vector of weights, each component referring to the estimated
  inverse variance of the allele effect at the corresponding variants.
  The length of the vector is the number of cohorts.

- n_pc:

  The number of PCs.

## Value

Output the data frame containing information of fitting the regression
model: the estimated coefficients, standard errors, the deviance of the
different constrained models and the corresponding degrees of freedom.

## Author

Siru Wang
