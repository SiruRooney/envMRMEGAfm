# The environment-adjusted meta-analysis output from env-MR-MEGA

The environment-adjusted meta-analysis results of genetic variants
across all the genome. Considering the impact of environmental exposures
that differ across GWAS, the environment-adjusted meta-regression model
is build upon the MR-MEGA meta-regression framework by adding
study-level environmental covariates. This allows us to identify genetic
variants that associated with the disease or trait while adjusting for
differing environmentl exposures between cohort.

## Usage

``` r
envMRMEGAout
```

## Format

A dataframe with 301 genetic variants and 27 variables:

- MARKERNAME:

  unique marker identification across input files.

- beta0:

  effect of the intercept of meta-regression.

- se0:

  std error of the intercept of meta-regression.

- beta_i:

  effect of the i-th PC or i-th environment covariate of
  environment-adjusted meta-regression.

- se_i:

  std error of the effect of the i-th PC or i-th environment covariate
  of environment-adjusted meta-regression.

- chisq_association:

  chisq value of the association.

- ndf_association:

  the number of degrees of freedom of the association.

- pvalue_association:

  p-value of the association.

- chisq_heter:

  chisq value of the ancestral and environmental heterogeneity.

- ndf_heter:

  the number of degrees of freedom of ancestral and environmental
  heterogeneity.

- pvalue_heter:

  p-value of ancestral and environmental heterogeneity.

- chisq_residual:

  chisq value of the residual heterogeneity.

- ndf_residual:

  the number of degrees of freedom of the residual heterogeneity.

- pvalue_residual:

  p-value of the residual heterogeneity.

- logBF:

  log of Bayesian Factors.

- chisq_env:

  chisq value of environmental heterogeneity.

- ndf_env:

  the number of degrees of freedom of environmental heterogeneity.

- pvalue_env:

  p-value of environmental heterogeneity.

- chisq_PC:

  chisq value of ancestral heterogeneity.

- ndf_PC:

  the number of degrees of freedom of ancestral heterogeneity.

- pvalue_PC:

  p-value of ancestral heterogeneity.
