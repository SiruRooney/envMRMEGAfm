# The meta-analysis output from MR-MEGA

A meta-analysis results of genetic variants across all the genome. This
approach uses genome-wide metrics of diversity between populations to
derive axes of genetic variation via multi-dimensional scaling. Allelic
effects of a variant across GWAS, weighted by their corresponding
standard errors, can then be modelled in a linear regression framework,
including the axes of genetic variation as covariates.

## Usage

``` r
MRMEGAout
```

## Format

A dataframe with 301 genetic variants and 17 variables:

- MARKERNAME:

  unique marker identification across input files.

- beta0:

  effect of the intercept of meta-regression.

- se0:

  std error of the intercept of meta-regression.

- beta_i:

  effect of the i-th PC of meta-regression.

- se_i:

  std error of the effect of the i-th PC of meta-regression.

- chisq_association:

  chisq value of the association.

- ndf_association:

  the number of degrees of freedom of the association.

- pvalue_association:

  p-value of the association.

- chisq_anceheter:

  chisq value of the heterogeneity due to different ancestry.

- ndf_anceheter:

  the number of degrees of freedom of ancestral heterogeneity.

- pvalue_anceheter:

  p-value of the ancestral heterogeneity.

- chisq_residual:

  chisq value of the residual heterogeneity.

- ndf_residual:

  the number of degrees of freedom of the residual heterogeneity.

- pvalue_residual:

  p-value of the residual heterogeneity.

- logBF:

  log of Bayesian Factors.
