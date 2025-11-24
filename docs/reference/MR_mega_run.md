# environment-adjusted MR-MEGA approach or MR-MEGA approach.

Fit the (environment-adjusted) meta-regression model using the
normalized statistic data across all cohorts.

## Usage

``` r
MR_mega_run(beta_pop, invse2_pop, cohort_count_filt, env, PCs, ncores)
```

## Arguments

- beta_pop:

  The estimated allelic effects of the target SNP conditioned on the SNP
  set across all the cohorts.

- invse2_pop:

  The inverse of the standard errors of the target SNP conditioned on
  the SNP set across all the cohorts.

- cohort_count_filt:

  The number of presence for the target SNP across all the cohorts.

- env:

  The study-level environment factors across all the cohorts. Each row
  refers to one population and each column refers to one environment
  covariate. For MR-MEGA approach, env=NULL.

- PCs:

  The axes of genetic variation, which can also be called the principal
  components (PCs). Each row refers to one population. Note: Each `env`
  row and `PCs` row should correspond to same population.

- ncores:

  The number of cores which would be used for running in parallel.

## Value

Output a file containing names of genetic variants, estimated
coefficients, standard errors, chisq value of the association, the
number of degrees of freedom of the association, p-value of the
association, chisq value of the heterogeneity due to different ancestry,
ndf of the heterogeneity due to different ancestry, p-value of the
heterogeneity due to different ancestry, chisq value of the residual
heterogeneity, ndf of the residual heterogeneity, p-value of the
residual heterogeneity.

## Author

Siru Wang
