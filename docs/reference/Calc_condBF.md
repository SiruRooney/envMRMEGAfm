# Conditional Bayesian factors

Calculate the Bayesian factor for each target SNP conditioned on the SNP
set

## Usage

``` r
Calc_condBF(
  gwas.list,
  ld.list,
  which.ld,
  meta.snp,
  sel.set,
  n_cohort,
  env,
  PCs,
  out_loc,
  ncores,
  cred.thr,
  actual.geno,
  collinear
)
```

## Arguments

- gwas.list:

  A list of length K_g which contains the pre-processed GWAS files. Each
  component contains one GWAS file comprising these required columns:
  "MARKERNAME","CHROMOSOME","POSITION","EA","NEA","EAF","N","BETA","SE".

- ld.list:

  A list of length K_ld which contains the pre-processed LD structures.
  Each component contains one LD structure.

- which.ld:

  A character vector of length K_g. Each component of the vector
  corresponds to one LD structure. The length of which.ld should equal
  to the number of gwas files.

- meta.snp:

  A character vector containing all SNP names present in the output of
  MR-MEGA method or env-MR-MEGA method.

- sel.set:

  A set contains all potential associated SNPs.

- n_cohort:

  The number of cohorts.

- env:

  The study-level environment factors across all the cohorts. Each row
  refers to one population and each column refers to one environment
  covariate. For MR-MEGA approach, env=NULL.

- PCs:

  The axes of genetic variation, which can also be called the principal
  components (PCs). Each row refers to one population. Note: Each `env`
  row and `PCs` row should correspond to same population.

- out_loc:

  Path to save pre-processed GWAS files and LD structures. By default,
  out_loc=NULL.

- ncores:

  The the number of cores which would be used for running in parallel.

- cred.thr:

  Credible threshold for the credible set for each selected potential
  SNP. By default, cred.thr=0.99 refers to 99% credible sets.

- actual.geno:

  An indicator to specify whether the true cohort-level LD structure is
  applied.

- collinear:

  A threshold to filter out the target SNP in high LD with the SNP set.
  If the squared multiple correlation between the target SNP exceeds the
  threshold, such as 0.9, the target SNP is ignored.

## Value

Output the results of (env-)MR-MEGAfm conditioned on the subset of the
selected potential SNP set.

## Author

Siru Wang
