# (env-)MR-MEGA fine-mapping approach (env-MR-MEGAfm)

Fine-mapping genetic associations across multiple cohorts, allowing for
multiple causal variants.

## Usage

``` r
env_MR_MEGA_fm(
  gwas.list,
  ld.list,
  which.ld,
  meta.file,
  PCs,
  env,
  out_loc = NULL,
  ncores = 1,
  collinear = 0.9,
  pvalue_cutoff = 5e-08,
  cred.thr = 0.99,
  actual.geno = FALSE
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

- meta.file:

  A data frame obtained from the output of MR-MEGA method or env-MR-MEGA
  method.

- PCs:

  The axes of genetic variation, which can also be called the principal
  components (PCs). Each row refers to one population and each column
  refers to one axes of genetic variation. The order of matrix row in
  `PCs` should be consistent with the order of gwas names.

- env:

  The matrix which contains the study-level environmental factors. Each
  row refers to on population and each column refers to one environment
  covariate. The order of matrix row in `env` should be consistent with
  the order of gwas names. If MR-MEGAfm is utilized, env=NULL.

- out_loc:

  Path to save the output of (env-)MR-MEGAfm. By default, out_loc=NULL.

- ncores:

  The the number of cores which would be used for running in parallel.

- collinear:

  A threshold to filter out the target SNP in high LD with the SNP set.
  If the squared multiple correlation between the target SNP exceeds the
  threshold, such as 0.9, the target SNP is ignored.

- pvalue_cutoff:

  A cutoff p-value to identify the potential associated SNP with the
  smallest p-value exceeding the cutoff p-value. By default,
  pvalue_cutoff=5e-8.

- cred.thr:

  Credible threshold for the credible set for each selected potential
  SNP. By default, cred.thr=0.99 refers to 99% credible sets.

- actual.geno:

  An indicator to specify whether the true cohort-level LD structure is
  applied. If actual.geno is TRUE, the inputted LD structures are
  derived from the true individual-level genotype data. By default,
  actual.geno=FALSE is recommended in GCTA-COJO software.

## Value

A list of length 2. One component contains the selected potential
associated SNP. Another component contains their credible sets with the
specified thresholds. By default, the credible threshold is set to 0.99.

## Author

Siru Wang
