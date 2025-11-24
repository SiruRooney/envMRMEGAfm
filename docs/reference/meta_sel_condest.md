# Conditional analysis of the linear regression model

Estimate the allelic effect size and the standard errors of the target
SNP conditional on the SNP set.

## Usage

``` r
meta_sel_condest(
  add.gwas,
  sel.gwas,
  addsel.ld,
  sel.ld,
  Vp.med,
  collinear,
  actual.geno
)
```

## Arguments

- add.gwas:

  A dataframe which contains the GWAS information for the target SNP.
  GWAS information comprise
  "MARKERNAME","CHROMOSOME","POSITION","EA","NEA","EAF","N","BETA","SE".

- sel.gwas:

  A dataframe which contains the GWAS information for SNP set.

- addsel.ld:

  A LD matrix containing the correlation matrix between the target SNP
  and the SNP set.

- sel.ld:

  A LD matrix containing the correlation between the SNP set.

- Vp.med:

  The median of variance of phenotype. Refer to GCTA-COJO.

- collinear:

  A threshold to filter out the target SNP in high LD with the SNP set.
  If the squared multiple correlation between the target SNP exceeds the
  threshold, such as 0.9, the target SNP is ignored.

- actual.geno:

  An indicator to specify whether the true cohort-level LD structure is
  applied. If actual.geno is TRUE, the inputted LD structures are
  derived from the true individual-level genotype data. If actual.geno
  is FALSE, the the inputted LD structures are derived from the
  approximation LD strucure.

## Value

Output a data.frame containing the estimated allelic effect and the
inverse of squared standard error for the target SNP.

## Author

Siru Wang
