# Joint analysis of a linear regression model

Estimate the allelic effect size and the standard errors of the target
SNP in the single SNP linear model.

## Usage

``` r
meta_sel_joinest(sel.gwas, sel.ld, actual.geno)
```

## Arguments

- sel.gwas:

  A dataframe which contains the GWAS information for SNP set. GWAS
  information comprise
  "MARKERNAME","CHROMOSOME","POSITION","EA","NEA","EAF","N","BETA","SE".

- sel.ld:

  A LD matrix containing the correlation between the SNP set.

- actual.geno:

  An indicator to specify whether the true cohort-level LD structure is
  applied. If actual.geno is TRUE, the inputted LD structures are
  derived from the true individual-level genotype data.

## Value

Output a data.frame containing the estimated allelic effect, standard
error, the inverse of squared standard error and the p-value of
association test for the target SNP.

## Author

Siru Wang
