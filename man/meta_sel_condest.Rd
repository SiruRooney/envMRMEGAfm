% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LRanalysis.R
\name{meta_sel_condest}
\alias{meta_sel_condest}
\title{Conditional analysis of the linear regression model}
\usage{
meta_sel_condest(
  add.gwas,
  sel.gwas,
  addsel.ld,
  sel.ld,
  Vp.med,
  collinear,
  actual.geno
)
}
\arguments{
\item{add.gwas}{A dataframe which contains the GWAS information for the target SNP. GWAS information comprise "MARKERNAME","CHROMOSOME","POSITION","EA","NEA","EAF","N","BETA","SE".}

\item{sel.gwas}{A dataframe which contains the GWAS information for SNP set.}

\item{addsel.ld}{A LD matrix containing the correlation matrix between the target SNP and the SNP set.}

\item{sel.ld}{A LD matrix containing the correlation between the SNP set.}

\item{Vp.med}{The median of variance of phenotype. Refer to GCTA-COJO.}

\item{collinear}{A threshold to filter out the target SNP in high LD with the SNP set. If the squared multiple correlation between the target SNP exceeds the threshold,
such as 0.9, the target SNP is ignored.}

\item{actual.geno}{An indicator to specify whether the true cohort-level LD structure is applied. If actual.geno is TRUE, the inputted LD structures are derived from the true individual-level genotype data.
If actual.geno is FALSE, the the inputted LD structures are derived from the approximation LD strucure.}
}
\value{
Output a data.frame containing the estimated allelic effect and  the inverse of squared standard error for the target SNP.
}
\description{
Estimate the allelic effect size and the standard errors of the target SNP conditional on the SNP set.
}
\author{
Siru Wang
}
