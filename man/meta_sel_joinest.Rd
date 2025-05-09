% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LRanalysis.R
\name{meta_sel_joinest}
\alias{meta_sel_joinest}
\title{Joint analysis of a linear regression model}
\usage{
meta_sel_joinest(sel.gwas, sel.ld, actual.geno)
}
\arguments{
\item{sel.gwas}{A dataframe which contains the GWAS information for SNP set. GWAS information comprise "MARKERNAME","CHROMOSOME","POSITION","EA","NEA","EAF","N","BETA","SE".}

\item{sel.ld}{A LD matrix containing the correlation between the SNP set.}

\item{actual.geno}{An indicator to specify whether the true cohort-level LD structure is applied. If actual.geno is TRUE, the inputted LD structures are derived from the true individual-level genotype data.}
}
\value{
Output a data.frame containing the estimated allelic effect, standard error, the inverse of squared standard error and the p-value of association test for the target SNP.
}
\description{
Estimate the allelic effect size and the standard errors of the target SNP in the single SNP linear model.
}
\author{
Siru Wang
}
