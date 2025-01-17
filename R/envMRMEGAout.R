#' The environment-adjusted meta-analysis output from env-MR-MEGA
#'
#' @description The environment-adjusted meta-analysis results of genetic variants across all the genome. Considering the impact of environmental exposures that differ across GWAS,
#' the environment-adjusted meta-regression model is build upon the MR-MEGA meta-regression framework by adding study-level environmental covariates.
#' This allows us to identify genetic variants that associated with the disease or trait while adjusting for differing environmentl exposures between cohort.
#'
#'@format A dataframe with 301 genetic variants and 27 variables:
#'\describe{
#'\item{MARKERNAME}{unique marker identification across input files.}
#'\item{beta0}{effect of the intercept of meta-regression.}
#'\item{se0}{std error of the intercept of meta-regression.}
#'\item{beta_i}{effect of the i-th PC or i-th environment covariate of environment-adjusted  meta-regression.}
#'\item{se_i}{std error of the effect of the i-th PC or i-th environment covariate of environment-adjusted meta-regression.}
#'\item{chisq_association}{chisq value of the association.}
#'\item{ndf_association}{the number of degrees of freedom of the association.}
#'\item{pvalue_association}{p-value of the association.}
#'\item{chisq_heter}{chisq value of the ancestral and environmental heterogeneity.}
#'\item{ndf_heter}{the number of degrees of freedom of ancestral and environmental heterogeneity.}
#'\item{pvalue_heter}{p-value of ancestral and environmental heterogeneity.}
#'\item{chisq_residual}{chisq value of the residual heterogeneity.}
#'\item{ndf_residual}{the number of degrees of freedom of the residual heterogeneity.}
#'\item{pvalue_residual}{p-value of the residual heterogeneity.}
#'\item{logBF}{log of Bayesian Factors.}
#'\item{chisq_env}{chisq value of environmental heterogeneity.}
#'\item{ndf_env}{the number of degrees of freedom of environmental heterogeneity.}
#'\item{pvalue_env}{p-value of environmental heterogeneity.}
#'\item{chisq_PC}{chisq value of ancestral heterogeneity.}
#'\item{ndf_PC}{the number of degrees of freedom of ancestral heterogeneity.}
#'\item{pvalue_PC}{p-value of ancestral heterogeneity.}
#'}
#'
#'
"envMRMEGAout"
