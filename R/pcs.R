#' Two axes of genetic variation.
#'
#' @description
#'  The data frame contains two principal components (PCs). In this simulation example, we focus on 10 sex-stratified continental African populations, excluding the ACB and ASW populations.
#'
#' To derive the axes of genetic variation, we used a subset of genetic variants with minor allele frequency (MAF) >5\% across all populations.
#' These variants were separated by 1Mb in the genome. For each reference panel, we calculated the matrix of pairwise Euclidean distances
#' between the populations using the selected genetic variants. The two axes of genetic variation were derived through the multi-dimensional
#' scaling of the distance matrix in the corresponding reference panel.
#' @source 5 populations of African ancestry were collected from Phase 3 of 1000 Genomes \url{https://ctg.cncr.nl/software/MAGMA/ref_data/}.
"pcs"
