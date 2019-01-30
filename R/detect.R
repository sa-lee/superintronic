#' Intron Signal Aggregator
#'
#'
#' @param x a GRanges object from `merge_coverage()`
#' @param design a DataFrame containing a primary key to join to `x`
#' @param ... columns names to partiton x over
#'
#' @details
#'
#' BAM processing
#' 1. Over each gene compute the log2 trimmed mean for exon and intron counts
#' and flag whether the gene is expressed (>= 3 reads in both exon and intron)
#' 2. Compute the run length of intron ranges exceeding a nominal threshold
#'
#'
#' @return a GRanges object with the following metadata columns
#' * GeneId
#' * IsExpressed: a flag for whether gene is expressed or not
#' * ExonLength
#' * IntronLength
#' * MeanExonCoverage
#' * MeanIntronCoverage
#'
#'
#' @export
summarise_ir <- function(x, design, ...) {
  # on dots
  # no dots just summarise coverage over gene over type
  # dots same but grouped by dots, should function actually just
  # be a wrapper to group by 
  
}



weighted_coverage <- function(score, width) sum(score*width) / sum(width)

trimmed_coverage <- function() NULL