#' Intron Signal Aggregator
#'
#'
#' @param x a GRanges object from `merge_coverage()\merge_design()`
#' @param ... columns to partition x 
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
summarise_intronic <- function(x, ..., threshold = 0L) {
  # on dots
  dots <- rlang::enquos(...)
  
  default_cols <- 
    rlang::quos(
      total_coverage = sum(score * width),
      average_coverage = mean_coverage(score, width),
      n_bases_above_threshold = sum(width[score > !!threshold])
    )
  # no dots just summarise coverage as follows
  # coverage over each exon/intron in a gene
  # this is computed as weighted sum 
  if (length(dots) == 0L) {
    x %>% 
      plyranges::group_by(gene_id, feature_type, feature_strand) %>% 
      plyranges::reduce_ranges(
        average_coverage = mean_coverage(score, width),
        
        
      )
  }
  # dots same but add partitions
  
}





# helper fun for width 
mean_coverage <- function(score, width) sum(score*width) / sum(width)

