#' Squish coverage over intronic/exonic parts
#'
#'
#' @param x a GRanges object from `merge_coverage()\merge_design()`
#' @param ... columns to partition x 
#' @param threshold Filter bases covered to above some integer threshold (default = 0L).
#'
#' @details Here we use `plyranges::reduce_ranges()` to squish
#' the coverage GRanges to have coverage over each exon/intron within
#' a gene. If `...` are supplied additional grouping is added by metadata
#' columns in x, using `plyranges::group_by()`.
#'
#'
#'
#' @return a GRanges object, where each gene is split into exonic
#' and intronic parts, and the coverage is summarised as follows over each
#' exon / intron:
#' 
#' * `gene_id` the gene_id as from `prepare_annotation()`
#' * `feature_type` whether a feature is an intron or exon
#' * `feature_length` the total length of the intronic/exonic part
#' * `total_coverage`
#' * `average_coverage`
#' * `n_bases_above_threshold`
#' * `prop_bases_covered_at_threshold`
#'
#' @export
squish_coverage <- function(x, ..., threshold = 0L) {
  
  groups <- set_grouping(...)
  # coverage over each exon/intron in a gene
  # this is computed as weighted sum of n bases covered (width) 
  # multiplied by the score divided by the sum of the width
  # we also compute the total coverage
  # and the number of bases above some threshold
  # if we use `reduce_ranges` we get the exon/intron features
  # for each gene.
  # if there were no added groups we collapse over all genes
  # additional dots add partitions
  x %>% 
    plyranges::group_by(!!!groups) %>% 
    plyranges::reduce_ranges(
      !!!.signal_cols(threshold = threshold)
    ) %>% 
    plyranges::mutate(
      prop_bases_covered_at_threshold = n_bases_above_threshold / width
    )
}




#' Spread intron/exon coverage over a gene
#' 
#' @inheritParams squish_coverage
#'
#' @details This function summarises the *total* coverage within an exon/intron
#' for each gene, it then widens the result so there are seperate columns
#' for each exon/intron variable. 
#'
#'@export
spread_coverage <- function(x, ..., threshold  = 0L) {
  groups <- set_grouping(...)
  cols <- .signal_cols(threshold = threshold)
  res <- x %>% 
    plyranges::group_by(!!!groups) %>% 
    plyranges::summarise(
      !!!cols
    ) %>% 
    S4Vectors::transform(
      prop_bases_covered_at_threshold = n_bases_above_threshold / total_coverage
    )
  
  # spread on feature_type
  res <- S4Vectors::split(res, res[["feature_type"]])
  values <- BiocGenerics::do.call(S4Vectors::cbind, 
                                  res[, c(names(cols), "prop_bases_covered_at_threshold")])
  keys <- res[, sapply(groups, rlang::quo_name)][[1]]
  S4Vectors::cbind(keys, values)
}

# set default groupings
set_grouping <- function(...) {
  # on dots
  dots <- rlang::enquos(...)
  # splice with default groups
  rlang::quos(
    !!!dots,
    gene_id,
    feature_strand,
    feature_type
  ) 
}

# weighted coverage computation 
mean_coverage <- function(score, width) sum(score*width) / sum(width)

# default columns to compute in summary calculations
.signal_cols <- function(threshold = 0L) {
  rlang::quos(
    feature_length = sum(unique(feature_length)),
    total_coverage = sum(score * width),
    average_coverage = mean_coverage(score, width),
    n_bases_above_threshold = sum(width[score > !!threshold]), 
  )
}



