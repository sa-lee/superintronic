#' Recommended filters for analysis
#' 
#' @param .data A GRanges object obtained from `prepare_annotation()`
#' @param ... A set of logical expressions to filter `.data` if 
#' no expressions are provided the default rules are applied.
#' 
#' @details This runs the default set of filters before adding coverage
#' information, which are:
#' 
#' * Protein coding genes only
#' * Does not overlap any other genes
#' * Greater than 1 exonic feature
#' 
filter_rules <- function(.data, ...) {
  rules <- ...length()
  
  if (length(rules) == 0L) {
    filter(.data, 
           gene_type == "protein_coding", 
           n_self_olaps == 0L,
           lengths(simple_exonic) > 1)
  }
  filter(.data, ...)
}