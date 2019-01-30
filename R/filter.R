#' Apply filters to exon/intron annotations
#' 
#' @param .data A GRanges object obtained from `prepare_annotation()`
#' @param ... A set of logical expressions to pass to `plyranges::filter()`. If
#' no expressions are passed then default filters are applied (see filter)
#' 
#' @details If no expressions are passed to `...` then the following
#' filters are applied by default:
#'
#' * Restrict genes to have `gene_type == "protein_coding"`
#' * Genes should not overlap any other genes `n_olaps == 1L`
#' * Genes should have at least 2 exons `lengths(simple_exonic) > 1`
#' 
#' @return a GRanges object
#' 
#' 
#' @export
filter_rules <- function(.data, ...) {
  rules <- ...length()
  if (rules == 0L) {
    return(filter(.data, 
           gene_type == "protein_coding", 
           n_olaps == 1L,
           lengths(simple_exonic) > 1))
  }
  
  filter(.data, ...)
}