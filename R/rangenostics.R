#' Compute diagnostics over Ranges
#' 
#' @param x a Ranges object 
#' @param .var A meta data column in `x` to summarise over either a bare variable name
#' or a character vector of length 1.
#' @param .funs A list of functions to summarise
#' @param ... additional parameters to pass to `.funs`
#' 
#' 
#' @importClassesFrom plyranges GroupedGenomicRanges GroupedIntegerRanges
#' @importClassesFrom IRanges Ranges
#' 
#' @export
#' @rdname range-diagnostics
setGeneric("rangle", function(x, .var, .funs, ...) {
  standardGeneric("rangle")
})

make_views <- function(x, .var) {
  subject <- S4Vectors::Rle(mcols(x)[[.var]], lengths = BiocGenerics::width(x))
  rng <- reduce(ranges(subject))
  IRanges::Views(subject, start = rng)
} 

make_views_list <- function(x, .var) {
  as(lapply(x, make_views, .var = .var), "List")
}

#' @rdname range-diagnostics
#' @export
setMethod("rangle", "GenomicRanges", 
          function(x, .var, .funs, ...) {
            message("No grouping variable available, using seqnames.")
            .var <- rlang::enquo(.var)
            rangle(group_by(x, seqnames), !!.var, .funs, ...)
          })


#' @rdname range-diagnostics
#'@export
setMethod("rangle", "GroupedGenomicRanges", 
          function(x, .var, .funs, ...) {
     
            .var <- rlang::enquo(.var)
            .var_c <- as.character(rlang::quo_get_expr(.var))
            y <- split_ranges(x)
            views <- make_views_list(y, .var_c)
            # need to check name input for list
            res <- lapply(
              .funs, 
              function(.f) { viewApply(views, .f) }
            )
            res <- lapply(res, function(x) {
              ans <- unlist(x)
              if (is(ans, "list")) {
                return(as(ans, "List"))
              }
              return(ans)
            })
            
            res <- S4Vectors::DataFrame(res)
            regions <- summarise(x, start = min(start), end = max(end))
            res <- cbind(regions, res)
            GenomicRanges::makeGRangesFromDataFrame(res, keep.extra.columns = TRUE)
            
})
