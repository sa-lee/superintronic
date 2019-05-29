#' Compute diagnoistics over Ranges (rango)
#' 
#' @param x a Ranges object 
#' @param .var A meta data column in `x` to summarise over 
#' @param .funs A list of functions to summarise
#' @param ... additional parameters to pass to `.funs`
#' 
#' 
#' @importClassesFrom plyranges GroupedGenomicRanges GroupedIntegerRanges
#' @importClassesFrom IRanges Ranges
#' 
#' @export
#' @rdname range-diagnostics
setGeneric("rango", function(x, .var, .funs, ...) {
  standardGeneric("rango")
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
setMethod("rango", "GenomicRanges", 
          function(x, .var, .funs, ...) {
            message("No grouping variable available, using seqnames.")
            .var <- rlang::enquo(.var)
            rango(group_by(x, seqnames), !!.var, .funs, ...)
          })


#' @rdname range-diagnostics
#'@export
setMethod("rango", "GroupedGenomicRanges", 
          function(x, .var, .funs, ...) {
     
            .var <- rlang::enquo(.var)
            .var_c <- as.character(rlang::quo_get_expr(.var))
            sub <- plyranges::select(x, !!.var)
            y <- split_ranges(sub)
            views <- make_views_list(y, .var_c)
            # need to check name input for list
            res <- lapply(
              .funs, 
              function(.f) { viewApply(views, .f, ...) }
            )
            res <- lapply(res, function(x) {
              ans <- unlist(x)
              if (is(ans, "list")) {
                return(as(ans, "List"))
              }
              return(ans)
            })
            res <- DataFrame(res)
            groups <- mcols(x@inx)
            ans <- unlist(reduce(y, ignore.strand = TRUE))
            mcols(ans) <- cbind(groups, res)
            ans
            
})
