#'  Cast a Ranges to a RangesList 
#'  
#' @param x a Ranges object
#' @param ... grouping variables to split by (passed to `plyranges::group_by()`)
#' @return a RangesList object
#' 
#' @importClassesFrom plyranges GroupedGenomicRanges GroupedIntegerRanges
#' @importClassesFrom IRanges Ranges
#' 
#' @export
#' @rdname split_ranges 
setGeneric("split_ranges", 
           function(x, ...) {
             standardGeneric("split_ranges")
           }
)


#' @rdname split_ranges
#' @export
setMethod("split_ranges", "Ranges",  function(x, ...) {
  .split_ranges(plyranges::group_by(x, ...))
})


.split_ranges <- function(x, ...) {
  grps <- plyranges::groups(x)
  ans <- S4Vectors::split(x@delegate, 
                   plyranges::select(x, !!!grps, .drop_ranges = TRUE),
                   drop = TRUE)
  unname(ans)
}

#' @rdname split_ranges 
#' @export
setMethod("split_ranges", "GroupedGenomicRanges", function(x, ...) {
  .split_ranges(x, ...)
})

#' @rdname split_ranges 
#' @export
setMethod("split_ranges", "GroupedIntegerRanges", function(x, ...) {
  .split_ranges(x, ...)
})

