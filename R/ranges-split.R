#'  Cast a Ranges to a RangesList 
#'  
#' @param x a Ranges object
#' @param ... grouping variables to split by (passed to `plyranges::group_by()`)
#' @return a RangesList object
#' 
#' @importClassesFrom plyranges GroupedGenomicRanges GroupedIntegerRanges
#' @importClassesFrom IRanges Ranges
#' 
#' @details If `x` is already grouped, then passing `...` are ignored. If `x`
#' is not grouped and variable names are provided to `...` these are forwarded
#' to `[plyranges::group_by()]`, if no variables are provided then x is returned.
#' 
#' @examples 
#' suppressPackageStartupMessages(library("GenomicRanges"))
#' suppressPackageStartupMessages(library("plyranges"))
#' lvls <- paste0("chr",c(1:23))
#' gr <- GRanges(
#'   seqnames = factor(sample(lvls, 1000, replace = TRUE), levels = lvls),
#'   ranges = IRanges(start = rpois(1000, 10000), width = rpois(1000, 100)),
#'   grp = sample(letters[1:4], 1000, replace = TRUE),
#'   gc = runif(1000)
#' )
#' split_ranges(gr, seqnames)
#' 
#' gr_by_grp <- group_by(gr, grp)
#' split_ranges(gr_by_grp)   
#' 
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
  grps <- rlang::enquos(...)
  if (length(grps) == 0) return(x)
  
  .split_ranges(plyranges::group_by(x, !!!grps))
})


.split_ranges <- function(x, ...) {
  dots <-  rlang::enquos(...)
  if (length(dots) > 0) {
    warning(paste("Ignoring arguments in `...`"))
  }
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

