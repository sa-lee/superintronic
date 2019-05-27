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

setMethod("rango", "Ranges", 
          function(x, .var, .funs, ...) {
            .var <- rlang::ensym(.var)
            sub <- plyranges::select(x, !!.var)
            message("No grouping variable given, using seqnames instead.")
            x <- split_ranges(x, seqnames)
            views <- make_views_list(x, as.character(.var))
            lapply(.funs, function(.f) IRanges::viewApply(views, .f, ...))
            
          })


setMethod("rango", "GroupedGenomicRanges", 
          function(x, .var, .funs, ...) {
            .var <- rlang::ensym(.var)
            sub <- plyranges::select(x, !!.var)
            x <- split_ranges(x)
            views <- make_views_list(x, as.character(.var))
            lapply(.funs, function(.f) IRanges::viewApply(views, .f, ...))
  
})