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

make_rle <- function(x, .var) {
  S4Vectors::Rle(mcols(x)[[.var]], lengths = BiocGenerics::width(x))
} 

make_views <- function(x, .var) {
  # rle-ise the grouping
  rle_list <- as(lapply(x, make_rle, .var = .var), "List")
  
  # for each group we merge the overlapping ranges 
  ranges_list <- reduce(IRanges::ranges(x), ignore.strand = TRUE)
  
  IRanges::RleViewsList(rleList = rle_list, rangesList = ranges_list)
  
}

setMethod("rango", "Ranges", 
          function(x, .var, .funs, ...) {
            .var <- rlang::ensym(.var)
            sub <- plyranges::select(x, !!.var)
            message("No grouping variable given, using seqnames instead.")
            x <- split_ranges(x, seqnames)
            views <- make_views(x, as.character(.var))
            lapply(.funs, function(.f) IRanges::viewApply(views, .f, ...))
            
          })


setMethod("rango", "GroupedGenomicRanges", 
          function(x, .var, .funs, ...) {
            .var <- rlang::ensym(.var)
            sub <- plyranges::select(x, !!.var)
            x <- split_ranges(x)
            views <- make_views(x, as.character(.var))
            lapply(.funs, function(.f) IRanges::viewApply(views, .f, ...))
  
})