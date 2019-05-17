#' Compute diagnoistics over Ranges (rangnostics)
#' 
#' @param ranges a GRanges object 
#' @param .var A column in `ranges` defining an index
#' @param .funs A list of functions to summarise
#' @param ... additional parameters to pass to ...
#' 
#' @export
setGeneric("rangenostics", function(ranges, .var, .funs, ...) {
  standardGeneric("rangenostics")
})

setMethod("rangenostics", "GenomicRanges", 
          function(ranges, .var, .funs, ...) {
            sub <- plyranges::select(ranges, !!.var)
            gr_l <- split(sub, seqnames(sub), drop = TRUE)
            # split by seqnames
            rle_list <- as(lapply(gr_l, 
                                 function(.y) {
                                   S4Vectors::Rle(mcols(.y)[[.var]], 
                                                  lengths = width(.y))
                                 }),
                             "List")
            ranges_list <- reduce(ranges(gr_l), 
                                  ignore.strand = TRUE)
            
            views_list <- RleViewsList(rleList = rle_list, 
                                       rangesList = ranges_list)

            lapply(.funs,
                   function(.f) IRanges::viewApply(views_list, .f)
                   )
            
          })