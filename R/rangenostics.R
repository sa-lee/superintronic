#' Compute diagnostics over Ranges
#' 
#' @param x a Ranges object 
#' @param .var A metadata column in `x` to summarise over either a bare variable name
#' or a character vector of length 1.
#' @param .index Bare variable names that index each Range
#' @param .funs A named-list of functions that compute diagnostics.
#' 
#' 
#' @importClassesFrom plyranges GroupedGenomicRanges GroupedIntegerRanges
#' @importClassesFrom IRanges Ranges
#' 
#' @examples
#' gr <- GRanges(
#' seqnames = factor(sample(lvls, 1000, replace = TRUE), levels = lvls),
#' ranges = IRanges(start = rpois(1000, 10000), width = rpois(1000, 100)),
#' grp = sample(letters[1:4], 1000, replace = TRUE),
#' score = rnbinom(1000, size = 5, mu = 25)
#' )
#' 
#' rangle(gr, score, seqnames, list(mean = mean))
#' 
#' 
#' @export
#' @rdname range-diagnostics
setGeneric("rangle", function(x, .var, .index,  .funs) {
  standardGeneric("rangle")
})


#' @rdname range-diagnostics
#' @export
setMethod("rangle", "GenomicRanges", 
          function(x, .var, .index, .funs) {
            .var <- rlang::enquo(.var)
            .index <- rlang::enquo(.index)
            rangle(group_by(x, !!.index), !!.var, !!.index, .funs)
          })


#' @rdname range-diagnostics
#'@export
setMethod("rangle", "GroupedGenomicRanges", 
          function(x, .var, .index, .funs) {
            .var <- rlang::enquo(.var)
            .var_c <- as.character(rlang::quo_get_expr(.var))
            .regroups <- set_regroups(x, !!rlang::enquo(.index))
            x <- group_by(ungroup(x), !!!.regroups)
            y <- split_ranges(x)
            
            views <- make_views_list(y, .var_c)
            # need to check name input for list
            res <- lapply(
              .funs, 
              function(.f) { viewApply(views, .f, simplify = FALSE) }
            )
            res <- lapply(res, function(x) {
              ans <- unlist(x)
              if (is(ans, "list")) {
                return(as(ans, "List"))
              }
              return(ans)
            })
            
            res <- S4Vectors::DataFrame(res)
            regions <- summarise(x,
                                 seqnames_x = unlist(BiocGenerics::unique(seqnames)),
                                 strand_x = unlist(BiocGenerics::unique(strand)),
                                 start_x = min(start),
                                 end_x= max(end)
                                 )
            res <- cbind(regions, res)
            
            plyranges::as_granges(res, 
                                  seqnames = seqnames_x, 
                                  strand = strand_x,
                                  start = start_x,
                                  end = end_x)
            
})



rng_vars <- function(...) {
  rlang::quos(...)
}



make_views <- function(x, .var) {
  subject <- S4Vectors::Rle(mcols(x)[[.var]], lengths = BiocGenerics::width(x))
  rng <- reduce(ranges(subject))
  IRanges::Views(subject, start = rng)
} 

make_views_list <- function(x, .var) {
  as(lapply(x, make_views, .var = .var), "List")
}

set_regroups <- function(x, .index) {
  .groups <- plyranges::group_vars(x)
  .index <- rlang::enquo(.index)
  .index_c <- as.character(rlang::quo_get_expr(.index))
  .regroups <- union(.groups, .index_c)
  rlang::syms(.regroups)
}
