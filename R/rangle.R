#' Rangle new features out of .ranges
#' 
#' @param x a Ranges object 
#' @param .var A metadata column in `x` to summarise over either a bare variable name
#' or a character vector of length 1.
#' @param .index A list of columns generated with `rng_vars()` 
#' that index the ranges in `x`.
#' @param .funs A named-list of functions that compute `rangenostics()`
#' 
#' 
#' @importClassesFrom plyranges GroupedGenomicRanges GroupedIntegerRanges
#' @importClassesFrom IRanges Ranges
#' 
#' @examples
#' suppressPackageStartupMessages(library("GenomicRanges"))
#' lvls <- paste0("chr", 1:23)
#' N <- 1e5
#' gr <- GRanges(
#' seqnames = factor(sample(lvls, N, replace = TRUE), levels = lvls),
#' ranges = IRanges(start = rpois(N, 10000), width = rpois(N, 100)),
#' grp = sample(letters[1:4], N, replace = TRUE),
#' score = rnbinom(N, size = 5, mu = 25)
#' )
#' 
#' rangle(gr, score, rng_vars(seqnames), list(mean = mean))
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
            if (!rlang::is_quosures(.index)) {
              .index <- rlang::enquos(.index)
            }
            
            rangle(group_by(x, !!!.index), 
                   !!.var, 
                   .index, 
                   .funs)
            
          })


#' @rdname range-diagnostics
#'@export
setMethod("rangle", "GroupedGenomicRanges", 
          function(x, .var, .index, .funs) {
            .var <- rlang::enquo(.var)
            .var_c <- as.character(rlang::quo_get_expr(.var))
            if (!rlang::is_quosures(.index)) {
              .index <- rlang::enquos(.index)
            }
            .regroups <- set_regroups(x, .index)
            # this would be faster if plyranges had an add = TRUE on groups
            x <- group_by(ungroup(x), !!!.regroups)
            y <- split_ranges(x)
            rle_list <- make_rlelist(y, .var_c)
            .funs <- check_funs(.funs, .var_c)
            
            res <- lapply(
              .funs, 
              function(.f) { mapper(rle_list, .f)}
            )
            
            res <- lapply(res, function(x) {
              if(is(x, "list")) {
                return(dplyr::bind_rows(x))
              }
              x
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


check_funs <- function(.funs, .var_c) {
  stopifnot(is(.funs, "list"))
  stopifnot(!any(names(.funs) == ""))
  stopifnot(length(unique(names(.funs))) == length(.funs))
  names(.funs) <- paste0(.var_c, "_", names(.funs))
  .funs
}

#' Select range variables
#' 
#' @param ... a set of variables 
#' 
#' @description This helper function is used to provide semantics
#' for selecting variables. It is useful for specifying the index 
#' in `rangle()` or passing to scoped variants in `dplyr` like 
#' `dplyr::mutate_at()`.
#' 
#' @return A list of quosures
#' @export
rng_vars <- function(...) {
  rlang::quos(...)
}

make_rlelist <- function(x, .var) {
  # ignoring self overlaps
  inx <- BiocGenerics::order(x)
  var <- as(lapply(mcols(x, level = "within"), 
                   function(.) .[[.var]]), "List")[inx]
  bases <- width(x)[inx]
  as(mapply(S4Vectors::Rle, 
            var, 
            bases), 
     "RleList")
}

set_regroups <- function(x, .index) {
  .groups <- plyranges::group_vars(x)
  if (!rlang::is_quosures(.index)) {
    .index <- rlang::enquos(.index)
  } 
  .index_c <- vapply(.index, 
                     function(q) as.character(rlang::quo_get_expr(q)),
                     character(1)
  )
  .regroups <- union(.groups, .index_c)
  rlang::syms(.regroups)
}

mapper <- function(x, .fun) {
  .fun <- rlang::as_function(.fun)
  
  res <- try(.fun(x), silent = TRUE)
  
  if (is(res, "try-error")) {
    res <- lapply(x, function(.) .fun(.))
  }
  return(res)
}
