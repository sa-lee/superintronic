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
            rango(group_by(x, seqnames), .var, .funs, ...)
          })


#' @rdname range-diagnostics
#'@export
setMethod("rango", "GroupedGenomicRanges", 
          function(x, .var, .funs, ...) {
            .var <- rlang::ensym(.var)
            sub <- plyranges::select(x, !!.var)
            y <- split_ranges(x)
            views <- make_views_list(y, as.character(.var))
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



# not look in to dplyr:::summarise_at ; manip_at and as_fun_list

tile_view <- function(x, width) {
  if (is(x, "Views")) {
    ln <- length(IRanges::subject(x))
  } else {
    ln <- length(x)
  }
  trim(successiveViews(x, rep.int(width, ln %/% width + 1)))
}

roll_view <- function(x, width, step) {
  
  if (is(x, "Views")) {
    ln <- length(IRanges::subject(x))
  } else {
    ln <- length(x)
  }
  
  rng <- IRanges(start = seq.int(1, ln, by = step), width = width)
  
  trim(Views(x, start = rng))
  
}

roll_map <- function(x, ...) {
  
}

stretch_view <- function(x, width, step) {
  if (is(x, "Views")) {
    ln <- length(IRanges::subject(x))
  } else {
    ln <- length(x)
  }
  
  rng <- IRanges(start = 1, width =  seq.int(width, ln, by = step))
  
  trim(Views(x, start = rng))
  
}





