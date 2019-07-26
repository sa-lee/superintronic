#' Map a function over an Rle in windows
#' 
#' @param .x An atomic vector or Rle object.
#' @param .size The (integer) size of the window.
#' @param .step The (integer) amount to shift the start of the window by.
#' @param .fun A function
#' @param ... Additional arguments passed on to the mapped function
#'
#' @details 
#' The map functions apply a function over non-overlapping windows `[tile_rle()]`,
#' overlapping windows `[roll_rle()]`, and windows with a fixed start but
#' increasing width [`stretch_rle()`]
#' 
#' @examples 
#' x <- Rle(1:10, 1:10)
#' tile_rle(x, .size = 2, mean)
#' roll_rle(1:5, .size = 2, .step = 1, mean)
#' stretch_map(1:5, .size = 1, .size = 2, mean)
#' 
#' 
#' @importFrom methods as 
#' @export
#' @rdname windows
tile_rle <- function(.x, .size = 1L, .fun, ...) {
  .check_tsibble()
  stopifnot(length(.size) == 1 | length(.size) == length(.x))
  if (is(.x, "List")) {
    res <- as(mapply(function(.y, .size) .tile(.y, .size, .fun = .fun, ...),
                     .x, .size, SIMPLIFY = FALSE),
              "RleList")
              
  } else {
    res <- .tile(.x, .size = .size, .fun = .fun, ...)
  }
  res
}

#' @importFrom methods is
.autosize <- function(x) {
  if(is(x, "List")) return(pmax(as.integer(lengths(x) / 30), 1))
  pmax(as.integer(length(x) / 30), 1)
}

.tile <- function(.x, .size, .fun, ...) {
  as(unlist(tsibble::tile(.x, .f = .fun, ..., .size = .size)), "Rle")
}
#' @rdname windows
#' @export
stretch_rle <- function(.x, .size = 1L, .step = 1L, .fun, ...) {
  .check_tsibble()
  if (is(.x, "List")) {
    res <- as(lapply(.x, 
                     function(.y) .stretch(.y, .size = .size, ,step = .step, .fun = .fun, ...)),
              "RleList")
    
  } else {
    res <- .stretch(.x, .size = .size, .step = .step, .fun = .fun, ...)
  }
  res
}

.stretch <- function(.x, .size, .step, .fun, ...) {
  as(unlist(tsibble::stretch(.x, .f = .fun, ..., .init = .size, .step = .step)),
     "Rle")
}

#' @rdname windows
#' @export
roll_rle <- function(.x, .size = 1L, .step = 1L, .fun, ...) {
  .check_tsibble()
  if (is(.x, "List")) {
    res <- as(lapply(.x, 
                     function(.y) .roll(.y, .size = .size, ,step = .step, .fun = .fun, ...)),
              "RleList")
    
  } else {
    res <- .roll(.x, .size = .size, .step = .step, .fun = .fun, ...)
  }
  res
}

.roll <- function(.x, .size, .step, .fun, ...) {
  as(unlist(tsibble::slide(.x, .f = .fun, ..., .step = .step, .size = .size)), 
     "Rle")
}

.check_tsibble <- function() {
  if (!requireNamespace("tsibble", quietly = TRUE)) {
    stop("Package: tsibble required to run windowing functions. 
         Please install it.")
  } 
}