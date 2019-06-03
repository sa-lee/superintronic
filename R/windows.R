#' Map a function over windows
#' 
#' @param .x An atomic vector or Rle object.
#' @param .size The (integer) size of the window.
#' @param .step The (integer) amount to shift the start of the window by.
#' @param .start The (integer) starting position of the window
#' @param .fun A function
#' @param ... Additional arguments passed on to the mapped function
#'
#' @details 
#' The map functions apply a function over non-overlapping windows `[tile_map()]`,
#' overlapping windows `[roll_map()]`, and windows with a fixed start but
#' increasing width [`stretch_map()`]
#' 
#' @examples 
#' tile_map(1:5, .size = 2, mean)
#' roll_map(1:5, .size = 2, .step = 1, mean)
#' stretch_map(1:5, .size = 1, .size = 2, mean)
#' 
#' @export
#' @rdname windows
tile_map <- function(.x, .size, .fun, ...) {
  viewApply(tile_view(.x, .size), .fun, ...)
}

#' @rdname windows
#' @export
roll_map <- function(x, .size, .step, .fun, ...) {
  viewApply(roll_view(x, .size, .step), .fun)
}

#' @rdname windows
#' @export
stretch_map <- function(.x, .start, .step, .fun, ...) {
  viewApply(stretch_view(.x, .start, .step), .fun, ...)
}

set_ln <- function(x) {
  if (is(x, "Views")) return(length(IRanges::subject(x)))
  return(length(x))
}

tile_view <- function(x, width) {
  l <- set_ln(x)
  r <- l %% width
  w <- c(rep.int(width, l %/% width), r)
  w <- w[w != 0]
  trim(successiveViews(as(x, "Rle"), w))
}

roll_view <- function(x, width, step) {
  rng <- IRanges(start = seq.int(1, set_ln(x), by = step), width = width)
  Views(as(x, "Rle"), start = rng)
}

stretch_view <- function(x, start, step) {
  rng <- IRanges(start = start, width =  seq.int(start, set_ln(x), by = step))
  trim(Views(as(x, "Rle"), start = rng))
}


mapper <- function(x, .fun) {
  .fun <- rlang::as_function(.fun)
  
  res <- try(.fun(x), silent = TRUE)
  
  if (is(res, "try-error")) {
    res <- viewApply(x, function(.) .fun(.))
  }
  return(res)
}