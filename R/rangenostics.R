#' Rangewise diagnostics
#'
#' @param x a numeric vector or Rle object
#' @param .size the integer size of the window
#'  
#' @details These are inspired by the concept of cognostics by John Tukey, and related work in the
#' area by Lee Wilkinson (graph theoretic scagnostics) and in time series analysis by Rob Hyndman.
#' 
#' The 
#' 

lumpy <- function(x, .size) {
  var(tile_map(x,  .size, function(.) mean(., na.rm = TRUE)), na.rm = TRUE)
}

stable <- function(x, .size) {
  var(tile_map(x,  .size, function(.) var(., na.rm = TRUE)), na.rm = TRUE)
}

sum_above <- function(x, .threshold) {
  sum(x > .threshold)
}

flat_spots <- function(x, .threshold = 0, .tol = .Machine$double.eps^0.5) {
  sum(abs(x - .threshold) < .tol)
}

count_grams <- function(x, .size) {
  bins <- tile_map(x, .size, function(.) sum(., na.rm = TRUE))
  c(count_min = min(bins),
    count_max = max(bins),
    count_q1 = unname(quantile(bins, 0.25)),
    count_q3 = unname(quantile(bins, 0.75)),
    count_mean = mean(bins),
    count_median = median(bins),
    count_var = var(bins)
  )
  
}

gammy <- function(x) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package: mgcv required to run gammy diagnostics. Please install it.")
  } 
}