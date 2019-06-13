#' Rangewise diagnostics
#'
#' @param x a numeric vector or Rle object
#' @param .size the integer size of the window
#'  
#' @details These are inspired by the concept of cognostics by John Tukey, and related work in the
#' area by Lee Wilkinson (graph theoretic scagnostics) and in time series analysis by Rob Hyndman.
#' 
#' The 
#' @rdname rangenostics
#' @export
lumpy <- function(x, .size) {
  var(tile_rle(x, .size, function(.) mean(., na.rm = TRUE)), na.rm = TRUE)
}

#' @rdname rangenostics
#' @export
bumpy <- function(x, .size) {
  var(tile_rle(x, .size, function(.) var(., na.rm = TRUE)), na.rm = TRUE)
}

#' @rdname rangenostics
#' @export
sum_above <- function(x, .threshold) {
  sum(x > .threshold)
}

#' @rdname rangenostics
#' @export
flat_spots <- function(x, .threshold = 0, .tol = .Machine$double.eps^0.5) {
  sum(abs(x - .threshold) < .tol)
}

#' @rdname rangenostics
#' @export
count_grams <- function(x, .size) {
  bins <- tile_rle(x, .size, function(.) sum(., na.rm = TRUE))
  list(count_min = min(bins),
    count_max = max(bins),
    count_q1 = unname(quantile(bins, 0.25)),
    count_q3 = unname(quantile(bins, 0.75)),
    count_mean = mean(bins),
    count_median = median(bins),
    count_var = var(bins)
  )
}

#' @rdname rangenostics
#' @export
max_mean_shift <- function(x, .size, .step) {
  mns <- roll_rle(x, .size, .step, mean, na.rm = TRUE)
  diffs <- abs(diff(mns, lag = .step))
  max_mean <- max(diffs, na.rm = TRUE)
  max_inx <- which.max(diffs) + 1L
  
  list(max_mean = max_mean, 
       max_inx = max_inx)
  
}

#' @rdname rangenostics
#' @export
max_var_shift <- function(x, .size, .step) {
  vars <- roll_rle(x, .size, .step, var, na.rm = TRUE)
  vars <- var(views, na.rm = TRUE)
  diffs <- abs(diff(vars, lag = .step))
  
  max_var <- max(diffs, na.rm = TRUE)
  max_inx <- which.max(diffs) + 1L
  
  c(max_var = max_var, max_inx = max_inx)
  
}

#' @rdname rangenostics
#' @export
gammy <- function(x) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package: mgcv required to run gammy diagnostics. Please install it.")
  } 
}