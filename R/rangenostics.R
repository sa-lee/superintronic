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
  c(lumpy = var(tile_map(x,  .size, function(.) mean(., na.rm = TRUE)), na.rm = TRUE))
}

stable <- function(x, .size) {
  c(stable = var(tile_map(x,  .size, function(.) var(., na.rm = TRUE)), na.rm = TRUE))
}

sum_above <- function(x, .threshold) {
  c(sum_above = sum(x > .threshold))
}

flat_spots <- function(x, .threshold = 0, .tol = .Machine$double.eps^0.5) {
  c(flat_spots = sum(abs(x - .threshold) < .tol))
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


max_mean_shift <- function(x, .size, .step) {
  views <- roll_view(x, .size, .step)
  mns <- mean(views)
  diffs <- abs(diff(mns))
  
  max_mean <- max(diffs, na.rm = TRUE)
  max_inx <- which.max(diffs) + 1L
  
  c(max_mean = max_mean, max_inx = max_inx)
  
}

max_var_shift <- function(x, .size, .step) {
  views <- roll_view(x, .size, .step)
  vars <- var(views, na.rm = TRUE)
  diffs <- abs(diff(vars))
  
  max_var <- max(diffs, na.rm = TRUE)
  max_inx <- which.max(diffs) + 1L
  
  c(max_var = max_var, max_inx = max_inx)
  
}

gammy <- function(x) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package: mgcv required to run gammy diagnostics. Please install it.")
  } 
}