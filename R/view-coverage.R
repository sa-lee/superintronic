#' View coverage as an area chart oriented from 5' to 3'
#' 
#' @param data a GRanges object 
#' @param score an expression or column in data that represents the coverage score
#' @param colour a column in the data to highlight ranges by  
#' @param facets a set of metadata or range variables in `data` quoted by `rng_vars()` that is passed 
#' down to `facet_wrap()`
#'
#' @importFrom ggplot2 ggplot aes facet_wrap labs theme geom_segment geom_rect element_blank geom_hline
#' @importFrom patchwork wrap_plots
#' @export
view_coverage <- function(data, score, colour = NULL, facets = rng_vars()) {
  
  score <- rlang::enquo(score) 
  colour <- rlang::enquo(colour) 

  no_col <- rlang::quo_is_null(colour)
  if (length(facets) > 0L) {
    data <- plyranges::group_by(data, !!!facets)
  } 
  
  # disjoin coverage scores
  if (!no_col) {
    coverage_view <- plyranges::disjoin_ranges_directed(
      data,
      score = mean(!!score),
      feature = unlist(unique(!!colour))
    ) 
  } else {
    coverage_view <- plyranges::disjoin_ranges_directed(
      data,
      score = mean(!!score)
    ) 
  }

  
  plot_tbl <- as.data.frame(coverage_view)
  
  cvg_hist <- ggplot(plot_tbl, 
                     aes(x = start, xend = end, y = 0, yend = score)) +
    rescale_by_width(data)
  
  if (!no_col)  {
    cvg_hist <- cvg_hist + geom_segment(aes(colour = feature))
  } else {
    cvg_hist <- cvg_hist + geom_segment()
  }
    
  if (length(facets) > 0L) {
    cvg_hist <- cvg_hist + facet_wrap(facets, ncol = 1L)
  }
  
  # cosmetic fixes 
  cvg_hist <- cvg_hist +
    labs(y = rlang::quo_get_expr(score), x = NULL)
  
  
  return(cvg_hist)
}

