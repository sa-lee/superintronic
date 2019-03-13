#' Generate a scatter (xy) plot over a Ranges or DataFrame
#' 
#' @param data a Ranges or DataFrame object obtained
#' @param x,y bare variable names for the x,y plot aesthetics 
#' @param facets a set of variables quoted by `ggplot2::vars()` that is passed 
#' down to `facet_wrap()`
#' @param geom produce a 'point' scatter or a 'hex' scatter or c('hex', 'point')
#' for both
#'
#' @return a ggplot object
#' 
#' @importFrom rlang enquo 
#' @importFrom ggplot2 ggplot geom_point geom_hex facet_wrap
#' @export
view_xy <- function(data, x, y, facets = ggplot2::vars(), geom = "point") {
  
  x <- rlang::enquo(x)
  y <- rlang::enquo(y)
  geom <- match.arg(geom, c("point", "hex"))
  
  p <- ggplot2::ggplot(as.data.frame(data), ggplot2::aes(!!x, !!y))
  
  if (geom %in% "point") {
    p <- p + ggplot2::geom_point() 
  }
  
  if (geom %in% "hex") {
    p  <-  p  + ggplot2::geom_hex()
  }
  
  if (length(facets) > 0) {
    p <- p + ggplot2::facet_wrap(facets)
  }
  
  return(p)
  
} 