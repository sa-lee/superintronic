#' An opinionated way of plotting intron/exon features
#' 
#' 
#' @param data a GRanges object from `collect_parts()`
#' @param colour an optional expression or bare variable in `data` that represents a feature 
#' 
#' @importFrom ggplot2 geom_segment scale_x_reverse aes scale_y_continuous scale_colour_brewer theme_bw theme guides scale_x_continuous expansion
#' @return a ggplot object
#' @export
view_segments <- function(data, colour) {
  
  colour <- enquo(colour)
  
  seg1 <- as.data.frame(filter(data, !!colour == "exon"))
  seg2 <- as.data.frame(filter(data, !!colour == "intron"))
  
  annotation_tracks <-  ggplot()  +
    geom_segment(data = seg1,
                 aes(x = start, xend = end, y = 0.5, yend = 0.5, 
                     colour = !!colour),
                 size = 10) +
    geom_segment(data = seg2,
                 aes(x = start, xend = end, y = 0.5, yend = 0.5, 
                     colour = !!colour),
                 lineend = "butt",
                 linejoin = "round",
                 size = 2) +
    scale_y_continuous(expand = c(0,0))+
    scale_colour_brewer(palette = "Dark2") +
    guides(colour = FALSE) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    )
  
  annotation_tracks + rescale_by_width(data)
  
}


layer_segment <- function(data, y, ...) {
  stopifnot(is(data, "GRanges"))
  y <- enquo(y)
  plot_tbl <- as.data.frame(data)
  
  geom_segment(data = plot_tbl, aes(y = !!y, yend = !!y), ...)
  
}

rescale_by_width <- function(data, expand = ggplot2::expansion(mult = 0.02)){
  
  if (S4Vectors::runValue(BiocGenerics::strand(data)) == "-") {
    return(
      scale_x_reverse(
        expand = expand,
        labels = genome_labeller
      )
    )
  }
  
  scale_x_continuous(
    expand = expand,
    labels = genome_labeller
  )
  
}


genome_labeller <- function(breaks) {
  n_bases <- max(breaks, na.rm = TRUE)
  if (n_bases < 10000L) {
    return(paste(breaks, "bp"))
  } else if (n_bases >= 10000L & n_bases < 1e5L) {
    return(paste(breaks / 1000, "kb"))
  } else {
    return(paste(breaks / 1e6, "Mb"))
  }
}