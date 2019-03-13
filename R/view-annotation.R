#' An opinionated way of plotting intron/exon features
#' 
#' 
#' @param data a GRanges object
#' @param color a column in data that represents a feature 
#' 
#' 
#' @importFrom ggplot2 geom_segment aes scale_y_continuous scale_colour_brewer theme_bw theme guides scale_x_continuous expand_scale
#' @return a ggplot object
#' @export
view_segments <- function(data, color) {
  
  color <- enquo(color)
  
  seg1 <- as.data.frame(filter(data, !!color == "exon"))
  seg2 <- as.data.frame(filter(data, !!color == "intron"))
  
  annotation_tracks <-  ggplot()  +
    geom_segment(data = seg1,
                 aes(x = start, xend = end, y = 0.5, yend = 0.5, 
                     colour = !!color),
                 size = 10) +
    geom_segment(data = seg2,
                 aes(x = start, xend = end, y = 0.5, yend = 0.5, 
                     colour = !!color),
                 lineend = "butt",
                 linejoin = "round",
                 size = 2) +
    scale_y_continuous(expand = c(0,0))+
    scale_colour_brewer(palette = "Dark2") +
    guides(colour = FALSE) +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    )
  
  annotation_tracks + rescale_by_width(data)
  
}

rescale_by_width <- function(data, expand = ggplot2::expand_scale(mult = 0.02)){
  
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