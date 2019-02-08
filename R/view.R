#' Visualisation Helpers
#' 
#' @param cvg a GRanges object obtained from `merge_coverage()`
#' @param target a GRanges object from `prepare_annotation()` containing a gene_id column
#' 
#' @importFrom ggplot2 scale_color_brewer aes facet_grid labs theme geom_segment geom_rect scale_x_reverse element_blank
#' @export
view_coverage_within_gene <- function(cvg, target, facets = ggplot2::vars()) {
  
  
  
  annotation <- bind_ranges(exon = unlist(target$simple_exonic, use.names = FALSE),
                            intron = unlist(target$simple_intronic, use.names = FALSE),
                            .id = "feature")
  
  annotation_tracks <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = annotation %>% 
                filter(feature == "exon") %>% 
                as.data.frame(),
                aes(xmin = start, xmax = end, ymin = 0, ymax = 1)
              ) +
    ggplot2::geom_segment(data = annotation %>%
                            filter(feature == "intron") %>% 
                            as.data.frame(),
                          aes(x = start, xend = end, y = 0.5, yend = 0.5)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
  

  
  coverage_view <- cvg %>% 
    plyranges::filter(gene_id == target$gene_id) %>% 
    plyranges::mutate(strand = feature_strand) 
  
  if (length(facets) > 0L) {
    coverage_view <- coverage_view %>% 
      group_by(!!!facets)
  } 
  
  # disjoin coverage scores
  coverage_view <- coverage_view %>%
    disjoin_ranges_directed(
      score = mean(log2(score) + 1)
    ) %>% 
    mutate(midpoint = IRanges::mid(.))
  
 
  
 
  cvg_hist <- ggplot2::ggplot(as.data.frame(coverage_view),
                     ggplot2::aes(x = midpoint, y= score)) + 
    ggplot2::geom_area() +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    labs(title = paste("Coverage profile:", target$gene_name))
  
  if (S4Vectors::runValue(BiocGenerics::strand(coverage_view)) == "-") {
    annotation_tracks <- annotation_tracks + ggplot2::scale_x_reverse()
    cvg_hist <- cvg_hist + ggplot2::scale_x_reverse()
  }
  
  if (length(facets) > 0L) {
    cvg_hist <- cvg_hist + facet_wrap(facets, ncol = 1L)
  }
  
  patchwork::wrap_plots(cvg_hist, 
                        annotation_tracks, 
                        nrow = 2L, 
                        heights = c(4,1))
  
}

#' Scatter exon-intron signal 
#' 
#' @param data a DataFrame object obtained from `spread_coverage()`
#' @param geom produce a 'point' scatter or a 'hex' scatter or c('hex', 'point')
#' for both
#' @param x,y bare variable names for the x,y plot aesthetics 
#' @param facets a set of variables quoted by `ggplot2::vars()` that is passed 
#' down to `facet_wrap()`
#'
#'    
#' 
#' @importFrom rlang enquo 
#' @importFrom ggplot2 ggplot geom_point geom_hex facet_wrap
#' @export
#' 
view_exin <- function(data, x, y, facets = ggplot2::vars(), geom = "point") {
  
  x <- rlang::enquo(x)
  y <- rlang::enquo(y)
  geom <- match.arg(geom, c("point", "hex"))
  
 
  p <- ggplot2::ggplot(as.data.frame(data), ggplot2::aes(!!x, !!y))
  
  
  if (geom == "point") {
     p <- p + ggplot2::geom_point() 
  }
  
  if (geom == "hex") {
    p <- p  + ggplot2::geom_hex()
  }
  
  if (length(facets) > 0) {
    p <- p + ggplot2::facet_wrap(facets)
  }
  
  return(p)
  
} 