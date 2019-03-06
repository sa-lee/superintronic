#' Visualisation Helpers
#' 
#' @param cvg a GRanges object obtained from `merge_coverage()`
#' @param target a GRanges object from `prepare_annotation()` containing a gene_id column
#' @param facets a set of variables quoted by `ggplot2::vars()` that is passed 
#' down to `facet_wrap()`
#'
#' 
#' @importFrom ggplot2 scale_color_brewer aes facet_grid labs theme geom_segment geom_rect scale_x_reverse element_blank
#' @export
view_coverage_within_gene <- function(cvg, target, hline = NULL, facets = ggplot2::vars()) {
  
  
  
  annotation <- bind_ranges(exon = unlist(target$simple_exonic, use.names = FALSE),
                            intron = unlist(target$simple_intronic, use.names = FALSE),
                            .id = "feature")
  
  annotation_tracks <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = annotation %>% 
                filter(feature == "exon") %>% 
                as.data.frame(),
                aes(x = start, xend = end, 
                    y = 0.5, yend = 0.5, colour = feature),
                size = 10
              ) +
    ggplot2::geom_segment(data = annotation %>%
                            filter(feature == "intron") %>% 
                            as.data.frame(),
                          aes(x = start, xend = end, 
                              y = 0.5, yend = 0.5, colour = feature),
                          lineend = "butt",
                          linejoin = "round",
                          size = 2) +
    ggplot2::scale_y_continuous(expand = c(0,0))+
    ggplot2::scale_colour_brewer(palette = "Dark2") +
    ggplot2::guides(colour = FALSE) +
    ggplot2::theme_bw() +
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
      score = mean(log2(score + 1)),
      feature_type = unlist(unique(feature_type))
    ) 
  
  
  plot_tbl <- as.data.frame(coverage_view)
 
  cvg_hist <- ggplot2::ggplot(plot_tbl, 
                              mapping = ggplot2::aes(x = start, xend = end, y = 0, yend = score)) + 
    ggplot2::geom_segment(aes(colour = feature_type)) 
  
  if (!is.null(hline)) {
    cvg_hist <- cvg_hist + ggplot2::geom_hline(data = hline, aes(yintercept = log2(E +1)))
  }
    
  expand <- ggplot2::expand_scale(mult = 0.02)
  
  if (S4Vectors::runValue(BiocGenerics::strand(coverage_view)) == "-") {
    annotation_tracks <- annotation_tracks + 
      ggplot2::scale_x_reverse(expand = expand, labels = ggbio:::trans_seq_format("kb"))
    cvg_hist <- cvg_hist + 
      ggplot2::scale_x_reverse(expand = expand, labels = ggbio:::trans_seq_format("kb"))
  } else {
    annotation_tracks <- annotation_tracks + 
      ggplot2::scale_x_continuous(expand = expand, labels = ggbio:::trans_seq_format("kb"))
    cvg_hist <- cvg_hist + 
      ggplot2::scale_x_continuous(expand = expand, labels = ggbio:::trans_seq_format("kb"))
  }
  
  if (length(facets) > 0L) {
    cvg_hist <- cvg_hist + facet_wrap(facets, ncol = 1L)
  }
  
  # cosmetic fixes 
  cvg_hist <- cvg_hist +
    ggplot2::theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    ggplot2::scale_colour_brewer(palette = "Dark2") +
    ggplot2::guides(colour = FALSE) +
    labs(title = paste("Coverage profile:", target$gene_name),
         y  =  "coverage score (log2)") 
  
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