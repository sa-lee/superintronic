#' View coverage as a stacked area chart
#' 
#' @param cvg a GRanges object obtained from `merge_coverage()`
#' @param target a GRanges object from `prepare_annotation()` containing a gene_id column
#' @param hline an optional data.frame giving marker points
#' @param facets a set of variables quoted by `ggplot2::vars()` that is passed 
#' down to `facet_wrap()`
#'
#' 
#' @importFrom ggplot2 ggplot aes facet_wrap labs theme geom_segment geom_rect element_blank geom_hline
#' @importFrom patchwork wrap_plots
#' @export
view_coverage <- function(cvg, target, hline = NULL, facets = ggplot2::vars()) {
  
  annotation <- bind_ranges(
    exon = unlist(target$simple_exonic, use.names = FALSE),
    intron = unlist(target$simple_intronic, use.names = FALSE),
    .id = "feature"
  )
  
  annotation_tracks <- view_segments(annotation, color = feature)
  
  coverage_view <- plyranges::filter(cvg, gene_id == target$gene_id)
  coverage_view <- plyranges::mutate(coverage_view, strand == feature_strand)

  if (length(facets) > 0L) {
    coverage_view <- plyranges::group_by(coverage_view, !!!facets)
  } 
  
  # disjoin coverage scores
  coverage_view <- plyranges::disjoin_ranges_directed(
    coverage_view,
    score = mean(log2(score + 1)),
    feature = unlist(unique(feature_type))
  ) 
  
  plot_tbl <- as.data.frame(coverage_view)
 
  cvg_hist <- ggplot(plot_tbl, 
                     aes(x = start, xend = end, y = 0, yend = score)) +
    geom_segment(aes(colour = feature)) +
    rescale_by_width(coverage_view)
  
  if (!is.null(hline)) {
    cvg_hist <- cvg_hist + geom_hline(data = hline, aes(yintercept = log2(E +1)))
  }
  
  if (length(facets) > 0L) {
    cvg_hist <- cvg_hist + facet_wrap(facets, ncol = 1L)
  }
  
  # cosmetic fixes 
  cvg_hist <- cvg_hist +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    scale_colour_brewer(palette = "Dark2") +
    guides(colour = FALSE) +
    labs(title = paste("Coverage profile:", target$gene_name),
         y  =  "coverage score (log2)") 
  
  patchwork::wrap_plots(cvg_hist, 
                        annotation_tracks, 
                        nrow = 2L, 
                        heights = c(4,1))
}

