#' Visualisation Helpers
#' 
#' @param x a GRanges object obtained from `merge_coverage()`
#' @param gene a GRanges object from `prepare_annotation()` containing a gene_id column
#' 
#' @importFrom ggbio autoplot   
#' @importFrom ggplot2 scale_color_brewer aes facet_grid labs theme
#' @export
plot_cvg_over_gene <- function(x, gene) {
  gene_id <- gene$gene_id
  
  t1 <- x %>% 
    plyranges::filter(gene_id == !!gene_id) %>% 
    ggbio::autoplot(., 
                    aes(
                      y = log2(score+1), 
                      group = Replicate, 
                      colour = feature_type
                    ), 
                    geom = "line") + 
    facet_grid(Kit ~ .) +
    scale_color_brewer(palette = "Dark2") +
    theme(axis.text.x = element_text(angle = 60, vjust = 0.5)) +
    labs(title = mcols(gene)[["gene_name"]])
  
  t1
}