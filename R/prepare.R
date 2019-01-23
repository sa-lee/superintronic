
#' Prepare GTF/GFF file for coverage analysis by
#' splitting into intron/exon features
#' 
#' @param .data a GRanges object obtained from `rtracklayer::import()` or
#' `plyranges::read_gff()`.
#' @param ...
prepare_annotation <- function(.data) {
  stopifnot(is(.data, "GRanges"))
  .data <- mutate(.data, gene_id = as(gene_id, "Rle"))
  
  pc_genes <- .data %>% 
    filter(type == "gene")  %>% 
    select(gene_id, gene_name, gene_type) %>% 
    arrange(gene_id)
  
  # prepare list columns for simple intronic/exonic features 
  pc_exons <- .data %>% 
    filter(type == "exon") %>% 
    select(gene_id) %>% 
    arrange(gene_id) %>% 
    S4Vectors::split(., .$gene_id) %>% 
    reduce()

  pc_introns <- pc_genes %setdiff% pc_exons

  # add columns
  md <- DataFrame(n_self_olaps = count_overlaps(pc_genes, pc_genes),
                  simple_exonic = pc_exons,
                  simple_intronic = pc_introns)
  mcols(pc_genes) <- cbind(mcols(pc_genes), md)
  
  
  pc_genes
  
}
