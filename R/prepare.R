#' GTF processing:
#'
#' 1. Import the GTF file and extract genes and exon features. Retain any relevant annotation
#' information.
#' 2. Mark genes if they self-overlap.
#' 3. Split each gene into exon/intron coordinates. An intron is any region within
#' a gene not annotated as an exon.
#' 4. Report total exon and intron lengths in each gene, where exons are the
#' union set of all transcripts from the same gene.


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
    IRanges::reduce()

  pc_introns <- pc_genes %setdiff% pc_exons

  # add columns
  md <- S4Vectors::DataFrame(n_self_olaps = count_overlaps(pc_genes, pc_genes),
                  simple_exonic = pc_exons,
                  simple_intronic = pc_introns)
  S4Vectors::mcols(pc_genes) <- cbind(S4Vectors::mcols(pc_genes), md)
  
  
  pc_genes
  
}
