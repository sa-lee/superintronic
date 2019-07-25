#' @importFrom BiocGenerics width strand
flatten <- function(x, type = "exon") {
  type <- match.arg(type, choices =c("exon", "intron"))
  n_features_in_genes <- lengths(x)
  gene_id <- S4Vectors::Rle(names(x), n_features_in_genes)
  feature_type <- S4Vectors::Rle(
    factor(type,levels = c("exon", "intron")), 
    sum(n_features_in_genes)
  )
  feature_rank <- unlist(lapply(n_features_in_genes, seq_len), 
                         use.names = FALSE)
  feature_length <- unlist(width(x), use.names = FALSE)
  feature_strand <- unlist(strand(x), use.names = FALSE)
  S4Vectors::DataFrame(
    gene_id, feature_type, feature_rank, feature_strand, feature_length
  )
}

#' Expand the exonic/intronic parts of a gene into a long GRanges object
#' 
#' @param x a GRanges object from `collect_parts()`
#' @details This function restricts the expands/unnests, the ranges from
#' `collect_parts` the result is a GRanges object with the following columns: 
#' 
#' * `gene_id`, the gene_id from `features` corresponding to an exon/intron.
#' * `feature_type`, whether the range corresponds to an exon or intro feature.
#' * `feature_rank`, the rank of the exon/intron feature within a gene.
#' * `feature_length`, the width of the exon/intron
#' * `feature_strand`, the strand of the exon/intron
#' 
#' @return a GRanges object
#' 
#' @importFrom S4Vectors mcols<- mcols
#' @export
unnest_parts <- function(x) {
  
  has_exin_cols <- sum(names(mcols(x)) %in% c("exonic_parts", "intronic_parts")) == 2
  if (!has_exin_cols) stop("Cannot find ")
  
  exonic <- mcols(x)[["exonic_parts"]]
  features_exons <- unlist(exonic, use.names = FALSE) 
  mcols(features_exons) <- flatten(exonic)
  
  intronic <- mcols(x)[["intronic_parts"]]
  features_introns <- unlist(intronic, use.names = FALSE)
  mcols(features_introns) <- flatten(intronic, "intron")
  
  BiocGenerics::sort(plyranges::bind_ranges(features_exons, features_introns)) 
    
}


#' Expand exonic/intronic parts and intersect with ranges
#' 
#' @param x a GRanges object
#' @param parts an annotation GRanges from `collect_parts()`
#' 
#' @details This is a wrapper to `unnest_parts()` and 
#' `plyranges::join_overlap_intersect()`. It restricts the GRanges object `x` 
#' to the ranges that intersect the `exonic_parts` and `intronic_parts` columns
#' from `parts`, the result is a GRanges object with additional columns: 
#' 
#' * `gene_id`, the gene_id from `features` corresponding to an exon/intron.
#' * `feature_type`, whether the range corresponds to an exon or intro feature.
#' * `feature_rank`, the rank of the exon/intron feature within a gene.
#' * `feature_length`, the width of the exon/intron
#' * `feature_strand`, the strand of the exon/intron
#' 
#' @return a GRanges object
#' 
#' @export
join_parts <- function(x, parts) {
  f <- unnest_parts(parts)
  plyranges::join_overlap_intersect(x, f)
}
