#' @importFrom BiocGenerics width strand
flatten <- function(x, type = "exon") {
  type <- match.arg(type, choices =c("exon", "intron"))
  n_features_in_genes <- BiocGenerics::lengths(x)
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

#' @importFrom S4Vectors mcols<- mcols
combine_exin <- function(x) {
  exonic <- mcols(x)[["simple_exonic"]]
  features_exons <- unlist(exonic, use.names = FALSE) 
  mcols(features_exons) <- flatten(exonic)
  
  intronic <- mcols(x)[["simple_intronic"]]
  features_introns <- unlist(intronic, use.names = FALSE)
  mcols(features_introns) <- flatten(intronic, "intron")
  
  plyranges::bind_ranges(features_exons, features_introns) %>% 
    BiocGenerics::sort()
}


#' Restrict coverage over features
#' 
#' 
#' @param cvg a GRanges object from `gather_coverage()`
#' @param features a GRanges object from `prepare_annotation()`
#' 
#' @details This function restricts the coverage GRanges, `cvg` to the ranges
#' that intersect the `simple_exonic` and `simple_intronic` ranges, the result
#' is a GRanges object with additional columns: 
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
merge_coverage <- function(cvg, features) {
  f <- combine_exin(features)
  join_overlap_intersect(cvg, f)
}