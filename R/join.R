flatten <- function(x, type = "exon") {
  type <- match.arg(type, choices =c("exon", "intron"))
  n_features_in_genes <- lengths(x)
  S4Vectors::DataFrame(
    gene_id = S4Vectors::Rle(names(x), lengths = n_features_in_genes),
    feature_type = S4Vectors::Rle(factor(type, levels = c("exon", "intron")), 
                       lengths = sum(n_features_in_genes)),
    feature_rank = unlist(lapply(n_features_in_genes, seq_len), use.names = FALSE),
    feature_length = unlist(BiocGenerics::width(x), use.names = FALSE),
    feature_strand = unlist(BiocGenerics::strand(x), use.names = FALSE)
  )
}

#' @importFrom S4Vectors `mcols<-`
combine <- function(x) {
  exonic <- x[["simple_exonic"]]
  features_exons <- unlist(exonic, use.names = FALSE) 
  mcols(features_exons) <- flatten(exonic)
  
  intronic <- x[["simple_intronic"]]
  features_introns <- unlist(intronic, use.names = FALSE)
  mcols(features_introns) <- flatten(intronic, "intron")
  
  plyranges::bind_ranges(features_exons, features_introns) %>% 
    BiocGenerics::sort()
}


#' Coverage over features
#' 
#' @param cvg a GRanges object from `gather_coverage()`
#' @param features a GRanges object from `prepare_annotation()`

merge_coverage <- function(cvg, features) {
  f <- combine(features)
  join_overlap_intersect(cvg, f)
}