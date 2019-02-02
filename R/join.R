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
#' * `feature_prop_covered`, the proportion of bases that overlap a feature
#' * `feature_strand`, the strand of the exon/intron
#' 
#' @return a GRanges object
#' 
#' @export
merge_coverage <- function(cvg, features) {
  f <- combine_exin(features)
  join_overlap_intersect(cvg, f) %>% 
    mutate(
      feature_prop_covered = width / feature_length
    )
}


merge_design <- function(cvg, design, to = NULL) {
  # this is an inner join
  # for some reason though merge is really slow
  # so have gone for a kludgy match approach
  stopifnot(any(names(mcols(cvg)) %in% "source"))
  if (is.null(to)) {
    to <- "source"
  }
  stopifnot(any(names(design) %in% to))
  
  diff <- BiocGenerics::setdiff(
    S4Vectors::runValue(mcols(cvg)[["source"]]),
    design[[to]]
  )
  if (length(diff) > 0) {
    cvg <- plyranges::filter(cvg, !(source %in% !!diff))
  }
  
  cvg <- cvg[BiocGenerics::order(mcols(cvg)[["source"]]), ]

  nr <- seq_len(S4Vectors::nrun(mcols(cvg)[["source"]]))
  rl <- S4Vectors::runLength(mcols(cvg)[["source"]])

  design <- design[BiocGenerics::order(design[[to]]), -c(to)]
  design <- design[BiocGenerics::rep.int(nr, rl), ]
  
  mcols(cvg) <- cbind(
    mcols(cvg),
    design
  )
  cvg
}


common_mcols <- function(x, y, by = NULL) {
  x_names <- names(mcols(x))
  y_names <- names(mcols(y))
  if (is.null(by)) {
    common_mcols <- intersect(x_names,  y_names)
    if (length(common_mcols) == 0) {
      stop("No common columns between x & y", call. = FALSE)
    }
    return(common_mcols)
  } else {
    named_by <- names(by)
    if (length(named_by) > 0) {
      stopifnot(named_by %in% x_names || by %in% y_names)
      by
      
    } else {
      stopifnot(by %in% x_names || by %in% y_names)
      by
    }
  }
}