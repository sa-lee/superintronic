#' Compute long-form coverage as a GRanges object
#' 
#' @param spec a `data.frame` containing an experimental design
#' @param source a column in the data identifying the name of BAM files
#' @param .target an optional GRanges object (default = NULL) to compute coverage
#' over a region (requires the BAM file to be indexed).
#' @param .genome_info a GRanges object containing reference annotation (default = NULL)
#' @param .drop_empty Filter ranges if they have zero coverage over an entire contig/chromosome (default TRUE)
#' @param .parallel a BiocParallel object (default = `BiocParallel::bpparam()`)
#' 
#' @details This function computes coverage as a GRanges object, with
#' a `source` column containing the name of the input BAM file(s), and a
#' `score` column containing the coverage score. 
#' 
#' The `.genome_info` argument takes a GRanges object containing 
#' contig/chromosome information. If supplied the resulting GRanges will be 
#' properly annotated with a `seqinfo` slot. This is important for ensure the
#' integrity of any downstream overlap operations. 
#' 
#' @return a GRanges object
#' 
#' @importFrom BiocParallel bpparam bplapply
#' @importFrom Rsamtools BamFile BamFileList
#' @importClassesFrom Rsamtools BamFileList BamFile
#' @importFrom GenomeInfoDb seqinfo seqinfo<- keepSeqlevels
#' @export
methods::setGeneric("compute_coverage_long",
                    signature = c("spec", "source"),
                    function(spec, source, ...) {
                      methods::standardGeneric("compute_coverage_long")
                    })


#'@export 
methods::setMethod("compute_coverage_long", 
                   signature = c("character", "missing"),
                   function(spec, source, .target = NULL, .genome_info = NULL, .drop_empty = TRUE, .parallel = BiocParallel::bpparam()) {
                     bfl <- BamFileList(spec)
                     .bamlist_coverage(bfl, .target, .genome_info, .parallel)
                   })


methods::setMethod("compute_coverage_long",
                   signature = c("DataFrame", "character"),
                   function(spec, source, .target = NULL, .genome_info = NULL, .drop_empty = TRUE, .parallel = BiocParallel::bpparam()) {
                     bfl <- BamFileList(spec[[source]])
                     spec[[source]] <- as(spec[[source]], "Rle")
                     S4Vectors::mcols(bfl) <- spec
                     ans <- .bamlist_coverage(bfl, 
                                              .target, 
                                              .genome_info, 
                                              .drop_empty, 
                                              .parallel)
                     S4Vectors::mcols(ans)[[source]] <-basename(S4Vectors::mcols(ans)[[source]])
                     ans
                   })

methods::setMethod("compute_coverage_long",
                   signature = c("data.frame", "character"), 
                   function(spec, source, .target = NULL, .genome_info = NULL, .drop_empty = TRUE, .parallel = BiocParallel::bpparam()) {
                     spec <- as(spec, "DataFrame")
                     compute_coverage_long(spec, 
                                           source, 
                                           .target, 
                                           .genome_info, 
                                           .drop_empty, 
                                           .parallel)
                   }
)
.bamlist_coverage <- function(bfl, .target, .genome_info, .drop_empty, .parallel) {
  
  sb <- Rsamtools::ScanBamParam()
  
  if (!is.null(.target)) {
    sb <- Rsamtools::ScanBamParam(which = .target)
  }
  
  cvg  <- BiocParallel::bplapply(
    bfl,
    FUN = function(x) {
        compute_coverage(x, param = sb)
      },
    BPPARAM = .parallel
    ) 
  
  cvg <- GenomicRanges::GRangesList(cvg)
  
  if (!is.null(S4Vectors::mcols(bfl))) {
    md  <- S4Vectors::mcols(bfl)
  } else {
    md <- S4Vectors::DataFrame(source = S4Vectors::Rle(names(cvg)))
  }
  
  if (!is.null(.genome_info)) {
    cvg <- GenomeInfoDb::keepSeqlevels(cvg, 
                                       GenomeInfoDb::seqnames(.genome_info), 
                                       "tidy")
    seqinfo(cvg) <- seqinfo(.genome_info)
  }
  
  if (.drop_empty) {
    .drop_rng <- as(seqinfo(cvg), "GRanges")
    cvg <- S4Vectors::endoapply(cvg,
                  function(x) {
                    IRanges::subsetByOverlaps(x, .drop_rng, 
                                              type = "equal", invert = TRUE)
                  }
    )
  }
  
  inx <- S4Vectors::rep.int(seq_along(cvg), S4Vectors::elementNROWS(cvg))
  cvg <- unlist(cvg, use.names = FALSE)
  md <- md[inx,, drop = FALSE]
  mcols(cvg) <- cbind(md, mcols(cvg))
  cvg
}

#' Merge a design matrix onto a GRanges object
#' 
#' @param ranges a GRanges object 
#' @param design a DataFrame object with 
#' @param by the key column either a character vector of the common columns
#' between `ranges` and `design` or a named character vector (default = NULL).
#' 
#' @details 
#' To summarise ranges features over variables in the experimental `design`,
#' we can join the `design`` matrix to a `ranges` object. 
#' This function computes an inner join on the common key between the two 
#' objects and sorts by the key.
#' 
#' @seealso `dplyr::inner_join()`
#'
#' @export
join_design <- function(ranges, design, by = NULL) {
  # for some reason though base::merge is really slow
  # so have gone for a kludgy indexing approach
  stopifnot(any(names(mcols(ranges)) %in% "source"))
  stopifnot(is(design, "data.frame") | is(design, "DataFrame"))
  
  by <- common_mcols(ranges, design, by = by)
  if (length(by) > 1) stop("Multiple keys not supported.")
  
  key_x <- by
  key_y <- by
  is_named_by <- length(names(by)) == 1
  
  col_x <- as(mcols(ranges)[[key_x]], "Rle")
  col_y <- design[[key_y]]
  
  if (is_named_by) key_x <- names(by)
  diff <- BiocGenerics::setdiff(
    S4Vectors::runValue(col_x),
    col_y
  )
  
  if (length(diff) > 0) {
    ranges <- plyranges::filter(ranges, !(!!key_x %in% !!diff))
  }
  
  ranges <- ranges[BiocGenerics::order(col_x), ]
  
  nr <- seq_len(S4Vectors::nrun(col_x))
  rl <- S4Vectors::runLength(col_x)
  
  design <- design[BiocGenerics::order(col_y), 
                   !key_y,
                   drop = FALSE]
  
  design <- design[BiocGenerics::rep.int(nr, rl), , drop = FALSE]
  
  
  mcols(ranges) <- BiocGenerics::cbind(
    mcols(ranges),
    design
  )
  ranges
}


# Helper function for keys
common_mcols <- function(x, y, by = NULL) {
  x_names <- names(mcols(x))
  y_names <- names(y)
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