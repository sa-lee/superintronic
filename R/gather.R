
#' Compute long-form coverage
#' 
#' 
#' @param bams a character vector containing a paths to indexed BAM files
#' @param genome_info a GRanges object containing reference annotation (default = NULL)
#' @param BPPARAM a BiocParallel object (default = `BiocParallel::bpparam()`)
#' 
#' @details This function computes coverage as a GRanges object, with
#' a `source` column containing the name of the input BAM file, and a
#' `score` column containing the coverage score. The `genome_info` argument
#' takes a GRanges object containing contig/chromosome information. If
#' supplied coverage will be computed only over those contigs, and the
#' annotation information will be set on the resulting GRanges to that availble
#' in `genome_info`.
#' 
#' @return a GRanges object
#' 
#' @importFrom BiocParallel bpparam bplapply
#' @importFrom Rsamtools BamFile BamFileList
#' @importFrom GenomeInfoDb seqinfo seqinfo<- keepSeqlevels
#' @export
gather_coverage <- function(bams, genome_info = NULL, BPPARAM = BiocParallel::bpparam()) {
  
  bfl <- Rsamtools::BamFileList(bams)

  sb <- Rsamtools::ScanBamParam()
  
  if (!is.null(genome_info)) {
    sb <- Rsamtools::ScanBamParam(which = genome_info)
  }
  
  cvg  <- BiocParallel::bplapply(
    bfl,
    FUN = function(x) {
        compute_coverage(x, param = sb)
      },
    BPPARAM = BPPARAM
    ) 
  
  cvg <- GenomicRanges::GRangesList(cvg)
  
  if (!is.null(genome_info)) {
    cvg <- GenomeInfoDb::keepSeqlevels(cvg, 
                                       GenomeInfoDb::seqnames(genome_info), 
                                       "tidy")
    seqinfo(cvg) <- seqinfo(genome_info)
  }
  
  md <- S4Vectors::DataFrame(source =  Rle(names(cvg), lengths(cvg)))
  
  cvg <- unlist(cvg, use.names = FALSE)
  mcols(cvg) <- cbind(md, mcols(cvg))
  
  cvg
}