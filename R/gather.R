
#' Generate long-form coverage as coverage
#' @param bams a character vector containing a list of paths to indexed BAM files
#' @param genome_info a GRanges object containing reference annotation
#' @param BPPARAM a BiocParallel object (default = `bpparam()`)
#' 
#' @importFrom BiocParallel bpparam bplapply
#' @importFrom Rsamtools BamFile BamFileList
#' @export
gather_coverage <- function(bams, genome_info = NULL, BPPARAM = BiocParallel::bpparam()) {
  
  bfl <- Rsamtools::BamFileList(bams)

  cvg  <- BiocParallel::bplapply(
    bfl,
    FUN = compute_coverage,
    BPPARAM = BPPARAM
    ) 
  
  cvg <- GenomicRanges::GRangesList(cvg)
  
  cvg <- unlist(cvg, use.names = FALSE)
  
  if (!is.null(genome_info)) {
    return(cvg %>% filter_by_overlaps(genome_info))
  }
  cvg
}