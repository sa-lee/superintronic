#' A wrapper to compute coverage but only on 'standard' chromosomes
standardise_coverage <- function(x, genome_info) {
  plyranges::compute_coverage(x) %>%
    plyranges::filter_by_overlaps(genome_info)
}

standardise_reference <- function(x, genome_info) {
  x %>%
    plyranges::get_genome_info() %>% 
    plyranges::filter_by_overlaps(genome_info)
}

#' Generate long-form coverage as coverage
#' @param bams a character vector containing a list of paths to indexed BAM files
#' @param genome_info a GRanges object containing reference annotation
#' @param BPARAM a BiocParallel object (default = `bpparam()`)
#' 
#' @importFrom BiocParallel bpparam bplapply
#' @importFrom Rsamtools BamFile BamFileList
#' @export
gather_coverage <- function(bams, genome_info, BPARAM = BiocParallel::bpparam()) {
  
  bfl <- lapply(bams,
                function(.f) Rsamtools::BamFile(.f)
                ) %>% 
    Rsamtools::BamFileList()
  
  genome_info <- standardise_reference(bfl, genome_info)
  
  cvg  <- BiocParallel::bplapply(
    bfl,
    FUN = function(x) standardise_coverage(x, genome_info = genome_info),
    BPARAM = BPARAM
  ) %>%
    GenomicRanges::GRangesList()
  
  unlist(cvg, use.names = FALSE)
}