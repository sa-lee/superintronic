#' Intron signal finder
#'
#' We make the function generic with dual dispatch,
#' the main signature can be BamFile or character + TxDb object which
#' we can use to represent a GFF file (and we can make a GTFFile from)
#' it can basically all be pushed down to that combo, so we really only
#' need one engine to do the work
#'
#'
#' @param reads an object representing a set of alignments such as an
#' [Rsamtools::BamFile]
#' @param annoationDb an object representing a database of reference annotations
#' such as GTFFile or a [TxDb](GenomicFeatures::`TxDb-class`) object.
#' @param ... additional arguments passed
#'
#' @details
#'
#' BAM processing
#' 1. For a given BAM file compute the coverage as an RleList
#' 2. Given the output of GTF above, group by gene id, and then reduce over
#' the ranges and aggregate the coverage over each group. (could do an overlap join here instead?)
#' 3. Over each gene compute the log2 trimmed mean for exon and intron counts
#' and flag whether the gene is expressed (>= 3 reads in both exon and intron)
#' 4. Compute the run length of intron ranges exceeding a nominal threshold
#'
#'
#' @return a GRanges object with the following metadata columns
#' * GeneId
#' * IsExpressed: a flag for whether gene is expressed or not
#' * ExonLength
#' * IntronLength
#' * MeanExonCoverage
#' * MeanIntronCoverage
#'
#'
#' @export
#' @name findIR
setGeneric("findIR", function(reads, annotationDb, ...) {
  standardGeneric("findIR")
})


#' with exons, is gene expressed?
#' using intron coverage, find genes with low average coverage out of those
#' are there any that have spikes?
#' runs above some threshold then finding max peak

#' @export
#' @name findIR
setMethod("findIR", signature = c("BamFile", "TxDb"),
          function(reads, annotationDb, threshold) {
            cvg <- plyranges::compute_coverage(reads)
            db <- processAnnoDb(annotationDb)
            join_overlap_inner(cvg, db)
          })

