
#' Internal generic to process a GFF/GTF file
#'
#'
#' GTF processing:
#'
#' 1. Import the GTF file and extract genes and exon features. Retain any relevant annotation
#' information.
#' 2. Mark genes if they self-overlap.
#' 3. Split each gene into exon/intron coordinates. An intron is any region within
#' a gene not annotated as an exon.
#' 4. Report total exon and intron lengths in each gene, where exons are the
#' union set of all transcripts from the same gene.
#'
setGeneric("processAnnoDb", function(db) standardGeneric("processAnnoDb"))

setMethod("processAnnoDb", "TxDb",
           function(db) {
             # extract all exonic and intronic parts
             # note these could be linked to multiple genes
             # i.e. a many to many relationship
             # we don't actually need to do any summaries at this point
             exonic_parts <- GenomicFeatures::exonicParts(db)
             # reduced from exonic parts, disjoin, linked to a single gene, no transcript information
             intronic_parts <- GenomicFeatures::intronicParts(db)
             all <- plyranges::bind_ranges(exonic_parts, intronic_parts)
             all <- plyranges::unnest(all, gene_id)
             all <- plyranges::mutate(all,
                                      gene_id = S4Vectors::Rle(gene_id),
                                      is_exonic = S4Vectors::Rle(all(is.na(exon_id))))

             ans <- plyranges::group_by(all, gene_id, is_exonic)
             ans <- plyranges::reduce_ranges_directed(ans)
             return(ans)
           })


setMethod("processAnnoDb", "GFF3File",
          function(db) {
            db <- makeTxDbFromGFF(db)
            callNextMethod()
          })

setMethod("processAnnoDb", "GTFFile",
          function(db) {
            db <- makeTxDbFromGFF(db)
            callNextMethod()
          })



rangeCoverage <- function(gr, cvg) {
  gr <- IRanges::reduce(gr)
  ans <- unlist(cvg[gr], use.names=FALSE)
  if (S4Vectors::runValue(BiocGenerics::strand(gr))[[1L]] == "-")
    ans <- rev(ans)
  ans
}
