#' Collect and combine intronic/exonic parts by gene 
#' 
#' @param annotation an object representing an annotation such as an
#' `GenomicFeatures::TxDb`, a `GenomicRanges::GRanges` object obtained 
#' from `rtracklayer::import()` or `plyranges::read_gff()` or the name of
#' GFF/GTF file.
#' @param ... optional arguments to be passed to downstream methods
#' 
#' @details This function takes an annotation object and returns
#' a GRanges object with number of rows equal to the number of gene 
#' features contained. The resulting GRanges object has the following 
#' additional metadata columns:
#' 
#' * `exonic_parts`: a GRangesList column containing the exonic regions 
#' for each gene, obtained from `type == 'exon'` in GFF/GTF file. In the 
#' case of a TxDb object, this is obtained by `GenomicFeatures::exonBy`.
#' * `intronic_parts`: a GRangesList column containing the intronic regions
#' for each gene. This is simply defined as the parallel set difference betweeen
#' the gene ranges and the exonic ranges.
#' * `n_olaps`: the integer count of the number times a gene overlaps
#' any other gene in the annotation.
#' 
#' @examples 
#' gtf <- system.file("extdata", "Homo_sapiens.GRCh37.75_subset.gtf", package = "airway")
#' collect_parts(gtf)
#' 
#' @return a GRanges object
#' @importClassesFrom GenomicFeatures TxDb
#' @importClassesFrom rtracklayer GFFFile GTFFile
#'  
#' @export
#' @rdname collect_parts   
setGeneric("collect_parts", 
                    function(annotation, ...) {
                      standardGeneric("collect_parts")
                    }
)

#' @export
#' @rdname collect_parts  
setMethod("collect_parts", 
                   signature = "GRanges",
                   function(annotation, ...) {
                     .collect_parts_default(annotation)
                   })


#' @export
#' @rdname collect_parts  
setMethod("collect_parts", 
                   signature = "character",
                   function(annotation, ...) {
                     is_gff_or_gtf <- grepl(".(gff)$|(gtf)$", 
                                            tolower(basename(annotation)))
                     if (!is_gff_or_gtf) {
                       stop(
                         sprintf("File '%s' is not a GFF/GTF file.", annotation),
                         call. = FALSE
                       )
                     }
                     stopifnot(file.exists(annotation))
                     annotation <- rtracklayer::import(annotation, ...)
                     collect_parts(annotation)
                   })

#' @export
#' @rdname collect_parts  
setMethod("collect_parts", 
                   signature = "GFFFile",
                   function(annotation, ...) {
                     collect_parts(rtracklayer::import(annotation,...))
                     })


#' @export
#' @rdname collect_parts  
setMethod("collect_parts", 
                   signature = "TxDb", 
                   function(annotation) {
                     .collect_parts_txdb(annotation)
                   } )

.collect_parts_txdb <- function(annotation) {
  
  # returns single strand only genes
  pc_genes <- GenomicFeatures::genes(annotation)
  # returns all genes
  pc_exons <- GenomicFeatures::exonsBy(annotation, by = "gene")
  pc_exons <- pc_exons[names(pc_exons) %in% pc_genes$gene_id]
  pc_introns <- IRanges::psetdiff(pc_genes, pc_exons)
  
  .nest_parts(pc_genes, pc_exons, pc_introns) 
}

.collect_parts_default <- function(annotation) {
  # required inputs are gene_id and type
  stopifnot(sum(names(mcols(annotation)) %in% c("gene_id", "type")) == 2L)
  
  # cast gene_id as Rle
  annotation <- plyranges::mutate(annotation, 
                                  gene_id = methods::as(gene_id, "Rle"))
  
  # extract genes from annotation
  pc_genes <- plyranges::filter(annotation, type == "gene")   
  pc_genes <- plyranges::select(pc_genes, 
                                gene_id, gene_name, gene_biotype,  gene_source, 
                                type, source) 
  pc_genes <- plyranges::arrange(pc_genes, gene_id)
  
  # prepare list columns for simple intronic/exonic features 
  pc_exons <- plyranges::filter(annotation, type == "exon") 
  pc_exons <- plyranges::select(pc_exons, gene_id) 
  pc_exons <- plyranges::arrange(pc_exons, gene_id)
  
  pc_exons_grl <- S4Vectors::split(pc_exons, 
                                   S4Vectors::mcols(pc_exons)[["gene_id"]]) 
  pc_exons_grl <- IRanges::reduce(pc_exons_grl)
  
  pc_introns_grl <- IRanges::psetdiff(pc_genes, pc_exons_grl)
  
  .nest_parts(pc_genes, pc_exons_grl, pc_introns_grl)
}

.nest_parts <- function(genes, exons, introns) {
  # add nested columns containing intronic/exonic parts
  md <- S4Vectors::DataFrame(n_olaps = plyranges::count_overlaps(genes, genes),
                             exonic_parts = exons,
                             intronic_parts = introns)
  S4Vectors::mcols(genes) <- cbind(S4Vectors::mcols(genes), md)
  genes
}