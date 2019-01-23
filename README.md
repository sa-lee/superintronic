
<!-- README.md is generated from README.Rmd. Please edit that file -->

# superintronic

Explore intron signal in high-throughput (RNA) sequencing data.

The interface is built on top of `plyranges` and is built around the
GRanges object.

The basic workflow centers around exploring coverage over exons and
intron regions.

## Setting up annotation GRanges

We assume you are starting from a GTF/GFF file for your given organism.
The reading of the GFF file is done external to this package, and can be
done via `rtracklayer::import` or `plyranges::read_gff`. The resulting
GRanges is passed to `prepare_annotation()`, then default filters are
applied using `filter_rules()`. Note that your own filters can be
supplied to `filter_rules()` which will override the default options.

The result is a GRanges object with number of rows equal to genes, and
columns containing the intronic and exonic features (as list columns).

``` r
library(superintronic)
library(magrittr)

gff <- "my.gff"

features <- gff %>% 
  plyranges::read_gff() %>% 
  prepare_annotation() %>% 
  filter_rules()
```

## Computing coverage over BAM(s)

Here we computer a GRanges containing coverage scores in parallel over a
set of BAM files. The (optional) parallel computation is done with
`BiocParallel`. The result is a long GRanges containing columns called
`score` and `file`.

``` r
bams <- dir(pattern = "*.bam")
cvg <- gather_coverage(bams, BPARAM = BiocParallel::bpparam(...))
```

## Summarising Overlaps

Then we can summarise the coverage scores that overlap our exon/introns
for each sample via `summarise_coverage(features, cvg)`, where we get
average exon coverage score, average exon score, number of bases above
some coverage threshold.

## Visualising coverage scores

After identifying some gene of interest - you can visualise the results
per gene via `plot_gene_coverage()`
