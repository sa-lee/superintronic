
<!-- README.md is generated from README.Rmd. Please edit that file -->

# superintronic

[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

Exploring coverage signal in high-throughput (RNA) sequencing data via
coverage estimation.

*superintronic* centers around exploring coverage over exonic and
intronic regions via computing simple summary statistics and
visualisations. The aim is to provide an extremely modular worklfow via
an interface built on top of the
[*plyranges*](https://sa-lee.github.io/plyranges/index.html) package.
This means that you can modify any of the steps provided with the
*plyranges* grammar or just use our defaults.

The basic workflow consists of three steps

1.  Setting up annotations
2.  Computing coverage
3.  Visualising coverage results

## Setting up annotation GRanges

We assume you are starting from a GTF/GFF file for your given organism.
The reading of the GFF file is done external to this package, and can be
done via `rtracklayer::import` or `plyranges::read_gff`. The resulting
GRanges is passed to `collect_parts()`, then you can use
`plyranges::filter` to select the criteria for genes you’re interested
in.

Note that `collect_parts()` should be a generic and dispatch on the
following annotation objects:

  - character
  - GFF/GTFFile
  - TxDb
  - EnsDb
  - Others (??)

The result is a GRanges object with number of rows equal to genes, and
columns containing the intronic and exonic features (as list columns).

``` r
library(superintronic)
library(magrittr)

gff <- "my.gff"

features <- gff %>% 
  plyranges::read_gff() %>% 
  collect_parts() 
```

## Computing coverage over BAM(s)

### Making coverage scores into a tidy GRanges

Here we compute a long form GRanges containing coverage scores in
parallel over a set of BAM files. We have structured it so you specify a
data frame, (like you would get from targets in `limma`), that contains
a column identified by `source`, that specifies the BAM file names.
Other options include specifying a GRanges that represents the genome
build, an optional target GRanges for restricting coverage, and
`BPPARAM` obejct for perfoming parallel computations.

``` r
# is generic - should be able to dispatch on just the names of a BAM file too
compute_coverage_long(design,
                      source,
                      .target = GenomicRanges::GRanges(),
                      .genome_info,
                      .parallel = BiocParallel::bpparam()
                      )
```

This function automatically propagates, the metadata associated with a
design onto the resulting GRanges

As an example, let’s use BAM files from the
[`airway`](https://bioconductor.org/packages/release/data/experiment/html/airway.html)
data package.

``` r
design <- read.csv(system.file("extdata", 
                               "sample_table.csv", 
                               package = "airway"),
                   row.names = FALSE)
design$bam <- dir(system.file("extdata", package = "airway"), 
                  pattern = "*.bam",
                  recursive = TRUE)

cvg <- compute_coverage_long(design, source = "bam")
cvg
```

Once the coverage is in long form, we can then merge overlapping genomic
features over the experimental design via an intersect overlap join and
nesting over the union ranges of an index index (usually gene\_id).

``` r
cvg_over_features <- nest_by_targets(cvg,
                                  target, 
                                  key = sample,
                                  range_index = gene_id
                                  )
```

### Computing coverage wide form

Returns a `RaggedExperiment` object (useful for doing things like PCA
etc.)

``` r
compute_coverage_wide(design,
                      source,
                      .target = GenomicRanges::GRanges(),
                      .genome_info,
                      .parallel = BiocParallel::bpparam()
                      )
```

## Coverage diagnostics

We can then compute a bunch of `covnostics`, over a given index accross
a column representing the coverage score.

  - mean
  - min, median, max
  - quantile
  - variance
  - variance\_shift
  - bases\_above\_score
  - near\_feature\_covnostic

This is conceptually similar to `dplyr::summarise_at`

``` r
covnostics(cvg_over_features, y = score, .funs, ...)
```

## Visualising coverage scores

Options for visualising coverage over a given range

``` r
view_coverage(cvg, target)
```

This returns a regular old ggplot object, so segments can be overlaid
with by adding to the plot object.
