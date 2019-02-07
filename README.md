
<!-- README.md is generated from README.Rmd. Please edit that file -->

# superintronic

Explore intron signal in high-throughput (RNA) sequencing data via
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
set of BAM files. The parallel computation is achieved with
`BiocParallel`. The result is a long GRanges containing columns called
`score` and `source`. An optional GRanges may be passed to this function
to restrict coverage to compute over selected contigs/chromosomes.

``` r
bams <- dir(pattern = "*.bam")
cvg <- gather_coverage(bams)
```

The coverage scores are then restricted by merging the filtered features

``` r
# maybe better name would be restrict_coverage
cvg_over_features <- merge_coverage(cvg, features)
```

Experimental design variables can be then added with

``` r
cvg_over_features <- merge_design(cvg_over_features, design)
```

This requires an additional `data.frame` that contains a column that
matches the BAM file input.

## Summarising Overlaps

There are two modes of summarisation that aggregate over the union set
of introns/exons (called `spread_coverage()`) or over each individual
intron/exon within a gene (called `squish_coverage()`). Both of these
functions are able to compute over variables in experimental design or
within samples.

``` r
# over all samples
spread_coverage(cvg_over_features)
# within each BAM file
spread_coverage(cvg_over_features, source)
# within a design variable
spread_coverage(cvg_over_features, var)

# over all samples
squish_coverage(cvg_over_features)
# within each BAM file
squish_coverage(cvg_over_features, source)
# within a design variable
squish_coverage(cvg_over_features, var)
```

## Visualising coverage scores

The output of `spread_coverage()` can be visualised with `view_exin()`
and coverage plots over genes with `view_coverage_within_gene()`.
