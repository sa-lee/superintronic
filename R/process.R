
# this is equivalent to join_overlap_intersect()
rangeCoverage <- function(gr, cvg) {
  gr <- IRanges::reduce(gr)
  ans <- unlist(cvg[gr], use.names=FALSE)
  if (S4Vectors::runValue(BiocGenerics::strand(gr))[[1L]] == "-")
    ans <- rev(ans)
  ans
}
