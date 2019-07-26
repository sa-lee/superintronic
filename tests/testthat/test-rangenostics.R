context("Computing rangewise diagnostics")

lvls <- paste0("chr",c(1:23))

gr <- GRanges(
  seqnames = factor(sample(lvls, 1000, replace = TRUE), levels = lvls),
  ranges = IRanges(start = rpois(1000, 10000), width = rpois(1000, 100)),
  grp = sample(letters[1:4], 1000, replace = TRUE),
  score = rnbinom(1000, size = 5, mu = 25)
)

test_that("Returns a GRanges", {
  res <- rangle(gr, score, rng_vars(seqnames), list(mean = mean))
  expect_s4_class(res, "GRanges")
  correct <- summarise(group_by(gr, seqnames), score_mean = sum(score*width) / sum(width))
  expect_equal(res$score_mean, correct$score_mean)
})
