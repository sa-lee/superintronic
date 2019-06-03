context("group by + split")

lvls <- paste0("chr",c(1:23))
gr <- GRanges(
  seqnames = factor(sample(lvls, 1000, replace = TRUE), levels = lvls),
  ranges = IRanges(start = rpois(1000, 10000), width = rpois(1000, 100)),
  grp = sample(letters[1:4], 1000, replace = TRUE),
  gc = runif(1000)
)

groupings <- unique(select(gr, grp, seqnames, .drop_ranges = TRUE))

test_that("returns a GRangesList with desired properties", {
  # no split just returns
  expect_identical(split_ranges(gr), gr)
  
  # splitting ungrouped produces an unnamed GRangesList
  ans <- split_ranges(gr, seqnames)
  expect_true(is(ans, "GRangesList"))
  expect_equal(length(ans), length(lvls))
  expect_null(names(ans))
  
  # using grouped produces the same answer
  ans2 <- split_ranges(group_by(gr, seqnames))
  expect_identical(ans, ans2)
  expect_true(is(ans2, "GRangesList"))
  expect_equal(length(ans2), length(lvls))
  expect_null(names(ans2))
  
  # grouping again produces a warning
  expect_warning(split_ranges(group_by(gr, seqnames), strand))
  
  # multiple groups
  ans <- split_ranges(gr, grp, seqnames)
  expect_true(is(ans, "GRangesList"))
  expect_equal(length(ans), nrow(groupings))
  expect_null(names(ans))
  
  ans2 <-split_ranges(group_by(gr, grp, seqnames))
  expect_identical(ans, ans2)
  expect_equal(length(ans2), nrow(groupings))
  expect_null(names(ans2))
  
})
