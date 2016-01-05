context("Reading and plotting vplots")

# Read in data needed for tests ------------------------------------------------

test_vmat <- system.file("extdata", "test.VMat", package = "NucleoATACR")

# Distances between calls ------------------------------------------------------

test_that("can read in VMat",{
  v <- read_vplot(test_vmat)
  expect_is(v, "data.frame")
})

test_that("can plot v",{
  v <- read_vplot(test_vmat)
  expect_is(plotV(v), "ggplot")
})
