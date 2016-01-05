context("Finding relative positions")

# Read in data needed for tests ------------------------------------------------

test_nucpos <- system.file("extdata", "test.nucpos.bed.gz", package = "NucleoATACR")
nucs <- readNucs(test_nucpos)

# Make up some positions -------------------------------------------------------

test_pos <- GenomicRanges::GRanges("chrII", IRanges::IRanges(start = c(707200,707500), width = 1), strand = c("+","-")) 

# Distances between calls ------------------------------------------------------

test_that("get_dist_between_calls works",{
    dists <- get_dist_between_calls(nucs[1:5])
    expect_equal(dists[[195]] , 1)  
    expect_equal(dists[[163]] , 1)
    expect_equal(dists[[413]] , 1)
    expect_equal(sum(dists[-c(195,163,413)]),0)
    expect_equal(length(dists),1000)
    })

# Test distRanges --------------------------------------------------------------

test_that("distRanges works", {
  dists <- distRanges(test_pos, nucs)
  expect_equal(dists[1],14)
  expect_equal(dists[2],-44)
})

# Test +1 and -1 ---------------------------------------------------------------

test_that("+1 finding works", {
  p1 <- get_p1_nuc(nucs.ranges = nucs, sites = test_pos)
  expect_equal(BiocGenerics::start(p1)[1],707349)
  expect_equal(BiocGenerics::start(p1)[2],707349)
  expect_equal(p1$dist[1], 149)
  expect_equal(p1$dist[2],151)
})

test_that("-1 finding works", {
  m1 <- get_m1_nuc(nucs.ranges = nucs, sites = test_pos)
  expect_equal(BiocGenerics::start(m1)[1],707186)
  expect_equal(BiocGenerics::start(m1)[2],707544)
  expect_equal(m1$dist[1], -14)
  expect_equal(m1$dist[2],-44)
})



