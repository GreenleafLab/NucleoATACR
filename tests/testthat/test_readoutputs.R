context("Reading NucleoATAC outputs")

# Read in data needed for tests ------------------------------------------------

test_nucpos <- system.file("extdata", "test.nucpos.bed.gz", package = "NucleoATACR")
test_occpeaks <- system.file("extdata", "test.occpeaks.bed.gz", package = "NucleoATACR")
test_nucmap <- system.file("extdata", "test.nucmap_combined.bed.gz", package = "NucleoATACR")
test_nfr <- system.file("extdata", "test.nfrpos.bed.gz", package = "NucleoATACR")
test_nucleoatac_signal <- system.file("extdata", "test.nucleoatac_signal.bedgraph.gz", package = "NucleoATACR")


# Reading in nucleosome positions

test_that("can read nucpos file", {
  nucs <- readNucs(test_nucpos)
  expect_is(nucs, "GenomicRanges")
  expect_equal_to_reference(nucs, "nucpos.rds")
})

test_that("can read occpeaks file", {
  nucs <- readNucs(test_occpeaks)
  expect_is(nucs, "GenomicRanges")
  expect_equal_to_reference(nucs, "occpeaks.rds")
})

test_that("can read nucmap_combined file", {
  nucs <- readNucs(test_nucmap)
  expect_is(nucs, "GenomicRanges")
  expect_equal_to_reference(nucs, "nucmap_combined.rds")
})

# Reading in NFR positions -----------------------------------------------------

test_that("can read nfr file", {
  nfrs <- readNFRs(test_nfr)
  expect_is(nfrs, "GenomicRanges")
  expect_equal_to_reference(nfrs, "nfrs.rds")
})

# Reading in bedgraph file -----------------------------------------------------

test_that("can read tabix indexed bedgraph file", {
  sig <- readBedgraph(test_nucleoatac_signal, "chrII", 706551, 707705)
  expect_is(sig, "numeric")
  expect_true(is.na(sig[1]))
  expect_equal(sig[2], -0.413026576715)
  expect_equal(sig[11], 0.548198534934)
  expect_equal(length(sig), 707705 - 706551)
  expect_equal_to_reference(sig, "nucleoatac_signal.rds")
})






