library(geneAttribution)
context("Output of geneAttribtion functions is as expected")

test_that("geneModels() output is a GenomicRanges object with the required fields", {
	      gm <- geneModels()
	      expect_that(gm, is_a("GRanges"))
	      expect_that(gm$symbol, is_a("character"))
})

test_that("loadBed() output is a list of Genomic Ranges object with required fields", {
	      fileName <- system.file("extdata", "hiCRegions.b38.bed", package="geneAttribution")
	      bed <- loadBed(c(fileName))
	      expect_that(bed, is_a("list"))
	      expect_that(bed[[1]], is_a("GRanges"))
	      expect_that(bed[[1]]$name, is_a("character"))
})

test_that("normP() output is a vector that adds up to 1", {
	      np <- sum(normP(c(5, 1, 1, 1, 1, 1, 0.1)))
	      expect_that(np, is_a("numeric"))
	      expect_equal(sum(np), 1)
})

test_that("geneAttribution output is a named vector that adds up to 1", {
	      geneLocs <- geneModels()
	      fileName <- system.file("extdata", "eqtlHaplotypeBlocks.b38.bed", package="geneAttribution")
	      empirical <- loadBed(fileName)
	      ga <- geneAttribution("chr2", 127156000, geneLocs, empirical)
	      expect_that(ga, is_a("numeric"))
	      expect_that(names(ga), is_a("character"))
	      expect_equal(sum(ga), 1)
})
