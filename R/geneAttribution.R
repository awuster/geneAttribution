#' geneAttribution: Identification of candidate genes associated with noncoding genetic variation
#'
#' Identification of the most likely gene or genes through which variation at a given genomic
#' locus in the human genome acts. The most basic functionality assumes that the closer gene is
#' to the input locus, the more likely the gene is to be causative. Additionally, any empirical
#' data that links genomic regions to genes (e.g. eQTL or genome conformation data) can be used
#' if it is supplied in UCSC .bed file format.
#' A typical workflow requires loading gene models and empirical data, then running ga() on the
#' locus of interest
#' @docType package
#' @name geneAttribution
NULL

#' Load gene models
#'
#' Get gene models as a GenomicRanges object, with gene names in the symbol column
#' For hg19, you may want to use TxDb.Hsapiens.UCSC.hg19.knownGene and for GRCh38,
#' TxDb.Hsapiens.UCSC.hg38.knownGene (set as default)
#' @param txdb GenomicFeatures TxDb object containing genomic coordinates of genes
#' @param maxGeneLength An integer. Gene models that are longer than this are excluded. Optional
#' @param genesToInclude A character vector of gene symbols of genes to include (e.g. only protein coding genes). Optional
#' @param genesToExclude A character vector of gene symbols of genes to exclude. Optional
#' @return A GenomicRanges object containing human gene models
#' @examples
#' geneModels()
#' geneModels(genesToInclude = c("CYYR1", "ADAMTS1", "ADAMTS5", "N6AMT1", "LTN1"))
geneModels <- function(txdb=TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, maxGeneLength=1e+06, genesToInclude, genesToExclude) {

    # Get data
    geneLocs <- GenomicFeatures::genes(txdb)
    symbolDf <- BiocGenerics::as.data.frame(org.Hs.eg.db::org.Hs.egSYMBOL2EG)
    geneLocs$symbol <- symbolDf[match(geneLocs$gene_id, symbolDf$gene_id), "symbol"]

    # Filtering
    geneLocs <- geneLocs[GenomicRanges::width(geneLocs) < maxGeneLength]
    if (!missing(genesToInclude)) { geneLocs <- geneLocs[which(geneLocs$symbol %in% genesToInclude)] }
    if (!missing(genesToExclude)) { geneLocs <- geneLocs[which(! geneLocs$symbol %in% genesToExclude)] }

    # Report gene number and genome build
    message("Gene models for ", length(geneLocs), " genes returned\n")
    build <- GenomeInfoDb::genome(geneLocs)[1]
    message("Gene models are in build ", build, "\n")
    return(geneLocs)
}

#' Load UCSC *.BED files containing empirical data
#'
#' Required *.BED file format (tab-separated): chr start end name (optional column: score).
#' Sample files supplied with package are limited to chromosome 2.
#' @param files A character vector containing *.BED file names
#' @param weights An integer vector containing weighting for each bed file. Optional
#' @return A list of GenomicRanges objects containing the data from the *.BED files, with weightings in the score column
#' @examples
#' fileName1 <- system.file("extdata", "hiCRegions.b38.bed", package="geneAttribution")
#' fileName2 <- system.file("extdata", "eqtlHaplotypeBlocks.b38.bed", package="geneAttribution")
#' loadBed(c(fileName1, fileName2), c(2, 5))
loadBed <- function(files, weights) {

    # Make sure that file vector and (optional) weighting vector are of the same length
    if (!missing(weights)) {
	if (length(files) != length(weights)) { stop("Input vectors specifying files and weights must have same length.\n") }
    }

    # Go through each BED file in input and add to list
    empOut <- list()
    for (ii in 1:length(files)) {

	empF <- files[ii]
	message("\nNow loading BED file ", empF, "\n")
	empG <- rtracklayer::import.bed(empF) # Read BED file as GRanges object

	# Add weighting if not already in score column in BED file
	# Report back on which weighting was used
	if ("score" %in% names(GenomicRanges::mcols(empG))) {
	    if (missing(weights)) {
		message("Using weighting specified in score column of BED file\n")
	    } else {
		empG$score <- weights[ii]
		message("Overriding score column in BED file and assigning weighting of ", weights[ii], "\n")
	    }
	 } else {
	    if (missing(weights)) {
		empG$score <- 2
		message("Setting weighting to 2 (default)\n")
	    } else {
		empG$score <- weights[ii]
		message("Setting weighting to ", weights[ii], "\n")
	    }
	}

	empOut <- c(empOut, empG) # Add to output list
    }

    # Return
    return(empOut)

}

#' Normalize likelihoods and return probabilities
#' 
#' @param pVector A numeric vector of pre-normalization likelihoods
#' @param minPP A float. Genes with a posterior probability < minPP are lumped as "Other". Optional
#' @return A sorted, numeric vector of normalized probabilities
#' @examples
#' normP(c(5, 1, 1, 1, 1, 1, 0.1))
#' normP(c(5, 1, 1, 1, 1, 1, 0.1), minPP=0.1)
normP <- function(pVector, minPP=0) {
    
    pVector <- pVector / sum(pVector) # Normalization
    pVector <- sort(pVector, decreasing=TRUE) # Sort by decreasing p

    # If there are gene with p < minPP, summarize as "Other"
    otherIdx <- which(pVector < minPP)
    if (length(otherIdx) > 0) {
	otherP <- sum(pVector[otherIdx])
	names(otherP) <- "Other"
	pVector <- c(pVector[-otherIdx], otherP)
    }

    return(pVector)
}

#' Given genomic coordinate, return normalized probability for each gene
#' 
#' @param chr A character string representing a chromosome (e.g. "chr2")
#' @param pos An integer representing a genomic position in the same genome build that gene models
#' @param geneCoordinates A GenomicRanges object, as generated by geneModels()
#' @param empiricalData A list of GenomicRanges objects, as generated by loadBed(). Optional
#' @param lambda Float. Variable for exponential distribution. Default based on empirical eQTL data from multiple tissues. Optional
#' @param maxDist Integer. Only genes within this distance of query locus are considered. Optional
#' @param minPP Float. Genes with a posterior probability < minPP are lumped as "Other". Can be set to 0 when all genes should be reported. Optional
#' @return A sorted, numeric vector of normalized gene probabilities
#' @examples
#' geneLocs <- geneModels()
#' fileName <- system.file("extdata", "eqtlHaplotypeBlocks.b38.bed", package="geneAttribution")
#' empirical <- loadBed(fileName)
#' ga("chr2", 127156000, geneLocs, empirical)
ga <- function(chr, pos, geneCoordinates, empiricalData, lambda=7.61e-06, maxDist=1e+06, minPP=0.01) {
    
    query <- GenomicRanges::GRanges(chr, IRanges::IRanges(pos, pos))

    # Distance to query
    geneDistancesAll <- GenomicRanges::distance(query, geneCoordinates, ignore.strand=TRUE) # Distance between input and all genes
    geneDistanceIdx <- which(!is.na(geneDistancesAll) & geneDistancesAll < maxDist) # Get indices for genes with defined distance (on same chromosome, within maxDist) 
    geneDistances <- geneDistancesAll[geneDistanceIdx] # Get distances
    names(geneDistances) <- geneCoordinates[geneDistanceIdx]$symbol # Get associated gene names

    # Gene likelihood based on distance to gene
    likelihood <- stats::dexp(geneDistances, lambda)

    # Prior, based on empirical data
    prior <- rep(1, length(likelihood))
    names(prior) <- names(likelihood)

    # Add empirical data
    if (!missing(empiricalData)) {
	for (iiEmp in empiricalData) {
	    
	    # Which empirical data ranges overlap query?
	    iiOverlapIdx <- GenomicRanges::findOverlaps(query, iiEmp)
	    empSubset <- iiEmp[as.data.frame(iiOverlapIdx)[[2]]]
	    empWeights <- empSubset$score
	    names(empWeights) <- empSubset$name

	    # Only keep weightings for genes that are within maxDist
	    empWeights <- empWeights[names(prior)]
	    empWeights[which(is.na(empWeights))] <- 0

	    # Add weights to prior
	    prior <- prior + empWeights

	}
	prior <- prior / sum(prior) # Normalization
    }

    # (Likelihood * Prior) / sum(Likelihood * Prior)
    return(normP(likelihood * prior, minPP=minPP))

}
