#' Create X Y linked SNP's allelic expression
#'
#' @param gatk_ai allelic imbalance from \code{readGATK.AllelicImbalance}
#' @param countsToCall how many counts to the alternative allele to call
#' @param total use tota as denominator in the allele ratios
#' @param callheteroThr the threshold to call hetorozygous
#' @param species either hsapiens or mmusculus. Choose the species to remove pseudoAutosomal regions
#' @param removeNonDBsnp should I remove the non DBsnp SNPs
#'
#' @return a list with
#'   \item{genes}{genes related object for visualization}
#'   \item{snps}{snps related object for visualization}
#'   \item{geneSummary}{tabular summary of the genes}
#'   \item{snpSummary}{tabular summary of the SNPs}
#'
#' @rdname processData
#'
#' @export
#'
processSexualData <- function(gatk_ai, countsToCall, total, callheteroThr, species=c("hsapien","mmusculus"), removeNonDBsnp = TRUE) {
  species <- species[1]
  removePseudoAutosomalRegions <- switch(species,
                                         hsapiens = removeHg38PsudosomalRegion,
                                         mmusculus = removeM38PsudosomalRegion
         )

  keep = gatk_ai$alleleAnnot$chr == "X" | gatk_ai$alleleAnnot$chr == "chrX" |
  gatk_ai$alleleAnnot$chr == "Y" | gatk_ai$alleleAnnot$chr == "chrY"
  chrxy_ai <- filter_rows(gatk_ai, keep)

  chrxy_ai <- removePseudoAutosomalRegions(chrxy_ai)

  chrXY <- computeAlleleRatios(chrxy_ai, countsToCall = countsToCall, total = TRUE)
  chrXY <- callHetero(chrXY, thr = callheteroThr)

  if (removeNonDBsnp) {
    keepSNPDB = chrXY$annot$id != "."
    chrXY <- filter_rows(chrXY, keepSNPDB)
  }

  chrXYsnpSummary <- collapseObject(chrXY, names(chrXY)[c(2,3,1,5)])
  ChrXYgenes <- callHeteroSNPperGene(chrXY) ## Rtn

  chrXYgenesSummary <- collapseObject(ChrXYgenes, names(ChrXYgenes)[c(2,1,3)])
  list(genes=ChrXYgenes, snps=chrXY, geneSummary=chrXYgenesSummary, snpSummary=chrXYsnpSummary)
}

#' Create Non Sexual linked SNP's allelic expression
#' @inheritParams processSexualData
#'
#' @rdname processData
#'
#' @export
#'
processNonSexualData <- function(gatk_ai, countsToCall,  total, callheteroThr, species=c("hsapien","mmusculus"), removeNonDBsnp = TRUE){
  species <- species[1]
  name <- switch(species,
         hsapiens = data("Hg38_gene_annotation"),
         mmusculus = data("M38_gene_annotation")
  )
  gene_annot <- get(name)

  keep = gatk_ai$alleleAnnot$chr != "X" & gatk_ai$alleleAnnot$chr != "chrX" &
    gatk_ai$alleleAnnot$chr != "Y" & gatk_ai$alleleAnnot$chr != "chrY"

  gatk_ai_nox <- filter_rows(gatk_ai, keep)

  aiGenes <- computeAlleleRatios(gatk_ai_nox, countsToCall = countsToCall, total = total)
  aiGenes <- callHetero(aiGenes, thr = callheteroThr)

  if (removeNonDBsnp) {
    keepSNPDB = aiGenes$annot$id != "."
    aiGenes <- filter_rows(aiGenes, keepSNPDB)
  }
  LOIsnpSummary <- collapseObject(aiGenes, names(aiGenes)[c(2,3,1,5)])
  aiGenesGenes <- callHeteroSNPperGene(aiGenes)
  geneWithNoSNPdetected <- nonUsedGenes(heteroGene=aiGenesGenes, gene_annot$gene)

  fakeData <- createArtificialHeteroGene(geneWithNoSNPdetected, colnames(aiGenesGenes$geneLOI), gene_annot = gene_annot)
  aiGenesGenes <- mergeHeteroGeneObjs(aiGenesGenes, fakeData)
  aiGenesGenes <- sortByGenes(aiGenesGenes, gene_annot$gene)

  LOIgenesSummary <- collapseObject(aiGenesGenes, names(aiGenesGenes)[c(2,1,3)])

  list(genes=aiGenesGenes, snps=aiGenes, geneSummary=LOIgenesSummary, snpSummary=LOIsnpSummary)
}

#' @export
createArtificialHeteroGene <- function(gene_list, samples, gene_annot){
  sel <- gene_annot$gene %in% gene_list
  if (sum(sel)==0)
    return(NULL)

  geneAnnot <- gene_annot[sel, c("chr", "pos", "gene"),drop=F]
  if (any(duplicated(geneAnnot$gene)))
    stop("Invalid: gene duplication")

  row.names(geneAnnot) <- geneAnnot$gene

  geneLOI <- matrix(0, ncol=length(samples), nrow=length(geneAnnot$gene), dimnames = list(geneAnnot$gene, samples))
  geneRatio <- matrix(0, ncol=length(samples), nrow=length(geneAnnot$gene), dimnames = list(geneAnnot$gene, samples))
  return(list(geneLOI=geneLOI, geneAnnot=geneAnnot, geneRatio=geneRatio))
}
