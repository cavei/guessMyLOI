guessSomatic <- function(file, annotCol=7, depth = 10, countsToCall = 5, total = TRUE, callheteroThr=0.2, plot=T) {
  gatk_ai <- readGATK.AllelicImbalance(file=file, annotCol = annotCol)
  gatk_ai <- filter_by_overall_depth(gatk_ai, depth = depth)
  keep = gatk_ai$alleleAnnot$chr != "chrX" & gatk_ai$alleleAnnot$chr != "chrY"
  gatk_ai_nox <- filter_rows(gatk_ai, keep)
  aiGenes <- computeAlleleRatios(gatk_ai_nox,
                                 countsToCall = countsToCall,
                                 total = total)

  aiGenes <- callHetero(aiGenes, thr = callheteroThr)
  keepSNPDB = aiGenes$annot$id != "."
  aiGenes <- filter_rows(aiGenes, keepSNPDB)

  igenes <- callHeteroSNPperGene(aiGenes)
  plot_guessedLOI(igenes, plot=plot)
}

createDatasets <- function(file, annotCol=7, depth = 10, countsToCall = 5, total = TRUE, callheteroThr=0.2) {
  gatk_ai <- readGATK.AllelicImbalance(file=file, annotCol = annotCol)
  gatk_ai <- filter_by_overall_depth(gatk_ai, depth = depth)

  ## Sexual
  keep = gatk_ai$alleleAnnot$chr == "X" | gatk_ai$alleleAnnot$chr == "chrX" |
    gatk_ai$alleleAnnot$chr == "Y" | gatk_ai$alleleAnnot$chr == "chrY"

  chrxy_ai <- filter_rows(gatk_ai, keep)
  chrxy_ai <- removeM38PsudosomalRegion(chrxy_ai)
  chrXY <- computeAlleleRatios(chrxy_ai, countsToCall = countsToCall, total = TRUE)
  chrXY <- callHetero(chrXY, thr = callheteroThr)

  keepSNPDB = chrXY$annot$id != "."
  chrXY <- filter_rows(chrXY, keepSNPDB)

  chrXYsnpSummary <- collapseObject(chrXY, names(chrXY)[c(2,3,1,5)])
  ChrXYgenes <- callHeteroSNPperGene(chrXY) ## Rtn

  chrXYgenesSummary <- collapseObject(ChrXYgenes, names(ChrXYgenes)[c(2,1,3)])

  # xchrData <- plot_guessedLOI(sortGeneByPosition(ChrXYgenes))

  ### NON sexual part
  keep = gatk_ai$alleleAnnot$chr != "X" & gatk_ai$alleleAnnot$chr != "chrX" &
    gatk_ai$alleleAnnot$chr != "Y" & gatk_ai$alleleAnnot$chr != "chrY"

  gatk_ai_nox <- filter_rows(gatk_ai, keep)

  aiGenes <- computeAlleleRatios(gatk_ai_nox, countsToCall = countsToCall, total = TRUE)
  aiGenes <- callHetero(aiGenes, thr = callheteroThr)

  keepSNPDB = aiGenes$annot$id != "."
  aiGenes <- filter_rows(aiGenes, keepSNPDB)
  LOIsnpSummary <- collapseObject(aiGenes, names(aiGenes)[c(2,3,1,5)])

  aiGenesGenes <- callHeteroSNPperGene(aiGenes)
  LOIgenesSummary <- collapseObject(aiGenesGenes, names(aiGenesGenes)[c(2,1,3)])
  return(list(X=list(genes=ChrXYgenes, snps=chrXY), S=list(genes=aiGenesGenes, snps=aiGenes)))
}

