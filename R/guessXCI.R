#' Guess X Chromosome Inactivation
#'
#' @export
#'
guessXCI <- function(file, annotCol=7, depth = 10, countsToCall = 5, total = TRUE, callheteroThr=0.2, plot=T) {
  gatk_ai <- readGATK.AllelicImbalance(file=file, annotCol = annotCol)
  gatk_ai <- filter_by_overall_depth(gatk_ai, depth = depth)
  keep = gatk_ai$alleleAnnot$chr == "X" | gatk_ai$alleleAnnot$chr == "chrX"

  chrx_ai <- filter_rows(gatk_ai, keep)
  chrx_ai <- removePsudosomalRegion(chrx_ai)

  chrX <- computeAlleleRatios(chrx_ai, countsToCall = countsToCall, total = TRUE)
  chrX <- callHetero(chrX, thr = callheteroThr)

  keepSNPDB = chrX$annot$id != "."
  chrX <- filter_rows(chrX, keepSNPDB)

  plots <- lapply(colnames(chrX$ratios), function(name) {
    plotAllelicRatios(chrX, samples = name, na.rm=T)$plot
  })
  names(plots) <- colnames(chrX$ratios)
  return(plots)
  # small <- plotAllelicRatios(chrX, samples=c("IMR90_F", "HPD08"))
}
