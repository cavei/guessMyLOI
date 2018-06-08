#' Make a plot of the LOI genes
#'
#' @importFrom pheatmap pheatmap
#' @export
#'
plotAllelicImbalance <- function(heteroGene, geneOrder=NULL, minToSee=1,flatTo=5,
                                 gaps_row=NULL, gaps_col=NULL, addThisMono=NULL,
                                 removeME=NULL, plot=T) {
  require(RColorBrewer)

  if(!is.null(addThisMono)) {
    heteroGene$geneLOI <- rbind(heteroGene$geneLOI, addThisMono$geneLOI)
    heteroGene$geneAnnot <- rbind(heteroGene$geneAnnot, addThisMono$geneAnnot)
  }

  if (!is.null(geneOrder)) {
    heteroGene <- sortByGenes(heteroGene, geneOrder)
  }

  geneSnpCount <- heteroGene$geneLOI
  flatteredGeneSnpCount <- clipToValues(geneSnpCount, 0, flatTo)

  if (FALSE) {
    color = colorRampPalette(c("#000000","#FF6600"))(6)
    bc="grey60"
  } else {
    color = brewer.pal(n = 9, name = "PuBu")[-c(2:4)]
    bc="black"
  }

  color[seq_len(minToSee)]=color[1]

  if (!is.null(removeME)) {
    good = !row.names(flatteredGeneSnpCount) %in% removeME
    flatteredGeneSnpCount <- flatteredGeneSnpCount[good, , drop=F]
    row.names(flatteredGeneSnpCount) <- sub('\\.1', replacement = "", row.names(flatteredGeneSnpCount))
  }

  height=7
  if(NROW(flatteredGeneSnpCount) ==1) {
    flatteredGeneSnpCount <- t(flatteredGeneSnpCount)
    height=3
  }

  if(plot) {
    pheatmap(flatteredGeneSnpCount, cluster_cols=F, cluster_rows=F,
             color= color,
             border_col=bc,
             na_col="#E0E0E0",
             gaps_col=gaps_col,
             gaps_row=gaps_row,
             width=10, height=height)
  }
  return(heteroGene)
}
