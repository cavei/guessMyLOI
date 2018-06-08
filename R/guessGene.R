guessGene <- function(file, gene, annotCol=7, depth = 10, countsToCall = 5, total = TRUE, callheteroThr=0.2) {
  gatk_ai <- readGATK.AllelicImbalance(file=file, annotCol = annotCol)

  gatk_ai <- filter_by_overall_depth(gatk_ai, depth = depth)

  if (!(gene %in% gatk_ai$alleleAnnot$gene))
    return("Gene not found.")

  keep = gatk_ai$alleleAnnot$gene == gene
  gene_ai <- filter_rows(gatk_ai, keep)

  g <- computeAlleleRatios(gene_ai, countsToCall = 5, total = TRUE)
  g <- callHetero(g, thr = 0.2)

  keepSNPDB = g$annot$id != "."
  g <- filter_rows(g, keepSNPDB)
  plotAllelicRatios(g, na.rm = F)
}


