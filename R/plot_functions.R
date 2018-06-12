#' Make a plot of the LOI genes
#'
#' @importFrom pheatmap pheatmap
#' @export
#'
plot_guessedLOI <- function(heteroGene, geneOrder=NULL, minToSee=1, flatTo=5, gaps_row=NULL, gaps_col=NULL, addThisMono=NULL, removeME=NULL, plot=T) {
  # library(pheatmap, lib.loc="/Users/paolo/bioinfotree/prj/guessMyLoi-r-package/loi-analysis/guess-LOI/rlib")
  require(RColorBrewer)

  if(!is.null(addThisMono)) {
    heteroGene$geneLOI <- rbind(heteroGene$geneLOI, addThisMono$geneLOI)
    heteroGene$geneAnnot <- rbind(heteroGene$geneAnnot, addThisMono$geneAnnot)
  }

  if (!is.null(geneOrder)) {
    heteroGene <- sortByGenes(heteroGene, geneOrder)
  } else {
    # heteroGene <- sortGeneByPosition(heteroGene)
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

  if (all(flatteredGeneSnpCount==0 | is.na(flatteredGeneSnpCount)))
    return("No SNP found in X Chromosomes.")

  height=80
  width=7
  if(NROW(flatteredGeneSnpCount) ==1) {
    flatteredGeneSnpCount <- t(flatteredGeneSnpCount)
    width=3
  }

  if(plot) {
    pheatmap(flatteredGeneSnpCount, cluster_cols=F, cluster_rows=F,
           color= color,
           border_col=bc,
           # filename=paste0(opts$'pdf-out-suffix', "gene-with-LOI.pdf"),
           na_col="#E0E0E0",
           gaps_col=gaps_col,
           gaps_row=gaps_row,
           # cellwidth = 12,
           # cellheight = 12,
           # width=width, height=height,
           fontsize = 6
           )

  }
  return(heteroGene)
}

# plotChrSample <- function(obj, na.rm=T) {
#   s=colnames(obj$ratios)
#   if (length(s)!=2)
#     stop("Anomalus imput data.")
#
#   values <- lapply(colnames(obj$ratios), function(s) {
#     cbind(obj$ratios[,s]*100, 100-obj$ratios[,s]*100)
#   })
#
#   if (na.rm) {
#     keep <- !apply(is.na(values[[1]]), 1, any)
#     values <- lapply(values, function(v) v[keep, , drop=F])
#     annot <- obj$annot[keep, , drop=F]
#   }
#
#   par(mfrow=c(2,1))
#   barplot(t(values[[1]]), col=c("orange","blue"),
#           names.arg=annot$id, main=s[1],
#           las=3, cex.axis=0.6, cex.names=0.4, ylab="% allelic Expression")
#   barplot(t(values[[2]]), col=c("orange","blue"),
#           names.arg=annot$id, main=s[2],
#           las=3, cex.axis=0.6, cex.names=0.4, ylab="% allelic Expression")
# }

#' Make a plot of the allelic ratios
#'
#' @export
#'
plotAllelicRatios <- function(obj, samples=NULL, na.rm=F, removeSpecificIds=NULL) {
  require(ggpubr)
  o <- obj
  if (!is.null(samples)) {
    o <- filter_columns(o, samples)
  }

  kp <- apply(o$isHetero==1, 1, any, na.rm=T)
  if (na.rm) {
    naKp <- !apply(is.na(o$isHetero), 1, any)
    kp = kp & naKp
  }

  o <- filter_rows(o, kp)
  if (nrow(o$ratios)==0) {
    return(list(data=data.frame(), plot="No SNP available for genes."))
  }

  data = data.frame(o$ratios, stringsAsFactors = F, check.names = F)
  inv <- 1-data

  data$allele <- rep("alt", NROW(data))
  data$id <- factor(o$annot$id, levels=o$annot$id[order(o$annot$pos)])
  data$pos <- o$annot$pos[order(o$annot$pos)]

  inv$allele <- rep("ref", NROW(inv))
  inv$id <- factor(o$annot$id, levels=o$annot$id[order(o$annot$pos)])
  inv$pos <- o$annot$pos[order(o$annot$pos)]

  data <- rbind(data,inv)

  idx <- which(colnames(data) %in% c("allele","id","pos"))

  filter <- !data$id %in% removeSpecificIds
  data <- data[filter, , drop=F]

  y <- colnames(data)[-idx]
  if (length(y)==1) {
    yidx <- match(y, colnames(data))
    colnames(data)[yidx] <- "sample"
    y <- "sample"
  }
  p = ggbarplot(data, x="id",
            y = y,
            palette=c("red","blue"),
            combine = TRUE,
            x.text.angle = 90,
            fill="allele",
            color = "white")
  return(list(data=data, plot=p))

}

#' Make a plot of CHR and SNPs
#'
#' @export
#'
plotChrXwithSNPs <- function(snps, karyogramm, useRepel=TRUE) {
  chrX <- karyogramm[seqnames(karyogramm) == "chrX"]
  p <- ggplot() + layout_karyogram(chrX, cytoband = TRUE, geom = NULL)

  thik <- data.frame(x1=snps$pos, x2=snps$pos+1, y1=0, y2=10.5, seqnames = "chrX")
  p <- p + ggplot2::geom_rect(data = thik,
                              do.call(aes, list(xmin = substitute(x1), xmax = substitute(x2),
                                                ymin = substitute(y1), ymax = substitute(y2))),
                              color = "red", fill = "red", size = 0.5, alpha = 0.2)

  lbs.df <- data.frame(x1=snps$pos, y2=11, ids=snps$id)

  if(useRepel) {
    p <- p + geom_text_repel(data = lbs.df, aes(x1, y2+0.5, label=ids),
                             angle = 90, size = 2,
                             direction = "x", segment.color = "white")
  } else {
    p <- p + ggplot2::geom_text(data = lbs.df, aes(x1, y2, label=ids),
                                angle = 90, size = 3,
                                check_overlap =TRUE,
                                vjust=0.5, hjust=0, nudge_x = 5)
  }

  p <- p + theme_alignment(grid = FALSE, ylabel = TRUE, border = FALSE) +
    scale_y_continuous(breaks = 5, labels = "chrX", limits=c(0,14)) +
    theme(strip.background = element_rect(colour = "NA", fill = "NA")) +
    theme(strip.text.y = element_text(colour = "white")) +
    theme(legend.position = "none") + ggplot2::xlab("")

  p <- p + theme(axis.ticks.y = element_blank())

  return(p)
}
