#' Compute Allelic Rations
#'
#' @export
#'
computeAlleleRatios <- function(obj, countsToCall, total = FALSE) {
  # countsToCall is the minimum value for the minor allele [HARD check]
  ratio <- tapply(seq_len(NROW(obj$alleleCounts)), obj$positionUID, function(x) {
    computeHeteroScore(obj$alleleCounts, x, counts=countsToCall, total = total) # MISS
  })

  if (NCOL(ratio)==1 & !is.list(ratio)) {
    ratiosMatrix <- data.matrix(ratio)
    colnames(ratiosMatrix) <- colnames(obj$alleleCounts)
  } else {
    ratiosMatrix <- do.call(rbind, ratio)
    ratiosMatrix <- suppressWarnings(data.matrix(ratiosMatrix))
  }

  ratioString <- tapply(seq_len(NROW(obj$alleleCounts)), obj$positionUID, function(x) {
    ratioInString(obj$alleleCounts, x, total = total) # MISS
  })

  if (NCOL(ratioString)==1 & !is.list(ratioString)) {
    ratioStringsMatrix <- data.matrix(ratioString)
    colnames(ratioStringsMatrix) <- colnames(obj$alleleCounts)
  } else {
    ratioStringsMatrix <- do.call(rbind, ratioString)
    ratioStringsMatrix <- suppressWarnings(data.matrix(ratioStringsMatrix))
  }

  snpAnnot <- obj$alleleAnnot
  snpAnnot$variant <- NULL
  snpAnnot <- unique(snpAnnot)
  snpAnnotPositionUID <- paste(snpAnnot$chr, snpAnnot$pos,sep="_")
  row.names(snpAnnot) <- snpAnnotPositionUID

  ratios=ratiosMatrix[snpAnnotPositionUID, ,drop=F]
  ratioStringsMatrix=ratioStringsMatrix[snpAnnotPositionUID, , drop=F]
  colnames(ratios) <- colnames(obj$alleleCounts)
  annot=snpAnnot
  posUID=snpAnnotPositionUID
  list(ratios=ratios, annot=annot, ratioStringsMatrix=ratioStringsMatrix, posUID=posUID)
}

#' Compute Hetero Score
#'
#' @export
#'
computeHeteroScore <- function(alleles, idx, counts, total=FALSE) {
  tot <- apply(alleles[idx, ,drop=F], 2, function(r) {
    if (all(is.na(r)))
      return(NA)

    if (any(is.na(r))) {
      return("NN")
    }

    if (min(r) >= counts) {
      if (total) {
        min(r)/(max(r)+min(r))
      } else {
        min(r)/max(r)
      }
    } else if (max(r) < counts) {
      NA
    } else {
      0
    }
  })
  if(any(tot=="NN", na.rm=T)) {
    save(tot, file="tot.RData")
    warning(paste0("NA vs Number at lines ", paste(idx, collapse=" ")))
  }
  suppressWarnings(as.numeric(tot))
}

#' Compute Call if Hetero
#'
#' @export
#'
callHetero <- function(objr, thr=0.2) {
  heteroMatrix <- objr$ratios
  heteroMatrix[heteroMatrix > thr] <- 1
  heteroMatrix[heteroMatrix <= thr] <- 0
  objr$isHetero <- heteroMatrix
  return(objr)
}

#' Call Hetero per gene
#'
#' @export
#'
callHeteroSNPperGene <- function(objr) {
  guessLOIgenes <- tapply(seq_len(NROW(objr$isHetero)), objr$annot$gene, function(x) {
    if (length(x) < 2) {
      objr$isHetero[x,]
    } else {
      colSums(objr$isHetero[x,,drop=F], na.rm=T)
    }
  })

  guessLOIratio <- tapply(seq_len(NROW(objr$ratios)), objr$annot$gene, function(x) {
    if (length(x) < 2) {
      objr$ratios[x,]
    } else {
      colMeans(objr$isHetero[x,,drop=F], na.rm=T)
    }
  })

  geneAnnot <- extractFirstSnpPositionByGene(objr$annot)

  if (NCOL(guessLOIgenes)==1 & (!is.list(guessLOIgenes))) {
    guessLOIgenes <- data.matrix(guessLOIgenes)
    colnames(guessLOIgenes) <- colnames(objr$isHetero)
  } else {
    guessLOIgenes <- do.call(rbind, guessLOIgenes)
  }

  if (NCOL(guessLOIratio)==1 & (!is.list(guessLOIratio))) {
    guessLOIratio <- data.matrix(guessLOIratio)
    colnames(guessLOIratio) <- colnames(objr$isHetero)
  } else {
    guessLOIratio <- do.call(rbind, guessLOIratio)
  }

  guessLOIratio[is.na(guessLOIratio)] <- NA

  if (!(NROW(guessLOIgenes)==NROW(geneAnnot)))
    stop("Matrix dimension differs.")

  if (!(NROW(guessLOIratio)==NROW(geneAnnot)))
    stop("Matrix dimension differs.")

  if (!identical(row.names(guessLOIgenes), row.names(geneAnnot)))
    guessLOIgenes <- guessLOIgenes[row.names(geneAnnot), , drop=F]

  if (!identical(row.names(guessLOIratio), row.names(geneAnnot)))
    guessLOIratio <- guessLOIratio[row.names(geneAnnot), , drop=F]

  list(geneLOI=guessLOIgenes, geneAnnot=geneAnnot, geneRatio=guessLOIratio)
}

#' Extract First Snp Position By Gene
#'
#' @export
#'
extractFirstSnpPositionByGene <- function(snpAnnot, return_matrix=TRUE) {
  r <- tapply(seq_len(NROW(snpAnnot)), snpAnnot$gene, function(idx) {
    if (length(idx) == 1) {
      snpAnnot[idx,c("chr","pos","gene")]
    } else {
      c(unique(snpAnnot[idx,"chr"]), min(snpAnnot[idx,"pos"]), unique(snpAnnot[idx,"gene"]))
    }
  })

  if (return_matrix)
    r <- do.call(rbind, r)
  colnames(r) <- c("chr","pos","gene")
  r
}

#' Clip to values
#'
#' @export
#'
clipToValues <- function(m, cmin, cmax) {
  t(apply(m, 1, function(x) {
    y <-ifelse(x <= cmin, cmin, x)
    y <-ifelse(y >cmax, cmax, y)
  }))
}

#' Sort By Gene
#'
#' @export
#'
sortByGenes <- function(heteroGene, geneOrder) {
  ordine <- sapply(geneOrder, function(g){
    which(sapply(strsplit(heteroGene$geneAnnot$gene, ";"), function(c) any(c %in% g)))
  })
  ordine <- unlist(ordine)

  geneLOI=heteroGene$geneLOI[ordine, , drop=F]
  geneAnnot=heteroGene$geneAnnot[ordine, ,drop=F]
  geneRatio=heteroGene$geneRatio[ordine, , drop=F]

  sel <- !(duplicated(geneAnnot$gene))
  list(geneLOI=geneLOI[sel, ,drop=F], geneAnnot=geneAnnot[sel, ,drop=F], geneRatio=geneRatio[sel, ,drop=F])
}

#' Get Genes With No SNPs
#'
#' @export
#'
nonUsedGenes <- function(heteroGene, geneOrder) {
  setdiff(geneOrder,unlist(strsplit(heteroGene$geneAnnot$gene, ";")))
}

#' Merge HeteroGene Objects
#'
#' @export
#'
mergeHeteroGeneObjs <- function(obj1, obj2) {
  if (length(intersect(row.names(obj1$geneLOI), row.names(obj2$geneLOI))))
    stop("Shared Row Names In geneLOI")
  geneLOI <- rbind(obj1$geneLOI, obj2$geneLOI)

  if (length(intersect(row.names(obj1$geneAnnot), row.names(obj2$geneAnnot))))
    stop("Shared Row Names In geneAnnot")
  geneAnnot <- rbind(obj1$geneAnnot, obj2$geneAnnot)

  if (length(intersect(row.names(obj1$geneRatio), row.names(obj2$geneRatio))))
    stop("Shared Row Names In geneRatio")
  geneRatio <- rbind(obj1$geneRatio, obj2$geneRatio)

  list(geneLOI=geneLOI,
       geneAnnot=geneAnnot,
       geneRatio=geneRatio)
}

#' Sort Gene By Position
#'
#' @export
#'
sortGeneByPosition <- function(heteroGene) {
  ord <- order(as.numeric(heteroGene$geneAnnot$pos))
  list(geneLOI=heteroGene$geneLOI[ord, , drop=F],
       geneAnnot=heteroGene$geneAnnot[ord, ,drop=F],
       geneRatio=heteroGene$geneRatio[ord, , drop=F])
}


#' set NA to FALSE
#'
#' @export
#'
setNA2FALSE <- function(x) {
  res <- x
  res[is.na(x)] <- FALSE
  return(res)
}

#' Set NA 2 TRUE
#'
#' @export
#'
setNA2TRUE <- function(x) {
  res <- x
  res[is.na(x)] <- TRUE
  return(res)
}

#' Create a string representing the ration in string
#'
#' @export
ratioInString <- function(alleles, idx, total=FALSE) {
  tot <- apply(alleles[idx, ,drop=F], 2, function(r) {
    if (all(is.na(r)))
      return(NA)

    if (any(is.na(r))) {
      return("NN")
    }

    if (total) {
      paste0(min(r),"/",(max(r)+min(r)))
    } else {
      paste0(min(r),"/",max(r))
    }
  })
}

#' Collapse object
#'
#' @export
collapseObject <- function(object, collapse) {
  unknow <- setdiff(collapse, names(object))
  if (length(unknow)>0)
    stop("Choose who to collapse properly")

  if (length(collapse)==1)
    return(object[[collapse]])

  if (!(is.matrix(object[[collapse[1]]])| is.data.frame(object[[collapse[1]]]))) {
    export <- matrix(object[[collapse[1]]], ncol=1)
  } else {
    export <- object[[collapse[1]]]
  }

  for (cllp in collapse[2:length(collapse)]) {
    if (is.matrix(object[[cllp]])) {
      if(!identical(row.names(export), row.names(object[[cllp]])))
        stop(paste0(cllp, " have different row.names."))
      export <- cbind(export, object[[cllp]])
    } else {
      if(!identical(row.names(export), names(object[[cllp]])))
        stop(paste0(cllp, " have different row.names."))
      export <- cbind(export, object[[cllp]])
    }
  }
  export
}
