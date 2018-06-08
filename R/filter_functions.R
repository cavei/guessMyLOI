#' Filter Rows
#'
#' @export
#'
filter_rows <- function(obj, keep) {
  lapply(obj, function(o) {
    if (is.matrix(o) || is.data.frame(o)) {
      return(o[keep, , drop=F])
    }
    if (is.character(o)) {
      return(o[keep])
    }
    warning("Don't recognize.")
  })
}

#' Filter Columns
#'
#' @export
#'
filter_columns <- function(obj, keep) {
  lapply(obj, function(o) {
    if (is.matrix(o) || is.data.frame(o)) {
      if (!("chr" %in% colnames(o))) {
        notFound <- setdiff(keep, colnames(o))
        if (length(notFound) != 0) {
          stop(paste0("The following samples are not present: ",
                      paste(notFound, collapse = ", ")))
        }
        return(o[, keep, drop=F])
      }
    }
    return(o)
  })
}

#' Filter By overal Depth
#'
#' @export
#'
filter_by_overall_depth <- function(obj, depth){
  goodSNP <- tapply(seq_len(NROW(obj$alleleCounts)), obj$positionUID, function(x) {
    any(colSums(obj$alleleCounts[x,,drop=F], na.rm = TRUE) >= depth)
  })
  keep <- names(which(goodSNP))
  keep <- obj$positionUID %in% keep
  filter_rows(obj, keep)
}

#' Filter By gene
#'
#' @export
#'
filter_by_gene <- function(obj, gene) {

  geneIdx <- which(obj$alleleAnnot$gene==gene)
  if(length(geneIdx)==0){
    warning("Gene not found")
    return(NULL)
  }
  list(alleleCounts=obj$alleleCounts[geneIdx, , drop=F],
       alleleAnnot=obj$alleleAnnot[geneIdx, , drop=F],
       positionUID=obj$positionUID[geneIdx])
}

#' Filter By gene
#'
#' @export
#'
filter_by_gene.ratios <- function(obj, gene) {
  # "ratios"   "annot"    "posUID"   "isHetero"
  if (!all(c("ratios", "annot", "posUID") %in% names(obj)))
    stop("obj must contain ratios, annot and posUID")

  geneIdx <- grep(gene, obj$annot$gene, fixed = T)
  # geneIdx <- match(gene, obj$annot$gene)
  if(length(geneIdx)==0){
    warning("Gene not found")
    return(NULL)
  }
  out <- list(ratios=obj$ratios[geneIdx, , drop=F],
              annot=obj$annot[geneIdx, , drop=F],
              posUID=obj$posUID[geneIdx])
  if ("isHetero" %in% names(obj))
    out$isHetero <- obj$isHetero[geneIdx, , drop=F]
  out
}

#' Filter Matrix By gene
#'
#' @export
#'
filter_matrix_by_genes <- function(mat, genes) {
  match <- sapply(mat[,1], function(geneD) {
    any(unlist(strsplit(geneD ,split = ";")) %in% genes)
  })
  mat[match, ,drop=F]
}

#' Filter By genes
#'
#' @export
#'
filter_by_genes <- function(heteroGene, genes) {
  match <- sapply(heteroGene$geneAnnot$gene, function(geneD) {
    any(unlist(strsplit(geneD ,split = ";")) %in% genes)
  })

  geneIdx = which(match)

  list(geneLOI=heteroGene$geneLOI[geneIdx, , drop=F],
       geneAnnot=heteroGene$geneAnnot[geneIdx, , drop=F])
}

.filterSamples <- function(samples, obj, cmpSamples=NULL, considerThisSNP=NULL) {
  if (is.null(obj$isHetero))
    stop("You need to call hetero SNPs")

  if (!all(samples %in% colnames(obj$ratios)))
    stop("Your sample was not found in ratios.")

  uninformative <- rep(FALSE, NROW(obj$annot))

  if (!is.null(cmpSamples) & length(samples) != length(cmpSamples))
    stop("Samples must have paired cmpSamples")

  selectHeteroSNP <- considerThisSNP
  if (is.null(considerThisSNP))
    selectHeteroSNP <- names(which(apply(obj$ratios[, samples, drop=F] > 0, 1, any, na.rm=T)))

  res <- lapply(seq_along(samples), function(si) {
    if (is.null(cmpSamples)) {
      filterBySample(samples[si], obj, cmpSamples, considerThisSNP=selectHeteroSNP)
    } else {
      filterBySample(samples[si], obj, cmpSamples[si], considerThisSNP=selectHeteroSNP)
    }
  })

  shrinkAnnot <- res[[1]]$a

  for (i in seq_len(length(samples)-1)){
    if (!identical(shrinkAnnot, res[[i]]$a))
      stop("annotation confusion")
  }

  samplesSnp <- lapply(res, function(r) {
    rbind(ref=1-r$r, alt=r$r)
  })
  names(samplesSnp) <- samples

  gidx <- match(unique(shrinkAnnot$gene),shrinkAnnot$gene)
  glbs <- rep("", length(shrinkAnnot$gene))
  glbs[gidx] <- shrinkAnnot$gene[gidx]

  list(snpsValue = samplesSnp, geneLbs=glbs, annotations=shrinkAnnot)
}

#' Filter By Sample
#'
#' @export
#'
filterBySample <- function(sample, obj, cmpSample=NULL, considerThisSNP=NULL, byGene=TRUE) {
  if (is.null(obj$isHetero))
    stop("You need to call hetero SNPs")

  if (!all(sample %in% colnames(obj$ratios)))
    stop("Your sample was not found in ratios.")

  uninformative <- rep(FALSE, NROW(obj$annot))
  if (!is.null(cmpSample)) {
    if (!cmpSample %in% colnames(obj$ratios))
      stop("Your cmpSample was not found in ratios.")

    selectUninformative <- setNA2FALSE(obj$isHetero[,cmpSample] > 0)

    uninformative <- selectUninformative
    if (byGene) {
      uninformativeGenes <- obj$annot[selectUninformative, "gene"]
      uninformative <- obj$annot[, "gene"] %in% uninformativeGenes
    }
  }

  if (!is.null(considerThisSNP)) {
    selectSNP <- row.names(obj$annot) %in% considerThisSNP
  } else {
    selectSNP <- setNA2FALSE(obj$isHetero[,sample] > 0)
  }

  uninformative <- uninformative[selectSNP]

  true_ratios <- obj$ratios[selectSNP, c(sample, cmpSample), drop=F]
  true_ratios[uninformative, sample] <- NA

  return(list(ratios=true_ratios, annot = obj$annot[selectSNP, , drop=F]))
}

#' Filter Samples
#'
#' @export
#'
filterSamples <- function(samples, ratios, heteroScores, annotations, cmpSamples=NULL, considerThisSNP=NULL) {
  if (!identical(colnames(ratios), colnames(heteroScores)))
    stop("Disorder found in colnames for ratios and heteroscore.")

  if (!identical(row.names(ratios), row.names(heteroScores)))
    stop("Disorder found in row.names for ratios and heteroscore.")

  if (!identical(row.names(ratios), row.names(annotations)))
    stop("Mismatched rownames between annotations and scores.")

  if (!all(samples %in% colnames(ratios)))
    stop("Your sample was not found in ratios.")

  uninformative <- rep(FALSE, NROW(annotations))

  if (!is.null(cmpSamples) & length(samples) != length(cmpSamples))
    stop("Samples must have paired cmpSamples")

  selectHeteroSNP <- considerThisSNP
  if (is.null(considerThisSNP))
    selectHeteroSNP <- names(which(apply(heteroScores[, samples, drop=F] > 0, 1, any)))

  res <- lapply(seq_along(samples), function(si) {
    if (is.null(cmpSamples)) {
      filterBySample(samples[si], ratios, heteroScores, annotations, cmpSamples, considerThisSNP=selectHeteroSNP)
    } else {
      filterBySample(samples[si], ratios, heteroScores, annotations, cmpSamples[si], considerThisSNP=selectHeteroSNP)
    }
  })

  shrinkAnnot <- res[[1]]$a

  for (i in seq_len(length(samples)-1)){
    if (!identical(shrinkAnnot, res[[i]]$a))
      stop("annotation confusion")
  }

  samplesSnp <- lapply(res, function(r) {
    rbind(ref=1-r$r, alt=r$r)
  })
  names(samplesSnp) <- samples

  gidx <- match(unique(shrinkAnnot$gene),shrinkAnnot$gene)
  glbs <- rep("", length(shrinkAnnot$gene))
  glbs[gidx] <- shrinkAnnot$gene[gidx]

  list(snpsValue = samplesSnp, geneLbs=glbs, annotations=shrinkAnnot)
}

#' Remove Mouse Pseudosomal Regions
#'
#' @export
#'
removeM38PsudosomalRegion <- function(obj, withChr="FALSE") {
  annotations <- obj$alleleAnnot
  regions_x <- matrix(c(169969756, 170931299), nrow=1)
  regions_y <- matrix(c(90745845, 91644698), nrow=1)

  if (!("chr" %in% colnames(annotations)))
    stop("chr column not found in annotations")

  if (!("pos" %in% colnames(annotations)))
    stop("pos column not found in annotations")

  X <- ifelse(withChr, "chrX", "X")
  Y <- ifelse(withChr, "chrY", "Y")

  removeMeX <- annotations$chr == X &
    (annotations$pos >= regions_x[1, 1] )

  removeMeY <- annotations$chr == Y &
    (annotations$pos >= regions_y[1, 1] )

  removeMe <- ! (removeMeX | removeMeY)
  list(alleleCounts=obj$alleleCounts[removeMe,],
       alleleAnnot=obj$alleleAnnot[removeMe,],
       positionUID=obj$positionUID[removeMe])

}

#' Remove Human pseudosomal regions
#'
#' @export
#'
removeHg38PsudosomalRegion <- function(obj, withChr="FALSE") {
  annotations <- obj$alleleAnnot
  hg38_regions_x <- rbind(par1=c(10000,2781479), par2=c(155701383,156030895))
  hg38_regions_y <- rbind(par1=c(10000,2781479), par2=c(56887903,57217415))

  if (!("chr" %in% colnames(annotations)))
    stop("chr column not found in annotations")

  if (!("pos" %in% colnames(annotations)))
    stop("pos column not found in annotations")

  X <- ifelse(withChr, "chrX", "X")
  Y <- ifelse(withChr, "chrY", "Y")
  removeMeX <- annotations$chr == X &
    ( annotations$pos <= hg38_regions_x["par1",2] | annotations$pos >= hg38_regions_x["par2", 1] )

  removeMeY <- annotations$chr == Y &
    ( annotations$pos <= hg38_regions_y["par1",2] | annotations$pos >= hg38_regions_y["par2", 1] )

  removeMe <- ! (removeMeX | removeMeY)
  annotations[removeMe, ]
  list(alleleCounts=obj$alleleCounts[removeMe,],
       alleleAnnot=obj$alleleAnnot[removeMe,],
       positionUID=obj$positionUID[removeMe])

}
