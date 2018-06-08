#' Read allelic expression file
#'
#' @param file a file name
#' @param annotCol the column number where the annotations ends
#'
#' @return a list with
#'   \item{alleleCounts}{each row is an allele (SNP) expression value}
#'   \item{alleleAnnot}{allele annotations}
#'   \item{positionUID}{uinque id for snp allele}
#' @rdname input_function
#' @export
#'
readGATK.AllelicImbalance <- function(file, annotCol=7) {
  data <- read.table(file=file, header=T, sep="\t", quote="\"",
                     check.names = FALSE,
                     comment.char = "",
                     stringsAsFactors = FALSE)

  colnames(data) <- checkNames(colnames(data))

  alleleAnnot <- data[, c(1:annotCol)]
  alleleCounts <- data[, -c(1:annotCol), drop=F]
  fi <- sapply(alleleCounts, is.factor)
  alleleCounts[fi] <- lapply(alleleCounts[fi], as.character)
  alleleCounts <- suppressWarnings(data.matrix(alleleCounts))

  if (!all(c("chr", "pos") %in% colnames(alleleAnnot)))
    stop(paste0("No columns chr or pos found in ", file))

  positionUID <- paste(alleleAnnot$chr,alleleAnnot$pos,sep="_")
  list(alleleCounts=alleleCounts, alleleAnnot=alleleAnnot, positionUID=positionUID)
}

#' Read allelic expression from VCF
#'
#' @param file.vcf a vcf file name
#'
#' @inheritParams readGATK.AllelicImbalance return
#'
#' @rdname input_function
#' @export
#'
readFromVCF <- function(file.vcf) {
  require(VariantAnnotation)

  vcf <- readVcf(file=file.vcf)
  ad <- readGeno(file=file.vcf, x="AD")
  gt <- readGT(file=file.vcf, nucleotide=T)

  vcfdf <- as.data.frame(rowRanges(vcf))
  vcfdf$ALT <- sapply(vcfdf$ALT, function(x) toString(x[[1]]))

  vcfdf$paramRangeID <- NULL
  vcfdf$strand <- NULL
  vcfdf$width <- NULL
  vcfdf$QUAL <- NULL
  vcfdf$FILTER <- NULL

  counts <- t(sapply(seq_len(NROW(ad)), function(i){
    sapply(ad[i, ], paste, collapse = "/")
  }))
  row.names(counts) <- row.names(vcf)

  keep <- apply(counts, 1, function(x) {!all(x=="NA/NA")})

  annoCounts <- cbind(ids=row.names(vcfdf)[keep], vcfdf[keep, , drop=F], counts[keep,, drop=FALSE])
  gt <- cbind(ids=row.names(gt)[keep], gt[keep, , drop=F])

  return(list(annoCounts=annoCounts, genotypes=gt))

}

checkNames <- function(nm) {
  nm <- gsub("-", "_", nm, fixed=T)
  nm <- gsub("*", ".", nm, fixed=T)
  nm <- gsub("/", ".", nm, fixed=T)
  gsub("+", "_", nm, fixed=T)
}
