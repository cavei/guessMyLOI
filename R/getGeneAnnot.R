#' @export
#'
collectGeneAnnotation <- function(file, colnames=c("gene", "chr","pos","DMR-Origin")) {
  df <- read.table(file, header=F, sep="\t", stringsAsFactors = F, check.names = F)
  if (NCOL(df) == length(colnames)) {
    colnames(df) <- colnames
  } else {
    warning("Columns number and colnames length differs. No colnames added", call. = FALSE)
  }
  df <- df[!duplicated(df$gene), , drop=F]
  row.names(df) <- NULL
  df
}
