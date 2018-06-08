drawChromosome <- function(data) {

  require(biovizBase)
  require(GenomicRanges)
  # library(ggbio, lib.loc="../guess-LOI/rlib")

  hg38 <- getIdeogram(genome = "hg38")

  bar.df <- data.frame(x1 = data$pos, x2 = data$pos+1, y1=-5, y2=12, seqnames = "chrX")
  lbs.df <- data.frame(x1 = data$pos, x2 = data$pos+1, y1=-5, y2=12, ids = as.character(data$id))

  # source("../loi-analysis/guess-LOI/teVedo.R")
  p = myIdeogram(hg38, lbs.df =lbs.df, bar.df=bar.df, subchr="chrX", xlabel = TRUE, aspect.ratio=1/25)

  ggsave(plot=p, filename="chrx-mapped-genes-inHouse-female.pdf", path=".", width = 15, height = 5)



  p = ggbio::Ideogram(hg38, subchr="chrX", xlabel = TRUE, aspect.ratio=1/25)

  chrX <- hg38[seqnames(hg38) == "chrX"]
  p <- ggplot() + ggbio::layout_karyogram(chrX, cytoband = TRUE, geom = NULL)

}
