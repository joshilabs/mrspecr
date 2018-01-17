#' @export
batch_alignment <- function(filepath,showplot=F)
{
  batchmatrix <- lcmodel_extract_batch(filepath)
  
  #getting euclidean mean and aligned mean
  betamean <- rowMeans(batchmatrix)
  time<-data.matrix(lcmodel_extract_coord_x(readLines(filepath, n=1)))
  
  align <- fdasrvf::align_fPCA(batchmatrix,time,showplot=F)
  mfn <- fdasrvf::srsf_to_f(align$mqn, time)
  align$mfn <- mfn
  #plots
  if (showplot == T)
  {
    plot(time,mfn, type='l', col = "green")
    lines(time,betamean, type='l', col = "blue" )  
  }

  return (align)
}

