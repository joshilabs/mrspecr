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

#' @export
metab_peak <- function(filepath,outfile=NULL)
{
  peaks<-matrix(ncol=2)
  
  ycoords <- lcmodel_extract_batch(filepath)
  for (i in 1:ncol(ycoords))
  {
    naa <- max(ycoords[,i])
    cre <- max(ycoords[60:190,i])
    peaks <- rbind(peaks,c(naa,cre))
    
  }
  colnames(peaks) <- c("NAA_peak","Cre_peak")
  
  return(peaks[2:nrow(peaks),])
}

#' @export
metab_concentrations <- function(filepaths,outfile=NULL)
{
  metabs<-lcmodel_extract_metab_batch(filepaths)
  concentrations<-cbind(metab_pull$NAA,metab_pull$Cr.PCr)
  colnames(concentrations)<-c("NAA_conc","Cre_conc")
  return(concentrations)
}

#' @export
correlation <- function(coordpath,metab_path)
{
  cor<-matrix(ncol=2)
  peaks<-metab_peak(coordpath)
  conc<-metab_concentrations(metab_path)
  for (i in 1:nrow(peaks))
  {
    cor <- rbind(cor,c(peaks[i,1]/peaks[i,2],conc[i,1]/conc[i,2]))
  }
  colnames(cor)<- c("peak_ratio", "conc_ratio")
  return(cor[2:nrow(cor),])
}

