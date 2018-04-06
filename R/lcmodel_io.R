#' @export
#extracts X and Y
lcmodel_extract_coord <- function(coordfile,outfile=NULL)
{
  c<-1
  coordtext=list(0)
  for (i in readLines(coordfile))
  {
    coordtext[c]<-trimws(i)
    c<-c+1
    #reading table line by line and remove leading and trailing white spaces
  }
  for (s in 1:length(coordtext)) 
  {
    if (coordtext[s] == "480 points on ppm-axis = NY")
    {
      idx1<-s
      #set idx1 to the point right before data pts begin
    }
    if (coordtext[s] == "NY phased data points follow")
    {
      idx2<-s
      #set idx2 to the point right after data pts end
    }
  }
  ppmlist<-""
  for (i in coordtext[(idx1+1):(idx2-1)])
  {
    ppmlist<-paste(ppmlist,i," ")
    #causing the first item in ppmlist to be "" when it should just be 4.0000
  }
  splitlist<-strsplit(ppmlist," +")
  ppmx<-as.double(splitlist[[1]][2:length(splitlist[[1]])])
  
  for (s in 1:length(coordtext)) 
  {
    if (coordtext[s] == "NY points of the fit to the data follow")
    {
      idx1<-s
      #set idx1 to the point right before data pts begin
    }
    if (coordtext[s] == "NY background values follow")
    {
      idx2<-s
      #set idx2 to the point right after data pts end
    }
  }
  yfitlist<-""
  for (i in coordtext[(idx1+1):(idx2-1)])
  {
    yfitlist<-paste(yfitlist,i," ")
    #causing the first item in yfitlist to be "" when it should just be 4.0000
  }
  ysplitlist<-strsplit(yfitlist," +")
  yfit<-as.double(ysplitlist[[1]][2:length(ysplitlist[[1]])])
  
  X <- data.frame(ppmx,yfit)
  if (!is.null(outfile))
  {
    utils::write.csv(X, outfile)
  }
  
  return (X)
}

#' @export
#extracts ppmx coords
lcmodel_extract_coord_x <- function(coordfile,outfile=NULL)
{
  c<-1
  coordtext=list(0)
  for (i in readLines(coordfile))
  {
    coordtext[c]<-trimws(i)
    c<-c+1
    #reading table line by line and remove leading and trailing white spaces
  }
  for (s in 1:length(coordtext)) 
  {
    if (coordtext[s] == "480 points on ppm-axis = NY")
    {
      idx1<-s
      #set idx1 to the point right before data pts begin
    }
    if (coordtext[s] == "NY phased data points follow")
    {
      idx2<-s
      #set idx2 to the point right after data pts end
    }
  }
  ppmlist<-""
  for (i in coordtext[(idx1+1):(idx2-1)])
  {
    ppmlist<-paste(ppmlist,i," ")
    #causing the first item in ppmlist to be "" when it should just be 4.0000
  }
  splitlist<-strsplit(ppmlist," +")
  ppmx<-as.double(splitlist[[1]][2:length(splitlist[[1]])])
  
  X <- data.frame(ppmx)
  if (!is.null(outfile))
  {
    utils::write.csv(X, outfile)
  }
  
  return (X)
}

#' @export
#extracts yfit coords
lcmodel_extract_coord_y <- function(coordfile,outfile=NULL)
{
  c<-1
  coordtext=list(0)
  for (i in readLines(coordfile))
  {
    coordtext[c]<-trimws(i)
    c<-c+1
    #reading table line by line and remove leading and trailing white spaces
  }
for (s in 1:length(coordtext)) 
{
  if (coordtext[s] == "NY points of the fit to the data follow")
  {
    idx1<-s
    #set idx1 to the point right before data pts begin
  }
  if (coordtext[s] == "NY background values follow")
  {
    idx2<-s
    #set idx2 to the point right after data pts end
  }
}
yfitlist<-""
for (i in coordtext[(idx1+1):(idx2-1)])
{
  yfitlist<-paste(yfitlist,i," ")
  #causing the first item in yfitlist to be "" when it should just be 4.0000
}
ysplitlist<-strsplit(yfitlist," +")
yfit<-as.double(ysplitlist[[1]][2:length(ysplitlist[[1]])])

X <- data.frame(yfit)
if (!is.null(outfile))
{
  utils::write.csv(X, outfile)
}

return (X)
}

#' @export
lcmodel_extract_batch <- function(coordfile,outfile=NULL)
{
coordslist=list()
#create empty list for dataframes to be pushed into
  
# for (x in readLines(coordfile))
# {
#   coordslist<-c(coordslist,lcmodel_extract_coord_y(x))
# }
coordslist <- lapply(readLines(coordfile),lcmodel_extract_coord_y)
coordmatrix <- matrix(unlist(coordslist), nrow = 480)
if (!is.null(outfile))
{
  utils::write.csv(coordmatrix, outfile)
}
return (coordmatrix)
}

#extracts metabolive values and S/N
#' @export
lcmodel_extract_metab <- function(lcmodel_outdir,outfile=NULL)
{
#customizing filepath to read to metabolite csv and sn table
csvfile <- file.path(lcmodel_outdir, 'spreadsheet.csv')
tblfile <- file.path(lcmodel_outdir, 'table')
  
#pull metabolites
met_levels<-read.csv(csvfile, nrows=2)
  
#remove whitespace from table
c<-1
tbl=list(0)
for (i in readLines(tblfile))
{
  tbl[c]<-trimws(i)
  c<- c+1
}
  
#pull S/N value
sn_val<-as.numeric(stringr::str_sub(tbl[45],-2,-1))
met_levels$SN<- sn_val
  
if (!is.null(outfile))
{
  utils::write.csv(met_levels, outfile)
}
return (met_levels)
}

#extracts batch met levels and S/N
#' @export
lcmodel_extract_metab_batch <- function(filepaths,outfile=NULL)
{
  
  metlist=list()
  metlist <- lapply(readLines(filepaths), lcmodel_extract_metab)
  #turn list of df into single df
  metlist<-plyr::ldply(metlist, data.frame)
  
  if (!is.null(outfile))
  {
    utils::write.csv(metlist, outfile)
  }
  return (metlist)
}

#create metabolite data matrix with subj_id and classification
#' @export
metab_matrix <- function(coordfile,sub_id,outfile=NULL)
{
  subj <- as.matrix(read.csv(sub_id,sep=",")) 
  mat <- as.matrix(lcmodel_extract_metab_batch(coordfile))
  metab_mat <- cbind(subj,mat[,-c(1,2)])
  
  if (!is.null(outfile))
  {
    utils::write.csv(metab_mat, outfile)
  }
  return (metab_mat)
}

#create feature fxn data matrix with subj_id and classification
#' @export
fxn_matrix <- function(coordfile,sub_id,outfile=NULL)
{
  subj <- as.matrix(read.csv(sub_id,sep=",")) 
  mat <- t(lcmodel_extract_batch(coordfile))
  fxn_mat <- cbind(subj,mat)
  
  if (!is.null(outfile))
  {
    utils::write.csv(fxn_mat, outfile)
  }
  return (fxn_mat)
}




