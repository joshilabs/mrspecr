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
