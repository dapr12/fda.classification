#' fderiv- Compute the numerical derivative for a Functional Data Object 
#' @param fdaobj fdClass Object
#' @param nderiv Number of derivative to be compute
#' @return A fda objects  
#' @references  Delaigle, A. and Hall, P. (2010). Defining probability density for a distribution of random functions. Ann. Statist., 38, 1171-1193.
#' @examples \dontrun{
#'  FWd<-fdaclass(nonworking$data)
#'  Fderiv<-fderiv(FWd,1)
#'  matplot(Fderiv$data, type="l")
#' }
#' @export
fderiv<-function(fdaobj,nderiv=1,...) 
  {
  
  if (!(inherits(fdaobj, "FdaClass"))) 
    stop("Argument FD  not a functional data object.")
  
  nas1<-is.na(fdaobj$data)
  
  if (any(nas1))  stop("Functional Data Object contains ",sum(nas1)," curves with some NA value \n")
  
  DATA<-fdaobj$data
  
  tt<- fdaobj$argvals 
  
  rtt<- fdaobj$rangevals
  
  ndist<-ncol(DATA)
  
  res<-matrix(NA,nrow=nrow(DATA),ncol=ncol(DATA))
  
  for (i in 1:nrow(DATA)) {
    
      a=diff(DATA[i,],differences=nderiv)/(tt[2:ndist]-tt[1:(ndist-1)])
      
      ab=matrix(NA,ncol=ndist,nrow=2)
      
      ab[1,2:ndist]=a
      
      ab[2,1:(ndist-1)]=a
      
      res[i,]=colMeans(ab,na.rm=TRUE)
    }
        
    res<-fdaclass(res,tt,rtt)
  
    return(res) 
}