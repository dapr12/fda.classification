#' Medianfd - Computes the functional spatial median of the functional data 
#' @param time An object of the FDA Class.
#' @param xmat  n object of the FDA Class.
#' @return A component list with: 
#'\itemize{
#' \item{"median"}{Spatial Median (1xm) vector}
#' \item{"weights"}{Weights (nx1) vector.} 
#' \item{"norms"}{Distances between each curve and the median (nx1) vector. Note: The median minimizes the sum(norms)}
#'}
#' @references  Gervini, D. (2008). Robust functional estimation using the spatial median and spherical principal components. Biometrika 95 587-600.
#' @examples \dontrun{
#' #Computing the Median
#' matplot(time,CVScores$YSmooth, type="l")
#' lines(time, Median$median, col="red",lwd = 4) 
#' }
#' @export
medianfd<-function( fdaobj )
{
 
  if (!(inherits(fdaobj, "FdaSmoothClass"))) 
    stop("Argument FD not a functional smooth data object.")
  
  xmat<- t(fdaobj$YSmooth)
  
  time<- fdaobj$argvals
  
  n<- dim(xmat)[1]
  
  m<- dim(xmat)[2]
  
  A <- 0.5 * ( xmat[, 1:m-1] %*% diag(time[2:m]-time[1:m-1]) %*% t(xmat[,1:m-1] )
               + xmat[, 2:m] %*% diag(time[2:m]-time[1:m-1]) %*% t(xmat[,2:m]) )
 
  w<-rep(1, n)/n
 
  norms<- sqrt( diag(A) + t(w) %*% A %*% w - 2 * A %*% w)
  
  f<- sum(norms)
  
  err<-1
 
  iter<- 0
  
  while( err > 1e-5 && iter < 50 )
  {
    iter <- iter +1
    
    f0 <- f
    
    if ( any(norms< .Machine$double.eps ))
      
    {i0<- which(norms<eps)
     
     w<-rep(0, length(n))
     
     w[i0]<- 1 /length(i0)}
    
    else {
      
      w<- 1/norms
      
      w<- w/sum(w) }
    
    norms<-sqrt(diag(A) + t(w) %*% A %*% w - 2 * A %*% w)
    
    f<-sum(norms)
    
    err<-abs(f/f0-1)
  }
  
  med<-t(w) %*% xmat
  
  list(median = as.vector(med),
       
       weights = w,
       
       norms = norms
       
  )
}