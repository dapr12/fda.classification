#' CVSelect - Select the Cross-Validation Bandwith described in (Foster, and ) for the Median of the PSE funcion based on Functional Data
#' @import KernSmooth
#' @param bandwith 
#' @param x Location of the discretization points. THis discretization points must be uniform and missing values are not accepted.  
#' @param y Typically a matrix or data frame which contains a set of curves stored in rows. Missing values are not accepte. 
#' @param degree Degree of the local Polynomial to be used. If Degree is missing takes by default degree = 1.
#' @return A bandwith that minimizes the Median of the Median PSE for the functional data set. 
#' @references  Foster and Stehpen. PhD Thesis. Manchester University
#' @examples \dontrun{
#' Mat<- fdaobjMale$data
#' h<- cv.select(c(0,10), 1:31,t(Mat),1)
#' }
#' @export
medianPSE<- function (bandwidth, x, y, degree) 
{
  
  y<-as.matrix(y)
  
  spacing <- diff(x)
  
  if(bandwidth <0)
    stop("'bandwithd' must be positive")
  
  if (any(spacing < 0)) 
    stop("'x' must be increasing")
  
  if (nrow(y) < 2) 
    stop("'y' must have at least two rows")
  
  if (length(x) != ncol(y)) 
    stop("length(x) and ncol(y) must be equal")
  
  n <- nrow(y)
  
  N <- ncol(y)
  
  y.hat <- apply(y, 1, function(z) locpoly(x = x, y = z, bandwidth = bandwidth, 
                                           gridsize = N, degree = degree)$y)
  
  mu.hat <- rowMeans(y.hat)
  
  residuals <- (n/(n - 1)) * mu.hat - y.hat/(n - 1) - t(y)
  
  PSEMedian<- apply(residuals^2, 2, median)
  
  return(median(PSEMedian))
    
}

