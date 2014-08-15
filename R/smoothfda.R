#' Smoothfda - Computes the functional spatial median of the functional data 
#' @import KernSmooth
#' @param fdaobj An object of the FDA Class.
#' @param bandwitdh  Bandwidth 
#' @param degree  Degree of the polynomial 
#' @return A component list with: 
#'\itemize{
#'\item{"Data"}{}
#'\item{"argvals"}{}
#'\item{"rangevals"}{}
#'\item{YSmooth}{}
#'\item{"CVScore"}{}
#'\item{Mean}{} 
#'}
#' @references  Ramsay, James O., and Silverman, Bernard W. (2002), Applied Functional Data Analysis, Springer,New York.
#' @examples \dontrun{
#' #Computing the Median
#' matplot(time,CVScores$YSmooth, type="l")
#' lines(time, Median$median, col="red",lwd = 4) 
#' }
#' @export
smoothfda <- function (fdaobj, bandwidth, degree = 1) 
{
  
  if (!(inherits(fdaobj, "FdaClass"))) 
    stop("Argument FD  not a functional data object.")
  
  x<- fdaobj$argvals
  
  rangevals<- fdaobj$rangeval
  
  y<-as.matrix(t(fdaobj$data))
  
  spacing <- diff(x)
  
  if ( bandwidth < 0) 
    stop("'bandwith' must be positive")
  
  if (any(spacing < 0)) 
    stop("'x' must be increasing")
  
  if (nrow(y) < 2) 
    stop("'y' must have at least two rows")
  
  #if (length(x) != ncol(y)) 
  #  stop("length(x) and ncol(y) must be equal")
  
  n <- nrow(y)
  
  N <- ncol(y)
  
  y.hat <- apply(y, 1, function(z) locpoly(x = x, y = z, bandwidth = bandwidth, 
                                           gridsize = N, degree = degree)$y)
  
  mu.hat <- rowMeans(y.hat)
  
  residuals <- (n/(n - 1)) * mu.hat - y.hat/(n - 1) - t(y)
  
  PSEMean<- apply(residuals^2, 2, mean, trim = .2)
  
  PSEMedian<- apply(residuals^2, 2, median)
  
  EmpiricalMean<-apply(y.hat, 1, mean)
  
  result<- list( Data = fdaobj$data, argvals = fdaobj$argvals, 
                 rangevals = fdaobj$rangeval , YSmooth= y.hat,
                 CVScore=mean(residuals^2),  Mean=EmpiricalMean ) 
  
  class(result) <- "FdaSmoothClass"
  
  return(result)
  
}