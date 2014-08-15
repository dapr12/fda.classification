#' CVSelect - Select the Cross-Validation Bandwith described in (Foster, and ) for the Median of the PSE funcion based on Functional Data
#' @param fdaobj  FddaObject  
#' @param degree Degree of the local Polynomial to be used. If Degree is missing takes by default degree = 1.
#' @param Interva Interval to search the smoothing parameter 
#' @return A bandwith that minimizes the Median of the Median PSE. 
#' @references  Peter Foster PhD Thesis. University of Manchester 
#' @examples \dontrun{
#' Int<-c(1,2)
#' h<- hselect(NoxNoOutliers, 1, Int)
#' SmoothNox<- smoothfda(fdaNox, bandwidth= h, degree=1)
#' matplot(SmoothNox$YSmooth, type="l",ylab="x(t)",main="Observations")
#' lines(SmoothNox$Mean,col="black",lwd = 4)  
#' }
#' @export
hselect <- function( fdaobj, degree, interval = NULL, ...) 
{     
  y<-fdaobj$data
  y<-t(y)
  x<-fdaobj$argvals
  if (is.null(interval)) {
    rangex   <- diff(range(x))
    meshx    <- rangex / (length(x) - 1)
    interval <- c( ifelse(degree < 2, meshx / 2, meshx), rangex / 2)
  }
  
  optimize(medianPSE, interval, x = x, y = y, degree = degree, ...)$minimum 
}
