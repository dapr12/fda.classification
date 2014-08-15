#' densityScores Compute the Density Scores
#' @import ks
#' @param pcaobj PcaClass Object
#' @param bivariate If bivariate = TRUE compute a comtour plot base on the two Scores 
#' @return The component list is: 
#'\itemize{
#'  \item{Estimate}{ Density Estimate}
#'  \item{h}{ Scalar bandwithd }
#'}
#'
#'   
#' @references  Delaigle, A. and Hall, P. (2010). Defining probability density for a distribution of random functions. Ann. Statist., 38, 1171-1193.
#' @examples \dontrun{
#' fdvector<-fdensity( fdaobjMale, pcaobj, bandwith=800)
#' plot(fdvector, type="p")
#' }
#' @export
densityScores<-function ( pcaobj, bivariate= FALSE){
  
  if (!(inherits(pcaobj, "PcaClass"))) 
    stop("Argument FD  not a functional data object.")
  
  Scores<- pcaobj$Scores
  
  n<- dim(Scores)[1]
  
  m<- dim(Scores)[2]
  
  fhat<- kde( Scores[,1] , positive=TRUE)
  
  plot(fhat, title="")
    
    if( bivariate == TRUE)
  {
    plot( kde(Scores[,1:2], positive=TRUE), title="Density Estimation", xlab="1st Score", ylab="2nd Score")
  
  }
  
  return( list(Estimate= fhat$estimate, h=fhat$h))

  
}



