#' pcafd - Functional Principal components analysis aims to display types of variation across a sample of
#' functions. Principal components analysis is an exploratory data analysis that tends to be an early part of many projects.
#' @import fda
#' @param fdaobject A functional Data Object 
#' @param nharm the number of harmonics or principal components to compute
#' @return A list containg the following components: 
#'\itemize{
#'\item{Harmonics}{ A functional data object for the harmonics or eigenfunctions.}
#'\item{values}{ The complete set of eigenvalues.}
#'\item{scores}{ A matrix of scores on the principal components or harmonics.}
#'\item{varprop}{A vector giving the proportion of variance explained by each eigenfunction.}
#' \item{cumvarprop}{A vector giving the cumulative proportion of variance explained by each eigenfunction.}
#'}
#' @references  Ramsay, James O., and Silverman, Bernard W. (2002), Applied Functional Data Analysis, Springer, New York.
#' @examples \dontrun{
#' FdaMaleHarm <- pcafd( fdaobjMale, nharm=3, Plot= TRUE)
#' plot(FdaMaleHarm)
#' }
#' @export
pcafd<- function ( fdaobj, nharm, Plot=FALSE ) 
{
  require(fda)
  
  if (!(inherits(fdaobj, "FdaClass"))) 
    stop("Argument FD  not a functional data object.")
  
  data<- fdaobj$data
  
  argvals<- fdaobj$argvals
  
  rangevals<- fdaobj$rangeval
  
  Data2FD <- Data2fd( argvals, data)
  
  PCA <- pca.fd(Data2FD, nharm, centerfns =FALSE)
  
  Harmnonics<-PCA$harmonics$coefs
  
  Scores<- PCA$scores
  
  values<- PCA$values[1:nharm]
  
  varprop<- PCA$varprop
  
  cumvarprop<- cumsum(PCA$varprop)
  
  if (Plot == TRUE){ plot (PCA) }
  
  pcaFD<-list(Harmonics = Harmnonics, 
              Scores= PCA$scores, 
              Values= values, 
              VarProp=PCA$varprop,
              CumVarProp=cumvarprop)
  
  class(pcaFD) <- "PcaClass"
  return(pcaFD)
  
  
}
