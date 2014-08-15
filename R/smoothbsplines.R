#' Bsplines- Evaluates a set of B-spline basis functions, or a derivative of these functions, at a set of arguments
#' @import splines 
#' @import Matrix
#' @import graphics
#' @param fdaobj An object of the FDA Class.
#' @param norder The order of the B-spline basis functions. The order less one is the degree of the piece-wise polynomials that make up any B-spline function.
#' @param lambda A nonnegative real number specifying the amount of smoothing to be applied to the estimated functional parameter.
#' @param Lf either a nonnegative integer or a linear differential operator object
#' @return A component list with: 
#'\itemize{
#'\item{"argvals "}{Range of values t_{1}, ..., t_{n}}
#'\item{"fdSmooth"}{Smooth Data Matrix}
#'\item{"df"}{Degrees of freedom}
#'\item{"SSE"}{Sum of the Square for Error}
#'\item{"PenMat"}{A roughness penalty matrix}
#'\item{"GCV"}{The generalized cross validation score}
#'}
#' @references  Ramsay, James O., and Silverman, Bernard W. (2006), Functional Data Analysis, 2nd ed., Springer, New York.
#' @examples \dontrun{
#' norder<-6
#' lambda <- 0.01
#' Lf<-4
#' fdaobj<-fdaobjMale
#' BSplinesGrowth <-bsplines( fdaobj, norder, lambda, Lf )
#' matplot(BSplinesGrowth$fdSmooth$coefs, type="l",main="Smooth Data using B-Splines",xlab="Age (years)", ylab="Heigth (cm)")
#' }
#' @export
smoothbsplines <- function ( fdaobj, norder, lambda, Lf )
{
  
  require('fda')
  
  if (!(inherits(fdaobj, "FdaClass"))) 
    stop("Argument FD  not a functional data object.")
  
  argvals<- fdaobj$argvals
    
  rangevals<- fdaobj$rangeval
  
  data<-as.matrix((fdaobj$data))
  
  norder<- norder
  
  nbasis <- length(argvals) + norder - 2
  
  lambda<- lambda
  
  basis <- create.bspline.basis(rangevals, nbasis, norder, argvals)
  
  fdPar  <- fdPar(basis, Lf, lambda)
  
  fdSmooth <- smooth.basis(argvals, data, fdPar)
  
  result<- list( argvals = fdaobj$argvals,  fdSmooth= fdSmooth$fd, df = fdSmooth$df,
                 SSE= fdSmooth$SSE, PenMat= fdSmooth$penmat, GCV= fdSmooth$gcv) 
  
  return(result)
  
}

