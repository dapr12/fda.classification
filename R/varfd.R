#' Varfd - Computes the nearest positive definite variance - covariance matrix of functional data matrix
#' @param mat Typically a matrix or data frame which contains a set of curves stored in rows. Missing values are not accepte. 
#' @return A list containg the following components: 
#'\itemize{
#' \item{"PositDefMat"}{A Functional Variance Covariance Matrix }
#'}
#' @references  NJ Higham (1988, Linear Algebra Appl. 103:103-118).
#' @examples \dontrun{
#' Mat<- fdaobjMale$data
#' time<- fdaobjMale$argvals
#' CVScores<-cv.score( 0.2, time, t(Mat))
#' }
#' @export
varfd<-function(mat) 
{
  nmat<-var(mat)
  if (!is.matrix(nmat)) 
    nmat = as.matrix(nmat)
  d = dim(nmat)[1]
  if (dim(nmat)[2] != d) 
    stop("Input matrix is not square!")
  es = eigen(nmat, symmetric = TRUE)
  esv = es$values
  tol = d * max(abs(esv)) * .Machine$double.eps
  delta = 2 * tol
  tau = pmax(0, delta - esv)
  dm = es$vectors %*% diag(tau, d) %*% t(es$vectors)
  return( PositDefMat= (nmat + dm) )
}