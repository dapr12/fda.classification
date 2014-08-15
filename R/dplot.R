#' dpout Graphical Bayesian selection of functional PCs.
#' @param lambda  Eigenvalues of covariance or correlation matrix. (Must be in decreasing order and strictly positive -- drop the zeros)
#' @return A list containg the following components: 
#'\itemize{
#'\item{fhat}{Function "F-hat"}
#'\item{ijump}{Indices of THMAX where d-function has a step.}
#'\item{bstck}{ Broken stick}
#'}
#' @references Auer, P. and Gervini, D. (2008). Choosing principal components: a new graphical method based on Bayesian model selection.  Communications in Statistics - Simulation and Computation 37 962-977.
#' @references NJ Higham (1988, Linear Algebra Appl. 103:103-118).
#' @examples \dontrun{
#' Var<-varfd((Smoothfda$YSmooth))
#' Eg<-eigen(Var)
#' Dplt<- dpout(Eg$val, plotting=TRUE)
#' }
#' @export 
dpout<-function(lambda, plotting=TRUE,...)
{

  lambda <- lambda[ lambda > 1e-5 ] 
    
  m<-length(lambda)
  
  d<-seq(0,m-1,1)
  
  fhat<-rep(0,m)
  
  for( i in 1:m)
    { 
        fhat[i]<-sum(log(lambda[i:m])) - (m-i+1)*log(mean(lambda[i:m]))
  }
  
  imaj<-chull(d, fhat)

  thmin<-rep(0,m)
  
  for( i in 2:m-1)
    {
      thmin[i]<-max((fhat[(i+1):m]- fhat[i])/(d[(i+1):m]-d[i]))
  }

  thmin[m]<-0

  thmax<-rep(0,m)

  thmax[1]<-Inf

  for( i in 2:m )
    {
    
    thmax[i]<-min((fhat[i]-fhat[1:(i-1)])/(d[i]-d[1:(i-1)]))
    
  }
  
  ijump<-which(thmin<=thmax)
  
  bstck<-rep(0,m)

  for(i in 1:m)
    {
    bstck[i]<-(1/m)*sum(1/(i:m))
          
  }

  bstck<-sum(lambda)*bstck
  
  if (plotting==TRUE) {  
    
   plot(thmax[ijump],d[ijump], main="D-Plot", xlab="Theta", ylab="d(Theta) - Harmonics", type="S")
   
          
  }   
  return( list( d=d, thmin=thmin,  thmax= thmax, fhat = fhat, ijump = ijump, bstck = bstck  ) ) 
  
}



