#' simulatefda - Performs simulation of Functional Data 
#' @import MASS
#' @param nsamples{ Number of Sample Curves}
#' @param ndrawn{ Number of Drawn in each sample}
#' @param rangeval{ Range of simulation} 
#' @param mean{ Constant Mean}
#' @param sigma{ Variance }
#' @return A component list with: 
#'\itemize{
#'\item{"RangeTime"}{ Simulation Range of Time }
#'\item{"Simulation"}{ Simulate Curves}
#'}
#' @references  Ramsay, James O., and Silverman, Bernard W. (2006), Functional Data Analysis, 2nd ed., Springer, New York.
#' @examples \dontrun{
#' simulation1<- simulatefda( nsamples=10, ndrawn=100, rangeval= c(-10,10), mean=0, sigma=1 )
#' simulation2<- simulatefda( nsamples=50, ndrawn=100, rangeval= c(-5,5), mean=3, sigma=1 )
#' n<-150
#' hues = seq(15, 375, length=n+1)
#' color<-hcl(h=hues, l=65, c=100)[1:n]
#' par(mfrow = c(2,1))
#' matplot(simulation1$rangetime,simulation1$simulation,type="l",col =color,xlab="",ylab="x(t)",main="Observations",lty=1)
#' matplot(simulation2$rangetime,simulation2$simulation,type="l",col =color,xlab="",ylab="x(t)",main="Observations",lty=1)
#' }
#' @export
simulatefda<-function ( nsamples, ndrawn, rangeval , mean , sigma )
  {
  
  samples<- nsamples
  
  nval<-4 ## Need to change
  
  sn<- sigma
  
  mn<- mean
  
  valmax<-rangeval[2]
  
  valmin<-rangeval[1]
    
  x.star<-seq(valmin,valmax,l=ndrawn)
  
  y<-rnorm(nval,mn,0)
  
  f <-data.frame(x.star,y)
  
  k.xx <- Sigma(f$x,f$x)
  
  k.xxs <- Sigma(f$x,x.star)
  
  k.xsx <- Sigma(x.star,f$x)
  
  k.xsxs <- Sigma(x.star,x.star)
  
  f.bar.star <-k.xsx%*%solve(k.xx+sn^2*diag(1,ncol(k.xx)))%*%f$y
  
  cov.f.star <-k.xsxs-k.xsx%*%solve(k.xx+sn^2*diag(1,ncol(k.xx)))%*%k.xxs
  
  values<-matrix(rep(0,length(x.star)*samples),ncol=samples)
  
  for (i in 1:samples){ 
        
    values[,i]<-mvrnorm(1,f.bar.star,cov.f.star)
    
  }
  
  list( rangetime= x.star, simulation = values  ) 
  
}


Sigma<-function( x, y ){
  
  Sigma<-matrix(rep(0,length(x)*length(y)),nrow=length(x))
  
  for(i in 1:nrow(Sigma)){
    
    for (j in 1:ncol(Sigma)){ 
      
      Sigma[i,j]<-exp(-1/2*(abs(x[i]-y[j]))^2)
    }
    
  }
  
  return(Sigma)
}
