#' fdadensity - Compute the Functional Data Density Estimation base on the Harmonics 
#' @param fdaobj fdClass Object
#' @param pcaobj PcaClass Object
#' @param bandwithd Bandwith for Density Estimation
#' @return A vector with the estimate density values  
#' @references  Delaigle, A. and Hall, P. (2010). Defining probability density for a distribution of random functions. Ann. Statist., 38, 1171-1193.
#' @examples \dontrun{
#' fdvector<-fdensity( fdaobjMale, pcaobj, bandwith=800)
#' plot(fdvector, type="p")
#' }
#' @export
fdensity<-function ( fdaobj, pcaobj,  bandwith){
  
  if (!(inherits(fdaobj, "FdaClass"))) 
    stop("Argument FD  not a functional data object.")
  
  if (!(inherits(pcaobj, "PcaClass"))) 
    stop("Argument FD  not a functional data object.")
  
  y<-as.matrix(t(fdaobj$data))
  
  Scores<- pcaobj$Scores
  
  n<- dim(y)[1] #Row
  
  m<- dim(y)[2] #Col
  
  dist<-rep(0,n)
  
  dist<-distance(pcaobj$Scores[,1],pcaobj$Scores)
   
  denest<-rep(0,n)
    
  for(i in 1:m){    
    
    denest[i]<- kern(dist[i]/bandwith)
  }
  
  return(denest)

}

kern <-function(u)
  
{
  
  if((u>0) && (u<1))
    
    (kern<-(1-u^2))
  
  else  (kern<-0)
  
    return(kern)
  
}


distance <-function ( x, xk ){
  
  nk<-dim(xk)[1]
  
  m<-dim(xk)[2]
  
  dist<-numeric(nk) 
  
  dist<-rep(0,nk) 
  
  for( i in 1: nk) {
    
    ssqd <- 0
    
    for( j in 1:m){
      
      ssqd<-ssqd+(x[j]-xk[i,j])^2}
    
    dist[i]<-sqrt(ssqd)
  }
  
  return(dist)
  
}

  