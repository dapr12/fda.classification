#' fdaclass - Create an functional data object
#' @param mdata Typically a matrix of (n,m) dimension which contains a set of n curvesdiscretized in m points or argvals. Missing values are not accepted
#' @param argval Locations of the discritization points. Missing values are not accepted.
#' @param rangeval Range of discretization points, by default range(argevals). Missing values are not accepted.
#' @return An object of class FDA. The component list is: 
#'\itemize{
#' \item{"data"}{The data}
#' \item{"argvals"}{Locations of the discretization points, T1=t1, ..., Tm=tm}
#' \item{"rangevals"}{Range of discretization points}
#'}
#' @examples \dontrun{
#' fdaobjMale<-fdaclass(Gwd$Male, Gwd$Age, c(1,18))
#' }
#' @export
fdaclass<- function(mdata, argval= NULL, rangeval = NULL) 
{
   
  
  if ( is.null(argval) && is.null(rangeval) )
  { 
    
    data<- mdata  
    
    m<-dim(data)[1] #Row
    
    n<- dim(data)[2] #Col
    
    argvals<- round(seq(1, m, length.out = m), 3)
    
    rangeval<- c(1,m)
    
    rownames(data) <- argvals
    
    colnames(data) <- paste("Obs",1:n, sep= " ")
    
    result<- list(data = data, 
         
         argvals = argvals, 
         
         rangevals = rangeval)
         
     
  }
  
  else {
    
    if( all(range(argval) == rangeval ) )
  { 
    
      
    data<- mdata  
    
    argvals<-argval
    
    m<-dim(data)[1] #Row
    
    n<- dim(data)[2] #Col
    
    argvalmax<-rangeval[2]
    
    argvalmin<-rangeval[1]
    
    rangeval<- range(argvals)
    
    rownames(data) <- argvals
    
    colnames(data) <- paste("Obs",1:n, sep= " ")
    
    result<- list(data = data, 
         
         argvals = argvals, 
         
         rangevals = rangeval
         
     )
    
    }
  
  else  {stop("Range should be equal")}
  
  }
  
  class(result) <- "FdaClass"
  return(result)
  
 }  
  



