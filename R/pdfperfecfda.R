#'Performs Near Perfect Classification for functional data 
#' @import sampling
#' @param Classes vector containing the categorical responses giving the group number for each 
# curve in the matrix CURVES (if nbclass is the number of groups, "Classes" contains numbers 1,2,...,nbclass)
#' @param "Data" matrix containing the curves dataset (row by row) used for the estimating stage
#' @param "test" matrix containing the data test.
#' @param "ind0" vector containing the index of the class 0
#' @param "ind1" vector containing the index of the class 1
#' @param "indt" vector containing the index of the test class
#' @return A component list with: 
#'\itemize{
#' \item{"Hatpsir"}{Data projection onto a space of dimension 1}
#' \item{"Missclasification"}{misclassification rate computed from estimated values and observed values}
#' \item{"MeanClass0"}{Mean for data in class 0}
#' \item{"MeanClass1"}{Mean for data in class 1}
#'}
#' @references  Delaigle, A. and Hall, P. (2012). Achieving near-perfect classification for functional data. J. Roy. Statist. Soc. Ser. B, 74, 267-286
#' @examples \dontrun{
#' fd3<-pdfclasf(XXdata, XXtest, indClass0, indClass1, indtest )
#' plot(fd3$Hatpsir, type="l")
#' fd3$Missclassification
#' }
#' @export
#source("Functions.R")
pdfclasf<-function(Data, test, ind0, ind1, indt )
{
  
  XX<-Data
  
  XXtest<-test
  
  indClass0<- ind0
  
  indClass1<- ind1
  
  indtest<- indt 
  
  Deltat<-1
  
  nall<-nrow(XX)
  
  B2<-5
  
  B<-5
  
  K<-5
  
  varXX<-var(as.vector(XX))
  
  muX<-mean(as.vector(XX))
  
  XX<-(XX-muX)/sqrt(varXX)
  
  XXtest<-(XXtest-muX)/sqrt(varXX)
  
  indClass1=1:length(indClass1)
  
  indClass0=1:length(indClass0)+length(indClass1)
  
  barXClass1=XX[indClass1[1],]
  
  for(i in indClass1[2:length(indClass1)])
    barXClass1=barXClass1+XX[i,]

  barXClass0=XX[indClass0[1],]
  for(i in indClass0[2:length(indClass0)])
    barXClass0=barXClass0+XX[i,]


  barXClass1=barXClass1/length(indClass1)
  barXClass0=barXClass0/length(indClass0)
  
  MeanClass1<-barXClass1
  MeanClass0<-barXClass0
  
  ZZ=XX[indClass0,]-outer(rep(1, length(indClass0)), barXClass0)
  ZZ=rbind(ZZ,XX[indClass1,]-outer(rep(1, length(indClass1)), barXClass1))
  compo=prcomp(ZZ)

  phi=(compo$rotation)
  phi=phi/sqrt(Deltat)
  lambda=(compo$sdev)
  lambda=lambda^2*Deltat
  
  cumul=cumsum(lambda)/sum(lambda)
  npsir=min(which(cumul>1-1/ntrain^2))
  npsir=min(npsir,ntrain/2)
  
  DiffbarX=barXClass1-barXClass0
  hatmuj=DiffbarX%*%phi[,]*Deltat
  
  KKCV=round(ntrain/K)
  ind5foldCV=matrix(0,nrow=B2,ncol=KKCV)

  for(b in 1:B2)
  {
    s=runif(KKCV,0,ntrain)
    s=ceiling(s)
    ind5foldCV[b,]=s
  }
  
  for(b in 1:B2)
  {
    
    i=ind5foldCV[b,]
    XXCV=ZZ[-i,]
    compoCV=prcomp(XXCV)
    phiCV=(compo$rotation)
    eval(parse(text=paste("phiCV",b,"=phiCV/sqrt(Deltat)", sep = "")))
    lambdaCV=(compoCV$sdev)
    eval(parse(text=paste("lambdaCV",b,"=lambdaCV^2*Deltat", sep = "")))
    
  }

  CV=rep(0,npsir)
  nbCV=0
  npsirA=1

  hatpsir=0*phi[,1]
  for(j in 1:npsirA)
    hatpsir=hatpsir+hatmuj[j]/lambda[j]*phi[,j]
  
  BBB=length(indtest)
  CV1=0
  CV2=0
  CVPLS=0
  
  for(b in 1:BBB)
  {
    X=XXtest[b,]
    Xproj=X%*%hatpsir
    XXproj=XX%*%hatpsir
    TofX=(Xproj-mean(XXproj[indClass1]))^2-(Xproj-mean(XXproj[indClass0]))^2
    
    if((TofX>=0)&(is.element(indtest[b], indClass1All)))
      CV1=CV1+1
    if((TofX<=0)&(is.element(indtest[b], indClass0All)))
      CV1=CV1+1
    
  }

return( list(Missclassification= CV1/nall, Hatpsir = -hatpsir, 
             MeanClass1= MeanClass1, MeanClass0= MeanClass0) )

}



CVpsi=function(Ddata,npsir2)
{#This function calculates the classification error when projecting on the finction psir
  
  CV=0
  longClass1=length(indClass1)
  longClass0=length(indClass0)
  
  for(b in 1:B)
  {
    
    eval(parse(text=paste("phiCV=phiCV",b, sep = "")))
    eval(parse(text=paste("lambdaCV=lambdaCV",b,sep = "")))
    
    subClass1=0*barXClass1
    subClass0=0*barXClass1
    sublongClass1=longClass1
    sublongClass0=longClass0
    
    for(i in ind5foldCV[b,])
    {
      
      if(is.element(i, indClass1))
      {
        subClass1=subClass1+Ddata[i,]
        sublongClass1=sublongClass1-1
      }
      if(is.element(i, indClass0))
      {
        subClass0=subClass0+Ddata[i,]
        sublongClass0=sublongClass0-1
      }
      
    }				
    
    DiffbarXCV=(barXClass1*longClass1-subClass1)/sublongClass1-(barXClass0*longClass0-subClass0)/sublongClass0
    hatmuj2CV=DiffbarXCV%*%phiCV[,]*Deltat
    
    psir=hatmuj2CV[1]/lambdaCV[1]*phiCV[,1]
    if(npsir2>1)
    {
      for(k in 2:npsir2)
        psir=psir+hatmuj2CV[k]/lambdaCV[k]*phiCV[,k]
    }
    
    
    #project each data value on projecftion
    XXproj=(Ddata)%*%psir
    XprojCV=XXproj
    indClass1B=setdiff(indClass1,ind5foldCV[b,])
    indClass0B=setdiff(indClass0,ind5foldCV[b,])
    
    for(i in ind5foldCV[b,])
    {
      
      Xproj=XXproj[i,]
      #Classify in 1st population iff T>0
      TofX=(Xproj-mean(XprojCV[indClass1B]))^2-(Xproj-mean(XprojCV[indClass0B]))^2
      
      #misclassify
      if((TofX>=0)&(is.element(i, indClass1)))
        CV=CV+1
      if((TofX<=0)&(is.element(i, indClass0)))
        CV=CV+1
      
    }#loop i
  }#loop b
  
  return(CV)
}

CVmplsr=function(Ddata,Y,m)
{#This function calculates the CV classification error when projecting  on the function psir using partial least squares
  #Ddata is the data matrix of training curves
  
  
  CV=0
  
  for(b in 1:B)
  {
    
    projecftion=mplsr(Ddata[-ind5foldCV[b,],],Y[-ind5foldCV[b,]],m)$COEF
    
    #project each data value on projecftion
    XXproj=Ddata%*%projecftion
    XprojCV=XXproj
    
    indClass1B=setdiff(indClass1,ind5foldCV[b,])
    indClass0B=setdiff(indClass0,ind5foldCV[b,])
    
    
    for(i in ind5foldCV[b,])
    {
      
      Xproj=XXproj[i]
      #Classify in 1st population iff T>0
      TofX=(Xproj-mean(XprojCV[indClass1B]))^2-(Xproj-mean(XprojCV[indClass0B]))^2
      
      if((TofX>=0)&(is.element(i, indClass1)))
        CV=CV+1
      
      if((TofX<=0)&(is.element(i, indClass0)))
        CV=CV+1
    }
  }
  
  
  return(CV)
}


mplsr<- function(X, Y, K = 5)
{
  # Copyright (c) October 1993, Mike Denham.
  # Comments and Complaints to: snsdenhm@reading.ac.uk
  #
  # Orthogonal Scores Algorithm for PLS (Martens and Naes, pp. 121--123)
  #
  # X: predictors (matrix)
  #
  # Y: multivariate response (matrix)
  #
  # K: The number of PLS factors in the model which must be less than or
  #    equal to the  rank of X.
  #
  # Returned Value is the vector of PLS regression coefficients
  #
  tol <- 1e-10
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  dx <- dim(X)
  nbclass <- ncol(Y)
  xbar <- apply(X, 2, sum)/dx[1]
  ybar <- apply(Y, 2, sum)/dx[1]
  X0 <- X - outer(rep(1, dx[1]), xbar)
  Y0 <- Y - outer(rep(1, dx[1]), ybar)
  W <- matrix(0, dx[2], K)
  P <- matrix(0, dx[2], K)
  Q <- matrix(0, nbclass, K)
  sumofsquaresY <- apply(Y0^2, 2, sum)
  u <- Y0[, order(sumofsquaresY)[nbclass]]
  tee <- 0
  for(i in 1:K) {
    test <- 1 + tol
    while(test > tol) {
      w <- crossprod(X0, u)
      w <- w/sqrt(crossprod(w)[1])
      W[, i] <- w
      teenew <- X0 %*% w
      test <- sum((tee - teenew)^2)
      tee <- teenew
      cee <- crossprod(tee)[1]
      p <- crossprod(X0, (tee/cee))
      P[, i] <- p
      q <- crossprod(Y0, tee)[, 1]/cee
      u <- Y0 %*% q
      u <- u/crossprod(q)[1]
    }
    Q[, i] <- q
    X0 <- X0 - tee %*% t(p)
    Y0 <- Y0 - tee %*% t(q)
  }
  COEF <- W %*% solve(crossprod(P, W)) %*% t(Q)
  b0 <- ybar - t(COEF) %*% xbar
  list(b0 = b0, COEF = COEF)
}

