#'Performs functional discrimination of a sample of curves when a categorical response is observed (supervised classification). 
#'A local bandwidth (i.e. local number of neighbours) is selected by a cross-validation procedure with a quadratic kernel. 
#' @param Classes vector containing the categorical responses giving the group number for each 
# curve in the matrix CURVES (if nbclass is the number of groups, "Classes" contains numbers 1,2,...,nbclass)
#' @param  Curve matrix containing the curves dataset (row by row) used for the estimating stage
#' @param  Pred matrix containing new curves stored row by row used for computing predictions
#' @param ..." arguments needed for the call of the function computing the semi-metric between curves
#' @return A component list with: 
#'\itemize{
#' \item{"Estimated.classnumber""}{Vector containing estimated class membership for each curve}
#' \item{"Predicted.classnumber"}{if PRED different from CURVES, this vector contains predicted class membership for each curve of PRED} 
#' \item{"Misclas"}{misclassification rate computed from estimated values and observed values}
#'}
#' @references  Ferraty, F. and Vieu, P. (2006). Nonparametric functional data analysis. Springer Series in Statistics, New York.
#' @examples \dontrun{
#' tr <- sample(1:100, 50)
#' AllGwdprotein<-cbind(Gwdprotein$Group1,Gwdprotein$Group2)
#' train<- AllGwdprotein[,tr]
#' test<- AllGwdprotein[,-tr]
#' Classlearn <- sort(rep(1:2,25)) 
#' Classifmplsr <- classfd(Classlearn, train, test)
#' }
#' @export
classfd <- function(Classes, CURVES, PRED)
{
  kernel = "quadratic"

  semimetric = "mplsr"
  
  Classes <- as.vector(Classes)
  
  CURVES<-t(CURVES)
  
  PRED<- t(PRED)
  
  if(is.vector(PRED)) PRED <- as.matrix(t(PRED))
  
  testfordim <- sum(dim(CURVES)==dim(PRED))==2
  
  twodatasets <- T
  
  if(testfordim) twodatasets <- sum(CURVES==PRED)!=prod(dim(CURVES))
  
  sm <- get(paste("semimetric.", semimetric, sep = ""))
  
  SEMIMETRIC1 <- sm(Classes, CURVES, CURVES, 1)
  
  kernel <- get(kernel)
  
  n1 <- ncol(SEMIMETRIC1)
  
  step <- ceiling(n1/100)
  
  if(step == 0)  {step <- 1}
  
  Knearest <- seq(from = 10, to = n1 %/% 2, by = step)
  
  kmax <- max(Knearest)
  
  Classes.estimated <- 0
  
  Bandwidth.opt <- 0
  
  nbclass <- max(Classes)
  
  BINARY <- matrix(0, n1, nbclass)
  
  for(g in 1:nbclass)  { BINARY[, g] <- as.numeric(Classes == g) }
  
  HAT.PROB <- matrix(0, nrow = nbclass, ncol = length(Knearest))
  
  Knn1 <- 0
  
  for(i in 1:n1) {
    
    Norm.diff <- SEMIMETRIC1[, i]
    
    Norm.order <- order(Norm.diff)
    
    zz <- sort(Norm.diff)[2:(kmax + 2)]
    
    Bandwidth <- 0.5 * (zz[-1] + zz[ - (kmax + 1)])
    
    z <- zz[ - (kmax + 1)]
    
    ZMAT <- matrix(rep(z, kmax), nrow = kmax, byrow = T)
    
    UMAT <- ZMAT/Bandwidth
    
    KMAT <- kernel(UMAT)
    
    KMAT[col(KMAT) > row(KMAT)] <- 0
    
    Ind.curves <- Norm.order[2:(kmax + 1)]
    
    for(g in 1:nbclass) {
      
      Ind.resp <- BINARY[Ind.curves, g]
      
      YMAT <- matrix(rep(Ind.resp, kmax), nrow = kmax, byrow = T)
      
      HAT.PROB[g,  ] <- apply(YMAT[Knearest,  ] * KMAT[ Knearest,  ], 1, sum)
    }
    
    Kmatsumbyrow <- apply(KMAT[Knearest,  ], 1, sum)
    
    HAT.PROB <- HAT.PROB/matrix(Kmatsumbyrow,nrow(HAT.PROB), ncol(HAT.PROB), byrow=T)
    
    Criterium <- t(rep(1, nbclass)) %*% (HAT.PROB - BINARY[i,  ])^2
    
    index <- order(as.vector(Criterium))[1]
    
    Knn1[i] <- Knearest[index]
    
    Classes.estimated[i] <- order(HAT.PROB[, index])[nbclass]
    
    Bandwidth.opt[i] <- Bandwidth[index]
  }
  
  Misclas.estimated <- sum(Classes.estimated != Classes)/n1
  
  SEMIMETRIC2 <- sm(Classes, CURVES, PRED, 1)
  
  Bandwidth2 <- 0
  
  n2 <- ncol(SEMIMETRIC2)
  
    for(k in 1:n2) {
      
      Sm2k <- SEMIMETRIC2[, k]
      
      Sm2k.ord <- order(SEMIMETRIC2[, k])
      
      knn <- Knn1[Sm2k.ord[1]]
      
      Bandwidth2[k] <- sum(sort(Sm2k)[knn:(knn+1)])*0.5
      
    }
  
    KERNEL <- kernel(t(t(SEMIMETRIC2)/Bandwidth2))
  
    KERNEL[KERNEL < 0] <- 0
  
    KERNEL[KERNEL > 1] <- 0
  
    Denom <- apply(as.matrix(KERNEL), 2, sum)
  
    PROB.PREDICTED <- matrix(0, nrow = n2, ncol = nbclass)
  
    for(g in 1:nbclass) {
      
      PROBKERNEL <- KERNEL * BINARY[, g]
      
      PROB.PREDICTED[, g] <- apply(as.matrix(PROBKERNEL), 2, sum)/Denom
      
    }
  
    Classes.predicted <- as.vector((PROB.PREDICTED == apply(PROB.PREDICTED, 1, max)) %*% (1:nbclass))
  
    return(list(Predicted.classnumber = Classes.predicted, 
                Misclas = Misclas.estimated, 
                MissClassifTable=table( Classes.estimated,  Classes.predicted) )
           
    )
}

semimetric.mplsr <- function(Classes1, DATA1, DATA2, q)
{
  if(is.vector(DATA1)) DATA1 <- as.matrix(t(DATA1))
  if(is.vector(DATA2)) DATA2 <- as.matrix(t(DATA2))
  testfordim <- sum(dim(DATA1)==dim(DATA2))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(DATA1==DATA2)!=prod(dim(DATA1))
  qmax <- ncol(DATA1)
  if(q > qmax) stop(paste("give a integer q smaller than ", qmax))
  n1 <- nrow(DATA1)
  nbclass <- max(Classes1)
  BINARY1 <- matrix(0, nrow = n1, ncol = nbclass)
  for(g in 1:nbclass) {
    BINARY1[, g] <- as.numeric(Classes1 == g)
  }
  mplsr.res <- mplsr(DATA1, BINARY1, q)
  COMPONENT1 <- DATA1 %*% mplsr.res$COEF
  COMPONENT1 <- outer(rep(1, n1), as.vector(mplsr.res$b0)) + COMPONENT1
  if(twodatasets) {
    n2 <- nrow(DATA2)
    COMPONENT2 <- DATA2 %*% mplsr.res$COEF
    COMPONENT2 <- outer(rep(1, n2), as.vector(mplsr.res$b0)) + 
      COMPONENT2
  }
  else {
    COMPONENT2 <- COMPONENT1
  }
  SEMIMETRIC <- 0
  for(g in 1:nbclass)
    SEMIMETRIC <- SEMIMETRIC + outer(COMPONENT1[, g], COMPONENT2[, 
                                                                 g], "-")^2
  return(sqrt(SEMIMETRIC))
}

mplsr <- function(X, Y, K = 5)
{
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

quadratic <- function(u)
{
  1 - (u)^2
}