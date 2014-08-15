#' @export
plot.FdaClass <- function(x, ...){
  
  matplot(x$argvals, x$data, type="l", main="Observations", xlab="Time", ylab="x(t)")
  
}
