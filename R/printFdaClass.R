#' @export
print.FdaClass <- function(x, ...){
  
  cat("Data:\n")
  print(x$data)
  
  cat("\nValues:\n")
  print(x$argvals)
  
  cat("\nRange:\n")
  print(x$rangevals)
  
}