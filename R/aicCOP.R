"aicCOP" <-
function(u,v, cop=NULL, para=NULL, ...) {
  if(is.null(cop)) {
    warning("must have copula argument specified, returning NULL")
    return(NULL)
  }
  theo.cop <-      cop(u,v, para=para,            ...)
  emp.cop  <- EMPIRcop(u,v, para=data.frame(u,v), ...)
  mse <- mean((theo.cop - emp.cop)^2)
  n <- length(emp.cop); m <- length(para)
  return(2*m + n*log(mse))
}
