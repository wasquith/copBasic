"EMPIRgridder" <-
function(empgrid=NULL, ...) {

  if(is.null(empgrid)) {
    warning("The gridded empirical copula (say from EMPIRgrid) is NULL")
    return(NULL)
  }
  if(! is.list(empgrid)) {
    warning("The gridded empirical copula is expected as a list from EMPIRgrid")
    return(NULL)
  }

  deluv <- empgrid$deluv
  empcop <- empgrid$empcop
  rc <- dim(empcop)

  n <- rc[1]
  if(n != rc[2]) {
     warning("grid is not square!")
     return(NA)
  }
  if(deluv != 1/(n-1)) {
     warning("concerns over value of deluv, not congruent with matrix size")
  }

  the.deriv <- matrix(nrow=n, ncol=n)
  for(i in 1:n) {
     section <- empcop[,i]
     diff.section  <- diff(section)
     derivative    <- c(0, diff.section/deluv)
     the.deriv[,i] <- derivative
  }

  for(i in 2:n) {
     #denominator <- the.deriv[i,n]
     #print(c(i,the.deriv[i,n]))
     #print(c(i,denominator))
     the.deriv[i,] <- the.deriv[i,]/the.deriv[i,n]
     #print(c(i,sum(the.deriv[i,])))
     #if(is.nan(the.deriv[i,1])) stop(sum(the.deriv[i,]))
  }
  the.deriv[1,] <- rep(NA, n)

  for(i in 2:n) {
    if(length(the.deriv[! is.finite(the.deriv[i,])]) > 0) {
      #print(empcop[,i]) # this would have been the second
      warning("found nonfinite values on row=",i," in grid derivative")
      next
    }
  }

  attributes(the.deriv) <- list(dim=dim(empcop),
                                rownames=empgrid$u,
                                colnames=empgrid$v,
                                message="use the rows!, wrt U")

  return(the.deriv)
}

