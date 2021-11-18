blomatrixCOPdec <-
function(cop=NULL, para=NULL, as.sample=FALSE, ...) {
  if(as.sample) {
    if(is.null(para)) {
      warning("Sample Blomqvist's Beta desired but para is NULL, returning NULL")
      return(NULL)
    }
    if(length(names(para)) != 2) {
      warning("para argument must be data.frame having only two columns, returning NULL")
      return(NULL)
    }
    u <- para[,1]; v <- para[,2]; n <- length(u)
    #A <- (1+n)/2
    #samBLOM <- (2/n)*(sum(as.numeric((rank(u) - A)*(rank(v) - A) >= 0))) - 1
    return("not yet implemented")
  } else {
    if(is.null(cop)) {
      warning("must have copula argument specified, returning NULL")
      return(NULL)
    }
    t <- c(0.10, 0.50, 0.90)
    A <- 1 / P(t, t[3])
    B <- 1 / P(t, t[2])
    C <- 1 / P(t, t[1])
    blom <- matrix(c(cop(t, t[3], para=para, ...) * A - 1,
                     cop(t, t[2], para=para, ...) * B - 1,
                     cop(t, t[1], para=para, ...) * C - 1), ncol=3)
    colnames(blom) <- paste0("U|V=", c("0.10", "0.50", "0.90"))
    rownames(blom) <- rev(colnames(blom))
    return(blom)
  }
}


blomatrixCOPiqr <-
function(cop=NULL, para=NULL, as.sample=FALSE, ...) {
  t <- c(0.25, 0.50, 0.75)
  if(as.sample) {
    if(is.null(para)) {
      warning("Sample Blomqvist's Beta desired but para is NULL, returning NULL")
      return(NULL)
    }
    if(length(names(para)) != 2) {
      warning("para argument must be data.frame having only two columns, returning NULL")
      return(NULL)
    }
    samBLOM <- matrix(EMPIRcop(t, t[3], para=para, ...),
                      EMPIRcop(t, t[2], para=para, ...),
                      EMPIRcop(t, t[1], para=para, ...),
    return("not yet implemented")
  } else {
    if(is.null(cop)) {
      warning("must have copula argument specified, returning NULL")
      return(NULL)
    }
    #A <- c(16/3, 16/6, 16/9)
    #B <- c(16/2, 16/4, 16/6)
    #C <- c(16/1, 16/2, 16/3)
    A <- 1 / P(t, t[3])
    B <- 1 / P(t, t[2])
    C <- 1 / P(t, t[1])
    blom <- matrix(c(cop(t, t[3], para=para, ...) * A - 1,
                     cop(t, t[2], para=para, ...) * B - 1,
                     cop(t, t[1], para=para, ...) * C - 1), ncol=3)
    colnames(blom) <- paste0("U|V=", c("0.25", "0.50", "0.75"))
    rownames(blom) <- rev(colnames(blom))
    return(blom)
  }
}
