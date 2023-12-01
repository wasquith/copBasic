"rCOP" <- function(n, cop=NULL, para=NULL, na.rm=TRUE, seed=NULL, resamv01=FALSE, ...) {
   if(! is.null(seed)) set.seed(seed)
   u <- runif(n)
   v <- simCOPmicro(u, cop=cop, para=para, ...)
   if(resamv01) {
     while(any(v <= 0)) {
       message("rCOP() has some v <= 0, resampling those")
       v[v <= 0] <- simCOPmicro(runif(length(v[v <= 0])), cop=cop, para=para, ...)
     }
     while(any(v >= 1)) {
       message("rCOP() has some v >= 1, resampling those")
       v[v >= 1] <- simCOPmicro(runif(length(v[v >= 1])), cop=cop, para=para, ...)
     }
   }
   data.frame(U=u, V=v)
}


"simCOP" <-
function(n=100, cop=NULL, para=NULL, na.rm=TRUE, seed=NULL, keept=FALSE,
                graphics=TRUE, ploton=TRUE, points=TRUE, snv=FALSE,
                infsnv.rm=TRUE, trapinfsnv=.Machine$double.eps,
                resamv01=FALSE, ...) {

  if(is.null(cop)) {
     warning("must have copula argument specified, returning NULL")
     return(NULL)
  }

  if(! graphics) {
     ploton <- FALSE
     points <- FALSE
  }

  if(! is.null(seed)) set.seed(seed)

  n <- as.integer(n)

  u <- runif(n); t <- runif(n)
  v <- sapply(seq_len(n), function(i) { derCOPinv(cop=cop, u[i], t[i], para=para, ...) })
  if(resamv01) {
    while(any(v <= 0)) {
      message("simCOP() has some v <= 0, resampling those")
      v[v <= 0] <- sapply(seq_len(length(v[v <= 0])),
                            function(i) { derCOPinv(cop=cop, u[i], t[i], para=para, ...) })
    }
    while(any(v >= 1)) {
      message("simCOP() has some v >= 1, resampling those")
      v[v >= 1] <- sapply(seq_len(length(v[v >= 1])),
                            function(i) { derCOPinv(cop=cop, u[i], t[i], para=para, ...) })
    }
  }

  dots <- list(...)
  ditches <- c("delu", "derdir", "trace")
  for(d in ditches) {
    if(d %in% names(dots)) dots <- dots[ - which(names(dots) == d)]
  }

  # Because z is a data.frame, it must be assigned within the ifelse()
  ifelse(keept, z <- data.frame(U=u, V=v, T=t), z <- data.frame(U=u, V=v))
  if(na.rm) {
     z <- z[complete.cases(z), ]
     m <- length(z[,1])
     if(m != n) {
        warning("user requested n=",n," simulations but only m=", m,
                " could be made without NA from derCOPinv (uniroot failure therein)")
        row.names(z) <- NULL # reset the rows to "1:m"
     }
  }
  if(snv) {
     if(infsnv.rm) {
        z <- z[z$U != 0, ]; z <- z[z$U != 0, ]
        z <- z[z$V != 1, ]; z <- z[z$V != 1, ]
        row.names(z) <- NULL
     } else if(trapinfsnv) {
        z$U[z$U == 0] <-   trapinfsnv; z$V[z$V == 0] <-   trapinfsnv
        z$U[z$U == 1] <- 1-trapinfsnv; z$V[z$V == 1] <- 1-trapinfsnv
     }
     z$U <- qnorm(z$U)
     z$V <- qnorm(z$V)
  }
  if(ploton) {
     if(snv) {
        plot(z$U, z$V, type="n",
             xlab="STANDARD NORMAL VARIATE OF U",
             ylab="STANDARD NORMAL VARIATE OF V")
     } else {
        plot(NA, NA, type="n", xlim=c(0,1), ylim=c(0,1),
             xlab="U, NONEXCEEDANCE PROBABILITY",
             ylab="V, NONEXCEEDANCE PROBABILITY")
     }
  }
  dots$x <- z$U; dots$y <- z$V
  if(points & ! is.null(dev.list())) do.call("points", dots)

  return(z)
}
