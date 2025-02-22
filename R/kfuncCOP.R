"kmeasCOP" <-
function(z, cop=NULL, para=NULL, wrtV=FALSE, as.sample=FALSE, verbose=FALSE,
            subdivisions=100L, rel.tol=.Machine$double.eps^0.25, abs.tol=rel.tol, ...) {
   kfuncCOP(z, cop=cop, para=para, wrtV=wrtV, as.sample=as.sample, verbose=verbose, ...)
}

"kfuncCOP" <-
function(z, cop=NULL, para=NULL, wrtV=FALSE, as.sample=FALSE, verbose=FALSE,
            subdivisions=100L, rel.tol=.Machine$double.eps^0.25, abs.tol=rel.tol, ...) {
   as.sample <- as.character(as.sample)
   if(as.sample != "FALSE") {
      if(length(para[1,]) != 2) {
         warning("ambiguous UV provided by 'para', need only two columns")
         return(NA)
      }

      if(as.sample=="genest" | as.sample=="Genest" ) {
         n <- nrow(para); ns <- seq_len(n)
         R <- rank(para[,1]); S <- rank(para[,2])
         "VIN" <- function(i) sum(as.numeric(R <= R[i] & S <= S[i])) / n
         FKin <- sapply(z, function(t) {
                sum(sapply(ns, function(j) as.numeric(VIN(j) <= t) )) / n })
         return(FKin)
      } else if(as.sample == "charpentier" | as.sample == "Charpentier") {
        # Charpentier, Arthur, 2012, Kendall's function for copulas:
        # http://freakonometrics.hypotheses.org/1126  (12/09/2012)
        n  <- nrow(para);               ns <- seq_len(n)
        i  <- rep(ns, each=n); j <- rep(ns,n)
        S  <- (para[i,1] > para[j,1]) & (para[i,2] > para[j,2])
        Z  <- tapply(S, i, sum) / (n-1)
        zz <- data.frame(t=sort(Z), Kc=(ns)/n)
        return(zz)
      } else {
         n     <- nrow(para)
         FKin  <- sort(EMPIRcop(para[,1], para[,2], para=para, ...))
         Zin   <- ( rank(FKin) - 0.5 ) / n
         suppressWarnings( empkc <- approx(c(0, FKin, 1), c(0, Zin, 1), xout=z)$y )
         return(empkc)
      }

   }

    if(is.null(cop)) {
      warning("must have copula argument specified, returning NULL")
      return(NULL)
   }

   "KCfunc1" <- function(the.z) {
       if(the.z == 0) return(0)
       myint <- NULL
       try( myint <- integrate(function(u) {
                     cvt <- sapply(1:length(u), function(i) {
                             COPinv(u=u[i], t=the.z, cop=cop, para=para, ...)})
                     derCOP(u, cvt, cop=cop, para=para, ...)  }, the.z, 1,
                     subdivisions=subdivisions, rel.tol=rel.tol, abs.tol=abs.tol),
            silent=! verbose)
       if(is.null(myint)) {
          if(verbose) warning("error on integration encountered ",
                              "(some copula singularity?)")
          return(0)
       }
       KC <- the.z + myint$value
       if(KC >  1) KC <- 1 # assume rounding errors just breaking through
       if(KC <  0) KC <- 0 # again for rounding errors
       return(KC)
   }
   "KCfunc2" <- function(the.z) {
       if(the.z == 0) return(0)
       myint <- NULL
       try( myint <- integrate(function(v) {
                     cut <- sapply(1:length(v), function(i) {
                             COPinv2(v=v[i],t=the.z, cop=cop, para=para, ...)})
                     derCOP2(u=cut, v=v, cop=cop, para=para, ...)  }, the.z, 1,
                     subdivisions=subdivisions, rel.tol=rel.tol, abs.tol=abs.tol),
            silent=! verbose)
       if(is.null(myint)) {
          if(verbose) warning("error on integration encountered ",
                              "(some copula singularity?)")
          return(0)
       }
       KC <- the.z + myint$value
       if(KC >  1) KC <- 1 # assume rounding errors just breaking through
       if(KC <  0) KC <- 0 # again for rounding errors
       return(KC)
   }
   KCfunc <- KCfunc1
   if(wrtV) KCfunc <- KCfunc2
   kc <- sapply(1:length(z), function(i) { KCfunc(z[i]) } )
   return(kc)
}

