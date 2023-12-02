"EuvCOP" <-
function(v=seq(0.01, 0.99, by=0.01), cop=NULL, para=NULL, asuv=FALSE, ...) {

   u <- sapply(v, function(t) {
          E <- NULL
          try(E <- integrate(function(k) { 1-derCOP2(cop=cop, para=para, k, t, ...) },
                    lower=0, upper=1)$value, silent=TRUE)
          if(is.null(E)) return(NA)
          return(E)
        })
   u[u <= 0] <- 0
   u[u >= 1] <- 1
   if(asuv) return(data.frame(U=u, V=v))
   return(u)
}
