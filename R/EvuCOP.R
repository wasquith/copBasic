"EvuCOP" <-
function(u=seq(0.01, 0.99, by=0.01), cop=NULL, para=NULL, asuv=FALSE, ...) {

   v <- sapply(u, function(t) {
          E <- NULL
          try(E <- integrate(function(k) { 1-derCOP(cop=cop, para=para, t, k, ...) },
                    lower=0, upper=1)$value, silent=TRUE)
          if(is.null(E)) return(NA)
          return(E)
        })
   v[v <= 0] <- 0
   v[v >= 1] <- 1
   if(asuv) return(data.frame(U=u, V=v))
   return(v)
}
