"MOcop" <- function(u,v, para=NULL, ...) { # Marshall-Olkin copula
   if(length(u) == 1) {
     u <- rep(u, length(v))
   } else if (length(v) == 1) {
     v <- rep(v, length(u))
   }
   alpha <- para[1]; beta <- para[2]
   f <- vector(mode="numeric", length(u))
   wnt <- u^alpha >= v^beta
   f[  wnt] <- u[  wnt]^(1-alpha) * v[  wnt]
   f[! wnt] <- u[! wnt]           * v[! wnt]^(1-beta)
   return(f)
}
#
# Dobrowolski, E., and Kumar, P., 2014, Some properties of the Marshall--Olkin and generalized Cuadras--Augé familes of copula: Australian Journal of Mathematical Analysis and Applications: v. 11, no. 1, art. 2, pp. 1--13, accessed on August 10, 2025, at https://ajmaa.org/searchroot/files/pdf/v11n1/v11i1p2.pdf.
