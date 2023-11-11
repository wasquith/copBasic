"ORDSUMcop" <-
function(u,v, cop=W, para=list(para=NULL, part=c(0,1)), ...) {

  if(length(u) > 1 & length(v) > 1 & length(u) != length(v)) {
    warning("length u = ",length(u), " and length v = ",length(v))
    warning("longer object length is not a multiple of shorter object length, no recycling in M()")
    return(NA)
  }
  # The extra hassle of vectorization made here is to handle situations
  # in which nested integrals are used where uneven vectors can be passed
  if(length(u) == 1) {
     u <- rep(u, length(v))
  } else if(length(v) == 1) {
     v <- rep(v, length(u))
  }

  if(! exists("part", para)) {
    warning("must have the part element of the para list populated")
    return(NULL)
  }
  J <- para$part
  num_partitions <- length(J)
  if(! exists("para", para)) {
    Q <- lapply(seq_len(num_partitions), function(k) NULL)
  } else if(is.list(  para$para) & length(para$para) == 1) {
    Q <- lapply(seq_len(num_partitions), function(k) para$para[[1]])
  } else if(is.vector(para$para) & length(para$para) == 1) {
    Q <- lapply(seq_len(num_partitions), function(k) para$para)
  } else {
    Q <- para$para
  }
  if(! exists("cop", para)) {
    C <- lapply(seq_len(num_partitions), function(k) W)
  } else if(length(para$cop) == 1) {
    if(is.list(para$cop)) {
      C <- lapply(seq_len(num_partitions), function(k) para$cop[[1]])
    } else {
      C <- lapply(seq_len(num_partitions), function(k) para$cop)
    }
  } else {
    C <- para$cop
  }

  if(length(J) != length(Q) | length(J) != length(C)) {
    warning("malformed para argument")
    return(NULL)
  }

  # Nelsen(2006,p63)
  zz <- sapply(seq_len(length(u)), function(i) {
           for(k in seq_len(length(J))) {
             a <- J[[k]][1]; b <- J[[k]][2]; ba <- b-a
             if(u[i] < a | u[i] > b | v[i] < a | v[i] > b) next # outside [a,b]^2
             return(a+ba*C[[k]]((u[i]-a)/ba, (v[i]-a)/ba, para=Q[[k]], ...))
           }; return(min(u[i], v[i])) }) # M(u,v) otherwise
  return(zz)
}


para = list(cop=c(CLcop, M, PLcop, GHcop),
            para=list(4, NA, 0.1, c(3,4)),
            part=list(c(0,0.25), c(0.25,0.35), c(0.35,0.85), c(0.85,1)))
UV <- simCOP(n=100, cop=ORDSUMcop, para=para, ploton=FALSE)
plot(c(0,1), c(0,1), xlab="U, NONEXCEEDANCE PROBABILITY",
                     ylab="V, NONEXCEEDANCE PROBABILITY")
for(k in seq_len(length(para$part))) {
  a <- para$part[[k]][1]; b <- para$part[[k]][2]
  px <- c(a, b, b, a, a); py <- c(a,a,b,b,a)
  polygon(px, py, lty=2, lwd=0.8, col="lightgreen")
  text((a+b)/2, (a+b)/2, k, cex=3, col="blue")
}
points(UV, pch=21, cex=0.8, col=grey(0.1), bg="white")


para = list(cop=c(GHcop),
            para=list(c(2,3)),
            part=list(c(0,0.2), c(0.2,0.3), c(0.3,0.5), c(0.5,0.7), c(0.7,1)))
UV <- simCOP(n=200, cop=ORDSUMcop, para=para, ploton=FALSE)
plot(c(0,1), c(0,1), xlab="U, NONEXCEEDANCE PROBABILITY",
                     ylab="V, NONEXCEEDANCE PROBABILITY")
for(k in seq_len(length(para$part))) {
  a <- para$part[[k]][1]; b <- para$part[[k]][2]
  px <- c(a, b, b, a, a); py <- c(a,a,b,b,a)
  polygon(px, py, lty=2, lwd=0.8, col="lightgreen")
  text((a+b)/2, (a+b)/2, k, cex=3, col="blue")
}
points(UV, pch=21, cex=0.8, col=grey(0.1), bg="white")


para = list(cop=P,
            para=list(NULL),
            part=list(c(0,0.257), c(0.257,0.358), c(0.358,1)))
DI <- diagCOP(cop=ORDSUMcop, para=para, delt=0.001)
if(sum(DI$diagcop == DI$t) >= 1) {
  message("The ORDSUMcop() operation is an ordinal sum if there exists\n",
          "a t=(0,1) exists such that C(t,t)=t by Nelsen (2006, theorem 3.2.1).")
} else {
  message("The ORDSUMcop() appears to not be an ordinal sum.")
}
abline(0,1, col="red")

