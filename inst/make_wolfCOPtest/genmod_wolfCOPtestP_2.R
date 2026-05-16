"PadeApproximant" <- # https://en.wikipedia.org/wiki/Pade_approximant
function(x, a=0, b=0) {
   # Defined for m >= 0 and n >= 1. Hence, there is an offset on the indexing
   # when we need to access A coefficients for the j indexing variable.
   j <- seq_len(length(a))-1 # number of A coefficients, but starting at m=0
   k <- seq_len(length(b))   # number of B coefficients with starting at n=1
   R <- vector(mode="numeric", length(x)) # The response R(x)
   for(i in seq_len(length(x))) { # for each of the values in x
     # nj <-     sapply(j, function(j) a[j+1]*x[i]^j)         # diagnostics
     # dk <-     sapply(k, function(k) b[k  ]*x[i]^k)         # diagnostics
     # print(nj); print(dk); nj <- sum(nj); dk <- 1 + sum(dk) # diagnostics
     nj <-     sum(sapply(j, function(j) a[j+1]*x[i]^j)) # j=0 to m
     dk <- 1 + sum(sapply(k, function(k) b[k  ]*x[i]^k)) # k=1 to n
     R[i] <- nj / dk
   }
   return(R)
}

Aenv <- new.env()
Benv <- new.env()

use_whole_sample  <- FALSE
use_weights <- FALSE
log10_increment   <- 0.02
max_nm_complexity <- 6


xo <- log10( Z$n )
wo <- length(Z$nsim)*sqrt(Z$nsim)/sum(sqrt(Z$nsim))

if(use_whole_sample) {
  x <- xo
  w <- wo
} else {
  x  <- seq(min(xo), max(xo), by=log10_increment)
  w  <- approx(xo, wo, xout=x, rule=2)$y
  w  <- length(x)*sqrt(w)/sum(sqrt(w))
}
if(! use_weights) w <- rep(1, length(w))


"ofunc" <- # objective function
function(par, m=NA, n=NA, x=NA, y=NA, w=1) {
  mix <- seq_len(m+1); nix <- m+1 + seq_len(n) # vectors of indices for A then B
  # print(par); print(c(m, n))   # diagnostics
  # print(c(mix, nix))           # diagnostics
  # a <- par[mix]; b <- par[nix] # diagnostics
  # print(a); print(b); stop()   # diagnostics
  yp  <- PadeApproximant(x, a=par[mix], b=par[nix])
  #plot(x, yp-y)
  #uv  <- data.frame(U=lmomco::pp(x,    sort=FALSE),
  #                  V=lmomco::pp(w*(y - yp)^2, sort=FALSE))
  #plot(uv)
  #print(nrow(uv))
  #print(wolfCOP(para=uv, as.sample=TRUE, nlarge=1, usefastgrid=TRUE), 16)
  #wlf <- wolfCOPtest2(uv[,1], uv[,2], asuv=TRUE, aslist=FALSE)[4]
  #print(wlf)
  #names(wlf) <- NULL

  #return(-wlf)
  err <- sum(w*(y - yp)^2, na.rm=TRUE) # accumulate a square error
  return(err)
}

strs <- c("logitmu", "logitlam2", "logittau3", "logittau4")
for(str in strs) {
  message("-------------------------------------------------------------")
  yo <- Z[,str] # response variable
  if(use_whole_sample) {
    y <- yo # revert to the whole sample
  } else {
    y  <- approx(xo, yo, xout=x, rule=2)$y # a bit suboptimal to interpolate but this way we get
    # equal parts of the sample size in logarithms spread evenly with some speed increase to study and the
    # interpolation of the weights already made
  }
  ERR <- +Inf; KEY <- ""; A <- B <- NA # characteristics of the optimal solution for str
  for(m in 0:max_nm_complexity) {
    for(n in 1:max_nm_complexity) {
      key <- paste0("m=", m, ":n=", n); message(key, appendLF=FALSE)
      for(o in 1:3) { # subloop to grab different initial conditions
        ai <- runif(m+1, min=-1, max=+1) # stress the system with random starting points
        bi <- runif(n,   min=-1, max=+1) # stress the system with random starting points
        init.par <- c(ai, bi) # initial vector of m+1, n concatenation of A,B coefficients
        rt <- optim(init.par, ofunc, m=m, n=n, x=x, y=y, w=w) # Nelder-Mead optimization
        a  <- rt$par[seq_len(m+1)]; b <- rt$par[m+1 + seq_len(n)] # optimized A and B coefficients
        yp <- PadeApproximant(x, a=a, b=b)   # make the final predictions
        err <- sum(w*(y - yp)^2, na.rm=TRUE) # accumulate a square error
        if(err < ERR) { # if we have found a new minimum error, let us preserve those results
          ERR <- err; KEY <- key; A <- a; B <- b
          assign(str, A, envir=Aenv); assign(str, B, envir=Benv)
          message("*", appendLF=FALSE)
          txt <- paste0(key, " err=", round(err, digits=6))
          plot( 10^x, y, log="x", pch=21, col=grey(0.7), bg="white",
                         xlab="Sample size", ylab=str)
          lines(10^x, yp); mtext(txt, line=1); #stop()
        }
      }
      message(ifelse(n != max_nm_complexity, ",", ""), appendLF=FALSE)
    }
    message("")
  }
}

pade_logitmu <- pade_logitlam2 <- pade_logittau3 <- pade_logittau4 <- NULL
Alst <- as.list(Aenv)
Blst <- as.list(Benv)
for(str in strs) {
  #a   <- get(str, envir=Aenv); b <- get(str, envir=Benv)
  a   <- get(str, Alst); b <- get(str, Blst)
  m   <- length(a)-1; n <- length(b)
  x   <- Z$n
  w   <- length(Z$nsim)*sqrt(Z$nsim)/sum(sqrt(Z$nsim))
  y   <- Z[,str]
  yp  <- PadeApproximant(log10(x), a=a, b=b)
  if(str == "logitmu"  ) pade_logitmu   <- yp
  if(str == "logitlam2") pade_logitlam2 <- yp
  if(str == "logittau3") pade_logittau3 <- yp
  if(str == "logittau4") pade_logittau4 <- yp
  err <- sum(w*(y - yp)^2)
  txt <- paste0("m=", m, " and n=", n, "\nerr=", round(err, digits=6))
  plot( x, y, pch=21, col=grey(0.7), bg="white", log="x", cex=w+.5,
        xlab="Sample size", ylab=str)
  lines(x, yp, lwd=3); mtext(txt, line=1)
}

# Now go and paste these into wolfCOPtest.R
txt <- paste0("myAlst <- list(\n",  paste(
 sapply(1:4, function(k) paste0("         ", names(Alst)[k],
                         paste0(" = c(", paste(get(names(Alst)[k], Alst), collapse=", ")), ")")),
                collapse=",\n"), ")")
message(txt)
eval(str2expression(txt))

txt <- paste0("myBlst <- list(\n",  paste(
 sapply(1:4, function(k) paste0("         ", names(Blst)[k],
                         paste0(" = c(", paste(get(names(Blst)[k], Blst), collapse=", ")), ")")),
                collapse=",\n"), ")")
message(txt)
eval(str2expression(txt))
