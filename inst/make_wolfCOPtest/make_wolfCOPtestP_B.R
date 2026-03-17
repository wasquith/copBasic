if(Sys.getenv("RSTUDIO") == "1") {
  # Automatic change directory to location of this script when RStudio is running.
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

library(copBasic)
library(lmomco  )
library(mgcv    )

# --------------------------------------------------------------------------------------------------
# Using some text files of the output of several runs (even if incomplete) of make_wolfCOPtestP_A.R
# with a different intermediate output file name, load some results.
#D <-          read.table("wolfPI.txt",   header=TRUE)              # from running make_wolfCOPtestP_A.R
#D <- rbind(D, read.table("wolfPIy.txt",  header=TRUE))             # from running make_wolfCOPtestP_A.R
#D <- rbind(D, read.table("wolfPIz.txt",  header=TRUE))             # from running make_wolfCOPtestP_A.R

D <- NULL
D <- rbind(D, read.table("montecarlo_wolfPI_0.txt",  header=TRUE)) # from running make_wolfCOPtestP_A.R
D <- rbind(D, read.table("montecarlo_wolfPI_1.txt",  header=TRUE)) # from running make_wolfCOPtestP_A.R
D <- rbind(D, read.table("montecarlo_wolfPI_2.txt",  header=TRUE)) # from running make_wolfCOPtestP_A.R
D <- rbind(D, read.table("montecarlo_wolfPI_3.txt",  header=TRUE)) # from running make_wolfCOPtestP_A.R
D <- rbind(D, read.table("montecarlo_wolfPI_4.txt",  header=TRUE)) # from running make_wolfCOPtestP_A.R
D <- rbind(D, read.table("montecarlo_wolfPI_5.txt",  header=TRUE)) # from running make_wolfCOPtestP_A.R
D <- rbind(D, read.table("montecarlo_wolfPI_6.txt",  header=TRUE)) # from running make_wolfCOPtestP_A.R
D <- rbind(D, read.table("montecarlo_wolfPI_7.txt",  header=TRUE)) # from running make_wolfCOPtestP_A.R
D <- rbind(D, read.table("montecarlo_wolfPI_8.txt",  header=TRUE)) # from running make_wolfCOPtestP_A.R
D <- rbind(D, read.table("montecarlo_wolfPI_9.txt",  header=TRUE)) # from running make_wolfCOPtestP_A.R
D <- rbind(D, read.table("montecarlo_wolfPI_10.txt", header=TRUE)) # from running make_wolfCOPtestP_A.R


# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
Z <- NULL # now compute weighted mean columns
for(k in sort(unique(D$n))) {
  y <- D[D$n == k,]
  z <- y[1,]
  for(i in 3:ncol(y)) {
    z[,i] <- weighted.mean(y[,i], y$nsim)
  }
  z[,1] <- sum(y$nsim)
  Z <- rbind(Z, z)
}
Z <- Z[is.finite(Z$logitmu),]
Z$wgts <- sqrt(Z$nsim); Z$wgts <- Z$wgts/sum(Z$wgts) * length(Z$wgts)
# --------------------------------------------------------------------------------------------------

plotlmrdia(lmrdia(), xlim=c(0.05, 0.3), ylim=c(0.1, 0.20), empty=TRUE)
points(Z$tau3,      Z$tau4,       pch=21, col="turquoise4", bg="turquoise1", lwd=0.6, cex=sqrt(Z$wgts))
points(Z$logittau3, Z$logittau4,  pch=21, col="salmon4",    bg="salmon1",    lwd=0.6, cex=sqrt(Z$wgts))
plotlmrdia(lmrdia(), add=TRUE, nopoints=TRUE, autolegend=TRUE, xleg="topleft")

# --------------------------------------------------------------------------------------------------
OPTIM.TYPE <- "press"
INIT.PWR <- -0.01 # initial parameter for the optimization scheme, being on the right side zero matters
BYLOGIT <- TRUE # run this script with both TRUE and FALSE and harvest the coefficients of the
# regressions and then emplace them in the wolfCOPtest function.


# --------------------------------------------------------------------------------------------------
#    88888888ba   88888888888  ,ad8888ba,   88888888ba   88888888888  ad88888ba    ad88888ba
#    88      "8b  88          d8"'    `"8b  88      "8b  88          d8"     "8b  d8"     "8b
#    88      ,8P  88         d8'            88      ,8P  88          Y8,          Y8,
#    88aaaaaa8P'  88aaaaa    88             88aaaaaa8P'  88aaaaa     `Y8aaaaa,    `Y8aaaaa,
#    88""""88'    88"""""    88      88888  88""""88'    88"""""       `"""""8b,    `"""""8b,
#    88    `8b    88         Y8,        88  88    `8b    88                  `8b          `8b
#    88     `8b   88          Y8a.    .a88  88     `8b   88          Y8a     a8P  Y8a     a8P
#    88      `8b  88888888888  `"Y88888P"   88      `8b  88888888888  "Y88888P"    "Y88888P"
# --------------------------------------------------------------------------------------------------
n  <- log10(Z$n)
if(BYLOGIT) {
  mu <- Z$logitmu; l2 <- Z$logitlam2; t3 <- Z$logittau3; t4 <- Z$logittau4
} else {
  mu <- Z$mu; l2 <- Z$lam2;t3 <- Z$tau3; t4 <- Z$tau4
}
x <- seq(log10(3), max(n), by=.01) # Vector of sample sizes for plotting the predictions of the
# regressions during the power optimization scheme.

"ofunc" <- function(pwr, type=c("press", "wolf"), v=NA, w=NA) {
  L  <- lm(v~n+I(n^pwr), weights=w)
  yl <- predict(L, newdata=data.frame(n=x))
  yp <- fitted(L); yr <- residuals(L)
  ht <- hatvalues(L); press <- sum(c(yr/(1-ht))^2)
  uv <- data.frame(u=lmomco::pp(yp, sort=FALSE), v=lmomco::pp(yr, sort=FALSE))
  wolf <- copBasic::wolfCOP(para=uv, as.sample=TRUE)
  print(c(pwr, press, wolf))
  lines(10^x, yl, col="lightgreen", lwd=0.5)
  ifelse(type == "press", return(press), return(wolf))
}
message("power   PRESS  WOLF")

y <- mu # MEAN LOGIT OF WOLFF SIGMA DISTRIBUTION
plot(10^n, y, log="x", col="grey50", xlab="Sample size", las=1,
     ylab="Mu of logit transform of Schweizer–Wolff Sigma distribution")
suppressWarnings( rtmu <- optim(INIT.PWR, ofunc, type=OPTIM.TYPE, v=y, w=Z$wgts) )
kmu <- rtmu$par
O   <- lm(y~n           , weights=Z$wgts); yp <- predict(O, newdata=data.frame(n=x))
lines(10^x, yp, col="salmon2", lwd=1)
Lmu <- lm(y~n+I(n^kmu), weights=Z$wgts);
yp  <- predict(Lmu, newdata=data.frame(n=x))
lines(10^x, yp, col="darkgreen", lwd=2)
G   <- gam(y~s(n)       , weights=Z$wgts); yp <- predict(G, newdata=data.frame(n=x))
lines(10^x, yp, col="orchid", lwd=3, lty=2)
print("-------------------------------------------")

y <- l2 # LAMBDA2 LOGIT OF WOLFF SIGMA DISTRIBUTION
plot(10^n, y, log="x", col="grey50", xlab="Sample size", las=1,
     ylab="Lscale of logit transform of Schweizer–Wolff Sigma distribution")
suppressWarnings( rtl2 <- optim(INIT.PWR, ofunc, type=OPTIM.TYPE, v=y, w=Z$wgts) )
kl2 <- rtl2$par
O   <- lm(y~n           , weights=Z$wgts); yp <- predict(O, newdata=data.frame(n=x))
lines(10^x, yp, col="salmon2", lwd=1)
Ll2 <- lm(y~n+I(n^kl2), weights=Z$wgts);
yp  <- predict(Ll2, newdata=data.frame(n=x))
lines(10^x, yp, col="darkgreen", lwd=2)
G   <- gam(y~s(n)       , weights=Z$wgts); yp <- predict(G, newdata=data.frame(n=x))
lines(10^x, yp, col="orchid", lwd=3, lty=2)
print("-------------------------------------------")

y <- t3 # TAU3 LOGIT OF WOLFF SIGMA DISTRIBUTION
plot(10^n, y, log="x", col="grey50", xlab="Sample size", las=1,
     ylab="Lskew of logit transform of Schweizer–Wolff Sigma distribution")
suppressWarnings( rtt3 <- optim(INIT.PWR, ofunc, type=OPTIM.TYPE, v=y, w=Z$wgts) )
kt3 <- rtt3$par
O   <- lm(y~n           , weights=Z$wgts); yp <- predict(O, newdata=data.frame(n=x))
lines(10^x, yp, col="salmon2", lwd=1)
Lt3 <- lm(y~n+I(n^kt3), weights=Z$wgts);
yp  <- predict(Lt3, newdata=data.frame(n=x))
lines(10^x, yp, col="darkgreen", lwd=2)
G   <- gam(y~s(n)       , weights=Z$wgts); yp <- predict(G, newdata=data.frame(n=x))
lines(10^x, yp, col="orchid", lwd=3, lty=2)
print("-------------------------------------------")

y <- t4 # TAU4 LOGIT OF WOLFF SIGMA DISTRIBUTION
plot(10^n, y, log="x", col="grey50", xlab="Sample size", las=1,
     ylab="Lskew of logit transform of Schweizer–Wolff Sigma distribution")
suppressWarnings( rtt4 <- optim(INIT.PWR, ofunc, type=OPTIM.TYPE, v=y, w=Z$wgts) )
kt4 <- rtt4$par
O   <- lm(y~n           , weights=Z$wgts); yp <- predict(O, newdata=data.frame(n=x))
lines(10^x, yp, col="salmon2", lwd=1)
Lt4 <- lm(y~n+I(n^kt4), weights=Z$wgts);
yp  <- predict(Lt4, newdata=data.frame(n=x))
lines(10^x, yp, col="darkgreen", lwd=2)
G   <- gam(y~s(n)       , weights=Z$wgts); yp <- predict(G, newdata=data.frame(n=x))
lines(10^x, yp, col="orchid", lwd=3, lty=2)
print("-------------------------------------------")


message("*****************************************************************************************")
message("BYLOGIT=", BYLOGIT); txt <- ifelse(BYLOGIT, "logit_", "      ")
message(  txt, "Lambda1 : ", paste(round(c(coefficients(Lmu), kmu), digits=8), collapse=", "), 8)
message(  txt, "Lambda2 : ", paste(round(c(coefficients(Ll2), kl2), digits=8), collapse=", "), 8)
message(  txt, "___Tau3 : ", paste(round(c(coefficients(Lt3), kt3), digits=8), collapse=", "), 8)
message(  txt, "___Tau4 : ", paste(round(c(coefficients(Lt4), kt4), digits=8), collapse=", "), 8)
message("*****************************************************************************************")
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
"wolfCOPtest" <-
function(x, y, asuv=FALSE, aslist=TRUE, bylogit=TRUE, dtype="gno",
               na.rm=TRUE, digits=6, probs=c(0.90, 0.95, 0.98, 0.99, 0.995), ...) {
  # The probs are quantile levels of the sigma to report, and these are useful to check against the
  # simulations but also to produce these as critical values should the user be interested in
  # these as well as the p-value.
  if(is.null(probs))      probs <- 0.95 # 95th percentile or rather the 5-percent critical value (upper tail).
  if( length(probs) == 0) probs <- 0.95 # 95th percentile or rather the 5-percent critical value (upper tail).

  if(length(x) == 1) { # If x is just one value, then it is treated as the Schweizer-Wolff Sigma
    rwolf <- x[1]; lwolf <- log(rwolf/(1-rwolf)); n <- ngiven <- y[1] # and the sample size is in y[1]
    if(n < 9) {
      warning("sample size is <9; so, will treat it as if it were 9 and preserve original sample by negation")
      ngiven <- -n; n <- 9
    }
  } else {
    if(length(x) != length(y)) {
      warning("length x != length y, returning NULL")
      return(NULL)
    }
    uv <- data.frame(u=x, v=y)
    if(na.rm) uv <- uv[complete.cases(uv),]
    n <- ngiven <- nrow(uv) # sample size
    if(n < 9) {
      warning("sample size is <9; so, will treat it as if it were 9 and preserve original sample by negation")
      ngiven <- -n; n <- 9
    }
    if(! asuv) { # if true, then the user has provided the paired observations of probability
      uv[,1] <- lmomco::pp(uv[,1], sort=FALSE, ...)
      uv[,2] <- lmomco::pp(uv[,2], sort=FALSE, ...)
    }

    rwolf <- copBasic::wolfCOP(para=uv, as.sample=TRUE) # Schweizer-Wolff Sigma : wolf in (0,1)
    lwolf <- log(rwolf / (1 - rwolf)) # logit transform of the Sigma
  }

  if(bylogit) {
    # Nonlinear regression coefficients computed PRESS minimization of residuals for the
    # exponent on log10(sample size) term. The regressions come from simulation of the Sigma
    # distribution (its logit) assuming the Independence copula.
    mucoe <- c(-0.14855834, -1.07395761, 1.18720595, -1.16858)
    l2coe <- c(0.11568377, 0.00223519, 0.10275754, -1.7678)
    t3coe <- c(0.05490461, 0.01096357, 0.179919, -1.688)
    t4coe <- c(0.13327184, -0.00291816, 0.0380223, -3.48)
  } else {
    # Nonlinear regression coefficients computed PRESS minimization of residuals for the
    # exponent on log10(sample size) term. The regressions come from simulation of the Sigma
    # distribution not using its logit assuming the Independence copula.
    mucoe <- c(-2.71311403, 0.14816878, 3.05730872, -0.252656258)
    l2coe <- c( 0.13046609, -0.19426932, 0.11774448, 1.238258)
    t3coe <- c( 0.2185038, 0.00103145, 0.00866943, -2.1848)
    t4coe <- c( 0.1653274, -0.00117758, -0.02007299, -2.4728)
  }

  # Apply the regressions using the hardwired coefficients herein
  mu    <- mucoe[1] + mucoe[2] * log10(n) + mucoe[3] * log10(n)^mucoe[4] # Mean    (Lambda1)
  l2    <- l2coe[1] + l2coe[2] * log10(n) + l2coe[3] * log10(n)^l2coe[4] # Lambda2 (L-scale)
  t3    <- t3coe[1] + t3coe[2] * log10(n) + t3coe[3] * log10(n)^t3coe[4] # Tau3    (L-skew)
  t4    <- t4coe[1] + t4coe[2] * log10(n) + t4coe[3] * log10(n)^t4coe[4] # Tau4    (L-kurtosis)
  lmrs  <- c(mu, l2, t3, t4) # Tidy list of the Lmoments of the logit(Sigma) distribution
  para  <-  lmomco::lmom2par(lmomco::vec2lmom(c(mu, l2, t3, t4)), type=dtype) # Lmoment ratio diagram shows
  # Really close adherence to a generalized normal distribution (3-parameter log-normal)
  # and hence that distribution is chosen here.
  if(bylogit) {
    neps  <-               lmomco::par2cdf(lwolf, para, paracheck=FALSE)   # CDF of logit(SIGMA)
    quans <- 1 / (1 + exp(-lmomco::par2qua(probs, para, paracheck=FALSE))) # QUA of inverse logit --> quans are Sigmas
  } else { # The neps are the nonexceedance probabilities.
    neps  <-               lmomco::par2cdf(rwolf, para, paracheck=FALSE)   # CDF of SIGMA
    quans <-               lmomco::par2qua(probs, para, paracheck=FALSE)   # QUA no retransformation --> are Sigmas
  }

  quans     <- round(quans,     digits=digits) # Estimated upper-tail quantiles of Sigma distribution
  lmrs      <- round(lmrs,      digits=digits) # Lmoments of the logit(Sigma) distribution
  para$para <- round(para$para, digits=digits) # parameters of the fitted distribution
  rwolf     <- round(rwolf,     digits=digits) # Sigma as the test statistic
  lwolf     <- round(lwolf,     digits=digits) # logit(Sigma) as needed for p-value lookup

  # Because of the potential confusion with the logit() transform involved, let us have long names.
  names(lmrs)      <- c("mu_TEXT_sigmas",    "lscale_TEXT_sigmas",  "tau3_TEXT_sigmas", "tau4_TEXT_sigmas" )
  names(para$para) <- paste0("para", seq_len(length(para$para)), "_TEXT_sigmas")
  names(lmrs)      <- gsub("_TEXT_", ifelse(bylogit, "logit", "_"), names(lmrs     ))
  names(para$para) <- gsub("_TEXT_", ifelse(bylogit, "logit", "_"), names(para$para))

  quatxt           <- paste0("fit_f", gsub("\\.", "p", as.character(100 * probs)))
  names(quans)     <- quatxt

  pval <- round(1 - neps, digits=8); names(pval) <- "p.value"

  zz <- c(n, ngiven, rwolf, lwolf, pval, para$para, lmrs, quans)
  names(zz) <- c("sample_size_in_comps", "sample_size_given", "sigma", "logit_sigma", "p.value",
                 names(para$para), names(lmrs), quatxt)
  names(zz) <- gsub("_TEXT_", ifelse(bylogit, "logit", "_"), names(zz))
  if(aslist) {
    wz <- c(rwolf, lwolf); names(wz) <- c("sigma", "logit_sigma")
    zz <- list(sample_size_in_comps=n, sample_size_given=ngiven,
               sigma=wz, p.value=pval, distpara_by_lmoms=para$para)
    if(bylogit) {
      zz$lmoms_logit_sigma <- lmrs # L-moments of the logit(SIGMAS) distribution
    } else {
      zz$lmoms_sigma       <- lmrs # L-moments of the       SIGMAS  distribution
    }
    zz$sigma_quantiles <- quans # Put these last because this length of vector is mutable, and it
    # visually makes these better on the right side of aslist=FALSE (vector return), in particular.
  }
  return(zz)
}


# --------------------------------------------------------------------------------------------------
uv <- simCOP(n=62, cop=PSP, graphics=FALSE)
wolfCOPtest(uv[,1], uv[,2], asuv=TRUE, bylogit=TRUE)
wolfCOPtest(uv[,1], uv[,2], asuv=TRUE, bylogit=TRUE, aslist=FALSE)

uv <- simCOP(n=62, cop=P,   graphics=FALSE)
wolfCOPtest(uv[,1], uv[,2], asuv=TRUE, bylogit=TRUE, aslist=FALSE)
# --------------------------------------------------------------------------------------------------
stop()

# --------------------------------------------------------------------------------------------------
neps <- c("f90", "f95", "f98", "f99", "f99p5")
for(nep in neps) {
  q <- sapply(Z$n, function(n) { k <- wolfCOPtest(0, n, )$sigma_quantiles
                                 k[grep(paste0(nep, "$"), names(k))] })
  plot(Z[,nep], q, log="xy", lwd=0.8, col="grey50", main=nep,
       xlab="Empirically computed sigma quantile during simulation",
       ylab="Estimated sigma quantile from L-moments")
  abline(0,1)
  print(summary(Z[,nep] - q))
}
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
alpha <- 0.1
nsim <- 500
H <- NULL
ns <- seq(10, 100, by=2)
for(n in ns) {
  wolf   <- replicate(nsim, { uv <- simCOP(n=n, cop=P, graphics=FALSE);
                                   wolfCOP(para=uv, as.sample=TRUE) })
  wolfpv_sigma <- sapply(wolf, function(w) { wolfCOPtest(w, n, asuv=TRUE, bylogit=FALSE, aslist=FALSE)[5] })
  wolfpv_logit <- sapply(wolf, function(w) { wolfCOPtest(w, n, asuv=TRUE, bylogit=TRUE,  aslist=FALSE)[5] })
  message("Schweizer–Wolff Sigma(n=", n, ") type I error rate = ", sum(wolfpv_sigma < alpha) / nsim)
  message("Schweizer–Wolff Sigma(n=", n, ") type I error rate = ", sum(wolfpv_logit < alpha) / nsim)

  rhopv <- replicate(nsim, { uv <- simCOP(n=n, cop=P, graphics=FALSE);
                        cor.test(uv[,1], uv[,2], method="spearman")$p.value })
  message("         Spearman Rho(n=", n, ") type I error rate = ", sum(rhopv        < alpha) / nsim)
  H <- rbind(H, data.frame(n=n, wolfpv_sigma=sum(wolfpv_sigma < alpha) / nsim,
                                wolfpv_logit=sum(wolfpv_logit < alpha) / nsim,
                                rhopv=       sum(rhopv        < alpha) / nsim))
}


ylim <- max( abs( range( c(H$wolfpv_sigma, H$wolfpv_logit, H$rhopv, alpha)) ) ) - alpha
ylim <- c(alpha - ylim, ylim + alpha)
plot( H$n, H$wolfpv_sigma, ylim=ylim, pch=21, col="turquoise4", bg="turquoise1", lwd=0.6, las=1,
     xlab="Sample size", ylab="Type I error rate (rejection of NULL but NULL is TRUE [P copula parent])")
lines( H$n, H$wolfpv_logit, col="salmon3", lwd=2)
points(H$n, H$wolfpv_logit, pch=21, col="salmon4", bg="salmon1", lwd=0.6)
lines( H$n, H$rhopv, col="orchid3"); xlim <- par()$usr[1:2]
lines(xlim, rep(               alpha, 2), lty=3                  )
lines(xlim, rep(mean(H$wolfpv_sigma), 2), lty=2, col="turquoise4")
lines(xlim, rep(mean(H$wolfpv_logit), 2), lty=6, col="salmon3"   )
lines(xlim, rep(mean(H$rhopv       ), 2), lty=4, col="orchid3"   )
# --------------------------------------------------------------------------------------------------
stop()


# --------------------------------------------------------------------------------------------------
alpha <- 0.05
n <- 80
nsim <- 1000
set.seed(1)
wolf   <- replicate(nsim, { uv <- simCOP(n=n, cop=CIRCcop, graphics=FALSE);
                                 wolfCOP(para=uv, as.sample=TRUE) })
  wolfpv_sigma <- sapply(wolf, function(w) { wolfCOPtest(w, n, asuv=TRUE, bylogit=FALSE, aslist=FALSE)[5] })
  wolfpv_logit <- sapply(wolf, function(w) { wolfCOPtest(w, n, asuv=TRUE, bylogit=TRUE,  aslist=FALSE)[5] })
  message("Schweizer–Wolff Sigma(n=", n, ") type II error rate = ", 1 - sum(wolfpv_sigma < alpha)/nsim)
  message("Schweizer–Wolff Sigma(n=", n, ") type II error rate = ", 1 - sum(wolfpv_logit < alpha)/nsim)

set.seed(1)
rhopv <- replicate(nsim, { uv <- simCOP(n=40, cop=CIRCcop, graphics=FALSE);
                           cor.test(uv[,1], uv[,2], method="spearman")$p.value })
message("           Spearman Rho(n=", n, ") type II error rate = ", 1 - sum(rhopv < alpha)/nsim)
# --------------------------------------------------------------------------------------------------

