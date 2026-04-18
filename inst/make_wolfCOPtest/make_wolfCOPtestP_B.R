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

D <- read.table("aa_lapsims.txt",  header=TRUE)
for(f in list.files(pattern="^mc_wolfPI")) {
  D <- rbind(D, read.table(f, header=TRUE))
}
D <- D[order(D$n),]
D <- D[complete.cases(D),]
write.table(D, file="aa_allsims.txt", sep="\t", row.names=FALSE, quote=FALSE)

h <- aggregate(D, by=list(D$n), length)
print(paste0(h$Group.1[h$nsim == 2], collapse=","))



# --------------------------------------------------------------------------------------------------
WEIGHT_BY_NSIM <- FALSE
Z <- NULL # now compute weighted mean columns
for(k in sort(unique(D$n))) {
  y <- D[D$n == k,]
  z <- y[1,]
  #print(nrow(y))
  if(WEIGHT_BY_NSIM) {
    for(i in 3:ncol(y)) z[,i] <- weighted.mean(y[,i], y$nsim)
  } else {
    strs <- c("f90", "f95", "f98", "f99", "f99p5")
    for(i in seq_len(length(strs))) z[,i] <- weighted.mean(y[,strs[i]], y$nsim)
    strs <- c("mu",   "var",  "lam2", "tau3", "tau4", "tau5")
    vtrs <- c("muse", "l2se", "l2se", "t3se", "t4se", "t5se")
    for(i in seq_len(length(strs))) z[,strs[i]] <- weighted.mean(y[,strs[i]], 1/y[,vtrs[i]])
    strs <- c("logitmu",   "logitvar",  "logitlam2", "logittau3", "logittau4", "logittau5")
    vtrs <- c("logitmuse", "logitl2se", "logitl2se", "logitt3se", "logitt4se", "logitt5se")
    for(i in seq_len(length(strs))) z[,strs[i]] <- weighted.mean(y[,strs[i]], 1/y[,vtrs[i]])
  }
  z[,1] <- sum(y$nsim)
  z[,2] <- y$n[1]
  Z <- rbind(Z, z)
}
Z <- Z[is.finite(Z$logitmu),]
Z$wgts <- sqrt(Z$nsim); Z$wgts <- Z$wgts/sum(Z$wgts) * length(Z$wgts)
row.names(Z) <- NULL
# --------------------------------------------------------------------------------------------------
#Z <- Z[Z$tau3 < 0.235,]
#Z <- Z[Z$tau4 < 0.18,]


plotlmrdia(lmrdia(), xlim=c(0.05, 0.3), ylim=c(0.1, 0.20), empty=TRUE)
points(Z$tau3,      Z$tau4,       pch=21, col="turquoise4", bg="turquoise1", lwd=0.6, cex=sqrt(Z$wgts))
points(Z$logittau3, Z$logittau4,  pch=21, col="salmon4",    bg="salmon1",    lwd=0.6, cex=sqrt(Z$wgts))
plotlmrdia(lmrdia(), add=TRUE, nopoints=TRUE, autolegend=TRUE, xleg="topleft")


# --------------------------------------------------------------------------------------------------
OPTIM.TYPE <- "press"
INIT.PWR <- -0.01 # initial parameter for the optimization scheme, being on the right side zero matters
BYLOGIT <- TRUE # run this script with both TRUE and FALSE and harvest the coefficients of the
# regressions and then emplace them in the wolfCOPtest function.

if(BYLOGIT) {
  #Z <- Z[Z$tau3 < 0.235,]
  #Z <- Z[Z$tau4 < 0.18, ]
}

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
x <- seq(log10(7), log10(10000), by=.05) # Vector of sample sizes for plotting the predictions of the
# regressions during the power optimization scheme.

"ofunc" <- function(pwr, type=c("press", "wolf"), v=NA, w=NA) {
  L  <- lm(v~n+I(n^pwr), weights=w)
  yl <- predict(L, newdata=data.frame(n=x))
  yp <- fitted(L); yr <- residuals(L)
  ht <- hatvalues(L); press <- sum(c(yr/(1-ht))^2)
  uv <- data.frame(u=lmomco::pp(yp, sort=FALSE), v=lmomco::pp(yr, sort=FALSE))
  #wolf <- copBasic::wolfCOP(para=uv, as.sample=TRUE)
  #print(c(pwr, press, wolf))
  print(c(pwr, press))
  lines(10^x, yl, col="lightgreen", lwd=0.5)
  ifelse(type == "press", return(press), return(wolf))
}
message("power   PRESS  WOLF")

ifelse(BYLOGIT, ylim <- c(-4,0.5), ylim <- c(0,0.7))
y <- mu # MEAN LOGIT OF WOLFF SIGMA DISTRIBUTION
plot(10^n, y, log="x", col="grey50", xlab="Sample size", las=1, xlim=10^range(x), ylim=ylim,
     ylab="Mu of logit transform of Schweizer–Wolff Sigma distribution")
suppressWarnings( rtmu <- optim(-1, ofunc, type=OPTIM.TYPE, v=y, w=Z$wgts) )
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
ifelse(BYLOGIT, ylim <- c(0.1,0.3), ylim <- c(0,0.1))
plot(10^n, y, log="x", col="grey50", xlab="Sample size", las=1, xlim=10^range(x), ylim=ylim,
     ylab="Lscale of logit transform of Schweizer–Wolff Sigma distribution")
suppressWarnings( rtl2 <- optim(-2, ofunc, type=OPTIM.TYPE, v=y, w=Z$wgts) )
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
ifelse(BYLOGIT, ylim <- c(0.05,0.3), ylim <- c(0.2,0.28))
plot(10^n, y, log="x", col="grey50", xlab="Sample size", las=1, xlim=10^range(x), ylim=ylim,
     ylab="Lskew of logit transform of Schweizer–Wolff Sigma distribution")
suppressWarnings( rtt3 <- optim(-1, ofunc, type=OPTIM.TYPE, v=y, w=Z$wgts) )
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
ifelse(BYLOGIT, ylim <- c(0.1,0.25), ylim <- c(0.12,0.2))
plot(10^n, y, log="x", col="grey50", xlab="Sample size", las=1, xlim=10^range(x), ylim=ylim,
     ylab="Lkurtosis of logit transform of Schweizer–Wolff Sigma distribution")
suppressWarnings( rtt4 <- optim(-2, ofunc, type=OPTIM.TYPE, v=y, w=Z$wgts) )
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
"wolfCOPtest2" <-
function(x, y, asuv=FALSE, aslist=TRUE, bylogit=TRUE, dtype="gno",
               na.rm=TRUE, digits=6, probs=c(0.90, 0.95, 0.98, 0.99, 0.995), ...) {
  # The probs are quantile levels of the sigma to report, and these are useful to check against the
  # simulations but also to produce these as critical values should the user be interested in
  # these as well as the p-value.
  if(is.null(probs))      probs <- 0.95 # 95th percentile or rather the 5-percent critical value (upper tail).
  if( length(probs) == 0) probs <- 0.95 # 95th percentile or rather the 5-percent critical value (upper tail).

  # -------------------------------
  dtype <- tolower(dtype)
  suppressWarnings( have_dtype <- lmomco::dist.list(dtype) )
  if(length(grep("does not match", as.character(have_dtype))) == 1) { # "The given type argument does not match a distribution") {
    warning("requested dtype='", dtype, "' does not exist in lmomco::dist.list() function")
    return(NULL)
  }
  if(have_dtype >= 5) {
    warning("function only supports up to 4-parameter distributions ")
  }
  if(have_dtype == "gld") {
    warning("function does support first four Lmoments, but a fifth would be needed for dtype='gld'")
  }
  # -------------------------------

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
    mucoe <- c(-0.01524626, -1.09592416, 1.07797559, -1.295410168)
    l2coe <- c(0.13651546, -0.00275936, 0.08849981, -2.1843758)
    t3coe <- c(0.07385001, 0.00678137, 0.16538229, -1.85781258)
    t4coe <- c(0.12225675, 0.00041096, 0.04507729, -2.63758)
  } else {
    # Nonlinear regression coefficients computed PRESS minimization of residuals for the
    # exponent on log10(sample size) term. The regressions come from simulation of the Sigma
    # distribution not using its logit assuming the Independence copula.
    mucoe <- c(-5.13011976, 0.17965821, 5.43973795, -0.14481258)
    l2coe <- c(0.1117815, -0.07903764, 0.02014877, 1.71258)
    t3coe <- c(0.22033934, 0.00035634, 0.00707766, -2.8248)
    t4coe <- c(0.16112678, -0.00011312, -0.01705013, -2.9168)
  }

 # if(n <= 50) dtype <- "gno"

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
  qT <- sapply(Z$n, function(n) { k <- wolfCOPtest(0, n, bylogit=TRUE, type="nor")$sigma_quantiles
                                 k[grep(paste0(nep, "$"), names(k))] })
  qF <- sapply(Z$n, function(n) { k <- wolfCOPtest(0, n, bylogit=FALSE, type="nor")$sigma_quantiles
                                 k[grep(paste0(nep, "$"), names(k))] })
  df <- data.frame(emp=c(Z[,nep], Z[,nep]), q=c(qT, qF),
                  bg=c(rep("salmon1", length(qT)), rep("turquoise1", length(qF))),
                  col=c(rep("salmon4", length(qT)), rep("turquoise4", length(qF)))
                  )
  df <- df[sample(seq_len(nrow(df)), nrow(df)),]
  plot(qnorm(df$emp), qnorm(df$q),
       pch=21, col=df$col, bg=df$bg, lwd=0.8, main=nep,
       xlab="Empirically computed sigma quantile during simulation",
       ylab="Estimated sigma quantile from L-moments")
  #points(df$emp, df$qF, pch=21, col=df$col, bg="turquoise1", lwd=0.6)
  abline(0,1)
}
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
alpha <- 0.05
nsim <- 20000
H <- NULL
ns <- seq(10, 200, by=1)
for(i in seq_len(length(ns))) {
  n <- ns[i]
  m <- as.integer(nsim/n^0.5)
  wolf   <- replicate(m, { wolfCOP(para=data.frame(matrix(runif(n*2), ncol=2)), as.sample=TRUE) })
  wolfpv_sigma <- sapply(wolf, function(w) { wolfCOPtest2(w, n, asuv=TRUE, bylogit=TRUE, aslist=FALSE, dtype="pe3")[5] })
  wolfpv_logit <- sapply(wolf, function(w) { wolfCOPtest2(w, n, asuv=TRUE, bylogit=TRUE,  aslist=FALSE, dtype="gno")[5] })
  message("Schweizer–Wolff Sigma(n=", n, "/", m, ") type I error rate = ",
   sum(wolfpv_sigma < alpha) / m)
  message("Schweizer–Wolff Sigma(n=", n, "/", m, ") type I error rate = ",
   sum(wolfpv_logit < alpha) / m)

  #rhopv <- replicate(nsim, { uv <- simCOP(n=n, cop=P, graphics=FALSE);
  #                      cor.test(uv[,1], uv[,2], method="spearman")$p.value })
  #message("         Spearman Rho(n=", n, ") type I error rate = ", sum(rhopv        < alpha) / m)
  H <- rbind(H, data.frame(n=n, wolfpv_sigma=sum(wolfpv_sigma < alpha) / m,
                                wolfpv_logit=sum(wolfpv_logit < alpha) / m))
                                #rhopv=       sum(rhopv        < alpha) / m))

#ylim <- max( abs( range( c(H$wolfpv_sigma, H$wolfpv_logit, alpha)) ) ) - alpha
ylim <- c(alpha - 0.03, alpha + 0.03)
plot( H$n, H$wolfpv_sigma, ylim=ylim, type="n", las=1,
     xlab="Sample size", ylab="Type I error rate (rejection of NULL but NULL is TRUE [P copula parent])")
#lines( H$n, H$wolfpv_sigma, col="turquoise3", lwd=1)
#lines( H$n, H$wolfpv_logit, col="salmon3",    lwd=1)
points(H$n, H$wolfpv_sigma, pch=21, col="turquoise4", bg="turquoise1", lwd=0.6, cex=2)
points(H$n, H$wolfpv_logit, pch=21, col="salmon4",    bg="salmon1",    lwd=0.6, cex=1)
#lines( H$n, H$rhopv, col="orchid3")
xlim <- par()$usr[1:2]
lines(xlim, rep(               alpha, 2), lty=3                  )
try(lines(xlim, rep(mean(H$wolfpv_sigma[i:(i-10)], na.rm=TRUE), 2), lty=2, col="turquoise4"), silent=TRUE)
try(lines(xlim, rep(mean(H$wolfpv_logit[i:(i-10)], na.rm=TRUE), 2), lty=6, col="salmon3"   ), silent=TRUE)
#lines(xlim, rep(mean(H$rhopv       ), 2), lty=4, col="orchid3"   )

lines(H$n, sapply(seq_len(length(H$wolfpv_sigma)),function(i) mean(H$wolfpv_sigma[1:i])), col="turquoise4", lwd=2)
lines(H$n, sapply(seq_len(length(H$wolfpv_logit)), function(i) mean(H$wolfpv_logit[1:i])), col="salmon4", lwd=2)
}


ylim <- max( abs( range( c(H$wolfpv_sigma, H$wolfpv_logit, H$rhopv, alpha)) ) ) - alpha
ylim <- c(alpha - 0.03, alpha + 0.03)
plot( H$n, H$wolfpv_sigma, ylim=ylim, type="n", las=1,
     xlab="Sample size", ylab="Type I error rate (rejection of NULL but NULL is TRUE [P copula parent])")
lines( H$n, H$wolfpv_sigma, col="turquoise3", lwd=2)
lines( H$n, H$wolfpv_logit, col="salmon3",    lwd=2)
points(H$n, H$wolfpv_sigma, pch=21, col="turquoise4", bg="turquoise1", lwd=0.6)
points(H$n, H$wolfpv_logit, pch=21, col="salmon4",    bg="salmon1",    lwd=0.6)
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

