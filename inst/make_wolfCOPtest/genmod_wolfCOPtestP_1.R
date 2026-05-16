if(Sys.getenv("RSTUDIO") == "1") {
  # Automatic change directory to location of this script when RStudio is running.
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

library(copBasic)
library(lmomco  )
library(mgcv    )

# --------------------------------------------------------------------------------------------------
D <- read.table("mc_wolfPI/aa_allsimsREF.txt",  header=TRUE)
# Built by make_wolfCOPtestP_A.R script.
#for(f in list.files(pattern="^mc_wolfPI")) {
#  D <- rbind(D, read.table(f, header=TRUE))
#}
D <- D[order(D$n),]
D <- D[complete.cases(D),]
# write.table(D, file="aa_allsims.txt", sep="\t", row.names=FALSE, quote=FALSE)

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

files <- list.files(path="wolfCOP_EMPIRgrid_fast/", pattern = "wolfCOP_EMPIRgrid_fast")
AUX <- NULL
for(file in files) {
  df  <- read.table(paste0("wolfCOP_EMPIRgrid_fast/", file), sep="\t", header=TRUE)
  AUX <- rbind(AUX, df)
}
AUX <- AUX[AUX$n >= 8,]



AUX$wgts <- sqrt(AUX$nsim); AUX$wgts <- AUX$wgts/sum(AUX$wgts) * length(AUX$wgts)
AUX <- AUX[order(AUX$n),]
LMR <- aggregate(AUX, by=list(AUX$n), mean); jjj <- aggregate(AUX, by=list(AUX$n), sum )
LMR$nsim <- jjj$nsim; rm(jjj); LMR$Group.1 <- NULL # compute means and total up sims by sample size
AUX <- LMR
AUX$wgts <- sqrt(AUX$nsim); AUX$wgts <- AUX$wgts/sum(AUX$wgts) * length(AUX$wgts)

Z <- Z[Z$n <= 1000,]
Z <- merge(Z, AUX, all=TRUE)
Z <- Z
Z <- Z[order(Z$n, decreasing=TRUE),]

Z <- AUX


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

ifelse(BYLOGIT, ylim <- c(-4.2,0.5), ylim <- c(0,0.7))
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
ifelse(BYLOGIT, ylim <- c(0.12,0.26), ylim <- c(0,0.1))
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
ifelse(BYLOGIT, ylim <- c(0.10,0.3), ylim <- c(0.2,0.28))
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
ifelse(BYLOGIT, ylim <- c(0.12,0.20), ylim <- c(0.12,0.2))
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

message("Suggestion, go and source genmod_wolfCOPtestP_2.R")
# source("genmod_wolfCOPtestP_2.R") # This line and the previous is to leave a paper trail that the
# Pade Approximant prediction structure is built by this next script.

