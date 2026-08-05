setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(copBasic)
library( copula )
library(parallel)

# Here is a cross check of sorts between wolfCOPtest and the Independence test
# available in the copula package. We use a Clayton copula with various bivariate
# associations and will see for these circumstances that wolfCOPtest out performs
# copula::indepTest by just a small amount by rejecting the NULL more often.
OUFILE <- "zz_chck_wolfCOPtestCLsim10k.txt"
NSIM <- as.integer( 25000 ); ALPHA <- 0.05
myCOPTXT <- "CLcop"
myCOP    <- CLcop
myTAUS   <- seq(0, 0.5, by=0.02)
myPARAf  <- function(tau=NA) (2*tau)/(1-tau)
myPARA   <- NULL
myRTdigits <- 6
myCORES    <- 4
sample_sizes <- c(seq(5, 20, by=1), seq(22, 50, by=2), seq(55, 100, by=5))

IWfunc <- function(n) {
  J <- copula::indepTestSim(n, 2, N=100000); WW <- II <- NULL
  for(i in seq_len(NSIM)) {
    uv <- copBasic::rCOP(n, cop=myCOP, para=myPARA)
    I <- W <- NULL
    try(I <-   copula::indepTest(  uv,         J             )$pvalues[1])
    try(W <- copBasic::wolfCOPtest(uv[,1], uv[,2], asuV=TRUE )$p.value[1])
    if(is.null(I)) I <- NA
    if(is.null(W)) W <- NA
    II <- c(II, I); WW <- c(WW, W)
  }
  II <- II[! is.na(II)]; WW <- WW[! is.na(WW)]
                         Irejrate <-                sum(II < ALPHA) / length(II)
                                        Wrejrate <- sum(WW < ALPHA) / length(WW)
  return( paste0(n, ":", Irejrate, ":", Wrejrate) )
}

D <- NULL; apnd <- FALSE
for(tau in myTAUS) {
  message(paste0(Sys.time(), "; Tau=", tau))
                                        myPARA <- myPARAf( tau )
  rhoi <- round( rhoCOP(cop=myCOP, para=myPARA), digits=myRTdigits )
  taui <- round( tauCOP(cop=myCOP, para=myPARA), digits=myRTdigits )
  # -----------------------------------------------------------------------
  cl <- parallel::makeCluster(getOption("cl.cores", myCORES))
        parallel::clusterExport( cl, c("ALPHA", "NSIM", "myCOP", "myPARA"))
   R <- parallel::parSapply(     cl, sample_sizes, IWfunc)
        parallel::stopCluster(   cl   )
  # -----------------------------------------------------------------------
  d <- NULL
  for(i in seq_len(length(R)) ) {
    r <- as.numeric( unlist( strsplit(R[i], ":") ) )
    d <- rbind(d, data.frame(copula=myCOPTXT, nsim=NSIM, n=as.integer(r[1]), alpha=ALPHA,
       rho_given=NA, tau_given=tau, para=myPARA, rho=rhoi, tau=taui, Irejrate=r[2], Wrejrate=r[3]))
  }
  D <- rbind( D, d); print( d )
  write.table(D, file=OUFILE, sep="\t", row.names=FALSE, quote=FALSE, append=apnd, col.names=! apnd)
  apnd <- TRUE
  print(Sys.time())
}
message("Done")


