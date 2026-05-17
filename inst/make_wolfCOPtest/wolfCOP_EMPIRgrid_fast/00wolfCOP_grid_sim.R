library(copBasic)
library(parallel)
library(Lmoments)

n <- 5
uv <- as.matrix(matrix(runif(n *2), ncol=2))
gg <- EMPIRgrid(para=as.data.frame(uv), deluv=0.2, gridonly=TRUE); image(gg)
empCOP <- EMPIRgrid_fast(uv, ctype="1/n", gridonly=TRUE); image(empCOP, main="1/n")
empCOP <- EMPIRgrid_fast(uv, ctype="bernstein", gridonly=TRUE); image(empCOP, main="Bernstein")

n <- 500
para <- list(cop=GHcop, para=c(12,15), breve=0.5)
uv <- as.matrix(copBasic::simCOP(n, cop=breveCOP, para=para, graphics=TRUE)) # datasetModelMatrix
empCOP <- EMPIRgrid_fast(uv, gridonly=TRUE)
image(empCOP)

stop()
# --------------------------------------------------------------------------------------------------

"wolfCOPvec" <-
function(n=10, m=1000, logit=FALSE, aslmoms=FALSE, nmom=5, silent=TRUE) {
  PiCOP <- matrix(data=0, nrow=n+1, ncol=n+1)
  ix <- 2:(n+1)
  for(i in ix) PiCOP[i, ix] <- (i-1)*(ix-1)
               PiCOP <- PiCOP / (n*n)

  f <- 12 / (n^2 - 1)
  w <- rep(NA, m)
  for(i in seq_len(m)) {
    itxt <- as.character(i)
    if(! silent & length(grep("[50]0$", itxt)) == 1) {
      message(paste0(itxt, "-"), appendLF=FALSE)
      if(length(grep("[05]00$", itxt)) == 1) message("")
    }

    w[i] <- f * sum(abs(EMPIRgrid_fast(matrix(runif(2*n), ncol=2), ctype="1/n", gridonly=TRUE) - PiCOP))
  }
  if(! silent ) message("done")
  if(! aslmoms) ifelse(logit, return(log(w/(1-w))), return(w))
  ifelse(logit, return(Lmoments::Lcoefs(log(w/(1-w)), rmax=nmom)),
                return(Lmoments::Lcoefs(w,            rmax=nmom)))

}

a <- wolfCOPvec(n=5, silent=FALSE, aslmoms=TRUE, nmom=5)

stop()

replicates <- 10
ncores <- 4
nsim  <- 1000
ns <- 2+ as.integer(10^seq(3,4, by=.1))[-1]
ns <- rep(10001, 10)
SMR <- NULL
for(n in ns) {
  tm1 <- as.numeric( Sys.time() )
  message("CLUSTER: n=", n, " @ ", date(), appendLF=FALSE)
  str <- c("n", "nsim", "wolfCOPvec", "EMPIRgrid_fast")
  cl <- parallel::makeCluster(getOption("cl.cores", ncores)) # make CPU cluster
        parallel::clusterExport(cl, str )                    # make variables visible to cluster
  ww <- parallel::parLapply(    cl=cl, seq_len(replicates),  # run the cluster with subreplication
          fun=function(i) { wolfCOPvec(n, nsim, logit=TRUE, aslmoms=TRUE) } ) # simulations on nsims
        parallel::stopCluster(  cl )                         # stop the cluster
  tm2 <- as.numeric( Sys.time() )
  message(" DONE in ", round((tm2-tm1)/60, digits=1), " minutes")

  lmr <- as.data.frame( t( sapply(ww, function(k) as.data.frame(k) ) ) ) # manipulations toward a
  for(i in seq_len(ncol(lmr))) lmr[,i] <- unlist(lmr[,i])    # simple data frame of the results
  lmr$nsim <- nsim; lmr$n <- n; lmr$var <- (sqrt(pi)*lmr$L2)^2  # with some additional columns
  SMR <- rbind(SMR, lmr)                                     # and build up table for sample size
}
SMR <- SMR[,c("nsim", "n", "L1", "var", "L2", "tau3", "tau4", "tau5")] # rearrange columns
nm <- names(SMR) # prepare for column renaming
nm[3:8] <- c("logitmu", "logitvar", "logitlam2", "logittau3", "logittau4", "logittau5")
names(SMR) <- nm; rm(nm) # rename the columns
SMR <- SMR[is.finite(SMR$logitmu),]
LMR <- aggregate(SMR, by=list(SMR$n), mean); jjj <- aggregate(SMR, by=list(SMR$n), sum )
LMR$nsim <- jjj$nsim; rm(jjj); LMR$Group.1 <- NULL # compute means and total up sims by sample size

write.table(LMR, file="wolfCOP_EMPIRgrid_fastA.txt", sep="\t", row.names=FALSE, quote=FALSE)

