"wolfCOP_C.n" <- function(uv) {
  uv <- as.matrix(uv)
  #ab <- cbind(rank(uv[,1])/nrow(uv), rank(uv[,2])/nrow(uv))
  12*sum(abs((copula::C.n(uv, uv, smoothing="none") - uv[,1]*uv[,2])))/nrow(uv)
  #12*sum(abs((copula::C.n(ab, uv, smoothing="none") - ab[,1]*ab[,2])))/nrow(uv)
}
#uv <- simCOP(1000, cop=P)
#wolfCOP_C.n(uv)
#wolfCOP(para=uv, as.sample=TRUE)


"wolfCOPlmr" <- function(n, m=1000, nmom=5, logit=FALSE) {
  w <- 12*sapply(seq_len(m), function(i) {
           uv <- matrix(runif(n*2), ncol=2)
           sum(abs((copula::C.n(uv, uv, smoothing="none") - uv[,1]*uv[,2])))/n })
  ifelse(logit, return(Lmoments::Lcoefs(log(w/(1-w)), rmax=nmom)),
                return(Lmoments::Lcoefs(w,            rmax=nmom)))
}

"wolfCOPvec" <- function(n, m=1000, logit=FALSE) {
  w <- 12*sapply(seq_len(m), function(i) {
           uv <- matrix(runif(n*2), ncol=2)
           sum(abs((copula::C.n(uv, uv, smoothing="none") - uv[,1]*uv[,2])))/n })
  ifelse(logit, return(log(w/(1-w))), return(w))
}

stop("SAFE STOP:")

nsim <- 1000; r <- 100; ncore <- 10 # number simulations, number replicates, number CPU cores
ns  <- ((1:10) * 1000 ) + 1
ns <- 20000
SMR <- NULL
for(n in ns) {
  tm1 <- as.numeric( Sys.time() )
  message("CLUSTER: n=", n, " @ ", date(), appendLF=FALSE)
  cl <- parallel::makeCluster(getOption("cl.cores", ncore))        # make CPU cluster
        parallel::clusterExport(cl, c("n", "nsim", "wolfCOPlmr") ) # make variables visible to cluster
  ww <- parallel::parLapply(    cl=cl, seq_len(r), fun=function(i) wolfCOPlmr(n,nsim,logit=TRUE) )
        parallel::stopCluster(  cl)                                # stop the cluster
  tm2 <- as.numeric( Sys.time() )
  message(" DONE in ", round((tm2-tm1)/60, digits=1), " minutes")

  lmr <- as.data.frame( t( sapply(ww, function(k) as.data.frame(k) ) ) ) # manipulations toward a
  for(i in seq_len(ncol(lmr))) lmr[,i] <- unlist(lmr[,i])    # simple data frame of the results
  lmr$nsim <- m; lmr$n <- n; lmr$var <- (sqrt(pi)*lmr$L2)^2  # with some additional columns
  SMR <- rbind(SMR, lmr)                                     # and build up table for sample size
}
SMR <- SMR[,c("nsim", "n", "L1", "var", "L2", "tau3", "tau4", "tau5")] # rearrange columns
nm <- names(SMR) # prepare for column renaming
nm[3:8] <- c("logitmu", "logitvar", "logitlam2", "logittau3", "logittau4", "logittau5")
names(SMR) <- nm; rm(nm) # rename the columns
LMR <- aggregate(SMR, by=list(SMR$n), mean); jjj <- aggregate(SMR, by=list(SMR$n), sum )
LMR$nsim <- jjj$nsim; rm(jjj); LMR$Group.1 <- NULL # compute means and total up sims by sample size

write.table(LMR, file="wolfCOPaux_D.txt", sep="\t", row.names=FALSE, quote=FALSE)

