if(Sys.getenv("RSTUDIO") == "1") {
  # Automatic change directory to location of this script when RStudio is running.
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
} else {
  print("RStudio is not running")
}

library(copBasic)
library(lmomco)

cols <- rev( hcl.colors(9, palette="Roma") )
d <- 16
probs <- seq(0.001, 0.999, by=0.001)
probs <- pnorm(seq(-4.26, 4.26, by=0.02))
FCORE <- "smlsam_mc_wolfPI"
stop("SAFE STOP")
H <- NULL
for(k in 9) {#seq_len(9)) {
  ALLUNI <- NULL

  ns  <- 3:40
  UNI <- NULL
  for(n in ns) {
    OUTFILE <- paste0(FCORE, "_", n, ".txt")
    D <- NULL
    nsim <- 100000*k
    #wolfP <- replicate(nsim, { copBasic::wolfCOP(para=simCOP(n=n, cop=P, graphics=FALSE),
    #                                               as.sample=TRUE) })
    wolfP <- rep(NA, nsim)
    for(i in seq_len(nsim)) {
      wolfP[i] <- copBasic::wolfCOP(para=as.data.frame(matrix(runif(n*2), ncol=2)), as.sample=TRUE)
      #print(c(nsim, n, i, wolfP[i]))
    }
    #print(summary(wolfP))
    wolfemp <- round(quantile(wolfP, probs=probs, type=4, digits=d), digits=d)
    names(wolfemp) <- NULL
    df <- data.frame(nsim=nsim, n=n, probs=probs, wolfemp=round(wolfemp, digits=d))
    UNI <- c( UNI,   length(unique(df$wolfemp)) )
    print(c(nsim, n, length(unique(df$wolfemp))))
    D <- NULL
    for(wolfemp in unique(df$wolfemp)) {
      tmp <- df[df$wolfemp == wolfemp,]
      if(nrow(tmp) == 1) {
        D <- rbind(D, tmp)
      } else {
        kmp <- tmp[c(1, nrow(tmp)),]
           kmp$wolfemp[1] <-   kmp$wolfemp[1] - .Machine$double.eps
        if(kmp$wolfemp[1] < 0) kmp$wolfemp[1] <- 0
           kmp$wolfemp[2] <-   kmp$wolfemp[2] + .Machine$double.eps
        if(kmp$wolfemp[2] > 1) kmp$wolfemp[2] <- 1
        D <- rbind(D, kmp)
      }
    }
    D <- rbind(D, data.frame(nsim=rep(nsim, 2), n=rep(n, 2), probs=c(0,1), wolfemp=c(0,1)))
    D <- D[order(D$probs),]
    D$wolfemp <- sort(D$wolfemp)
    row.names(D) <- NULL
    write.table(D, file=OUTFILE, row.names=FALSE, sep="\t", quote=FALSE)
    #plot(D$probs, D$emp, type="l", main=n, col="salmon3")
  }
  if(k == 1) {
    plot( ns, UNI, type="b", col=cols[k], lwd=(9-k+1)/1.5,
         xlab="Sample size of Schweizer-Wolff",
         ylab="No. unique values Schweizer-Wolff given Independence")
  } else {
    lines(ns, UNI, type="b", col=cols[k], lwd=(9-k+1)/1.5,)
  }
  ALLUNI <- c(ALLUNI, UNI)
}


files <- list.files(pattern=FCORE)
A <- NULL
for(file in files) {
  A <- rbind(A, read.table(file, sep="\t", header=TRUE))
}
A <- A[order(A$n, A$probs),]
A$nsim <- NULL
# print(A$wolfemp[A$n == 5], 16)

wolfCOPtest_data_smlsam <- A
save(wolfCOPtest_data_smlsam, file="wolfCOPtest_data_smlsam.RData")

plot(qnorm(A$probs), A$wolfemp, type="l",
     xlab="Standard normal variate", ylab="Schweizer-Wolff sigma")
for(k in unique(A$n)) {
  readline(prompt=paste0("Next sample size (", k, "): "))
  para <- pargno(vec2lmom(wolfCOPtest(0,k )$lmoms_logit_sigma))
  quans <- 1 / (1 + exp(-lmomco::par2qua(probs, para, paracheck=FALSE)))
  lines(qnorm(probs), quans, col="red")
  para <- parpe3(vec2lmom(wolfCOPtest(0,k )$lmoms_logit_sigma))
  quans <- 1 / (1 + exp(-lmomco::par2qua(probs, para, paracheck=FALSE)))
  lines(qnorm(probs), quans, col="darkgreen", lwd=2)
}
# We see that the PE3 has the better upper tail performance almost immediately by about sample 9
# than the GNO.

